function within_sz_severity
%% ======================= CONFIG =======================
outCsv      = '../data/SN_counts/spike_counts_summary.csv';                     % spike counts with Duration_sec
%reportFile  = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-14_1505.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-20_1418.csv';
pairsOut    = '../output/eeg_visit_pairs_ALLTYPES_hasSz.csv';                   % audit output

maxGapDays        = 365;     % visit must be within this many days of EEG
useLogTransform   = true;    % log1p(spike rate) for modeling/plots
removeConstBinary = false;   % drop patients with no within-patient variation in HasSz across pairs

% EXCLUDE these epilepsy-type categories entirely
badTypes = lower([
    "Non-Epileptic Seizure Disorder"
    "Uncertain if Epilepsy"
    "Unknown or MRN not found"
    "" % empty / missing
]);

% KEEP these for analysis (exact labels expected in your data)
wantTypes = ["Focal","General"];
%wantTypes = ["General"];

%% ================== LOAD & PREP SPIKE SUMMARY ==================
S = readtable(outCsv, 'TextType','string', 'VariableNamingRule','preserve');  
reqS = {'EEG_Name','Patient','Session','Total_Spikes','Left_Spikes','Right_Spikes','Duration_sec'};
assert(all(ismember(reqS, S.Properties.VariableNames)), 'Spike summary missing required columns.');
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end
S.SpikeRate_perHour = S.Total_Spikes ./ (S.Duration_sec/3600);

%% ================== LOAD & PREP REPORT SHEET ===================
R = readtable(reportFile, 'TextType','string', 'VariableNamingRule','preserve');
reqR = {'patient_id','session_number','start_time_deid','visit_dates_deid','epilepsy_type','visit_hasSz'};
assert(all(ismember(reqR, R.Properties.VariableNames)), ...
    'Report file missing required columns. Found: %s', strjoin(R.Properties.VariableNames, ', '));
if ~isnumeric(R.('patient_id')),     R.('patient_id')     = double(str2double(string(R.('patient_id')))); end
if ~isnumeric(R.('session_number')), R.('session_number') = double(str2double(string(R.('session_number')))); end

% Parse EEG datetime
rawStart = string(R.('start_time_deid'));
okStart  = ~ismissing(rawStart) & strlength(strtrim(rawStart)) > 0;
R.EEG_Date = NaT(height(R),1);
try
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart)), 'InputFormat','M/d/yy HH:mm');
catch
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart))); % autodetect mixed formats
end

% Parse visit dates & hasSz (binary) into cell arrays
rawVisits  = string(R.('visit_dates_deid'));
rawHasSz   = string(R.('visit_hasSz'));
R.VisitDates = cell(height(R),1);
R.HasSzVec   = cell(height(R),1);
for i = 1:height(R)
    % Visit dates
    vstr = strtrim(rawVisits(i));
    if strlength(vstr) > 0 && vstr ~= "[]"
        try
            vcell = jsondecode(vstr); 
            vdt   = datetime(vcell, 'InputFormat','yyyy-MM-dd');
        catch
            toks = regexp(vstr, '\d{4}-\d{2}-\d{2}', 'match');
            vdt  = datetime(toks, 'InputFormat','yyyy-MM-dd');
        end
        R.VisitDates{i} = vdt(:);
    else
        R.VisitDates{i} = NaT(0,1);
    end

    % hasSz (0/1 list)
    hstr = strtrim(rawHasSz(i));
    if strlength(hstr) > 0 && hstr ~= "[]"
        try
            hv = jsondecode(char(hstr));
            R.HasSzVec{i} = double(hv(:));
        catch
            nums = regexp(hstr, '\d+', 'match');
            R.HasSzVec{i} = double(str2double(string(nums)));
        end
    else
        R.HasSzVec{i} = nan(0,1);
    end
end

%% ===== Per-patient rescue flag: HasSzAllZero (ALL visit_hasSz values are 0) =====
upids = unique(R.patient_id(~isnan(R.patient_id)));
hasAllZeroMap = containers.Map('KeyType','double','ValueType','logical');

for k = 1:numel(upids)
    pid = upids(k);
    rr = R(R.patient_id==pid, :);

    % Collect all hasSz entries across this patient's rows
    cells = rr.HasSzVec;                                % cell array, each is numeric vector or empty
    nonEmpty = ~cellfun(@isempty, cells);
    if any(nonEmpty)
        vals = vertcat(cells{nonEmpty});                % concatenate all visit-level 0/1s
        vals = vals(:);
        vals = vals(isfinite(vals));                    % drop NaNs just in case
    else
        vals = [];
    end

    % ALL-zeros rescue criterion
    hasAllZeroMap(pid) = (~isempty(vals) && all(vals==0));
end

nAllZero = sum(cell2mat(values(hasAllZeroMap)));
fprintf('Patients with ALL visit_hasSz == 0: %d\n', nAllZero);


%% ===== Build per-patient epilepsy_type and FILTER to Focal/General only =====
epi_str = string(R.epilepsy_type);
isEmpty = ismissing(epi_str) | strlength(strtrim(epi_str))==0;
Rt = sortrows(R(~isEmpty, {'patient_id','epilepsy_type'}), 'patient_id');
[uniq_pid, ia] = unique(Rt.patient_id, 'stable');
PerPatType = table(uniq_pid, Rt.epilepsy_type(ia), 'VariableNames', {'Patient','EpilepsyType'});

epi_norm_pat = lower(strtrim(PerPatType.EpilepsyType));
isBad  = ismember(epi_norm_pat, badTypes);
isWant = ismember(string(PerPatType.EpilepsyType), wantTypes);
validPatients = PerPatType.Patient(~isBad & isWant);

fprintf('Total patients in report: %d\n', numel(unique(R.patient_id(~isnan(R.patient_id)))));
fprintf('Patients kept (Focal/General only): %d\n', numel(validPatients));
fprintf('Patients excluded by type filter: %d\n', sum(~ismember(PerPatType.Patient, validPatients)));

% Restrict R and S to valid patients
R = R(ismember(R.patient_id, validPatients), :);
S = S(ismember(S.Patient,    validPatients), :);

%% ====== JOIN spike summary with report to get EEG date & visits & hasSz ======
JR = innerjoin(S, R, ...
    'LeftKeys',  {'Patient','Session'}, ...
    'RightKeys', {'patient_id','session_number'}, ...
    'RightVariables', {'EEG_Date','VisitDates','HasSzVec'});

% Keep rows with a valid EEG date and at least one visit date
JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) && all(~isnat(v)), JR.VisitDates), :);

%% ====== BUILD EEG–VISIT PAIRS (closest visit within maxGapDays) with ALL-ZERO rescue ======
pairs = table('Size',[0 10], ...
    'VariableTypes', {'double','double','string','datetime','double','datetime','double','logical','logical','string'}, ...
    'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour','Visit_Date','HasSz','RescueApplied','HasSzAllZero','EpilepsyType'});

% For fast type lookup
typeMap = containers.Map(double(PerPatType.Patient), cellstr(PerPatType.EpilepsyType));

for i = 1:height(JR)
    eegDate = JR.EEG_Date(i);
    vdates  = JR.VisitDates{i};
    hasV    = JR.HasSzVec{i};   % numeric vector (0/1) aligned to vdates
    pid     = JR.Patient(i);

    if isempty(vdates) || isempty(hasV), continue; end
    nAlign = min(numel(vdates), numel(hasV));
    vdates = vdates(1:nAlign);
    hasV   = hasV(1:nAlign);

    [dmin, idx] = min(abs(days(vdates - eegDate)));
    if isempty(dmin) || dmin > maxGapDays, continue; end

    thisHas = hasV(idx);
    rescueApplied = false;

    % ALL-zero rescue: if chosen hasSz missing/NaN and patient all-zero across visits, impute 0
    if ~(isfinite(thisHas) && (thisHas==0 || thisHas==1))
        if isKey(hasAllZeroMap, pid) && hasAllZeroMap(pid)
            thisHas = 0;
            rescueApplied = true;
        else
            continue % still missing => skip this pair
        end
    end

    % Attach epilepsy type label
    etype = "";
    if isKey(typeMap, pid), etype = string(typeMap(pid)); end

    pairs = [pairs; {pid, JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perHour(i), vdates(idx), thisHas, rescueApplied, ...
                     (isKey(hasAllZeroMap,pid) && hasAllZeroMap(pid)), etype}]; %#ok<AGROW>
end

fprintf('Built %d EEG–visit pairs (<=%d days). Rescued pairs: %d\n', ...
    height(pairs), maxGapDays, sum(pairs.RescueApplied));
if isempty(pairs)
    warning('No EEG–visit pairs were constructed.'); return
end

% Require ≥2 pairs per patient
counts   = groupcounts(pairs, 'Patient');
eligible = counts.Patient(counts.GroupCount >= 2);
pairs    = pairs(ismember(pairs.Patient, eligible), :);
fprintf('Eligible patients with ≥2 pairs: %d\n', numel(eligible));
if isempty(pairs), warning('Nothing remains after ≥2-pair filter.'); return; end

% Remove patients whose EEGs all matched the same visit date
[G, pidg] = findgroups(pairs.Patient);
nUniqueVisits = splitapply(@(d) numel(unique(d)), pairs.Visit_Date, G);
keepPatients = pidg(nUniqueVisits >= 2);
pairs = pairs(ismember(pairs.Patient, keepPatients), :);
if isempty(pairs), warning('Nothing left after visit-date filter.'); return; end

% Keep only Focal/General labels (defensive)
pairs = pairs(ismember(pairs.EpilepsyType, wantTypes), :);

%% ==================== MODELING PREP ====================
tbl = pairs;
tbl.Patient      = categorical(tbl.Patient);
tbl.EpilepsyType = categorical(string(tbl.EpilepsyType), wantTypes); % enforce order Focal, General
tbl.HasSz        = double(tbl.HasSz); % ensure numeric 0/1

% Optional: drop patients with constant HasSz across their pairs (no within-person contrast)
if removeConstBinary
    [Gc, pids2] = findgroups(tbl.Patient);
    varHas = splitapply(@(v) nanvar(v), tbl.HasSz, Gc);
    keepP2 = pids2(varHas > 0);
    tbl = tbl(ismember(tbl.Patient, keepP2), :);
    if isempty(tbl), warning('Nothing left after constant-HasSz filter.'); return; end
end

% Response variable
if useLogTransform
    tbl.SpikeRate_log = real(double(log1p(tbl.SpikeRate_perHour)));
    respVar = 'SpikeRate_log';
    yLabel  = 'log(1 + spike rate, spikes/hour)';
else
    tbl.SpikeRate_perHour = double(tbl.SpikeRate_perHour);
    respVar = 'SpikeRate_perHour';
    yLabel  = 'Spike rate (spikes/hour)';
end

ok = isfinite(tbl.(respVar)) & isfinite(tbl.HasSz) & ~isundefined(tbl.Patient) & ~isundefined(tbl.EpilepsyType);
tbl = tbl(ok,:);

% Need ≥2 pairs per patient (defensive after drops)
gc  = groupcounts(tbl,'Patient');
keepP3 = gc.Patient(gc.GroupCount >= 2);
tbl = tbl(ismember(tbl.Patient, keepP3), :);
if height(tbl) < 3 || numel(unique(tbl.Patient)) < 2
    warning('Not enough data to fit a mixed model.'); return
end

% (Centering HasSz is optional; we keep raw 0/1 so coefficient is interpretable as difference)
hasVar = 'HasSz';

%% =============== MIXED-EFFECTS MODELS: HasSz * EpilepsyType =================
% Random intercept per patient
lme_intercept = fitlme(tbl, sprintf('%s ~ %s * EpilepsyType + (1|Patient)', respVar, hasVar), ...
    'DummyVarCoding','effects');
disp(lme_intercept);

% Add random slope for HasSz (try; can be singular if few 0/1 switches)
lme_slope = fitlme(tbl, sprintf('%s ~ %s * EpilepsyType + (%s|Patient)', respVar, hasVar, hasVar), ...
    'DummyVarCoding','effects');
disp(lme_slope);
fprintf('\n--- Likelihood-Ratio Test (intercept vs slope) ---\n');
compare(lme_intercept, lme_slope);


% Report key fixed effects
b = lme_intercept.Coefficients;
% Effect of HasSz in reference level (Focal is reference due to order above)
ixHas = strcmp(b.Name, hasVar);
fprintf('\n[ALL] Effect of HasSz within Focal: Δ = %.3f  SE = %.3f  p = %.3g\n', ...
    b.Estimate(ixHas), b.SE(ixHas), b.pValue(ixHas));

% Interaction term(s): HasSz:EpilepsyType_General
ixInt = strcmp(b.Name, sprintf('%s:EpilepsyType_General', hasVar));
if any(ixInt)
    fprintf('[ALL] General vs Focal difference in HasSz effect: Δ = %.3f  SE = %.3f  p = %.3g\n', ...
        b.Estimate(ixInt), b.SE(ixInt), b.pValue(ixInt));
else
    % If naming differs, print available names for debugging
    fprintf('[ALL] Interaction name not found. Available terms:\n');
    disp(b.Name);
end

%% ===================== VISUALIZATION =====================
% 1) Box (log scale) by HasSz overall and split by EpilepsyType
figure('Color','w');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% Overall
nexttile; 
boxchart(categorical(tbl.HasSz), (useLogTransform)*tbl.SpikeRate_log + (~useLogTransform)*tbl.SpikeRate_perHour);
xticklabels({'HasSz=0','HasSz=1'});
ylabel(yLabel);
title('All types');
grid on; box off;

% Focal
subF = tbl(tbl.EpilepsyType=="Focal",:);
nexttile;
if ~isempty(subF)
    boxchart(categorical(subF.HasSz), (useLogTransform)*subF.SpikeRate_log + (~useLogTransform)*subF.SpikeRate_perHour);
    xticklabels({'0','1'}); title('Focal'); grid on; box off;
else
    text(0.5,0.5,'Focal: none','HorizontalAlignment','center'); axis off;
end

% General
subG = tbl(tbl.EpilepsyType=="General",:);
nexttile;
if ~isempty(subG)
    boxchart(categorical(subG.HasSz), (useLogTransform)*subG.SpikeRate_log + (~useLogTransform)*subG.SpikeRate_perHour);
    xticklabels({'0','1'}); title('General'); grid on; box off;
else
    text(0.5,0.5,'General: none','HorizontalAlignment','center'); axis off;
end

% 2) Optional per-patient time trends of spike rate with marker shape by HasSz
try
    figure('Color','w'); hold on; grid on; box off;
    [Gp, ~] = findgroups(pairs.Patient);
    splitapply(@(vd,sp,hs,pid) plot(vd, log1p(sp), '-o', 'MarkerFaceColor','auto', ...
        'DisplayName', sprintf('P%d', pid(1))), ...
        pairs.Visit_Date, pairs.SpikeRate_perHour, pairs.HasSz, pairs.Patient, Gp);
    datetick('x','yyyy'); xlabel('Visit date'); ylabel('log(1 + spike rate, spikes/hour)');
    title('Temporal spike rate (marker indicates paired visits)');
    legend('Location','bestoutside');
catch
    % skip quietly if plotting fails
end

%% =============== SAVE MATCHED PAIRS FOR AUDIT =================
writetable(pairs, pairsOut);
fprintf('Saved EEG–visit pairs to: %s\n', pairsOut);




%% ===== Option C: Patient-level paired analysis (HasSz=1 vs HasSz=0) =====
% Requirements in `pairs`:
%   pairs.Patient (double/int), pairs.HasSz (0/1), pairs.SpikeRate_perHour (double)
%   pairs.EpilepsyType (string/categorical with 'Focal'/'General')
% Optional: pairs.RescueApplied (logical) to exclude rescued pairs

% --- Config
useLogTransform   = true;   % log1p(spikes) for robustness vs skew
excludeRescued    = false;  % set true to drop pairs with RescueApplied==1 from this analysis
wantTypes         = ["Focal","General"];  % for stratified analysis

% --- Prepare working table
T = pairs;
if excludeRescued && ismember('RescueApplied', T.Properties.VariableNames)
    T = T(~T.RescueApplied, :);
end
T = T(ismember(string(T.EpilepsyType), wantTypes), :);
T = T(isfinite(T.SpikeRate_perHour) & isfinite(T.HasSz), :);

% Response for aggregation
if useLogTransform
    T.SpikeResp = log1p(T.SpikeRate_perHour);
    yLabel = 'log(1 + spikes/hour)';
else
    T.SpikeResp = T.SpikeRate_perHour;
    yLabel = 'spikes/hour';
end

% --- Summarize within patient & HasSz
[G, pid, has] = findgroups(T.Patient, T.HasSz);      % has is 0/1 per group
mu = splitapply(@(y) mean(y,'omitnan'), T.SpikeResp, G);
n  = splitapply(@numel,                    T.SpikeResp, G);
PerPH = table(pid, has, mu, n, 'VariableNames', {'Patient','HasSz','MeanResp','Npairs'});

% --- Manual pivot (robust; no unstack quirks)
t0 = PerPH(PerPH.HasSz==0, {'Patient','MeanResp','Npairs'});
t1 = PerPH(PerPH.HasSz==1, {'Patient','MeanResp','Npairs'});
t0.Properties.VariableNames(2:3) = {'MeanResp_0','Npairs_0'};
t1.Properties.VariableNames(2:3) = {'MeanResp_1','Npairs_1'};

PerP = outerjoin(t0, t1, 'Keys','Patient', 'MergeKeys', true);

% --- Attach one epilepsy type per patient (first non-missing)
[Gp, pid_pt] = findgroups(T.Patient);
etype = splitapply(@(x) string(x(find(~ismissing(x),1,'first'))), T.EpilepsyType, Gp);
Ttype = table(pid_pt, etype, 'VariableNames', {'Patient','EpilepsyType'});

PerP = innerjoin(PerP, Ttype, 'Keys','Patient');

% --- Coverage flags
has0 = ismember('MeanResp_0', PerP.Properties.VariableNames) & isfinite(PerP.MeanResp_0);
has1 = ismember('MeanResp_1', PerP.Properties.VariableNames) & isfinite(PerP.MeanResp_1);
PerP.HasBoth = has0 & has1;


fprintf('\n=== Patient-level paired analysis (HasSz 1 vs 0) ===\n');
fprintf('Patients with >=1 pair & any HasSz=0: %d\n', sum(isfinite(PerP.MeanResp_0)));
fprintf('Patients with >=1 pair & any HasSz=1: %d\n', sum(isfinite(PerP.MeanResp_1)));
fprintf('Patients with BOTH conditions: %d\n', sum(PerP.HasBoth));

% Keep only patients with both conditions for paired test
PP = PerP(PerP.HasBoth, :);
if isempty(PP)
    warning('No patients have both HasSz=0 and HasSz=1 pairs.'); 
    return
end

% Compute within-patient delta (1 - 0)
PP.Delta = PP.MeanResp_1 - PP.MeanResp_0;

% --- Overall paired signed-rank test
[p_overall,~,stats_overall] = signrank(PP.MeanResp_1, PP.MeanResp_0);
fprintf('Overall paired signed-rank: median(Δ) = %.3f  (n=%d),  p = %.3g,  z = %.3g\n', ...
    median(PP.Delta,'omitnan'), height(PP), p_overall, stats_overall.zval);

% --- Stratified by EpilepsyType
for t = 1:numel(wantTypes)
    et = wantTypes(t);
    sub = PP(strcmpi(string(PP.EpilepsyType), et), :);
    if height(sub) < 3
        fprintf('%s: insufficient patients with both conditions (n=%d)\n', et, height(sub));
        continue
    end
    [p_t,~,stats_t] = signrank(sub.MeanResp_1, sub.MeanResp_0);
    fprintf('%s paired: median(Δ) = %.3f  (n=%d),  p = %.3g,  z = %.3g\n', ...
        et, median(sub.Delta,'omitnan'), height(sub), p_t, stats_t.zval);
end

% --- Compact summary table (useful to save)
Summary = PP(:, {'Patient','EpilepsyType','MeanResp_0','MeanResp_1','Delta'});
if useLogTransform
    Summary.Properties.VariableNames(3:5) = {'Mean_logSpikes_HasSz0','Mean_logSpikes_HasSz1','Delta_log'};
else
    Summary.Properties.VariableNames(3:5) = {'Mean_spikesHr_HasSz0','Mean_spikesHr_HasSz1','Delta_raw'};
end
disp(head(Summary, 10));

% --- Visualization 1: paired lines (per patient)
figure('Color','w'); hold on; grid on; box off;
plot([zeros(height(PP),1) ones(height(PP),1)]', [PP.MeanResp_0 PP.MeanResp_1]', '-o');
xlim([-0.2 1.2]); xticks([0 1]); xticklabels({'HasSz=0','HasSz=1'});
ylabel(yLabel);
title('Within-patient mean spike rate: HasSz=1 vs HasSz=0');
% overlay medians
med0 = median(PP.MeanResp_0,'omitnan'); med1 = median(PP.MeanResp_1,'omitnan');
plot([0 0], [med0 med0], 'k-', 'LineWidth', 2);
plot([1 1], [med1 med1], 'k-', 'LineWidth', 2);
legend('Patients','Group medians','Location','bestoutside');

% --- Visualization 2: distribution of per-patient deltas
figure('Color','w'); hold on; grid on; box off;
boxchart(categorical(ones(height(PP),1)), PP.Delta);
yline(0,'k--');
ylabel(sprintf('Δ %s (HasSz=1 minus HasSz=0)', yLabel));
title('Within-patient change in mean spike rate');
set(gca,'XTickLabel',{'All patients'});

% --- Optional: stratified delta plots by type
figure('Color','w'); hold on; grid on; box off;
cats = categorical(string(PP.EpilepsyType), cellstr(wantTypes));
boxchart(cats, PP.Delta);
yline(0,'k--');
ylabel(sprintf('Δ %s (HasSz=1 - HasSz=0)', yLabel));
title('Within-patient change by epilepsy type');

end