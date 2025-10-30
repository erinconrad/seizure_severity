function within_sz_severity
%% ======================= CONFIG =======================
outCsv      = '../data/SN_counts/spike_counts_summary.csv';                     % spike counts with Duration_sec
reportFile  = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-20_1418.csv';
pairsOut    = '../output/eeg_visit_pairs_ALLTYPES_hasSz.csv';                   % audit output

maxGapDays        = 365;     % visit must be within this many days of EEG
useLogTransform   = true;    % log1p(spike rate) for modeling/plots
removeConstBinary = false;   % drop patients with no within-patient variation in HasSz across pairs

% EXCLUDE these epilepsy-type categories entirely (case-insensitive match)
badTypes = lower([
    "Non-Epileptic Seizure Disorder"
    "Uncertain if Epilepsy"
    "Unknown or MRN not found"
    "" % empty / missing
]);

% Three-group covariate order (reference = "General" when using DummyVarCoding='reference')
wantTypes4 = ["General","Temporal","Frontal"];

%% ================== LOAD & PREP SPIKE SUMMARY ==================
S = readtable(outCsv, 'TextType','string', 'VariableNamingRule','preserve');
reqS = {'EEG_Name','Patient','Session','Total_Spikes','Left_Spikes','Right_Spikes','Duration_sec'};
assert(all(ismember(reqS, S.Properties.VariableNames)), 'Spike summary missing required columns.');
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end
S.SpikeRate_perHour = S.Total_Spikes ./ (S.Duration_sec/3600);

%% ================== LOAD & PREP REPORT SHEET ===================
R = readtable(reportFile, 'TextType','string', 'VariableNamingRule','preserve');
reqR = {'patient_id','session_number','start_time_deid','visit_dates_deid','visit_hasSz','epilepsy_type'};
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
    cells = rr.HasSzVec;
    nonEmpty = ~cellfun(@isempty, cells);
    if any(nonEmpty)
        vals = vertcat(cells{nonEmpty});
        vals = vals(:);
        vals = vals(isfinite(vals));
    else
        vals = [];
    end
    hasAllZeroMap(pid) = (~isempty(vals) && all(vals==0));
end
nAllZero = sum(cell2mat(values(hasAllZeroMap)));
fprintf('Patients with ALL visit_hasSz == 0: %d\n', nAllZero);

%% ===== Build per-patient EpilepsyType and EpilepsySpecific (if present) =====
epi_str = string(R.epilepsy_type);
isEmpty = ismissing(epi_str) | strlength(strtrim(epi_str))==0;
Rt = sortrows(R(~isEmpty, {'patient_id','epilepsy_type'}), 'patient_id');
[uniq_pid, ia] = unique(Rt.patient_id, 'stable');
PerPatType = table(uniq_pid, Rt.epilepsy_type(ia), 'VariableNames', {'Patient','EpilepsyType'});

% Attach per-patient epilepsy_specific (first non-missing) if available
hasSpec = ismember('epilepsy_specific', R.Properties.VariableNames);
if hasSpec
    if ~isstring(R.epilepsy_specific), R.epilepsy_specific = string(R.epilepsy_specific); end
    spec_ok = ~ismissing(R.epilepsy_specific) & strlength(strtrim(R.epilepsy_specific))>0;
    RtSpec = sortrows(R(spec_ok, {'patient_id','epilepsy_specific'}), 'patient_id');
    [uniq_pid_s, ia_s] = unique(RtSpec.patient_id, 'stable');
    PerPatSpec = table(uniq_pid_s, RtSpec.epilepsy_specific(ia_s), ...
        'VariableNames', {'Patient','EpilepsySpecific'});
else
    PerPatSpec = table([],[], 'VariableNames', {'Patient','EpilepsySpecific'});
end

% Merge specific onto type
PerPat = outerjoin(PerPatType, PerPatSpec, 'Keys','Patient', 'MergeKeys', true);
if ~isstring(PerPat.EpilepsyType),    PerPat.EpilepsyType    = string(PerPat.EpilepsyType);    end
if ~isstring(PerPat.EpilepsySpecific),PerPat.EpilepsySpecific= string(PerPat.EpilepsySpecific); end
PerPat.EpilepsyType     = strtrim(PerPat.EpilepsyType);
PerPat.EpilepsySpecific = strtrim(PerPat.EpilepsySpecific);

% Exclude global "bad" types
%{
epi_norm_pat = lower(PerPat.EpilepsyType);
isBad = ismember(epi_norm_pat, badTypes);
PerPat = PerPat(~isBad, :);
%}

% Derive 3-level covariate: General / Temporal / Frontal
etLower   = lower(PerPat.EpilepsyType);
esLower   = lower(PerPat.EpilepsySpecific);
isGeneral = contains(etLower, "general");
isTemporal = false(height(PerPat),1);
isFrontal  = false(height(PerPat),1);
if hasSpec && ~isempty(esLower)
    isTemporal = contains(esLower, "temporal");
    isFrontal  = contains(esLower, "frontal");
end

EpiType4 = strings(height(PerPat),1);
EpiType4(isGeneral)  = "General";
EpiType4(~isGeneral & isTemporal) = "Temporal";
EpiType4(~isGeneral & ~isTemporal & isFrontal) = "Frontal";
keep3 = ismember(EpiType4, wantTypes4);
PerPat = PerPat(keep3, :);
PerPat.EpiType4 = EpiType4(keep3);

fprintf('Per-patient counts by EpiType4:\n');
disp(groupsummary(PerPat,'EpiType4'));

% Restrict R and S to valid patients
validPatients = PerPat.Patient;
R = R(ismember(R.patient_id, validPatients), :);
S = S(ismember(S.Patient,    validPatients), :);

%% ====== JOIN spike summary with report to get EEG date & visits & hasSz ======
JR = innerjoin(S, R, ...
    'LeftKeys',  {'Patient','Session'}, ...
    'RightKeys', {'patient_id','session_number'}, ...
    'RightVariables', {'EEG_Date','VisitDates','HasSzVec'});

% Keep rows with a valid EEG date and at least one visit date
JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) & any(~isnat(v)), JR.VisitDates), :);

%% ====== BUILD EEG–VISIT PAIRS (closest visit within maxGapDays) with ALL-ZERO rescue ======
pairs = table('Size',[0 11], ...
    'VariableTypes', {'double','double','string','datetime','double','datetime','double','logical','logical','string','string'}, ...
    'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour','Visit_Date','HasSz','RescueApplied','HasSzAllZero','EpilepsyType','EpiType4'});

% fast maps
typeMap   = containers.Map(double(PerPat.Patient), cellstr(PerPat.EpilepsyType));
type4Map  = containers.Map(double(PerPat.Patient), cellstr(PerPat.EpiType4));

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

    % ALL-zero rescue
    if ~(isfinite(thisHas) && (thisHas==0 || thisHas==1))
        if isKey(hasAllZeroMap, pid) && hasAllZeroMap(pid)
            thisHas = 0;
            rescueApplied = true;
        else
            continue
        end
    end

    etype  = ""; if isKey(typeMap,  pid), etype  = string(typeMap(pid));  end
    etype4 = ""; if isKey(type4Map, pid), etype4 = string(type4Map(pid)); end
    if ~ismember(etype4, wantTypes4), continue; end

    pairs = [pairs; {pid, JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perHour(i), ...
                     vdates(idx), thisHas, rescueApplied, ...
                     (isKey(hasAllZeroMap,pid) && hasAllZeroMap(pid)), ...
                     etype, etype4}]; %#ok<AGROW>
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

% Keep only the 3 canonical labels
pairs = pairs(ismember(pairs.EpiType4, wantTypes4), :);

%% ==================== MODELING PREP ====================
tbl = pairs;
tbl.Patient  = categorical(tbl.Patient);
tbl.HasSz    = double(tbl.HasSz);
tbl.EpiType4 = categorical(string(tbl.EpiType4), wantTypes4); % enforce order

if useLogTransform
    tbl.SpikeRate_log = real(double(log1p(tbl.SpikeRate_perHour)));
    respVar = 'SpikeRate_log';
    yLabel  = 'log(1 + spike rate, spikes/hour)';
else
    tbl.SpikeRate_perHour = double(tbl.SpikeRate_perHour);
    respVar = 'SpikeRate_perHour';
    yLabel  = 'spike rate (spikes/hour)';
end

ok = isfinite(tbl.(respVar)) & isfinite(tbl.HasSz) & ~isundefined(tbl.Patient) & ~isundefined(tbl.EpiType4);
tbl = tbl(ok,:);

% Need ≥2 pairs per patient (defensive after drops)
gc  = groupcounts(tbl,'Patient');
keepP3 = gc.Patient(gc.GroupCount >= 2);
tbl = tbl(ismember(tbl.Patient, keepP3), :);
if height(tbl) < 3 || numel(unique(tbl.Patient)) < 2
    warning('Not enough data to fit a mixed model.'); return
end

% Optional: drop patients with constant HasSz across pairs
if removeConstBinary
    [Gc2, pids2] = findgroups(tbl.Patient);
    varHas = splitapply(@(v) nanvar(v), tbl.HasSz, Gc2);
    keepP2 = pids2(varHas > 0);
    tbl = tbl(ismember(tbl.Patient, keepP2), :);
    if isempty(tbl), warning('Nothing left after constant-HasSz filter.'); return; end
end

%% =============== MAIN MIXED-EFFECTS MODEL: random intercept only ===============
% Effects coding (grand-mean reference)
lme_intercept = fitlme(tbl, sprintf('%s ~ HasSz * EpiType4 + (1|Patient)', respVar), ...
    'DummyVarCoding','effects');
disp(lme_intercept);

% Reference coding with "General" as reference
lme_ref = fitlme(tbl, sprintf('%s ~ HasSz * EpiType4 + (1|Patient)', respVar), ...
    'DummyVarCoding','reference'); % reference = first level in EpiType4 (General)
disp(lme_ref);

% Quick console report
br = lme_ref.Coefficients;
ixHas_General = strcmp(br.Name, 'HasSz');
fprintf('\n[General] Effect of HasSz: Δ = %.3f  SE = %.3f  p = %.3g\n', ...
    br.Estimate(ixHas_General), br.SE(ixHas_General), br.pValue(ixHas_General));
for lab = ["Temporal","Frontal"]
    nm = sprintf('HasSz:EpiType4_%s', lab);
    ix = strcmp(br.Name, nm);
    if any(ix)
        fprintf('[%s vs General] Difference in HasSz effect: Δ = %.3f  SE = %.3f  p = %.3g\n', ...
            lab, br.Estimate(ix), br.SE(ix), br.pValue(ix));
    end
end

%% ===================== FIGURES =====================

% ---------- Panel A: Population-level marginal means (FIXED) ----------
types = categories(tbl.EpiType4);

% Build population-level prediction grid with required vars
newT = table();
for k = 1:numel(types)
    for hs = [0 1]
        newT = [newT; table(categorical(types(k), types), double(hs), ...
            'VariableNames', {'EpiType4','HasSz'})]; %#ok<AGROW>
    end
end

% Ensure categories/types match training data exactly
newT.EpiType4 = categorical(newT.EpiType4, categories(tbl.EpiType4));  % same order
newT.HasSz    = double(newT.HasSz);                                     % same type

% IMPORTANT: add grouping variable used in the model formula
% Use any existing patient ID (since Conditional=false, random effects are 0)
refPatient = tbl.Patient(1);
newT.Patient = repmat(refPatient, height(newT), 1);  % categorical with same level

% Predict using FIXED effects only
[yhatA, yCIA] = predict(lme_ref, newT, 'Conditional', false);

% Plot
figure('Color','w'); hold on; box off; grid on; set(gca,'FontSize',18)
xA = repelem(1:numel(types), 2);
xA = xA + repmat([-0.12 0.12], 1, numel(types));
for i = 1:height(newT)
    line([xA(i) xA(i)], yCIA(i,:), 'LineWidth',2);
end
for k = 1:numel(types)
    idx = (newT.EpiType4==types(k));
    plot(xA(idx), yhatA(idx), 'o-','LineWidth',2,'MarkerFaceColor',[.8 .8 .8])
end
set(gca,'XTick',1:numel(types),'XTickLabel',cellstr(types))
xlabel('Epilepsy type'); ylabel(yLabel);
title('Panel A: Population-level marginal means (HasSz=0 vs 1)')
legend({'95% CI','Prediction'},'Location','northwest'); legend boxoff

% ---------- Panel B: Within-patient paired “spaghetti” ----------
if useLogTransform
    ycol = 'SpikeRate_log';
else
    ycol = 'SpikeRate_perHour';
end
Tg = groupsummary(tbl, {'Patient','EpiType4','HasSz'}, 'mean', ycol);
Tg.Properties.VariableNames(end) = {'MeanY'};

figure('Color','w'); tiledlayout(1, numel(types), 'TileSpacing','compact','Padding','compact');
for k = 1:numel(types)
    nexttile; hold on; grid on; box off; set(gca,'FontSize',16)
    Tk = Tg(Tg.EpiType4==types(k),:);
    [Gp, pid] = findgroups(Tk.Patient);
    haveBoth = splitapply(@(v) numel(unique(v))==2, Tk.HasSz, Gp);
    keepP = pid(haveBoth);
    Tk = Tk(ismember(Tk.Patient, keepP),:);

    T0 = Tk(Tk.HasSz==0, {'Patient','MeanY'});
    T1 = Tk(Tk.HasSz==1, {'Patient','MeanY'});
    W  = innerjoin(T0, T1, 'Keys','Patient', ...
        'LeftVariables',{'Patient','MeanY'}, 'RightVariables','MeanY');
    W.Properties.VariableNames = {'Patient','Y0','Y1'};

    for i = 1:height(W)
        plot([0 1], [W.Y0(i) W.Y1(i)], '-', 'Color',[0.6 0.6 0.6])
    end
    scatter(zeros(height(W),1), W.Y0, 20, 'filled')
    scatter(ones(height(W),1),  W.Y1, 20, 'filled')

    xlim([-0.25 1.25]); set(gca,'XTick',[0 1],'XTickLabel',{'HasSz=0','HasSz=1'})
    ylabel(yLabel); title(sprintf('Panel B: %s', string(types(k))))

    medDelta = median(W.Y1 - W.Y0, 'omitnan');
    text(0.02, max([W.Y0;W.Y1]) - 0.05*range([W.Y0;W.Y1]), sprintf('Median Δ = %.2f', medDelta), 'FontSize',14)
end

% ---------- Panel C: Fixed-effects forest (reference coding) [FIXED] ----------
FE  = lme_ref.Coefficients;

% Collect the terms you care about; keep only those that exist
wantNames = {'HasSz','HasSz:EpiType4_Temporal','HasSz:EpiType4_Frontal'};
mask = ismember(FE.Name, wantNames);
FEs = FE(mask, :);

% Build readable labels in same order as FEs
labels = strings(height(FEs),1);
for i = 1:height(FEs)
    switch FEs.Name{i}
        case 'HasSz'
            labels(i) = "HasSz (General)";
        case 'HasSz:EpiType4_Temporal'
            labels(i) = "Temporal vs General (Δ HasSz)";
        case 'HasSz:EpiType4_Frontal'
            labels(i) = "Frontal  vs General (Δ HasSz)";
        otherwise
            labels(i) = FEs.Name{i};
    end
end

% y must be increasing
y = 1:height(FEs);

% 95% CI (normal approx)
ciL = FEs.Estimate - 1.96*FEs.SE;
ciU = FEs.Estimate + 1.96*FEs.SE;

figure('Color','w'); hold on; grid on; box off; set(gca,'FontSize',16)

% draw CI bars
for i = 1:numel(y)
    line([ciL(i) ciU(i)], [y(i) y(i)], 'LineWidth', 2);
end
% point estimates
plot(FEs.Estimate, y, 'o', 'LineWidth', 2, 'MarkerFaceColor', [.85 .85 .85])

% cosmetics
set(gca, 'YTick', y, 'YTickLabel', labels);
ylim([0.5, max(y)+0.5]);   % give a little padding
xlabel(sprintf('Effect on %s', yLabel));
title('Panel C: Fixed-effects (reference coding)');
xline(0, 'k:');

% (optional) show p-values at right
for i = 1:numel(y)
    text(ciU(i) + 0.02*range([ciL;ciU]), y(i), ...
        sprintf('p = %.3g', FEs.pValue(i)), 'FontSize', 12, 'VerticalAlignment','middle');
end

% ---------- Panel D: Sensitivity of HasSz effect by EEG–visit gap ----------
pairs.GapDays = days(pairs.Visit_Date - pairs.EEG_Date);
edges  = [0 30 90 180 365];
labels = ["0–30","31–90","91–180","181–365"];
binIdx = discretize(pairs.GapDays, edges, 'IncludedEdge','right');
pairs.GapBin = categorical(binIdx, 1:numel(labels), labels);

tblD = pairs(ismember(pairs.EpiType4, wantTypes4), :);
tblD.Patient  = categorical(tblD.Patient);
tblD.HasSz    = double(tblD.HasSz);
tblD.EpiType4 = categorical(string(tblD.EpiType4), wantTypes4);

if useLogTransform
    tblD.SpikeRate_log = real(double(log1p(tblD.SpikeRate_perHour)));
    respVarD = 'SpikeRate_log';
else
    tblD.SpikeRate_perHour = double(tblD.SpikeRate_perHour);
    respVarD = 'SpikeRate_perHour';
end

okD = isfinite(tblD.(respVarD)) & isfinite(tblD.HasSz) & ...
      ~isundefined(tblD.Patient) & ~isundefined(tblD.EpiType4) & ...
      ~isnan(tblD.GapDays) & ~isundefined(tblD.GapBin);
tblD = tblD(okD,:);

bins = categories(tblD.GapBin);
nb   = numel(bins);
est  = nan(nb,1); se = nan(nb,1); nobs = nan(nb,1); npid = nan(nb,1);

for i = 1:nb
    Ti = tblD(tblD.GapBin==bins{i}, :);

    % remove unused categories (important for rank)
    Ti.EpiType4 = removecats(Ti.EpiType4);
    Ti.Patient  = removecats(Ti.Patient);

    % Require ≥2 pairs per patient *within bin* (keeps within-patient contrast)
    gc_i   = groupcounts(Ti,'Patient');
    keepPi = gc_i.Patient(gc_i.GroupCount >= 2);
    Ti     = Ti(ismember(Ti.Patient, keepPi), :);

    % Basic sufficiency checks
    if height(Ti) < 5 || numel(categories(Ti.Patient)) < 2
        continue
    end
    if numel(unique(Ti.HasSz)) < 2
        % HasSz doesn’t vary in this bin → can’t estimate effect
        continue
    end

    % Choose formula depending on how many EpiType4 levels are present
    if numel(categories(Ti.EpiType4)) >= 2
        fml = sprintf('%s ~ HasSz + EpiType4 + (1|Patient)', respVarD);
    else
        % Only one type present → drop EpiType4 to avoid collinearity
        fml = sprintf('%s ~ HasSz + (1|Patient)', respVarD);
    end

    try
        mdl = fitlme(Ti, fml, 'DummyVarCoding','reference');
    catch
        % Final safeguard: if still rank-deficient, skip bin
        continue
    end

    b = mdl.Coefficients;
    ixHas = strcmp(b.Name, 'HasSz');
    if any(ixHas)
        est(i)  = b.Estimate(ixHas);
        se(i)   = b.SE(ixHas);
        nobs(i) = height(Ti);
        npid(i) = numel(categories(Ti.Patient));
    end
end


ciL = est - 1.96*se;
ciU = est + 1.96*se;

figure('Color','w'); hold on; grid on; box off; set(gca,'FontSize',18)
x = 1:nb;
for i = 1:nb
    if ~isnan(est(i))
        line([x(i) x(i)], [ciL(i) ciU(i)], 'LineWidth',2);
    end
end
plot(x, est, 'o-', 'LineWidth',2, 'MarkerFaceColor',[.85 .85 .85])

xticks(x); xticklabels(bins); xtickangle(0)
xlabel('EEG–visit gap (days)')
ylabel(sprintf('HasSz effect on %s', yLabel))
title('Panel D: Sensitivity of HasSz effect by EEG–visit gap')
yline(0,'k:');

for i = 1:nb
    if ~isnan(est(i))
        if useLogTransform
            pct = 100*(exp(est(i)) - 1);
            txt = sprintf('  N=%d (P=%d)  %+0.0f%%', nobs(i), npid(i), pct);
        else
            txt = sprintf('  N=%d (P=%d)', nobs(i), npid(i));
        end
        text(x(i), est(i), txt, 'VerticalAlignment','bottom','FontSize',12)
    else
        text(x(i), 0, 'Insufficient data', 'HorizontalAlignment','center','FontSize',12)
    end
end

%% =============== SAVE MATCHED PAIRS FOR AUDIT =================
writetable(pairs, pairsOut);
fprintf('Saved EEG–visit pairs to: %s\n', pairsOut);

end
