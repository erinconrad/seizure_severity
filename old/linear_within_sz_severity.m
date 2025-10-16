function within_sz_severity
%% ======================= CONFIG =======================
outCsv      = '../data/SN_counts/spike_counts_summary.csv';                     % spike counts with Duration_sec
reportFile  = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-14_1505.csv';
pairsOut    = '../output/eeg_visit_pairs_GENERAL.csv';                          % audit output

maxGapDays        = 365;     % visit must be within this many days of EEG
NEG_POLICY        = 'nan';   % 'nan' (drop negatives) or 'zero' (set to 0)
useLogTransform   = true;    % log1p transform for modeling/plots (recommended)
removeConstSzFreq = false;   % drop patients with no Sz_Freq variation across pairs

% Keep ONLY this epilepsy type
keepType = "General";

% EXCLUDE these epilepsy-type categories entirely
badTypes = lower([
    "Non-Epileptic Seizure Disorder"
    "Uncertain if Epilepsy"
    "Unknown or MRN not found"
    "" % empty / missing
]);

%% ================== LOAD & PREP SPIKE SUMMARY ==================
S = readtable(outCsv, 'TextType','string', 'VariableNamingRule','preserve');  
reqS = {'EEG_Name','Patient','Session','Total_Spikes','Left_Spikes','Right_Spikes','Duration_sec'};
assert(all(ismember(reqS, S.Properties.VariableNames)), 'Spike summary missing required columns.');
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end
S.SpikeRate_perHour = S.Total_Spikes ./ (S.Duration_sec/3600);

%% ================== LOAD & PREP REPORT SHEET ===================
R = readtable(reportFile, 'TextType','string', 'VariableNamingRule','preserve');
reqR = {'patient_id','session_number','start_time_deid','visit_dates_deid','sz_freqs','epilepsy_type'};
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

% Parse visit dates & sz freqs into cell arrays
rawVisits = string(R.('visit_dates_deid'));
rawSz     = string(R.('sz_freqs'));
R.VisitDates = cell(height(R),1);
R.SzFreqs    = cell(height(R),1);
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
    % Seizure freqs
    sstr = strtrim(rawSz(i));
    if strlength(sstr) > 0 && sstr ~= "[]"
        sfix = regexprep(sstr, 'null', 'NaN', 'ignorecase');
        try
            vv = jsondecode(char(sfix)); 
            R.SzFreqs{i} = vv(:);
        catch
            nums = regexp(sfix, '[-+]?\d+(\.\d+)?([eE][-+]?\d+)?', 'match');
            R.SzFreqs{i} = str2double(string(nums));
        end
    else
        R.SzFreqs{i} = nan(0,1);
    end
end

%% ===== Per-patient rescue flag: HasSzAllZero (ALL visit_hasSz values are 0) =====
upids = unique(R.patient_id(~isnan(R.patient_id)));
hasAllZeroMap = containers.Map('KeyType','double','ValueType','logical');
for k = 1:numel(upids)
    pid = upids(k);
    rr = R(R.patient_id==pid, :);
    vals = [];
    if ismember('visit_hasSz', R.Properties.VariableNames)
        for j = 1:height(rr)
            rawHas = strtrim(string(rr.visit_hasSz(j)));
            if strlength(rawHas)==0 || rawHas=="[]" || rawHas==""
                continue
            end
            try
                vHas = jsondecode(char(rawHas));   % numeric vector (0/1)
                vals = [vals; vHas(:)]; %#ok<AGROW>
            catch
                nums = regexp(rawHas, '\d+', 'match');
                if ~isempty(nums)
                    vals = [vals; str2double(string(nums))]; %#ok<AGROW>
                end
            end
        end
    end
    hasAllZeroMap(pid) = (~isempty(vals) && all(vals==0));
end
nAllZero = sum(cell2mat(values(hasAllZeroMap)));
fprintf('Patients with ALL visit_hasSz == 0: %d\n', nAllZero);

%% ===== Build per-patient epilepsy_type and FILTER: keep ONLY "General" =====
epi_str = string(R.epilepsy_type);
isEmpty = ismissing(epi_str) | strlength(strtrim(epi_str))==0;
Rt = sortrows(R(~isEmpty, {'patient_id','epilepsy_type'}), 'patient_id');
[uniq_pid, ia] = unique(Rt.patient_id, 'stable');
PerPatType = table(uniq_pid, Rt.epilepsy_type(ia), 'VariableNames', {'Patient','EpilepsyType'});

epi_norm_pat = lower(strtrim(PerPatType.EpilepsyType));
isBad  = ismember(epi_norm_pat, badTypes);
isKeep = strcmpi(PerPatType.EpilepsyType, keepType);
validPatients = PerPatType.Patient(~isBad & isKeep);

fprintf('Total patients in report: %d\n', numel(unique(R.patient_id(~isnan(R.patient_id)))));
fprintf('Patients kept (GENERAL only): %d\n', numel(validPatients));
fprintf('Patients excluded by type filter: %d\n', sum(~ismember(PerPatType.Patient, validPatients)));

% Restrict R and S to valid patients
R = R(ismember(R.patient_id, validPatients), :);
S = S(ismember(S.Patient,    validPatients), :);

%% ====== JOIN spike summary with report to get EEG date & visits ======
JR = innerjoin(S, R, ...
    'LeftKeys',  {'Patient','Session'}, ...
    'RightKeys', {'patient_id','session_number'}, ...
    'RightVariables', {'EEG_Date','VisitDates','SzFreqs'});

JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) && all(~isnat(v)), JR.VisitDates), :);

%% ====== BUILD EEG–VISIT PAIRS (closest visit within maxGapDays) with ALL-ZERO rescue ======
pairs = table('Size',[0 9], ...
    'VariableTypes', {'double','double','string','datetime','double','datetime','double','logical','logical'}, ...
    'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour','Visit_Date','Sz_Freq','RescueApplied','HasSzAllZero'});

for i = 1:height(JR)
    eegDate = JR.EEG_Date(i);
    vdates  = JR.VisitDates{i};
    szf     = JR.SzFreqs{i};
    pid     = JR.Patient(i);

    if isempty(vdates) || isempty(szf), continue; end
    nAlign = min(numel(vdates), numel(szf));
    vdates = vdates(1:nAlign);
    szf    = szf(1:nAlign);

    [dmin, idx] = min(abs(days(vdates - eegDate)));
    if isempty(dmin) || dmin > maxGapDays, continue; end

    thisSz = szf(idx);
    rescueApplied = false;

    % ALL-zero rescue: if the chosen visit's Sz_Freq is NaN and patient is all-zero, impute 0
    if (isnan(thisSz) || ~isfinite(thisSz))
        if isKey(hasAllZeroMap, pid) && hasAllZeroMap(pid)
            thisSz = 0;
            rescueApplied = true;
        else
            % still missing — can't use this pair
            continue
        end
    end

    pairs = [pairs; {pid, JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perHour(i), vdates(idx), thisSz, rescueApplied, ...
                     (isKey(hasAllZeroMap,pid) && hasAllZeroMap(pid))}]; %#ok<AGROW>
end

fprintf('Built %d EEG–visit pairs (within %d days). Rescued pairs: %d\n', ...
    height(pairs), maxGapDays, sum(pairs.RescueApplied));
if isempty(pairs)
    warning('No EEG–visit pairs were constructed. Check date formats and maxGapDays.');
    return
end

% Require ≥2 pairs per patient
counts   = groupcounts(pairs, 'Patient');
eligible = counts.Patient(counts.GroupCount >= 2);
pairs    = pairs(ismember(pairs.Patient, eligible), :);
fprintf('Eligible GENERAL patients with ≥2 pairs: %d\n', numel(eligible));
if isempty(pairs)
    warning('After requiring ≥2 pairs per patient, nothing remains.');
    return
end

% Remove patients whose EEGs all matched the same visit date
[G, pid] = findgroups(pairs.Patient);
nUniqueVisits = splitapply(@(d) numel(unique(d)), pairs.Visit_Date, G);
keepPatients = pid(nUniqueVisits >= 2);
pairs = pairs(ismember(pairs.Patient, keepPatients), :);
if isempty(pairs), warning('Nothing left after visit-date filter.'); return; end

%% ===== Clean Sz_Freq (negatives) & drop non-finite =====
pairs.SpikeRate_perHour = double(pairs.SpikeRate_perHour);
pairs.Sz_Freq           = double(pairs.Sz_Freq);
negMask = isfinite(pairs.Sz_Freq) & pairs.Sz_Freq < 0;
switch lower(NEG_POLICY)
    case 'nan',  pairs.Sz_Freq(negMask) = NaN;
    case 'zero', pairs.Sz_Freq(negMask) = 0;
    otherwise,   error('Unknown NEG_POLICY: %s', NEG_POLICY);
end
pairs = pairs(isfinite(pairs.SpikeRate_perHour) & isfinite(pairs.Sz_Freq), :);
if isempty(pairs)
    warning('All pairs dropped after Sz_Freq cleaning. Consider NEG_POLICY="zero" or inspect data.');
    return
end

%% ===== Optional: remove patients with constant Sz_Freq across their pairs =====
if removeConstSzFreq
    [Gc, pids2] = findgroups(pairs.Patient);
    varSz = splitapply(@(v) nanvar(v), pairs.Sz_Freq, Gc);
    keepP2 = pids2(varSz > 0);
    pairs = pairs(ismember(pairs.Patient, keepP2), :);
    if isempty(pairs), warning('Nothing left after constant-frequency filter.'); return; end
end

%% ==================== MODELING PREP ====================
tbl = pairs;
tbl.Patient = categorical(tbl.Patient);

if useLogTransform
    tbl.SpikeRate_log = real(double(log1p(tbl.SpikeRate_perHour)));
    tbl.SzFreq_log    = real(double(log1p(tbl.Sz_Freq)));
    respVar = 'SpikeRate_log';
    predVar = 'SzFreq_log';
    xLabel  = 'log(1 + seizure frequency)';
    yLabel  = 'log(1 + spike rate, spikes/hour)';
else
    tbl.SpikeRate_perHour = double(tbl.SpikeRate_perHour);
    tbl.Sz_Freq           = double(tbl.Sz_Freq);
    respVar = 'SpikeRate_perHour';
    predVar = 'Sz_Freq';
    xLabel  = 'Seizure frequency';
    yLabel  = 'Spike rate (spikes/hour)';
end

ok  = isfinite(tbl.(respVar)) & isfinite(tbl.(predVar)) & ~isundefined(tbl.Patient);
tbl = tbl(ok,:);

% Require ≥2 pairs per patient (post transforms)
gc  = groupcounts(tbl,'Patient');
keepP3 = gc.Patient(gc.GroupCount >= 2);
tbl = tbl(ismember(tbl.Patient, keepP3), :);
if height(tbl) < 3 || numel(unique(tbl.Patient)) < 2
    warning('Not enough GENERAL data to fit a mixed model after filtering.');
    return
end

% z-score predictor (across GENERAL rows) for stability
z = (tbl.(predVar) - mean(tbl.(predVar),'omitnan')) ./ std(tbl.(predVar),'omitnan');
z(~isfinite(z)) = 0;
predVarZ = [predVar '_z'];
tbl.(predVarZ) = z;

%% =============== MIXED-EFFECTS MODELS (GENERAL ONLY) =================
lme_intercept = fitlme(tbl, sprintf('%s ~ %s + (1|Patient)', respVar, predVarZ), ...
    'DummyVarCoding','effects');
disp(lme_intercept);

% Random slope (optional; wrapped in try in case of singular fit)

lme_slope = fitlme(tbl, sprintf('%s ~ %s + (%s|Patient)', respVar, predVarZ, predVarZ), ...
    'DummyVarCoding','effects');
disp(lme_slope);
fprintf('\n--- Likelihood-Ratio Test (intercept vs slope) ---\n');
compare(lme_intercept, lme_slope);


% Report fixed effect from intercept model
b = lme_intercept.Coefficients;
ix = strcmp(b.Name, predVarZ);
fprintf('\n[GENERAL ONLY] Fixed effect of seizure frequency: beta = %.3f  SE = %.3f  p = %.3g\n', ...
    b.Estimate(ix), b.SE(ix), b.pValue(ix));

%% ===================== VISUALIZATION =====================
% 1) Scatter by patient (z-scored X) + population fit
figure('Color','w'); hold on;
x = tbl.(predVarZ); 
y = tbl.(respVar);
gscatter(x, y, tbl.Patient, [], [], 15, 'off');
xlabel(sprintf('%s (z-scored)', xLabel));
ylabel(yLabel);
title('Within-patient EEG–visit pairs (GENERAL only)');
grid on; box off;

xx = linspace(min(x), max(x), 200)';
X = table(xx, 'VariableNames', {predVarZ});
X.Patient = repmat(tbl.Patient(1), height(X), 1);
X.Patient = categorical(X.Patient, categories(tbl.Patient));
yy = predict(lme_intercept, X, 'Conditional', false);
plot(xx, yy, 'k-', 'LineWidth', 2);
legend('Patients','Population fit','Location','bestoutside');

% 2) Time-aligned view (optional): spike rate vs visit date per patient
try
    figure('Color','w'); hold on; grid on; box off;
    [Gp, ~] = findgroups(pairs.Patient);
    splitapply(@(vd,sp,pid) plot(vd, log1p(sp), '-o', 'DisplayName', sprintf('P%d', pid(1))), ...
        pairs.Visit_Date, pairs.SpikeRate_perHour, pairs.Patient, Gp);
    datetick('x','yyyy'); xlabel('Visit date'); ylabel('log(1 + spike rate, spikes/hour)');
    title('Temporal trend of spike rate (GENERAL only)');
    legend('Location','bestoutside');
catch
    % silently skip if plotting fails
end

%% =============== SAVE MATCHED PAIRS FOR AUDIT =================
writetable(pairs, pairsOut);
fprintf('Saved GENERAL EEG–visit pairs to: %s\n', pairsOut);

end
