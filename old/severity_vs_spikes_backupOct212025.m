%% ===== Paths =====
outCsv     = '../data/SN_counts/spike_counts_summary.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-20_1418.csv';
% reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-14_1505.csv';
outPatient = '../output/patient_szfreq_vs_spikes.csv';

%% ===== Load spike counts with durations =====
S = readtable(outCsv, 'TextType','string');

% Add per-EEG spike rates
S.SpikeRate_Hz      = S.Total_Spikes ./ S.Duration_sec;
S.SpikeRate_perHour = S.SpikeRate_Hz * 3600;

% Per-patient average spike rate across EEGs
Sg = groupsummary(S, 'Patient', 'mean', {'Total_Spikes','SpikeRate_perHour'});
Sg.Properties.VariableNames = {'Patient','GroupCount','MeanTotalSpikes','MeanSpikeRate_perHour'};
Sg.GroupCount = [];

%% ===== Load report =====
R = readtable(reportFile, 'TextType','string');
if ~isnumeric(R.patient_id), R.patient_id = double(str2double(R.patient_id)); end
if ~isstring(R.epilepsy_type), R.epilepsy_type = string(R.epilepsy_type); end
if ~ismember('sz_freqs', R.Properties.VariableNames)
    error('Column "sz_freqs" not found in report file: %s', reportFile);
end

%% ===== Filter by epilepsy_type BEFORE computing sz_freqs =====
% Define disallowed categories (case-insensitive)
badTypes = lower([
    "Non-Epileptic Seizure Disorder"
    "Uncertain if Epilepsy"
    "Unknown or MRN not found"
    "" % empty
]);

% Normalize strings for comparison
epi_norm = lower(strtrim(R.epilepsy_type));
isEmpty  = ismissing(R.epilepsy_type) | strlength(strtrim(R.epilepsy_type))==0;

% Build a per-patient epilepsy_type table (take the first non-missing label if multiple)
Rt = sortrows(R(~isEmpty, {'patient_id','epilepsy_type'}), 'patient_id');  % deterministic
[uniq_pid, ia] = unique(Rt.patient_id, 'stable');
PerPatType = table(uniq_pid, Rt.epilepsy_type(ia), 'VariableNames', {'Patient','EpilepsyType'});

% Mark valid patients (not in badTypes and not empty)
epi_norm_pat = lower(strtrim(PerPatType.EpilepsyType));
isBad = ismember(epi_norm_pat, badTypes);
validPatients = PerPatType.Patient(~isBad);

% Optionally show counts
fprintf('Total patients in report: %d\n', numel(unique(R.patient_id(~isnan(R.patient_id)))));
fprintf('Patients with valid epilepsy_type kept: %d\n', numel(validPatients));
fprintf('Patients excluded by epilepsy_type: %d\n', sum(isBad));

% Keep only rows in R belonging to valid patients
R = R(ismember(R.patient_id, validPatients), :);

%% ===== Compute per-patient mean sz_freq with rescue from visit_hasSz =====
pids = unique(R.patient_id(~isnan(R.patient_id)));
MeanSzFreq = nan(numel(pids),1);

for k = 1:numel(pids)
    pid = pids(k);
    rows = R.patient_id == pid;
    vals = [];
    rr = R(rows,:);

    for j = 1:height(rr)
        raw = strtrim(rr.sz_freqs(j));
        if strlength(raw)==0 || raw=="[]" || raw==""
            continue
        end
        s = regexprep(raw, 'null', 'NaN', 'ignorecase');
        try
            v = jsondecode(char(s));
            vals = [vals; v(:)]; %#ok<AGROW>
        catch
            nums = regexp(s, '[-+]?\d+(\.\d+)?([eE][-+]?\d+)?', 'match');
            if ~isempty(nums)
                vals = [vals; str2double(string(nums))]; %#ok<AGROW>
            end
        end
    end

    % ---- Rescue rule: use visit_hasSz if sz_freq missing ----
    if (isempty(vals) || all(isnan(vals))) && ismember('visit_hasSz', R.Properties.VariableNames)
        rawHas = string(rr.visit_hasSz);
        allVals = [];
        for jj = 1:numel(rawHas)
            sHas = strtrim(rawHas(jj));
            if strlength(sHas)==0 || sHas=="[]" || sHas==""
                continue
            end
            try
                vHas = jsondecode(char(sHas));
                allVals = [allVals; vHas(:)]; %#ok<AGROW>
            catch
                nums = regexp(sHas, '\d+', 'match');
                if ~isempty(nums)
                    allVals = [allVals; str2double(string(nums))]; %#ok<AGROW>
                end
            end
        end
        % if we found any hasSz entries and *all* are 0, assume seizure frequency = 0
        if ~isempty(allVals) && all(allVals==0)
            vals = 0;  % rescue as 0 seizures
        end
    end

    MeanSzFreq(k) = mean(vals,'omitnan');
end

Rg = table(pids, MeanSzFreq, 'VariableNames', {'Patient','MeanSzFreq'});

%% ===== Attach EpilepsyType and filter =====
Rg = innerjoin(Rg, PerPatType(~isBad, :), 'Keys','Patient');

%% ===== Per-patient longest spike burst (max across EEGs) =====
if ~ismember('Longest_Spike_Duration_sec', S.Properties.VariableNames)
    error(['Column "Longest_Spike_Duration_sec" not found in %s. ', ...
           'Re-run the counting script that writes this column.'], outCsv);
end
Sb = groupsummary(S, 'Patient', 'max', 'Longest_Spike_Duration_sec');
Sb.Properties.VariableNames = {'Patient','GroupCount','MaxLongestBurst_sec'};
Sb.GroupCount = [];

%% ===== Merge per-patient spike metrics (rate + burst) =====
Spp = outerjoin(Sg, Sb, 'Keys','Patient', 'MergeKeys', true);
Spp.MaxLongestBurst_sec(~isfinite(Spp.MaxLongestBurst_sec)) = 0;

%% ===== Single clean join: Spike metrics + SzFreq + EpilepsyType =====
P = innerjoin(Spp, Rg, 'Keys','Patient');

% ---- Coalesce any EpilepsyType* columns into a single 'EpilepsyType' ----
etMask = startsWith(P.Properties.VariableNames, 'EpilepsyType', 'IgnoreCase', true);
etCols = P.Properties.VariableNames(etMask);
if isempty(etCols)
    error('EpilepsyType column not found after joins. Available columns: %s', ...
          strjoin(P.Properties.VariableNames, ', '));
end
if ~any(strcmp(etCols, 'EpilepsyType'))
    P = renamevars(P, etCols{1}, 'EpilepsyType');
end
% Drop any extra EpilepsyType-like columns
etCols = P.Properties.VariableNames(startsWith(P.Properties.VariableNames, 'EpilepsyType', 'IgnoreCase', true));
etCols = setdiff(etCols, 'EpilepsyType', 'stable');
if ~isempty(etCols), P = removevars(P, etCols); end

% Ensure types
if ~isstring(P.EpilepsyType), P.EpilepsyType = string(P.EpilepsyType); end
P.EpilepsyType = strtrim(P.EpilepsyType);

% ---- Normalize labels exactly to {'Focal','General'} and drop the rest ----
et = lower(strtrim(P.EpilepsyType));
isGeneral = contains(et, "general");
isFocal   = contains(et, "focal");
normType = strings(size(et));
normType(isGeneral) = "General";
normType(isFocal)   = "Focal";
P.EpilepsyType = normType;

% Keep rows with valid labels + predictors/outcome
P = P(~isnan(P.MeanSzFreq) & isfinite(P.MeanSpikeRate_perHour) & ...
      isfinite(P.MaxLongestBurst_sec) & strlength(P.EpilepsyType)>0, :);

fprintf('Patients with both measures (post-filter): %d\n', height(P));

%% ===== Save patient-level CSV (optional) =====
% writetable(P, outPatient);
% fprintf('Saved patient-level table to: %s\n', outPatient);

%% ===== Overall: Spike rate vs SzFreq =====
x = P.MeanSpikeRate_perHour;   % spikes/hour
y = P.MeanSzFreq;              % seizures/month (units from report)
mask = isfinite(x) & isfinite(y); x = x(mask); y = y(mask);

% Spearman
[Rho, Pval] = corr(x, y, 'Type','Spearman', 'Rows','complete');
fprintf('Spearman (raw scale): r = %.3f, p = %.3g (n=%d)\n', Rho, Pval, numel(x));

% Log–log scatter
minpos_x = min(x(x>0)); if isempty(minpos_x), minpos_x = 1e-6; end
minpos_y = min(y(y>0)); if isempty(minpos_y), minpos_y = 1e-6; end
x_log = log(x + (x<=0).*0.5*minpos_x);
y_log = log(y + (y<=0).*0.5*minpos_y);

figure('Color','w');
scatter(x_log, y_log, 36, 'filled'); grid on; box off;
xlabel('Mean Spike Rate (log spikes/hour)');
ylabel('Mean Seizure Frequency (log units)');
title(sprintf('Sz freq vs spike rate (Spearman r=%.3f, p=%.3g; n=%d)', Rho, Pval, numel(x)));

% Robust sensitivity: winsorize 1% tails
lo = 1; hi = 99;
x_w = min(max(x, prctile(x,lo)), prctile(x,hi));
y_w = min(max(y, prctile(y,lo)), prctile(y,hi));
[RhoW, PvalW] = corr(x_w, y_w, 'Type','Spearman', 'Rows','complete');
fprintf('Spearman after 1%% winsorization: r = %.3f, p = %.3g\n', RhoW, PvalW);

%% ===== Stratified: Spike rate vs SzFreq by EpilepsyType =====
wantTypes = ["Focal","General"];
results = table('Size',[0 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'Type','N','Spearman_r','p_value'});

figure('Color','w'); hold on; grid on; box off;
xlabel('Mean Spike Rate (log spikes/hour)');
ylabel('Mean Seizure Frequency (log units)');
title('Sz freq vs spike rate by epilepsy_type');

for t = 1:numel(wantTypes)
    tname = wantTypes(t);
    idxType = strcmpi(P.EpilepsyType, tname);
    x_t = P.MeanSpikeRate_perHour(idxType);
    y_t = P.MeanSzFreq(idxType);

    mask = isfinite(x_t) & isfinite(y_t);
    x_t = x_t(mask); y_t = y_t(mask);

    if numel(x_t) >= 3
        [rho_t, p_t] = corr(x_t, y_t, 'Type','Spearman', 'Rows','complete');

        minpos_x = min(x_t(x_t>0)); if isempty(minpos_x), minpos_x = 1e-6; end
        minpos_y = min(y_t(y_t>0)); if isempty(minpos_y), minpos_y = 1e-6; end
        eps_x = 0.5*minpos_x; eps_y = 0.5*minpos_y;
        xlog_t = log(x_t + (x_t<=0).*eps_x);
        ylog_t = log(y_t + (y_t<=0).*eps_y);

        scatter(xlog_t, ylog_t, 42, 'filled', 'DisplayName', ...
            sprintf('%s (r=%.3f, p=%.3g, n=%d)', tname, rho_t, p_t, numel(x_t)));

        results = [results; {tname, numel(x_t), rho_t, p_t}]; %#ok<AGROW>
    else
        scatter(nan, nan, 42, 'filled', 'DisplayName', sprintf('%s (insufficient data)', tname));
    end
end
legend('Location','best'); hold off;
disp('=== Stratified Spearman correlations by epilepsy_type ===');
disp(results);

fprintf('Total patients in report: %d\n', numel(unique(R.patient_id)));
fprintf('Patients in spike file: %d\n', numel(unique(S.Patient)));
nBoth = numel(intersect(unique(S.Patient), unique(R.patient_id)));
fprintf('Patients with both EEG and report entries: %d\n', nBoth);

%% ===== Overall: Longest burst vs SzFreq =====
xb = P.MaxLongestBurst_sec;
yb = P.MeanSzFreq;
mask = isfinite(xb) & isfinite(yb); xb = xb(mask); yb = yb(mask);
[r_burst, p_burst] = corr(xb, yb, 'Type','Spearman', 'Rows','complete');
fprintf('Overall (Longest burst vs SzFreq): Spearman r = %.3f, p = %.3g, n = %d\n', r_burst, p_burst, numel(xb));

minpos_xb = min(xb(xb>0)); if isempty(minpos_xb), minpos_xb = 1e-6; end
minpos_yb = min(yb(yb>0)); if isempty(minpos_yb), minpos_yb = 1e-6; end
eps_xb = 0.5*minpos_xb; eps_yb = 0.5*minpos_yb;

xb_log = log(xb + (xb<=0).*eps_xb);
yb_log = log(yb + (yb<=0).*eps_yb);

figure('Color','w');
scatter(xb_log, yb_log, 36, 'filled'); grid on; box off;
xlabel('Longest spike burst (log seconds)');
ylabel('Mean seizure frequency (log units)');
title(sprintf('Overall: Sz freq vs Longest burst (Spearman r=%.3f, p=%.3g, n=%d)', r_burst, p_burst, numel(xb)));

%% ===== Stratified: Longest burst vs SzFreq by EpilepsyType =====
wantTypes = ["Focal","General"];
results_burst = table('Size',[0 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'Type','N','Spearman_r','p_value'});

figure('Color','w'); hold on; grid on; box off;
xlabel('Longest spike burst (log seconds)');
ylabel('Mean seizure frequency (log units)');
title('Sz freq vs longest spike burst by epilepsy_type');

for t = 1:numel(wantTypes)
    tname = wantTypes(t);
    idxT = strcmpi(P.EpilepsyType, tname);

    x_t = P.MaxLongestBurst_sec(idxT);
    y_t = P.MeanSzFreq(idxT);

    mask = isfinite(x_t) & isfinite(y_t);
    x_t = x_t(mask); y_t = y_t(mask);

    if numel(x_t) >= 3
        [rho_t, p_t] = corr(x_t, y_t, 'Type','Spearman', 'Rows','complete');

        minpos_xt = min(x_t(x_t>0)); if isempty(minpos_xt), minpos_xt = 1e-6; end
        minpos_yt = min(y_t(y_t>0)); if isempty(minpos_yt), minpos_yt = 1e-6; end
        eps_xt = 0.5*minpos_xt; eps_yt = 0.5*minpos_yt;

        xlog_t = log(x_t + (x_t<=0).*eps_xt);
        ylog_t = log(y_t + (y_t<=0).*eps_yt);

        scatter(xlog_t, ylog_t, 42, 'filled', 'DisplayName', ...
            sprintf('%s (r=%.3f, p=%.3g, n=%d)', tname, rho_t, p_t, numel(x_t)));

        results_burst = [results_burst; {tname, numel(x_t), rho_t, p_t}]; %#ok<AGROW>
    else
        scatter(nan, nan, 42, 'filled', 'DisplayName', sprintf('%s (insufficient data)', tname));
    end
end
legend('Location','best'); hold off;
disp('=== Stratified Spearman: Longest burst vs SzFreq ===');
disp(results_burst);

%% ===== Attach and normalize per-patient epilepsy_specific =====
if ismember('epilepsy_specific', R.Properties.VariableNames)
    if ~isstring(R.epilepsy_specific), R.epilepsy_specific = string(R.epilepsy_specific); end
    spec_is_ok = ~ismissing(R.epilepsy_specific) & strlength(strtrim(R.epilepsy_specific))>0;
    RtSpec = sortrows(R(spec_is_ok, {'patient_id','epilepsy_specific'}), 'patient_id');
    [uniq_pid_s, ia_s] = unique(RtSpec.patient_id, 'stable');
    PerPatSpecific = table(uniq_pid_s, RtSpec.epilepsy_specific(ia_s), ...
        'VariableNames', {'Patient','EpilepsySpecific'});

    % Merge into P
    P = outerjoin(P, PerPatSpecific, 'Keys','Patient', 'MergeKeys', true);

    % Coalesce to single 'EpilepsySpecific'
    specMask = startsWith(P.Properties.VariableNames, 'EpilepsySpecific', 'IgnoreCase', true);
    specCols = P.Properties.VariableNames(specMask);
    if ~isempty(specCols)
        if ~any(strcmp(specCols,'EpilepsySpecific'))
            P = renamevars(P, specCols{1}, 'EpilepsySpecific');
        end
        baseSpec = 'EpilepsySpecific';
        othersSpec = setdiff(specCols, baseSpec, 'stable');
        if ~isstring(P.(baseSpec)), P.(baseSpec) = string(P.(baseSpec)); end
        P.(baseSpec) = strtrim(P.(baseSpec));
        for ii = 1:numel(othersSpec)
            v = P.(othersSpec{ii});
            if ~isstring(v)
                if iscellstr(v), v = string(v); else, v = string(v); end
            end
            v = strtrim(v);
            fillIdx = (ismissing(P.(baseSpec)) | strlength(P.(baseSpec))==0) & ...
                      (~ismissing(v) & strlength(v)>0);
            P.(baseSpec)(fillIdx) = v(fillIdx);
        end
        if ~isempty(othersSpec), P = removevars(P, othersSpec); end
    else
        P.EpilepsySpecific = strings(height(P),1);
    end

    % Normalize to canonical buckets
    spec_norm = lower(strtrim(P.EpilepsySpecific));
    isTemporal   = contains(spec_norm, "temporal");
    isFrontal    = contains(spec_norm, "frontal");
    hasFocalWord = contains(spec_norm, "focal");
    isUnlocalized = contains(spec_norm, "unlocalized") | ...
                    contains(spec_norm, "unknown") | ...
                    (hasFocalWord & ~isTemporal & ~isFrontal);
    SpecCanon = strings(size(spec_norm));
    SpecCanon(isTemporal)    = "Temporal Lobe";
    SpecCanon(isFrontal)     = "Frontal Lobe";
    SpecCanon(isUnlocalized) = "Unlocalized Focal";
    P.EpilepsySpecific = SpecCanon;

else
    warning('Column "epilepsy_specific" not found in report; skipping specific stratifications.');
    P.EpilepsySpecific = strings(height(P),1);
end

%% ===== Define strata for epilepsy_specific =====
mask_General = ismember(P.EpilepsyType,"General");

mask_Temporal = ismember(P.EpilepsySpecific,"Temporal Lobe");
mask_Frontal  = ismember(P.EpilepsySpecific,"Frontal Lobe");
mask_Unloc    = ismember(P.EpilepsySpecific,"Unlocalized Focal");

% (Optional) Require focal type for the three focal-specific groups:
% mask_Temporal = mask_Temporal & ismember(P.EpilepsyType,"Focal");
% mask_Frontal  = mask_Frontal  & ismember(P.EpilepsyType,"Focal");
% mask_Unloc    = mask_Unloc    & ismember(P.EpilepsyType,"Focal");

Strata = {
    'General (by EpilepsyType)', mask_General;
    'Temporal Lobe',              mask_Temporal;
    'Frontal Lobe',               mask_Frontal;
    'Unlocalized Focal',          mask_Unloc
};

%% ===== Correlations per stratum: spike rate vs sz freq =====
disp('=== Stratified: Spike rate vs Seizure frequency ===');
results_strata_sr = table('Size',[0 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'Group','N','Spearman_r','p_value'});

for i = 1:size(Strata,1)
    gname = Strata{i,1};
    m = Strata{i,2};

    x = P.MeanSpikeRate_perHour(m);
    y = P.MeanSzFreq(m);

    [r,p,n] = doSpearmanPlot(x, y, ...
        sprintf('%s: Spike rate vs SzFreq', gname), ...
        'Mean Spike Rate (log spikes/hour)', ...
        'Mean Seizure Frequency (log units)');

    if n >= 3
        results_strata_sr = [results_strata_sr; {string(gname), n, r, p}]; %#ok<AGROW>
    end
end
disp(results_strata_sr);

%% ===== Correlations per stratum: longest burst vs sz freq =====
disp('=== Stratified: Longest burst vs Seizure frequency ===');
results_strata_lb = table('Size',[0 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'Group','N','Spearman_r','p_value'});

for i = 1:size(Strata,1)
    gname = Strata{i,1};
    m = Strata{i,2};

    x = P.MaxLongestBurst_sec(m);
    y = P.MeanSzFreq(m);

    [r,p,n] = doSpearmanPlot(x, y, ...
        sprintf('%s: Longest burst vs SzFreq', gname), ...
        'Longest spike burst (log seconds)', ...
        'Mean Seizure Frequency (log units)');

    if n >= 3
        results_strata_lb = [results_strata_lb; {string(gname), n, r, p}]; %#ok<AGROW>
    end
end
disp(results_strata_lb);

%% ===== Local helper (place at end of script) =====
function [r, p, n] = doSpearmanPlot(x, y, ttl, xlab, ylab)
    % Keep complete, finite cases
    mask = isfinite(x) & isfinite(y);
    x = x(mask); 
    y = y(mask);
    n = numel(x);

    if n < 3
        r = NaN; p = NaN;
        fprintf('%s: insufficient data (n=%d)\n', ttl, n);
        return;
    end

    % Spearman on raw scale
    [r, p] = corr(x, y, 'Type','Spearman', 'Rows','complete');
    fprintf('%s: Spearman r=%.3f, p=%.3g, n=%d\n', ttl, r, p, n);

    % Safe log transforms for plotting (handle zeros)
    minpos_x = min(x(x>0)); if isempty(minpos_x), minpos_x = 1e-6; end
    minpos_y = min(y(y>0)); if isempty(minpos_y), minpos_y = 1e-6; end
    eps_x = 0.5 * minpos_x; 
    eps_y = 0.5 * minpos_y;

    x_log = log(x + (x<=0).*eps_x);
    y_log = log(y + (y<=0).*eps_y);

    % Plot
    figure('Color','w');
    scatter(x_log, y_log, 36, 'filled'); 
    grid on; box off;
    xlabel(xlab); 
    ylabel(ylab);
    title(sprintf('%s (r=%.3f, p=%.3g, n=%d)', ttl, r, p, n));
end

%% ===== Linear model & 3-panel figure (A/B/C) =====
% Model: log(SpikeRate) ~ z(log(SzFreq)) * EpilepsyType

% --- Prep analysis table (safe logs + z-scoring) ---
% Safe eps for logs (handle zeros)
minpos_sr = min(P.MeanSpikeRate_perHour(P.MeanSpikeRate_perHour>0)); if isempty(minpos_sr), minpos_sr = 1e-6; end
minpos_sf = min(P.MeanSzFreq(P.MeanSzFreq>0));                   if isempty(minpos_sf), minpos_sf = 1e-6; end
eps_sr = 0.5*minpos_sr; 
eps_sf = 0.5*minpos_sf;

T = table;
T.Patient       = P.Patient;
T.EpilepsyType  = categorical(P.EpilepsyType);  % expects "Focal" and "General"
T.logSpikeRate  = log(P.MeanSpikeRate_perHour + (P.MeanSpikeRate_perHour<=0).*eps_sr);
T.logSzFreq     = log(P.MeanSzFreq           + (P.MeanSzFreq<=0).*eps_sf);

% Standardize logSzFreq
[T.zLogSzFreq, mu_logSz, sd_logSz] = zscore(T.logSzFreq);

% Keep only valid rows
ok = isfinite(T.logSpikeRate) & isfinite(T.zLogSzFreq) & ~ismissing(T.EpilepsyType);
T  = T(ok,:);

% --- Fit the linear interaction model ---
mdl = fitlm(T, 'logSpikeRate ~ zLogSzFreq * EpilepsyType');

% === Figure layout ===
f = figure('Color','w','Name','SpikeRate ~ z(log(SzFreq)) * EpilepsyType','Position',[100 100 1200 700]);
tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');

%% Panel A: Raw data (log–log), faceted by EpilepsyType
nexttile(1);
hold on; grid on; box off;
types = categories(T.EpilepsyType);
clr = lines(numel(types));
for i = 1:numel(types)
    m = T.EpilepsyType == types{i};
    scatter(T.logSzFreq(m), T.logSpikeRate(m), 18, 'filled', 'MarkerFaceAlpha', 0.35, 'DisplayName', char(types{i}));
    % add within-type robust fit line on the log–log scale for context
    if nnz(m) >= 3
        X = [ones(nnz(m),1), T.logSzFreq(m)];
        b = X\T.logSpikeRate(m);
        xx = linspace(min(T.logSzFreq(m)), max(T.logSzFreq(m)), 100)';
        yy = b(1) + b(2)*xx;
        plot(xx, yy, 'LineWidth', 1.5, 'Color', clr(i,:));
    end
end
xlabel('log Seizure Frequency');
ylabel('log Spike Rate');
title('A. Raw data (log–log), faceted colors by type');
legend('Location','best');

%% Panel B: Model-predicted lines ± 95% CI (interaction)
nexttile(2);
hold on; grid on; box off;

% Create a zLogSzFreq grid spanning observed range
zmin = min(T.zLogSzFreq); zmax = max(T.zLogSzFreq);
zgrid = linspace(zmin, zmax, 200)';

% Prediction per type
for i = 1:numel(types)
    newT = table;
    newT.zLogSzFreq   = zgrid;
    newT.EpilepsyType = repmat(categorical({types{i}}, types), numel(zgrid),1);
    [yhat, yci] = predict(mdl, newT);  % CI for mean prediction

    % Plot line and 95% CI band on log scale (as modeled)
    plot(zgrid, yhat, 'LineWidth', 2, 'DisplayName', sprintf('%s (fit)', types{i}), 'Color', clr(i,:));
    % Shaded CI
    fill([zgrid; flipud(zgrid)], [yci(:,1); flipud(yci(:,2))], clr(i,:), ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility','off');
end
xlabel('z(log Seizure Frequency)');
ylabel('Predicted log Spike Rate');
title('B. Model predictions with 95% CI (interaction)');
legend('Location','best');

%% Panel C: Coefficient / CI (“forest”) plot
nexttile(3, [1 2]); % span full bottom row
hold on; grid on; box off;

coefNames = mdl.CoefficientNames;          % cellstr
beta      = mdl.Coefficients.Estimate;     % estimates
ci        = coefCI(mdl, 0.05);             % Nx2
se        = mdl.Coefficients.SE;           %#ok<NASGU> (available if you prefer)

% Order: intercept first, then main effects, then interaction (keep as-is)
N = numel(beta);
ypos = N:-1:1;  % top to bottom

% horizontal error bars
for k = 1:N
    line([ci(k,1), ci(k,2)], [ypos(k), ypos(k)], 'LineWidth', 3, 'Color', [0.3 0.3 0.3]);
end
plot(beta, ypos, 'o', 'MarkerFaceColor',[0 0.45 0.74], 'MarkerEdgeColor','k', 'MarkerSize',6);

% zero reference
xline(0,'--','Color',[0.5 0.5 0.5]);

% y-axis labels
set(gca, 'YTick', ypos, 'YTickLabel', coefNames, 'YLim',[0.5 N+0.5]);
xlabel('Coefficient estimate (log Spike Rate units)');
title('C. Model coefficients with 95% CIs');

% --- Annotation: brief model summary in the figure ---
txt = sprintf('n = %d | R^2 adj = %.3f', mdl.NumObservations, mdl.Rsquared.Adjusted);
annotation(f,'textbox',[0.68 0.05 0.3 0.06],'String',txt,'FitBoxToText','on','EdgeColor','none');

%% (Optional) Back-transform helper for reporting (not plotted):
% To report on the original SpikeRate scale:
%   If Δ on the y-axis is d (in log units), then multiplicative change ≈ exp(d).
% Example: a coefficient of 0.69 ~ 2× spike rate.
