%% ===== Paths =====
outCsv     = '../data/SN_counts/spike_counts_summary.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-14_1505.csv';
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

%% ===== Compute per-patient mean sz_freq among VALID patients only =====
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
    MeanSzFreq(k) = mean(vals,'omitnan');
end

Rg = table(pids, MeanSzFreq, 'VariableNames', {'Patient','MeanSzFreq'});

% Attach epilepsy_type to Rg for reference
Rg = innerjoin(Rg, PerPatType(~isBad, :), 'Keys','Patient');

%% ===== Join per-patient tables (now only valid epilepsy_type patients remain) =====
P = innerjoin(Sg, Rg, 'Keys','Patient');
P = P(~isnan(P.MeanSzFreq) & ~isnan(P.MeanSpikeRate_perHour), :);

%% ===== Prepare vectors (handle zeros before log) =====
x = P.MeanSpikeRate_perHour;   % spikes/hour
y = P.MeanSzFreq;              % seizures/month (or whatever units sz_freqs represent)

% Keep complete, finite cases
mask = isfinite(x) & isfinite(y);
x = x(mask); y = y(mask);

% Epsilon: half of the smallest positive value in each variable (fallback tiny)
minpos_x = min(x(x>0)); if isempty(minpos_x), minpos_x = 1e-6; end
minpos_y = min(y(y>0)); if isempty(minpos_y), minpos_y = 1e-6; end
eps_x = 0.5*minpos_x; eps_y = 0.5*minpos_y;

x_log = log(x + (x<=0).*eps_x);
y_log = log(y + (y<=0).*eps_y);

%% ===== Spearman correlation (rank-based) =====
[Rho, Pval] = corr(x, y, 'Type','Spearman', 'Rows','complete');
fprintf('Patients with both measures (post-filter): %d\n', numel(x));
fprintf('Spearman (raw scale): r = %.3f, p = %.3g\n', Rho, Pval);

%% ===== Log–log scatter (for visualization) =====
figure('Color','w');
scatter(x_log, y_log, 36, 'filled'); grid on; box off;
xlabel('Mean Spike Rate (log spikes/hour)');
ylabel('Mean Seizure Frequency (log units)');
title(sprintf('Sz freq vs spike rate (Spearman r=%.3f, p=%.3g; n=%d)', Rho, Pval, numel(x)));

%% ===== (Optional) Robust sensitivity: winsorize top/bottom 1%% and recompute =====
lo = 1; hi = 99;  % percentiles
x_w = x; y_w = y;
px = prctile(x,[lo hi]); py = prctile(y,[lo hi]);
x_w = min(max(x_w, px(1)), px(2));
y_w = min(max(y_w, py(1)), py(2));
[RhoW, PvalW] = corr(x_w, y_w, 'Type','Spearman', 'Rows','complete');
fprintf('Spearman after 1%% winsorization: r = %.3f, p = %.3g\n', RhoW, PvalW);

%% ===== Save patient-level CSV (includes EpilepsyType for future stratification) =====
%writetable(P, outPatient);
%fprintf('Saved patient-level table to: %s\n', outPatient);
%% ===== Stratified correlations: Focal vs General =====
% (Place this AFTER the table P is created and cleaned in your script.)

wantTypes = ["Focal","General"]; % labels must match your EpilepsyType values

results = table('Size',[0 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'Type','N','Spearman_r','p_value'});

figure('Color','w'); hold on; grid on; box off;
xlabel('Mean Spike Rate (log spikes/hour)');
ylabel('Mean Seizure Frequency (log units)');
title('Sz freq vs spike rate by epilepsy_type');

for t = 1:numel(wantTypes)
    tname = wantTypes(t);

    % Select patients of this epilepsy type
    idxType = strcmpi(P.EpilepsyType, tname);
    x_t = P.MeanSpikeRate_perHour(idxType);
    y_t = P.MeanSzFreq(idxType);

    % Keep complete, finite cases
    mask = isfinite(x_t) & isfinite(y_t);
    x_t = x_t(mask); y_t = y_t(mask);

    % Only compute if we have enough data
    if numel(x_t) >= 3
        % Spearman on raw scale
        [rho_t, p_t] = corr(x_t, y_t, 'Type','Spearman', 'Rows','complete');

        % Log transform for plotting with safe eps
        minpos_x = min(x_t(x_t>0)); if isempty(minpos_x), minpos_x = 1e-6; end
        minpos_y = min(y_t(y_t>0)); if isempty(minpos_y), minpos_y = 1e-6; end
        eps_x = 0.5*minpos_x; eps_y = 0.5*minpos_y;

        xlog_t = log(x_t + (x_t<=0).*eps_x);
        ylog_t = log(y_t + (y_t<=0).*eps_y);

        % Plot this group
        scatter(xlog_t, ylog_t, 42, 'filled', 'DisplayName', ...
            sprintf('%s (r=%.3f, p=%.3g, n=%d)', tname, rho_t, p_t, numel(x_t)));

        % Store results
        results = [results; {tname, numel(x_t), rho_t, p_t}]; %#ok<AGROW>
    else
        % Not enough data to run correlation
        scatter(nan, nan, 42, 'filled', 'DisplayName', ...
            sprintf('%s (insufficient data)', tname));
    end
end

legend('Location','best');
hold off;

% Show a compact summary in the console
disp('=== Stratified Spearman correlations by epilepsy_type ===');
disp(results);

fprintf('Total patients in report: %d\n', numel(unique(R.patient_id)));
fprintf('Patients in spike file: %d\n', numel(unique(S.Patient)));

nBoth = numel(intersect(unique(S.Patient), unique(R.patient_id)));
fprintf('Patients with both EEG and report entries: %d\n', nBoth);

noSzFreq = sum(isnan(Rg.MeanSzFreq));
noSpike  = sum(isnan(Sg.MeanSpikeRate_perHour));
fprintf('Patients missing sz freq: %d, missing spike rate: %d\n', noSzFreq, noSpike);
