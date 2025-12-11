function summarize_patients
% Summarize patients:
%   - N patients (unique in spike summary)
%   - N with epilepsy (exclude NEAD/uncertain/unknown)
%   - N General, N Temporal, N Frontal (by patient)
%   - Median and range of number of EEGs across patients
%
% Uses:
%   - outCsv: spike summary with Patient, Session (rows = EEGs)
%   - reportFile: patient-level epilepsy_type and epilepsy_specific

%% ============ CONFIG ============
outCsv     = '../data/SN_counts/spike_counts_summary.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-10_1443.csv';
saveCsv    = '../output/patient_summary.csv';   % written at end

% Definitions of "non-epilepsy" buckets (case-insensitive)
badTypes = lower([
    "Non-Epileptic Seizure Disorder"
    "Uncertain if Epilepsy"
    "Unknown or MRN not found"
    ""  % empty/missing
]);

%% ============ LOAD EEG SUMMARY ============
S = readtable(outCsv, 'TextType','string', 'VariableNamingRule','preserve');
assert(all(ismember({'Patient','Session'}, S.Properties.VariableNames)), ...
    'Spike summary must contain Patient and Session.');
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end

% Universe of patients = those that appear in the EEG summary
patients_all = unique(S.Patient(~isnan(S.Patient)));
n_patients   = numel(patients_all);

% EEGs per patient
[~, ~, gi] = unique(S.Patient);
eegs_per_patient = accumarray(gi, 1);
% but accumarray counts NaN group too if present; restrict to valid patients
eegs_per_patient = eegs_per_patient(~isnan(unique(S.Patient))); %#ok<NASGU>

% Safer: explicitly compute from valid list
eegs_per_patient = arrayfun(@(p) sum(S.Patient==p), patients_all);

med_eegs = median(eegs_per_patient);
min_eegs = min(eegs_per_patient);
max_eegs = max(eegs_per_patient);

%% ============ LOAD REPORT (types) ============
R = readtable(reportFile, 'TextType','string', 'VariableNamingRule','preserve');
assert(all(ismember({'patient_id','epilepsy_type','epilepsy_specific'}, R.Properties.VariableNames)), ...
    'Report sheet missing required columns.');

if ~isnumeric(R.patient_id), R.patient_id = double(str2double(string(R.patient_id))); end
etype = string(R.epilepsy_type);
espec = string(R.epilepsy_specific);

% Per-patient first non-missing epilepsy_type/specific
T = table(R.patient_id, etype, espec, 'VariableNames', {'Patient','EpilepsyType','EpilepsySpecific'});
T = T(~isnan(T.Patient), :);
T = sortrows(T, 'Patient');
[up, ia] = unique(T.Patient, 'stable');
PerPat = T(ia, :);

% Join onto EEG-patient universe (so counts match patients with EEGs)
PP = outerjoin(table(patients_all, 'VariableNames', {'Patient'}), PerPat, ...
    'Keys', 'Patient', 'MergeKeys', true);

% Normalize strings
PP.EpilepsyType     = string(PP.EpilepsyType);
PP.EpilepsySpecific = string(PP.EpilepsySpecific);
etLower = lower(strtrim(PP.EpilepsyType));
esLower = lower(strtrim(PP.EpilepsySpecific));

% "Has epilepsy" = NOT in badTypes list
is_epilepsy = ~ismember(etLower, badTypes);

% Bucket to General/Temporal/Frontal (like your other scripts)
is_general  = contains(etLower, 'general');
is_temporal = contains(esLower, 'temporal');
is_frontal  = contains(esLower, 'frontal');

% Only count subtype buckets among "with epilepsy"
n_general  = sum(is_epilepsy & is_general);
n_temporal = sum(is_epilepsy & ~is_general & is_temporal);
n_frontal  = sum(is_epilepsy & ~is_general & ~is_temporal & is_frontal);

n_with_epilepsy = sum(is_epilepsy);

%% ============ BUILD SUMMARY TABLE ============
Summary = table( ...
    n_patients, ...
    n_with_epilepsy, ...
    n_general, ...
    n_frontal, ...
    n_temporal, ...
    med_eegs, ...
    min_eegs, ...
    max_eegs, ...
    'VariableNames', { ...
        'N_patients', ...
        'N_with_epilepsy', ...
        'N_general', ...
        'N_frontal', ...
        'N_temporal', ...
        'Median_EEGs_per_patient', ...
        'Min_EEGs_per_patient', ...
        'Max_EEGs_per_patient' ...
    } ...
);

disp(Summary);

% Save to CSV
writetable(Summary, saveCsv);
fprintf('Wrote patient summary to: %s\n', saveCsv);

end
