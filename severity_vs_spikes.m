%% ===== Paths =====
outCsv     = '../data/SN_counts/spike_counts_summary.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-07_1611.csv';
outPatient = '../output/patient_szfreq_vs_spikes.csv';

%% ===== Load spike counts =====
S = readtable(outCsv, 'TextType','string');  % EEG_Name, Patient, Session, Total_Spikes, Left_Spikes, Right_Spikes
% Per-patient average of Total_Spikes across EEGs
Sg = groupsummary(S, 'Patient', 'mean', 'Total_Spikes');
Sg.Properties.VariableNames = {'Patient','GroupCount','MeanTotalSpikes'};
Sg.GroupCount = [];

% Ensure numeric Patient
if ~isnumeric(Sg.Patient)
    Sg.Patient = double(str2double(string(Sg.Patient)));
end

%% ===== Load report and compute per-patient mean sz_freq =====
R = readtable(reportFile, 'TextType','string');  % has patient_id and sz_freqs columns
% Coerce IDs to numeric
if ~isnumeric(R.patient_id), R.patient_id = double(str2double(R.patient_id)); end

% Parse sz_freqs per patient
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
        % Replace 'null' with NaN, then JSON-decode
        s = regexprep(raw, 'null', 'NaN', 'ignorecase');
        try
            v = jsondecode(char(s));      % numeric array with NaNs
            vals = [vals; v(:)]; %#ok<AGROW>
        catch
            % Fallback: extract numbers if JSON is malformed
            nums = regexp(s, '[-+]?\d+(\.\d+)?([eE][-+]?\d+)?', 'match');
            if ~isempty(nums)
                vals = [vals; str2double(string(nums))]; %#ok<AGROW>
            end
        end
    end
    MeanSzFreq(k) = mean(vals, 'omitnan');
end

Rg = table(pids, MeanSzFreq, 'VariableNames', {'Patient','MeanSzFreq'});

%% ===== Join per-patient tables =====
P = innerjoin(Sg, Rg, 'Keys','Patient');
P = P(~isnan(P.MeanSzFreq) & ~isnan(P.MeanTotalSpikes), :);

%% ===== Correlate and report =====
[Rho, Pval] = corr(P.MeanTotalSpikes, P.MeanSzFreq, 'Type','Pearson', 'Rows','complete');
fprintf('Patients with both measures: %d\n', height(P));
fprintf('Pearson correlation: r = %.3f, p = %.3g\n', Rho, Pval);

% Optional scatter
figure('Color','w');
scatter(P.MeanTotalSpikes, P.MeanSzFreq, 36, 'filled');
xlabel('Mean Total Spikes (per patient)');
ylabel('Mean Seizure Frequency (per patient)');
title(sprintf('Mean sz freq vs mean total spikes (r=%.3f, p=%.3g)', Rho, Pval));
grid on; box off;

%% ===== Save patient-level CSV =====
%writetable(P, outPatient);
%fprintf('Saved patient-level table to: %s\n', outPatient);
