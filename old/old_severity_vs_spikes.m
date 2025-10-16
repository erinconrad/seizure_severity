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

%% ===== Load report and compute per-patient mean sz_freq =====
R = readtable(reportFile, 'TextType','string');
if ~isnumeric(R.patient_id), R.patient_id = double(str2double(R.patient_id)); end

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

%% ===== Join per-patient tables =====
P = innerjoin(Sg, Rg, 'Keys','Patient');
P = P(~isnan(P.MeanSzFreq) & ~isnan(P.MeanSpikeRate_perHour), :);

%% ===== Prepare vectors (handle zeros before log) =====
x = P.MeanSpikeRate_perHour;   % spikes/hour
y = P.MeanSzFreq;              % e.g., seizures/month

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
fprintf('Patients with both measures: %d\n', numel(x));
fprintf('Spearman (raw scale): r = %.3f, p = %.3g\n', Rho, Pval);

%% ===== Log–log scatter (for visualization) =====
figure('Color','w');
scatter(x_log, y_log, 36, 'filled'); grid on; box off;
xlabel('Mean Spike Rate (log spikes/hour)');
ylabel('Mean Seizure Frequency (log units)');
title(sprintf('Sz freq vs spike rate (Spearman r=%.3f, p=%.3g; n=%d)', Rho, Pval, numel(x)));

% Optional: add LOWESS trend in log space for visual aid
hold on;
[sl, idx] = sort(x_log); yl = y_log(idx);
%plot(sl, smooth(yl, max(5,round(0.1*numel(yl)))), 'LineWidth', 1.5);
hold off;

%% ===== (Optional) Robust sensitivity: winsorize top/bottom 1%% and recompute =====
lo = 1; hi = 99;  % percentiles
x_w = x; y_w = y;
px = prctile(x,[lo hi]); py = prctile(y,[lo hi]);
x_w = min(max(x_w, px(1)), px(2));
y_w = min(max(y_w, py(1)), py(2));
[RhoW, PvalW] = corr(x_w, y_w, 'Type','Spearman', 'Rows','complete');
fprintf('Spearman after 1%% winsorization: r = %.3f, p = %.3g\n', RhoW, PvalW);


%% ===== Save patient-level CSV =====
%writetable(P, outPatient);
%fprintf('Saved patient-level table to: %s\n', outPatient);
