%% Parameters
only_ambulatories = 2; % 1 = only amb, 2 = only routine, 0 = all
which_runs = 1; % 1 = first hour, 2 = first 24 hours, 0 = all

%% ===== Paths you gave =====
outCsv     = '../data/SN_counts/spike_counts_summary.csv';   % from previous step
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-08_0557.csv';

% Column in the report that encodes present/absent (string or categorical)
reportColName = 'report_SPORADIC_EPILEPTIFORM_DISCHARGES';

%% ===== Load spike-count summary =====
S = readtable(outCsv, 'TextType','string');

if only_ambulatories == 1
    S(S.Duration_sec < 3600*12,:) = [];
elseif only_ambulatories == 2
    S(S.Duration_sec > 3600*12,:) = [];
end

% If multiple CSVs exist per (Patient, Session), aggregate to one row per pair
Sgrp = groupsummary(S, {'Patient','Session'}, 'sum', ...
    {'Total_Spikes','Left_Spikes','Right_Spikes'});
Sgrp.Properties.VariableNames = {'Patient','Session','GroupCount','Total_Spikes','Left_Spikes','Right_Spikes'};
Sgrp.GroupCount = [];

%% ===== Load report CSV and normalize keys =====
R = readtable(reportFile, 'TextType','string');

% Make sure patient_id and session_number exist and are numeric
assert(all(ismember({'patient_id','session_number'}, R.Properties.VariableNames)), ...
    'Expected columns patient_id and session_number in the report CSV.');

% Coerce to numeric (they may import as string)
if ~isnumeric(R.patient_id);     R.patient_id     = double(str2double(R.patient_id)); end
if ~isnumeric(R.session_number); R.session_number = double(str2double(R.session_number)); end

R.Patient = R.patient_id;
R.Session = R.session_number;

% Ensure the report present/absent column exists
assert(ismember(reportColName, R.Properties.VariableNames), ...
    'Report column "%s" not found in the report CSV.', reportColName);

% Normalize present/absent
raw = lower(strtrim(string(R.(reportColName))));
isPresent = ismember(raw, ["present","yes","y","1","true","pos","positive"]);
isAbsent  = ismember(raw, ["absent","no","n","0","false","neg","negative","none","normal"]);
repNorm = strings(height(R),1);
repNorm(isPresent) = "present";
repNorm(isAbsent)  = "absent";
repNorm(~(isPresent | isAbsent)) = missing;   % leave others as missing

R.ReportStatus = categorical(repNorm, ["absent","present"]);  % fixed order
R2 = R(~ismissing(R.ReportStatus), {'Patient','Session','ReportStatus'});

%% ===== Join spike counts ↔ report status on (Patient, Session) =====
J = innerjoin(Sgrp, R2, 'Keys', {'Patient','Session'});  % keeps only matched rows
fprintf('Matched %d (Patient, Session) rows.\n', height(J));

%% ===== Compare Total_Spikes between groups =====
absMask = (J.ReportStatus=="absent");
preMask = (J.ReportStatus=="present");

x_abs = J.Total_Spikes(absMask);
x_pre = J.Total_Spikes(preMask);

fprintf('N (absent)  = %d\n', numel(x_abs));
fprintf('N (present) = %d\n', numel(x_pre));

% Descriptive stats
stats = @(x) [mean(x,'omitnan'), median(x,'omitnan'), std(x,'omitnan'), iqr(x), prctile(x,25), prctile(x,75)];
sa = stats(x_abs); sp = stats(x_pre);
fprintf('\nTotal spikes — ABSENT:  mean=%.2f  median=%.2f  sd=%.2f  IQR=%.2f  Q1=%.2f  Q3=%.2f\n', sa);
fprintf('Total spikes — PRESENT: mean=%.2f  median=%.2f  sd=%.2f  IQR=%.2f  Q1=%.2f  Q3=%.2f\n', sp);

% Nonparametric test (Mann–Whitney U)
p_rs = ranksum(x_abs, x_pre);
fprintf('\nRank-sum test (absent vs present), Total_Spikes: p = %.3g\n', p_rs);

%% ===== Optional: also compare Left and Right separately =====
% Uncomment if desired
% p_rs_L = ranksum(J.Left_Spikes(absMask),  J.Left_Spikes(preMask));
% p_rs_R = ranksum(J.Right_Spikes(absMask), J.Right_Spikes(preMask));
% fprintf('Rank-sum Left_Spikes p=%.3g | Right_Spikes p=%.3g\n', p_rs_L, p_rs_R);

%% ===== Plot mean ± SEM for Total_Spikes =====
means = [mean(x_abs,'omitnan'), mean(x_pre,'omitnan')];
sems  = [std(x_abs,'omitnan')/sqrt(max(1,numel(x_abs))), std(x_pre,'omitnan')/sqrt(max(1,numel(x_pre)))];

figure('Color','w');
bar(means); hold on;
errorbar(1:2, means, sems, 'linestyle','none', 'LineWidth', 1.5);
set(gca,'XTick',1:2,'XTickLabel',{'Absent','Present'});
ylabel('Total spikes per (Patient, Session)');
title(sprintf('Spike counts by report status (p = %.3g, rank-sum)', p_rs));
grid on; box off;

%% ===== Save merged table for auditing =====
%mergedOut = '../data/SN_counts/spike_counts_with_report_status.csv';
%writetable(J, mergedOut);
%fprintf('Saved merged table to: %s\n', mergedOut);
