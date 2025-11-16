%% Parameters
only_ambulatories = 1; % 1 = only amb, 2 = only routine, 0 = all
which_runs = 1;        % 1 = first run (~first hour), 2 = first 24 runs (~24h), 0 = whole file

%% ===== Paths you gave =====
outCsv     = '../data/SN_counts/spike_counts_summary.csv';   % from previous step
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-08_0557.csv';

% Column in the report that encodes present/absent (string or categorical)
reportColName = 'report_SPORADIC_EPILEPTIFORM_DISCHARGES';

%% ===== Load spike-count summary =====
S = readtable(outCsv, 'TextType','string');

% Keep the original duration to classify amb vs routine (your previous rule)
if only_ambulatories == 1
    S(S.Duration_sec < 3600*12,:) = [];
elseif only_ambulatories == 2
    S(S.Duration_sec > 3600*12,:) = [];
end

% Choose which columns to analyze based on which_runs
switch which_runs
    case 1   % first run (~first hour)
        colCount = "FirstRun_Spikes";
        colLeft  = "FirstRun_Left_Spikes";
        colRight = "FirstRun_Right_Spikes";
        colDur   = "FirstRun_Duration_sec";
        segmentLabel = 'First run (~1h)';
    case 2   % first 24 runs (~24 hours)
        colCount = "First24Runs_Spikes";
        colLeft  = "First24Runs_Left_Spikes";
        colRight = "First24Runs_Right_Spikes";
        colDur   = "First24Runs_Duration_sec";
        segmentLabel = 'First 24 runs (~24h)';
    otherwise % 0 → whole file
        colCount = "Total_Spikes";
        colLeft  = "Left_Spikes";
        colRight = "Right_Spikes";
        colDur   = "Duration_sec";
        segmentLabel = 'Whole file';
end

% Sanity: ensure the selected columns exist and are numeric
needCols = [colCount colLeft colRight colDur];
assert(all(ismember(needCols, S.Properties.VariableNames)), ...
    'Selected columns not found in summary: %s', strjoin(needCols, ', '));
for nm = needCols
    if ~isnumeric(S.(nm))
        S.(nm) = double(S.(nm));
    end
end

%% ===== Aggregate to one row per (Patient, Session) without groupsummary =====
[grp, pKeys, sKeys] = findgroups(S.Patient, S.Session); %#ok<ASGLU>

nanSum = @(x) sum(x(~isnan(x)));
Sum_Total = splitapply(nanSum, S.(colCount), grp);
Sum_Left  = splitapply(nanSum, S.(colLeft),  grp);
Sum_Right = splitapply(nanSum, S.(colRight), grp);
Sum_Dur   = splitapply(nanSum, S.(colDur),   grp);

Sgrp = table(pKeys, sKeys, Sum_Total, Sum_Left, Sum_Right, Sum_Dur, ...
    'VariableNames', {'Patient','Session','Sum_Total','Sum_Left','Sum_Right','Sum_Dur'});

% spikes/min
Sgrp.SpikesPerMin = nan(height(Sgrp),1);
validDur = Sgrp.Sum_Dur > 0;
Sgrp.SpikesPerMin(validDur) = (Sgrp.Sum_Total(validDur) ./ Sgrp.Sum_Dur(validDur)) * 60;

%% ===== Load report CSV and normalize keys =====
R = readtable(reportFile, 'TextType','string');

% Make sure patient_id and session_number exist and are numeric
assert(all(ismember({'patient_id','session_number'}, R.Properties.VariableNames)), ...
    'Expected columns patient_id and session_number in the report CSV.');

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
repNorm(~(isPresent | isAbsent)) = missing;

R.ReportStatus = categorical(repNorm, ["absent","present"]);  % fixed order
R2 = R(~ismissing(R.ReportStatus), {'Patient','Session','ReportStatus'});

%% ===== Join spike rates ↔ report status on (Patient, Session) =====
J = innerjoin(Sgrp(:,{'Patient','Session','SpikesPerMin'}), R2, 'Keys', {'Patient','Session'});
fprintf('Matched %d (Patient, Session) rows for analysis (%s).\n', height(J), segmentLabel);

%% ===== Compare spikes/min between groups =====
absMask = (J.ReportStatus=="absent");
preMask = (J.ReportStatus=="present");

x_abs = J.SpikesPerMin(absMask);
x_pre = J.SpikesPerMin(preMask);

fprintf('N (absent)  = %d\n', numel(x_abs));
fprintf('N (present) = %d\n', numel(x_pre));

% Descriptive stats
stats = @(x) [mean(x,'omitnan'), median(x,'omitnan'), std(x,'omitnan'), iqr(x), prctile(x,25), prctile(x,75)];
sa = stats(x_abs); sp = stats(x_pre);
fprintf('\nSpikes/min — ABSENT:  mean=%.3f  median=%.3f  sd=%.3f  IQR=%.3f  Q1=%.3f  Q3=%.3f\n', sa);
fprintf('Spikes/min — PRESENT: mean=%.3f  median=%.3f  sd=%.3f  IQR=%.3f  Q1=%.3f  Q3=%.3f\n', sp);

% Nonparametric test (Mann–Whitney U) on rates
p_rs = ranksum(x_abs, x_pre);
fprintf('\nRank-sum test (absent vs present), Spikes/min — %s: p = %.3g\n', segmentLabel, p_rs);

%% ===== Plot mean ± SEM for Spikes/min =====
means = [mean(x_abs,'omitnan'), mean(x_pre,'omitnan')];
sems  = [std(x_abs,'omitnan')/sqrt(max(1,sum(~isnan(x_abs)))), ...
         std(x_pre,'omitnan')/sqrt(max(1,sum(~isnan(x_pre))))];

figure('Color','w');
bar(means); hold on;
errorbar(1:2, means, sems, 'linestyle','none', 'LineWidth', 1.5);
set(gca,'XTick',1:2,'XTickLabel',{'Absent','Present'});
ylabel(sprintf('Spikes/min (%s)', segmentLabel));
title(sprintf('Spikes/min by report status — %s (p = %.3g, rank-sum)', segmentLabel, p_rs));
grid on; box off;

%% ===== (Optional) Save merged table for auditing =====
% mergedOut = '../data/SN_counts/spike_rates_with_report_status.csv';
% writetable(J, mergedOut);
% fprintf('Saved merged table to: %s\n', mergedOut);
