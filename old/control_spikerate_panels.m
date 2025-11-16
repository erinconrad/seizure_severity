function control_spikerate_panels(which_runs)
% control_spikerate_panels(which_runs)
%   0 = whole file (default)
%   1 = first run (~1h)
%   2 = first 24 runs (~24h)
%
% One figure with 3 panels, plotted as log10(spikes/min), ORDERED:
%   C) Report Present vs Absent (reportColName)
%   A) Epilepsy (any, excluding badTypes) vs NESD
%   B) General vs Temporal vs Frontal
%
% Also:
%   - Global dotted zero line at EPS_PER_MIN (default 1e-3 spikes/min → y = -3)
%   - Fixed y-limits across all panels
%   - Significance bars + p-values for C and A (two-group comparisons)

if nargin < 1, which_runs = 0; end

%% ======================= CONFIG =======================
outCsv      = '../../data/SN_counts/spike_counts_summary.csv';
reportFile  = '../../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-10_1443.csv';

figOut      = '../../figures/spikerate_control_panels.png'; % single combined figure

only_amb    = 0;   % 1 = ONLY ambulatory (>=12h); 2 = ONLY routine (<=12h); 0 = all
only_outpt  = 1;   % 1 = keep EEGs with acquired_on containing 'Spe' or 'radnor'
reportColName = 'report_SPORADIC_EPILEPTIFORM_DISCHARGES'; % Panel C

badTypes = lower(["Uncertain if Epilepsy","Unknown or MRN not found",""]); % excluded from "Epilepsy"
NESD_LABEL = "Non-Epileptic Seizure Disorder";
canon3 = ["General","Temporal","Frontal"];  % Panel B mapping

% Unified plotting limits (log10(spikes/min))
Y_LIMS = [-3.2 2];
EPS_PER_MIN = 1e-3;                % global zero substitute → yZero = -3
Y_ZERO = log10(EPS_PER_MIN);       % location of dotted line

%% ====================== LOAD SPIKE SUMMARY ======================
S = readtable(outCsv,'TextType','string','VariableNamingRule','preserve');

switch which_runs
    case 1, colCount="FirstRun_Spikes";    colDur="FirstRun_Duration_sec";    segLabel='First run (~1h)';
    case 2, colCount="First24Runs_Spikes"; colDur="First24Runs_Duration_sec"; segLabel='First 24 runs (~24h)';
    otherwise, colCount="Total_Spikes";    colDur="Duration_sec";             segLabel='Whole file';
end
segLabel = string(segLabel);

needS = {'Patient','Session','Duration_sec', char(colCount), char(colDur)};
assert(all(ismember(needS, S.Properties.VariableNames)), 'Missing required columns in spike CSV.');
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end

% Ambulatory/routine filter based on FULL duration
if only_amb == 1
    S(S.Duration_sec < 12*3600,:) = [];
elseif only_amb == 2
    S(S.Duration_sec > 12*3600,:) = [];
end

% Compute spikes/hour for selected segment
if ~isnumeric(S.(colCount)), S.(colCount) = double(S.(colCount)); end
if ~isnumeric(S.(colDur)),   S.(colDur)   = double(S.(colDur));   end
S.SpikeRate_perHour = nan(height(S),1);
okDur = S.(colDur) > 0;
S.SpikeRate_perHour(okDur) = (S.(colCount)(okDur) ./ S.(colDur)(okDur)) * 3600;

%% ====================== LOAD REPORT ======================
R = readtable(reportFile,'TextType','string','VariableNamingRule','preserve');
assert(ismember('patient_id',R.Properties.VariableNames),'patient_id missing');
assert(ismember('session_number',R.Properties.VariableNames),'session_number missing');
assert(ismember('epilepsy_type',R.Properties.VariableNames),'epilepsy_type missing');

if ~isnumeric(R.patient_id),     R.patient_id     = double(str2double(string(R.patient_id))); end
if ~isnumeric(R.session_number), R.session_number = double(str2double(string(R.session_number))); end
if ~isstring(R.epilepsy_type),   R.epilepsy_type  = string(R.epilepsy_type); end

% Optional outpatient filter via acquired_on
if only_outpt == 1
    if ismember('acquired_on', R.Properties.VariableNames)
        acq = lower(strtrim(string(R.acquired_on)));
        keep_out = ~ismissing(acq) & (contains(acq,"spe") | contains(acq,"radnor"));
        R = R(keep_out,:);
    else
        warning('Report lacks "acquired_on"; only_outpt filter skipped.');
    end
end

% Per-patient epilepsy_type (stable)
Rt = sortrows(R(~(ismissing(R.epilepsy_type) | strlength(strtrim(R.epilepsy_type))==0), ...
                {'patient_id','epilepsy_type'}), 'patient_id');
[uniq_pid, ia] = unique(Rt.patient_id, 'stable');
PerPatType = table(uniq_pid, Rt.epilepsy_type(ia), 'VariableNames', {'Patient','EpilepsyType'});

% Per-patient epilepsy_specific, if available
hasSpec = ismember('epilepsy_specific', R.Properties.VariableNames);
if hasSpec
    if ~isstring(R.epilepsy_specific), R.epilepsy_specific = string(R.epilepsy_specific); end
    Rs = sortrows(R(~(ismissing(R.epilepsy_specific) | strlength(strtrim(R.epilepsy_specific))==0), ...
                    {'patient_id','epilepsy_specific'}), 'patient_id');
    [up_s, ia_s] = unique(Rs.patient_id,'stable');
    PerPatSpec = table(up_s, Rs.epilepsy_specific(ia_s), 'VariableNames', {'Patient','EpilepsySpecific'});
else
    PerPatSpec = table([],[],'VariableNames',{'Patient','EpilepsySpecific'});
end

PerPat = outerjoin(PerPatType, PerPatSpec, 'Keys','Patient','MergeKeys',true);
if ~isstring(PerPat.EpilepsyType),    PerPat.EpilepsyType    = string(PerPat.EpilepsyType);    end
if ~isstring(PerPat.EpilepsySpecific),PerPat.EpilepsySpecific= string(PerPat.EpilepsySpecific); end
PerPat.EpilepsyType     = strtrim(PerPat.EpilepsyType);
PerPat.EpilepsySpecific = strtrim(PerPat.EpilepsySpecific);

% Map to 3 canonical categories using type + specific
et = lower(PerPat.EpilepsyType);
es = lower(PerPat.EpilepsySpecific);
isGeneral  = contains(et, "general");
isTemporal = contains(es, "temporal");
isFrontal  = contains(es, "frontal");

EpiType4 = strings(height(PerPat),1);
EpiType4(isGeneral)                            = "General";
EpiType4(~isGeneral & isTemporal)              = "Temporal";
EpiType4(~isGeneral & ~isTemporal & isFrontal) = "Frontal";
PerPat.EpiType4 = EpiType4;

% --- Restrict S to the exact (Patient,Session) pairs retained in R --- 
% (after the outpatient filter has already been applied to R)
RSess = unique(R(:, {'patient_id','session_number'}));
RSess.Properties.VariableNames = {'Patient','Session'};
RSess.Patient = double(RSess.Patient);
RSess.Session = double(RSess.Session);

% If only_outpt==1, R already has only outpatient rows, so this keeps only outpatient EEGs.
% If only_outpt==0, this keeps all EEGs that have a matching report row.
S = innerjoin(S, RSess, 'Keys', {'Patient','Session'});

% If you still want to exclude patients without a stable epilepsy type:
S = S(ismember(S.Patient, PerPat.Patient), :);


%% ================== MERGE & PER-PATIENT MEAN SPIKE RATE ==================
S = S(ismember(S.Patient, PerPat.Patient), :);
[gp, pids] = findgroups(S.Patient);
MeanSpikeRate = splitapply(@(x) mean(x,'omitnan'), S.SpikeRate_perHour, gp); % takes mean spike rate across eegs within a patient
Sg = table(pids, MeanSpikeRate, 'VariableNames', {'Patient','MeanSpikeRate_perHour'});
Pg = innerjoin(Sg, PerPat, 'Keys','Patient');

etype_norm = lower(Pg.EpilepsyType);
isNESD     = etype_norm == lower(NESD_LABEL);
isBad      = ismember(etype_norm, badTypes);

Epilepsy_any = ~isNESD & ~isBad & strlength(Pg.EpilepsyType)>0;
NESD_only    = isNESD;
in3          = Epilepsy_any & ismember(Pg.EpiType4, canon3);

%% ================== PANEL A: Epilepsy vs NESD ==================
x_ep  = Pg.MeanSpikeRate_perHour(Epilepsy_any); % spikes/hour
% remaining types should be combined generalized and focal, focal, general,
% and unclassified or unspecified
x_nes = Pg.MeanSpikeRate_perHour(NESD_only);    % spikes/hour
n_ep  = nnz(isfinite(x_ep));
n_nes = nnz(isfinite(x_nes));

% Stats on raw scale (unchanged)
p_rank = NaN; 
if n_ep >= 3 && n_nes >= 3
    p_rank = ranksum(x_ep, x_nes, 'method','approx');
end
m_ep  = median(x_ep,'omitnan');   iqr_ep  = iqr(x_ep);
m_nes = median(x_nes,'omitnan');  iqr_nes = iqr(x_nes);

% Plot data: convert hour->min and log10 with global epsilon
Y_A = [to_log10_per_min(x_ep(:), EPS_PER_MIN); to_log10_per_min(x_nes(:), EPS_PER_MIN)];
G_A = [repmat("Epilepsy", n_ep, 1); repmat("NESD", n_nes, 1)];

%% ================== PANEL B: General vs Temporal vs Frontal ==================
x3 = Pg(in3, {'EpiType4','MeanSpikeRate_perHour'});
x3.EpiType4 = categorical(string(x3.EpiType4), canon3);

% Stats (raw)
[g3, cats3] = findgroups(x3.EpiType4);
medVals = splitapply(@(x) median(x,'omitnan'), x3.MeanSpikeRate_perHour, g3);
iqrVals = splitapply(@(x) prctile(x,75) - prctile(x,25), x3.MeanSpikeRate_perHour, g3);
nVals   = splitapply(@(x) sum(isfinite(x)), x3.MeanSpikeRate_perHour, g3);
stats3  = table(cats3, nVals, medVals, iqrVals, ...
                'VariableNames', {'EpiType4','GroupCount','Median','IQR'});

p_kw = NaN;
if height(x3) >= 6 && numel(categories(removecats(x3.EpiType4))) >= 2
    p_kw = kruskalwallis(x3.MeanSpikeRate_perHour, x3.EpiType4, 'off');
end

pairs = ["General","Temporal"; "General","Frontal"; "Temporal","Frontal"];
p_pw = NaN(3,1);
for i = 1:3
    A = pairs(i,1); B = pairs(i,2);
    xa = x3.MeanSpikeRate_perHour(x3.EpiType4==A);
    xb = x3.MeanSpikeRate_perHour(x3.EpiType4==B);
    if nnz(isfinite(xa)) >= 3 && nnz(isfinite(xb)) >= 3
        p_pw(i) = ranksum(xa, xb, 'method','approx');
    end
end
p_pw_bonf = min(p_pw * 3, 1);

% Plot data (log10 spikes/min with global epsilon)
Y_B = to_log10_per_min(x3.MeanSpikeRate_perHour, EPS_PER_MIN);

%% ================== PANEL C: Report Present vs Absent ==================
% Per-(Patient,Session) spikes/min for selected segment
[grpPS, pKeys, sKeys] = findgroups(S.Patient, S.Session); %#ok<ASGLU>
nanSum = @(x) sum(x(~isnan(x)));
Sum_Total = splitapply(nanSum, S.(colCount), grpPS);
Sum_Dur   = splitapply(nanSum, S.(colDur),   grpPS);
Sgrp = table(pKeys, sKeys, Sum_Total, Sum_Dur, ...
    'VariableNames', {'Patient','Session','Sum_Total','Sum_Dur'});

Sgrp.SpikesPerMin = nan(height(Sgrp),1);
validDur = Sgrp.Sum_Dur > 0;
Sgrp.SpikesPerMin(validDur) = (Sgrp.Sum_Total(validDur) ./ Sgrp.Sum_Dur(validDur)) * 60;

% Normalize Present/Absent
assert(ismember(reportColName, R.Properties.VariableNames), ...
    'Report column "%s" not found in the report CSV.', reportColName);
raw = lower(strtrim(string(R.(reportColName))));
isPresent = ismember(raw, ["present","yes","y","1","true","pos","positive"]);
isAbsent  = ismember(raw, ["absent","no","n","0","false","neg","negative","none","normal"]);
repNorm = strings(height(R),1);
repNorm(isPresent) = "present";
repNorm(isAbsent)  = "absent";
repNorm(~(isPresent | isAbsent)) = missing;
R.ReportStatus = categorical(repNorm, ["absent","present"]);

R2 = R(~ismissing(R.ReportStatus), {'patient_id','session_number','ReportStatus'});
R2.Properties.VariableNames = {'Patient','Session','ReportStatus'};
J = innerjoin(Sgrp(:,{'Patient','Session','SpikesPerMin'}), R2, 'Keys', {'Patient','Session'});

x_abs = J.SpikesPerMin(J.ReportStatus=="absent");
x_pre = J.SpikesPerMin(J.ReportStatus=="present");
p_rs  = NaN;
if nnz(isfinite(x_abs))>=3 && nnz(isfinite(x_pre))>=3
    p_rs = ranksum(x_abs, x_pre, 'method','approx');
end

% Plot data (already per-min)
Y_C = to_log10_from_per_min([x_abs(:); x_pre(:)], EPS_PER_MIN);
G_C = [repmat("Absent", numel(x_abs), 1); repmat("Present", numel(x_pre), 1)];

%% ================== MAKE ONE FIGURE WITH THREE PANELS (ORDER: C, A, B) ==================
f = figure('Color','w','Position',[60 60 1500 520]);
tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');

% -------- Panel C first: Report Present vs Absent --------
axC = nexttile; hold(axC,'on'); box(axC,'off'); grid(axC,'on');
if ~isempty(Y_C)
    boxchart(axC, categorical(G_C), Y_C, 'BoxFaceAlpha',0.25);
    try
        swarmchart(axC, categorical(G_C), Y_C, 18, 'filled','MarkerFaceAlpha',0.45);
    catch
        scatter(axC, double(categorical(G_C)) + 0.05*randn(size(Y_C)), Y_C, 14, 'filled','MarkerFaceAlpha',0.4);
    end
end
yline(axC, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axC, Y_LIMS);
ylabel(axC, 'log_{10}(spikes/min)');
title(axC, sprintf('A. Reported presence or absence of spikes'));
if isfinite(p_rs)
    % significance bar connecting the two groups at the top
    add_sigbar(axC, 1, 2, Y_LIMS(2) - 0.08*range(Y_LIMS), sprintf('p = %.3g', p_rs));
else
    text(axC, 0.5,0.97,'Insufficient N','Units','normalized','HorizontalAlignment','center','FontSize',20);
end
set(axC,'FontSize',20);

% -------- Panel A second: Epilepsy vs NESD --------
axA = nexttile; hold(axA,'on'); box(axA,'off'); grid(axA,'on');
if ~isempty(Y_A)
    boxchart(axA, categorical(G_A), Y_A, 'BoxFaceAlpha',0.25);
    try
        swarmchart(axA, categorical(G_A), Y_A, 18, 'filled','MarkerFaceAlpha',0.45);
    catch
        scatter(axA, double(categorical(G_A)) + 0.05*randn(size(Y_A)), Y_A, 14, 'filled','MarkerFaceAlpha',0.4);
    end
end
yline(axA, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axA, Y_LIMS);
ylabel(axA, 'log_{10}(spikes/min)');
title(axA, sprintf('B. Epilepsy versus PNEE'));
if isfinite(p_rank)
    add_sigbar(axA, 1, 2, Y_LIMS(2) - 0.08*range(Y_LIMS), sprintf('p = %.3g', p_rank));
else
    text(axA, 0.5,0.97,'Insufficient N','Units','normalized','HorizontalAlignment','center','FontSize',20);
end
set(axA,'FontSize',20);

% -------- Panel B third: General vs Temporal vs Frontal --------
axB = nexttile; hold(axB,'on'); box(axB,'off'); grid(axB,'on');
if ~isempty(Y_B)
    boxchart(axB, x3.EpiType4, Y_B, 'BoxFaceAlpha',0.25);
    try
        swarmchart(axB, x3.EpiType4, Y_B, 18, 'filled','MarkerFaceAlpha',0.45);
    catch
        jitterX = double(x3.EpiType4) + 0.05*randn(height(x3),1);
        scatter(axB, jitterX, Y_B, 14, 'filled','MarkerFaceAlpha',0.4);
    end
end
yline(axB, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axB, Y_LIMS);
ylabel(axB, 'log_{10}(spikes/min)');
title(axB, sprintf('C: Epilepsy subtype'));
set(axB,'FontSize',20);

% -------- KW text --------
%{
if isfinite(p_kw)
    text(axB, 0.5,0.97,sprintf('KW p=%.3g', p_kw), ...
        'Units','normalized','HorizontalAlignment','center','FontSize',20);
else
    text(axB, 0.5,0.97,'Insufficient N', ...
        'Units','normalized','HorizontalAlignment','center','FontSize',20);
end
%}

% -------- Pairwise bars for ALL 3 comparisons (Bonferroni-adjusted) --------
% pairs = ["General","Temporal"; "General","Frontal"; "Temporal","Frontal"];
% p_pw_bonf is 3x1

% y positions stacked downward from the top
yTop  = Y_LIMS(2);
yStep = 0.07 * range(Y_LIMS);
y0    = yTop - 0.05 * range(Y_LIMS);

cats = categorical(canon3);  % positions 1,2,3 -> General, Temporal, Frontal

for i = 1:3
    A = pairs(i,1); B = pairs(i,2);
    x1 = find(cats == categorical(A));
    x2 = find(cats == categorical(B));

    pval = p_pw_bonf(i);
    % scalar label selection (strings, no arithmetic with strings)
    if ~isfinite(pval)
        lab = "n/a";
    elseif pval < 1e-3
        lab = "***";
    elseif pval < 1e-2
        lab = "**";
    elseif pval < 5e-2
        lab = "*";
    else
        lab = "ns";
    end

    add_sigbar(axB, x1, x2, y0 - (i-1)*yStep, lab);
end

% Save
if ~exist(fileparts(figOut),'dir'), mkdir(fileparts(figOut)); end
exportgraphics(f, figOut, 'Resolution', 300);
fprintf('Saved combined control figure: %s\n', figOut);

%% ================== Console summary (raw-scale stats) ==================
fprintf('\n=== Panel A: Epilepsy vs NESD — %s ===\n', segLabel);
fprintf('Epilepsy_any:  n=%d  median=%.3g  IQR=%.3g (spikes/hour)\n', n_ep,  m_ep,  iqr_ep);
fprintf('NESD:          n=%d  median=%.3g  IQR=%.3g (spikes/hour)\n', n_nes, m_nes, iqr_nes);
if isfinite(p_rank), fprintf('rank-sum p=%.3g\n', p_rank); else, fprintf('Insufficient N\n'); end

fprintf('\n=== Panel B: General vs Temporal vs Frontal — %s ===\n', segLabel);
disp(stats3(:,{'EpiType4','Median','IQR','GroupCount'})); % spikes/hour
if isfinite(p_kw), fprintf('Kruskal–Wallis p=%.3g\n', p_kw); else, fprintf('Insufficient N for KW\n'); end
for i = 1:3
    fprintf('Pair %s vs %s: p=%.3g  (Bonferroni: %.3g)\n', pairs(i,1), pairs(i,2), p_pw(i), p_pw_bonf(i));
end

fprintf('\n=== Panel C: %s — Present vs Absent — %s ===\n', reportColName, segLabel);
fprintf('N absent=%d, N present=%d (per-session)\n', numel(x_abs), numel(x_pre));
if isfinite(p_rs), fprintf('rank-sum p=%.3g\n', p_rs); else, fprintf('Insufficient N\n'); end

end % function

%% ====================== HELPERS ======================
function ylog = to_log10_per_min(x_per_hour, eps_per_min)
% Convert spikes/hour -> spikes/min, then log10 with global epsilon
pm = double(x_per_hour) / 60;
pm(~isfinite(pm) | pm <= 0) = eps_per_min;
ylog = log10(pm);
end

function ylog = to_log10_from_per_min(x_per_min, eps_per_min)
% spikes/min -> log10 with global epsilon
pm = double(x_per_min);
pm(~isfinite(pm) | pm <= 0) = eps_per_min;
ylog = log10(pm);
end

function add_sigbar(ax, x1, x2, y, ptext)
% Draw a horizontal significance bar between x1 and x2 at height y,
% with small vertical ticks and the p-value centered above.
tick = 0.03 * diff(ax.YLim);
plot(ax, [x1 x1 x2 x2], [y- tick, y, y, y- tick], 'k-', 'LineWidth', 1.3);
text(ax, mean([x1 x2]), y + 0.02*diff(ax.YLim), ptext, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);
end
