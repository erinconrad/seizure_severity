function compare_spikerates_groups(which_runs)
% compare_spikerates_groups(which_runs)
%   0 = whole file (default)
%   1 = first run (~1h)
%   2 = first 24 runs (~24h)

if nargin < 1, which_runs = 0; end

%% ======================= CONFIG =======================
outCsv      = '../data/SN_counts/spike_counts_summary.csv';
reportFile  = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-10_1443.csv';

figOut1     = '../figures/spikerate_epilepsy_vs_nesd.png';
figOut2     = '../figures/spikerate_general_temporal_frontal.png';
resultsCsv  = '../output/spikerate_group_comparisons.csv';

only_amb    = 0;  % 1 = ONLY ambulatory (>=12h); 2 = ONLY routine (<=12h); 0 = all
only_outpt  = 0;  % 1 = keep EEGs with acquired_on containing 'Spe' or 'radnor'

badTypes = lower(["Uncertain if Epilepsy","Unknown or MRN not found",""]); % excluded from "Epilepsy"
NESD_LABEL = "Non-Epileptic Seizure Disorder";
canon3 = ["General","Temporal","Frontal"];

%% =============== LOAD SPIKE SUMMARY ===============
S = readtable(outCsv,'TextType','string','VariableNamingRule','preserve');
switch which_runs
    case 1, colCount="FirstRun_Spikes";    colDur="FirstRun_Duration_sec";     segLabel='Whole file'; segLabel='First run (~1h)';
    case 2, colCount="First24Runs_Spikes"; colDur="First24Runs_Duration_sec";  segLabel='First 24 runs (~24h)';
    otherwise, colCount="Total_Spikes";    colDur="Duration_sec";              segLabel='Whole file';
end
segLabel = string(segLabel);   % ensure string scalar for all downstream uses

needS = {'Patient','Session','Duration_sec', char(colCount), char(colDur)};
assert(all(ismember(needS, S.Properties.VariableNames)), 'Missing required columns in spike CSV.');
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end

if only_amb == 1
    S(S.Duration_sec < 12*3600,:) = [];
elseif only_amb == 2
    S(S.Duration_sec > 12*3600,:) = [];
end

if ~isnumeric(S.(colCount)), S.(colCount) = double(S.(colCount)); end
if ~isnumeric(S.(colDur)),   S.(colDur)   = double(S.(colDur));   end
S.SpikeRate_perHour = nan(height(S),1);
okDur = S.(colDur) > 0;
S.SpikeRate_perHour(okDur) = (S.(colCount)(okDur) ./ S.(colDur)(okDur)) * 3600;

%% =============== LOAD REPORT ===============
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

% Per-patient epilepsy_type (first non-missing, stable)
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

%% =============== MERGE & PER-PATIENT MEAN SPIKE RATE ===============
S = S(ismember(S.Patient, PerPat.Patient), :);
[gp, pids] = findgroups(S.Patient);
MeanSpikeRate = splitapply(@(x) mean(x,'omitnan'), S.SpikeRate_perHour, gp);
Sg = table(pids, MeanSpikeRate, 'VariableNames', {'Patient','MeanSpikeRate_perHour'});
Pg = innerjoin(Sg, PerPat, 'Keys','Patient');

etype_norm = lower(Pg.EpilepsyType);
isNESD     = etype_norm == lower(NESD_LABEL);
isBad      = ismember(etype_norm, badTypes);

Epilepsy_any = ~isNESD & ~isBad & strlength(Pg.EpilepsyType)>0;
NESD_only    = isNESD;
in3 = Epilepsy_any & ismember(Pg.EpiType4, canon3);

%% =============== (1) EPILEPSY vs NESD ===============
x_ep  = Pg.MeanSpikeRate_perHour(Epilepsy_any);
x_nes = Pg.MeanSpikeRate_perHour(NESD_only);
n_ep  = nnz(~isnan(x_ep));
n_nes = nnz(~isnan(x_nes));

p_rank = NaN; stattxt = '';
if n_ep >= 3 && n_nes >= 3
    p_rank = ranksum(x_ep, x_nes, 'method','approx');
    stattxt = sprintf('rank-sum p=%.3g', p_rank);
else
    stattxt = 'Insufficient N for test';
end

m_ep  = median(x_ep,'omitnan');   iqr_ep  = iqr(x_ep);
m_nes = median(x_nes,'omitnan');  iqr_nes = iqr(x_nes);

fprintf('\n=== Spike rate: Epilepsy (any) vs NESD — %s ===\n', segLabel);
fprintf('Epilepsy_any:  n=%d  median=%.3g  IQR=%.3g\n', n_ep,  m_ep,  iqr_ep);
fprintf('NESD:          n=%d  median=%.3g  IQR=%.3g\n', n_nes, m_nes, iqr_nes);
fprintf('%s\n', stattxt);

f1 = figure('Color','w','Position',[70 70 720 540]); hold on; box off; grid on;
G1 = [repmat("Epilepsy", n_ep, 1); repmat("NESD", n_nes, 1)];
Y1 = [x_ep(:); x_nes(:)];
boxchart(categorical(G1), Y1, 'BoxFaceAlpha',0.25);
try
    swarmchart(categorical(G1), Y1, 20, 'filled','MarkerFaceAlpha',0.5);
catch
    scatter(double(categorical(G1)) + 0.05*randn(size(Y1)), Y1, 18, 'filled','MarkerFaceAlpha',0.4);
end
ylabel('Spike rate (spikes/hour)');
title(sprintf('Epilepsy vs NESD — %s', segLabel));
text(0.5, 0.97, stattxt, 'Units','normalized','HorizontalAlignment','center','FontSize',12);
set(gca,'FontSize',14);
if ~exist(fileparts(figOut1),'dir'), mkdir(fileparts(figOut1)); end
exportgraphics(f1, figOut1, 'Resolution', 300);
fprintf('Saved figure: %s\n', figOut1);

%% =============== (2) GENERAL vs TEMPORAL vs FRONTAL (within epilepsy) ===============
x3 = Pg(in3, {'EpiType4','MeanSpikeRate_perHour'});
x3.EpiType4 = categorical(string(x3.EpiType4), canon3);

% Robust IQR function (works without Statistics Toolbox name in groupsummary)
iqrFun = @(x) prctile(x,75) - prctile(x,25);

% --- Robust per-group stats (no groupsummary naming surprises) ---
x3 = Pg(in3, {'EpiType4','MeanSpikeRate_perHour'});
x3.EpiType4 = categorical(string(x3.EpiType4), canon3);

[g3, cats3] = findgroups(x3.EpiType4);
medVals = splitapply(@(x) median(x,'omitnan'), x3.MeanSpikeRate_perHour, g3);
iqrVals = splitapply(@(x) prctile(x,75) - prctile(x,25), x3.MeanSpikeRate_perHour, g3);
nVals   = splitapply(@(x) sum(isfinite(x)), x3.MeanSpikeRate_perHour, g3);

stats3 = table(cats3, nVals, medVals, iqrVals, ...
               'VariableNames', {'EpiType4','GroupCount','Median','IQR'});

% Kruskal–Wallis (overall)
p_kw = NaN;
if height(x3) >= 6 && numel(categories(removecats(x3.EpiType4))) >= 2
    p_kw = kruskalwallis(x3.MeanSpikeRate_perHour, x3.EpiType4, 'off');
end

% Pairwise ranksum with Bonferroni
pairs = ["General","Temporal";
         "General","Frontal";
         "Temporal","Frontal"];
p_pw = NaN(3,1); nA = NaN(3,1); nB = NaN(3,1);
for i = 1:3
    A = pairs(i,1); B = pairs(i,2);
    xa = x3.MeanSpikeRate_perHour(x3.EpiType4==A);
    xb = x3.MeanSpikeRate_perHour(x3.EpiType4==B);
    nA(i) = nnz(isfinite(xa)); nB(i) = nnz(isfinite(xb));
    if nA(i) >= 3 && nB(i) >= 3
        p_pw(i) = ranksum(xa, xb, 'method','approx');
    end
end
p_pw_bonf = min(p_pw * 3, 1);

% Print table cleanly
fprintf('=== Spike rate: General vs Temporal vs Frontal — %s ===\n', segLabel);
disp(stats3(:,{'EpiType4','Median','IQR','GroupCount'}));
if ~isnan(p_kw)
    fprintf('Kruskal–Wallis p=%.3g\n', p_kw);
else
    fprintf('Kruskal–Wallis: insufficient data\n');
end
for i = 1:3
    fprintf('Pair %s vs %s: n=(%d,%d)  p=%.3g  p_bonf=%.3g\n', ...
        pairs(i,1), pairs(i,2), nA(i), nB(i), p_pw(i), p_pw_bonf(i));
end

if ~isnan(p_kw)
    fprintf('Kruskal–Wallis p=%.3g\n', p_kw);
else
    fprintf('Kruskal–Wallis: insufficient data\n');
end
for i = 1:3
    fprintf('Pair %s vs %s: n=(%d,%d)  p=%.3g  p_bonf=%.3g\n', ...
        pairs(i,1), pairs(i,2), nA(i), nB(i), p_pw(i), p_pw_bonf(i));
end

% Figure (box + swarm)
f2 = figure('Color','w','Position',[90 80 820 560]); hold on; box off; grid on;
boxchart(x3.EpiType4, x3.MeanSpikeRate_perHour, 'BoxFaceAlpha',0.25);
try
    swarmchart(x3.EpiType4, x3.MeanSpikeRate_perHour, 20, 'filled','MarkerFaceAlpha',0.5);
catch
    jitterX = double(x3.EpiType4) + 0.05*randn(height(x3),1);
    scatter(jitterX, x3.MeanSpikeRate_perHour, 18, 'filled','MarkerFaceAlpha',0.4);
end
ylabel('Spike rate (spikes/hour)');
title(sprintf('General vs Temporal vs Frontal — %s', segLabel));
if ~isnan(p_kw)
    text(0.5, 0.97, sprintf('Kruskal–Wallis p=%.3g', p_kw), ...
        'Units','normalized','HorizontalAlignment','center','FontSize',12);
end
set(gca,'FontSize',14);
if ~exist(fileparts(figOut2),'dir'), mkdir(fileparts(figOut2)); end
exportgraphics(f2, figOut2, 'Resolution', 300);
fprintf('Saved figure: %s\n', figOut2);

%% =============== SAVE SUMMARY CSV ===============
if ~exist(fileparts(resultsCsv),'dir'), mkdir(fileparts(resultsCsv)); end

% Safe extractor for stats3 table
getStat = @(grp,field) extractStat(stats3, grp, field);

% Row 1: Epilepsy_any vs NESD
row1 = cell2table({ ...
    "Epilepsy_any_vs_NESD", segLabel, ...
    double(n_ep), double(m_ep), double(iqr_ep), ...
    double(n_nes), double(m_nes), double(iqr_nes), ...
    double(p_rank), NaN }, ...
    'VariableNames', {'Comparison','Segment','n_A','median_A','IQR_A','n_B','median_B','IQR_B','p_value','p_bonf'});

% Row 2: KW 3-group summary (General vs Temporal vs Frontal)
row2 = cell2table({ ...
    "KW_General_Temporal_Frontal", segLabel, ...
    double(getStat("General","GroupCount")), double(getStat("General","Median")), double(getStat("General","IQR")), ...
    double(getStat("Temporal","GroupCount")), double(getStat("Temporal","Median")), double(getStat("Temporal","IQR")), ...
    double(p_kw), NaN }, ...
    'VariableNames', {'Comparison','Segment','n_A','median_A','IQR_A','n_B','median_B','IQR_B','p_value','p_bonf'});

% Pairwise rows with Bonferroni
pairs = ["General","Temporal"; "General","Frontal"; "Temporal","Frontal"];
rowsPW = cell(0,10);
for i = 1:3
    A = pairs(i,1); B = pairs(i,2);
    rowsPW(end+1,:) = { ...
        "Pair_" + A + "_vs_" + B, segLabel, ...
        double(nA(i)), NaN, NaN, ...
        double(nB(i)), NaN, NaN, ...
        double(p_pw(i)), double(p_pw_bonf(i)) }; %#ok<AGROW>
end
rowsPW = cell2table(rowsPW, ...
    'VariableNames', row1.Properties.VariableNames);

writetable([row1; row2; rowsPW], resultsCsv);
fprintf('Saved results CSV: %s\n', resultsCsv);


end % function

function val = extractStat(T, grp, field)
% Safe extractor for stats3 table: returns NaN if group/field missing
val = NaN;
if isempty(T), return; end
if ~ismember(field, T.Properties.VariableNames), return; end
ix = strcmp(string(T.EpiType4), string(grp));
if any(ix)
    v = T.(field)(ix);
    if ~isempty(v), val = v(1); end
end
end
