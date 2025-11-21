%% ======================= CONFIG =======================
% Paths
spikeSummaryCsv   = '../data/SN_counts/spike_counts_summary.csv';
reportCsv         = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-19_1356.csv';

show_intermodality_comparisons = 0; % show ambulatory vs routine

% Which segment of the EEG to use for spike rates
%   0 = whole file, 1 = first run (~1h), 2 = first 24 runs (~24h)
which_runs  = 0;

% Cohort filters
only_amb    = 2;   % 1 = ONLY ambulatory (>=12h); 2 = ONLY routine (<=12h); 0 = all
only_outpt  = 1;   % 1 = keep only outpatient EEGs (as defined above)

% Logic for has sz
hassz_replace = 1; % if 1, then make visits where has Sz == 0 be 0 sz frequency  - this definitely improves stats

% Label constants
NESD_LABEL = "Non-Epileptic Seizure Disorder";
badTypes   = lower(["Uncertain if Epilepsy","Unknown or MRN not found",""]); % excluded from 'Epilepsy'
canonical3 = ["General","Temporal","Frontal"];   % canonical 3 subtypes of epilepsy

% ERIN MOVE THIS DOWN BELOW
% Plotting params (shared)
EPS_PER_MIN = 1e-3;                 % → dotted "zero" line at y = -3 in log10(spikes/min)
Y_ZERO      = log10(EPS_PER_MIN);
Y_LIMS      = [-3.2 2];             % fixed y-lims for box/swarm plots

% ERIN MOVE THIS DOWN BELOW
% Spearman-figure axes (fixed across panels)
spearman_xLims = [-3.5, 4];         % log10(seizures/month)
spearman_yLims = [-3, 2];           % log10(spikes/min)

% Outputs
fig1_out   = '../figures/spikerate_control_panels.png';
fig2_out   = '../figures/spearman_spikerate_szfreq_log10.png';
resultsCsv = '../output/spearman_spikerate_szfreq_log10.csv';

%% ======================= CORE LOAD: REPORT  =======================
% ---------- Load raw report (for outpatient rule + everything else) ----------
ReportTable = readtable(reportCsv,'TextType','string','VariableNamingRule','preserve');
if ~isnumeric(ReportTable.patient_id)
    ReportTable.patient_id = double(str2double(string(ReportTable.patient_id)));
end
if ~isnumeric(ReportTable.session_number)
    ReportTable.session_number = double(str2double(string(ReportTable.session_number)));
end
if ~isstring(ReportTable.epilepsy_type)
    ReportTable.epilepsy_type = string(ReportTable.epilepsy_type);
end

% Required outpatient classification columns
reqCols = {'report_PATIENT_CLASS','jay_in_or_out'};
missingReq = setdiff(reqCols, ReportTable.Properties.VariableNames);
if ~isempty(missingReq)
    error('Report table lacks required outpatient columns: %s', strjoin(missingReq, ', '));
end


hasAcqCol = ismember('acquired_on', ReportTable.Properties.VariableNames);
if ~hasAcqCol
    error('Report lacks "acquired_on"; Spe/Radnor outpatient-by-site rule will be unavailable.');
end

acqStr_report = lower(strtrim(string(ReportTable.acquired_on)));


%% ======================= CORE LOAD: SPIKE SUMMARY =======================
SpikeSummaryTable = readtable(spikeSummaryCsv,'TextType','string','VariableNamingRule','preserve');

switch which_runs
    case 1
        countCol = "FirstRun_Spikes";
        durCol   = "FirstRun_Duration_sec";
        segLabel = 'First run (~1h)';
    case 2
        countCol = "First24Runs_Spikes";
        durCol   = "First24Runs_Duration_sec";
        segLabel = 'First 24 runs (~24h)';
    otherwise
        countCol = "Total_Spikes";
        durCol   = "Duration_sec";
        segLabel = 'Whole file';
end
segLabel = string(segLabel);

SpikeSummaryTable = ensure_spikerates(SpikeSummaryTable, countCol, durCol);   % adds SpikeRate_perHour/SpikeRate_perMin

% Standardize key types for joining
if ~isnumeric(SpikeSummaryTable.Patient)
    SpikeSummaryTable.Patient = double(str2double(string(SpikeSummaryTable.Patient)));
end
if ~isnumeric(SpikeSummaryTable.Session)
    SpikeSummaryTable.Session = double(str2double(string(SpikeSummaryTable.Session)));
end

%% ======================= OUTPATIENT FILTER =======================
if only_outpt == 1
    fprintf('\n=== BUILDING OUTPATIENT KEYS (site rule OR report_patient_class OR jay_in_or_out) ===\n');

    % ---------- 1) Outpatient by acquired_on (Spe/Radnor) ----------
    OutptBySite = table('Size',[0 2], ...
                         'VariableTypes',{'double','double'}, ...
                         'VariableNames',{'Patient','Session'});

    isOutpt_site = (contains(acqStr_report,"spe",'IgnoreCase',true) | contains(acqStr_report,"radnor",'IgnoreCase',true));

    if any(isOutpt_site)
        OutptBySite = unique(ReportTable(isOutpt_site, {'patient_id','session_number'}));
        OutptBySite.Properties.VariableNames = {'Patient','Session'};

        % Ensure numeric
        if ~isnumeric(OutptBySite.Patient)
            OutptBySite.Patient = double(str2double(string(OutptBySite.Patient)));
        end
        if ~isnumeric(OutptBySite.Session)
            OutptBySite.Session = double(str2double(string(OutptBySite.Session)));
        end
    end

    fprintf('[Outpatient rule] Site-based (Spe/Radnor) outpatients: %d sessions.\n', ...
        height(OutptBySite));
    

    % ---------- 2) Outpatient by report_patient_class == "Outpatient" ----------
    OutptByClass = table('Size',[0 2], ...
                         'VariableTypes',{'double','double'}, ...
                         'VariableNames',{'Patient','Session'});
    classStr = lower(strtrim(string(ReportTable.report_PATIENT_CLASS)));
    isOutpt_class = (classStr == "outpatient");

    if any(isOutpt_class)
        OutptByClass = unique(ReportTable(isOutpt_class, {'patient_id','session_number'}));
        OutptByClass.Properties.VariableNames = {'Patient','Session'};

        if ~isnumeric(OutptByClass.Patient)
            OutptByClass.Patient = double(str2double(string(OutptByClass.Patient)));
        end
        if ~isnumeric(OutptByClass.Session)
            OutptByClass.Session = double(str2double(string(OutptByClass.Session)));
        end
    end

    fprintf('[Outpatient rule] report_PATIENT_CLASS=="Outpatient": %d sessions.\n', ...
        height(OutptByClass));

    % ---------- 3) Outpatient by jay_in_or_out == "out" ----------
    OutptByJay = table('Size',[0 2], ...
                       'VariableTypes',{'double','double'}, ...
                       'VariableNames',{'Patient','Session'});
    jayStr = lower(strtrim(string(ReportTable.jay_in_or_out)));
    isOutpt_jay = (jayStr == "out");

    if any(isOutpt_jay)
        OutptByJay = unique(ReportTable(isOutpt_jay, {'patient_id','session_number'}));
        OutptByJay.Properties.VariableNames = {'Patient','Session'};

        if ~isnumeric(OutptByJay.Patient)
            OutptByJay.Patient = double(str2double(string(OutptByJay.Patient)));
        end
        if ~isnumeric(OutptByJay.Session)
            OutptByJay.Session = double(str2double(string(OutptByJay.Session)));
        end
    end

    fprintf('[Outpatient rule] jay_in_or_out=="out": %d sessions.\n', ...
        height(OutptByJay));

    % ---------- 4) Finalize outpatient key set: union of all three ----------
    OutptKeys = [OutptBySite; OutptByClass; OutptByJay];
    if ~isempty(OutptKeys)
        OutptKeys = unique(OutptKeys);
    else
        error(['No outpatient sessions identified from any rule:\n' ...
               '  Spe/Radnor site, report_PATIENT_CLASS=="Outpatient", or jay_in_or_out=="out".']);
    end

    % ---------- 5) Apply filter to spike summary + report ----------
    SpikeSummaryTable = innerjoin(SpikeSummaryTable, OutptKeys, 'Keys', {'Patient','Session'});
    ReportTable = innerjoin(ReportTable, OutptKeys, ...
        'LeftKeys', {'patient_id','session_number'}, ...
        'RightKeys', {'Patient','Session'});

    fprintf(['\n[Outpatient filter] Outpatients defined as ANY of:\n' ...
         '  - acquired_on containing "Spe" or "Radnor" (case-insensitive) OR\n' ...
         '  - report_PATIENT_CLASS == "Outpatient" OR\n' ...
         '  - jay_in_or_out == "out"\n']);
    fprintf('[Outpatient filter] Kept %d spike rows and %d report rows.\n', ...
        height(SpikeSummaryTable), height(ReportTable));
end  % only_outpt == 1


% Keep a FULL (already outpatient-restricted if chosen) copy for modality comparisons
SpikeSummaryTable_FULL = SpikeSummaryTable;


%% ======================= BUILD PATIENT-LEVEL METRICS ONCE =======================
[PatientTypingAll, SzFreqPerPatient] = build_patient_metrics_from_report(ReportTable, canonical3,hassz_replace); % looks good!

%% ======================= RESTRICT TO STUDY TYPE =======================
% View 1: FULL (no only_amb restriction) — for ambulatory vs routine comparisons
ViewsFull.Sessions      = SpikeSummaryTable_FULL;
ViewsFull.Report        = ReportTable;
ViewsFull.SzFreqPerPat  = SzFreqPerPatient;

% View 2: FILTERED (apply only_amb) — for Fig 1 + Fig 2 
Views = build_filtered_view(SpikeSummaryTable, ReportTable, PatientTypingAll, SzFreqPerPatient, ...
                            only_amb, countCol, durCol, NESD_LABEL, badTypes, canonical3);


%% ======================= QUICK MODALITY COMPARISONS  =======================
if show_intermodality_comparisons
    compare_spikerate_by_modality(ViewsFull.Sessions, EPS_PER_MIN);
    compare_szfreq_by_modality(ViewsFull.Sessions, ViewsFull.SzFreqPerPat);
end


%% ======================= FIGURE 1: CONTROL PANELS =======================
SessionLevelSpikeRates   = Views.SessionLevelSpikeRates;      % per-session spikes/min
ReportForKeptSessions    = Views.ReportForKeptSessions;       % report rows matching sessions used
PatientLevelSpikeRates   = Views.PatientLevelSpikeRates;      % per-patient mean spike rates
SubtypePairs             = Views.Canonical3_Pairs;            % pair labels for subtype comparisons
p_pair_bonf              = Views.PvalsPairwiseBonf;
p_pair_raw               = Views.PvalsPairwise;

isEpilepsyMask           = Views.IsEpilepsyMask;
isNESDMask               = Views.IsNESDMask;
SubtypeSubsetTable       = Views.Canonical3_SubsetTable;
SubtypeStatsTable        = Views.Canonical3_Stats;

% ---- Panel A: Report Present vs Absent (session-level, FILTERED cohort) ----

% Normalize all four columns to "present"/"absent"/(blank)
tokPresent = ["present","yes","y","1","true","pos","positive"];
tokAbsent  = ["absent","no","n","0","false","neg","negative","none","normal"];

% Main "any spikes" column
any_spikes_epic_report = ReportForKeptSessions.report_SPORADIC_EPILEPTIFORM_DISCHARGES;
isMainPresent = any_spikes_epic_report == "present";
isMainAbsent  = any_spikes_epic_report == "absent";

% jay_* columns
rawF = lower(strtrim(string(ReportForKeptSessions.jay_focal_epi)));
rawM = lower(strtrim(string(ReportForKeptSessions.jay_multifocal_epi)));
rawG = lower(strtrim(string(ReportForKeptSessions.jay_gen_epi)));

isF_P = rawF=="present";
isF_A = rawF=="absent";
isM_P = rawM == "present";
isM_A = rawM == "absent";
isG_P = rawG == "present";
isG_A = rawG == "absent";

% Convenience masks
presentJay_any  = isF_P | isM_P | isG_P;           % any jay present
allJay_absent   = isF_A & isM_A & isG_A;          % all three explicitly "absent"
allJay_present  = isF_P & isM_P & isG_P;          % all three explicitly "present"

blankMain   = ~(isMainPresent | isMainAbsent);    % nothing interpretable in main col
blankJayAll = ~(isF_P | isF_A) & ~(isM_P | isM_A) & ~(isG_P | isG_A);  % all jay blank

% ---------- Discordance checks ----------
% 1) All jay columns say ABSENT but main column says PRESENT -> error
disc1 = allJay_absent & isMainPresent;

% 2) Main column says ABSENT and all jay columns say PRESENT -> error
disc2 = isMainAbsent & allJay_present;

discordMask = disc1 | disc2;
if any(discordMask)
    error('Discordant spike presence between main and jay_* columns');
end

% ---------- Final combined spike-present status ----------
repCombined = strings(height(ReportForKeptSessions),1);

% Rule 1: if ANY of [main, jay_*] say PRESENT -> PRESENT
presentAny = isMainPresent | presentJay_any;
repCombined(presentAny) = "present";

% Rule 2: if all jay_* say ABSENT and main is BLANK -> ABSENT
mask_absent_case1 = allJay_absent & blankMain;
repCombined(mask_absent_case1) = "absent";

% Rule 3: if main says ABSENT and all jay_* are BLANK -> ABSENT
mask_absent_case2 = isMainAbsent & blankJayAll;
repCombined(mask_absent_case2) = "absent";

% Anything else stays missing (repCombined == "")
repCombined(repCombined=="") = missing;

ReportForKeptSessions.ReportStatus = categorical(repCombined, ["absent","present"]);

% Keep only rows with a resolved "present"/"absent" status
ReportSlim = ReportForKeptSessions(~ismissing(ReportForKeptSessions.ReportStatus), ...
                                   {'Patient','Session','ReportStatus'});

JoinA = innerjoin(SessionLevelSpikeRates(:,{'Patient','Session','SpikesPerMin'}), ...
                  ReportSlim, 'Keys', {'Patient','Session'});

x_abs = JoinA.SpikesPerMin(JoinA.ReportStatus=="absent");
x_pre = JoinA.SpikesPerMin(JoinA.ReportStatus=="present");
p_rankSum_A  = NaN;
if nnz(isfinite(x_abs))>=3 && nnz(isfinite(x_pre))>=3
    p_rankSum_A = ranksum(x_abs, x_pre, 'method','approx');
end
Y_A = to_log10_from_per_min([x_abs(:); x_pre(:)], EPS_PER_MIN);
G_A = [repmat("Absent", numel(x_abs), 1); repmat("Present", numel(x_pre), 1)];




% ---- Panel B calculations: Epilepsy (any) vs NESD (patient-level means) ----
x_ep  = PatientLevelSpikeRates.MeanSpikeRate_perHour(isEpilepsyMask); % combined gen + focal, focal, gen, unclassified (leaves out uncertain, unknown, PNEE)
x_nes = PatientLevelSpikeRates.MeanSpikeRate_perHour(isNESDMask); % just pnee
n_ep  = nnz(isfinite(x_ep)); % how many finite numerical values
n_nes = nnz(isfinite(x_nes));
p_rankSum_B = NaN;
if n_ep >= 3 && n_nes >= 3
    p_rankSum_B = ranksum(x_ep, x_nes, 'method','approx'); % rank sum test
end
m_ep  = median(x_ep,'omitnan');   iqr_ep  = iqr(x_ep); % median and iqr
m_nes = median(x_nes,'omitnan');  iqr_nes = iqr(x_nes);
Y_B = [to_log10_per_min(x_ep(:), EPS_PER_MIN); to_log10_per_min(x_nes(:), EPS_PER_MIN)];
G_B = [repmat("Epilepsy", n_ep, 1); repmat("NESD", n_nes, 1)];

% ---- Panel C calculations: General vs Temporal vs Frontal (patient-level) ----
Y_C = to_log10_per_min(SubtypeSubsetTable.MeanSpikeRate_perHour, EPS_PER_MIN);



% ---- Draw figure ----
f1 = figure('Color','w','Position',[60 60 1500 520]);
tiledlayout(f1,1,3,'TileSpacing','compact','Padding','compact');

% A) Report Present vs Absent
axA = nexttile; hold(axA,'on'); box(axA,'off'); grid(axA,'on');
boxchart(axA, categorical(G_A), Y_A, 'BoxFaceAlpha',0.25);
swarmchart(axA, categorical(G_A), Y_A, 18, 'filled','MarkerFaceAlpha',0.25);    
yline(axA, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axA, Y_LIMS);
ylabel(axA, 'log_{10}(spikes/min)');
title(axA, 'A. Reported presence or absence of spikes');
add_sigbar(axA, 1, 2, Y_LIMS(2) - 0.08*range(Y_LIMS), sprintf('p = %.3g', p_rankSum_A));

set(axA,'FontSize',20);

% B) Epilepsy vs NESD
axB = nexttile; hold(axB,'on'); box(axB,'off'); grid(axB,'on');
boxchart(axB, categorical(G_B), Y_B, 'BoxFaceAlpha',0.25);
swarmchart(axB, categorical(G_B), Y_B, 18, 'filled','MarkerFaceAlpha',0.25);    
yline(axB, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axB, Y_LIMS);
ylabel(axB, 'log_{10}(spikes/min)');
title(axB, 'B. Epilepsy versus PNEE');
add_sigbar(axB, 1, 2, Y_LIMS(2) - 0.08*range(Y_LIMS), sprintf('p = %.3g', p_rankSum_B));
set(axB,'FontSize',20);

% C) General vs Temporal vs Frontal
axC = nexttile; hold(axC,'on'); box(axC,'off'); grid(axC,'on');
boxchart(axC, SubtypeSubsetTable.EpiType4, Y_C, 'BoxFaceAlpha',0.25);
swarmchart(axC, SubtypeSubsetTable.EpiType4, Y_C, 18, 'filled','MarkerFaceAlpha',0.25);   
yline(axC, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axC, Y_LIMS);
ylabel(axC, 'log_{10}(spikes/min)');
title(axC, 'C. Epilepsy subtype');
set(axC,'FontSize',20);

% Pairwise bars with Bonferroni-coded labels
yTop  = Y_LIMS(2);
yStep = 0.07 * range(Y_LIMS);
y0    = yTop - 0.05 * range(Y_LIMS);
cats3 = categorical(canonical3);
for i = 1:3
    A = SubtypePairs(i,1); B = SubtypePairs(i,2);
    x1 = find(cats3 == categorical(A));
    x2 = find(cats3 == categorical(B));
    pval = p_pair_bonf(i);
    if pval < 1e-3, lab = "***";
    elseif pval < 1e-2, lab = "**";
    elseif pval < 5e-2, lab = "*";
    else, lab = "ns";
    end
    add_sigbar(axC, x1, x2, y0 - (i-1)*yStep, lab);
end

if ~exist(fileparts(fig1_out),'dir'), mkdir(fileparts(fig1_out)); end
exportgraphics(f1, fig1_out, 'Resolution', 300);
fprintf('Saved Fig 1 (controls): %s\n', fig1_out);

% Console summary
fprintf('\n=== Panel A: Spikes — Present vs Absent — %s ===\n', segLabel);
fprintf('N absent=%d, N present=%d (per-session)\n', numel(x_abs), numel(x_pre));
fprintf('rank-sum p=%.3g\n', p_rankSum_A);
fprintf('\n=== Panel B: Epilepsy vs NESD — %s ===\n', segLabel);
fprintf('Epilepsy_any:  n=%d  median=%.3g  IQR=%.3g (spikes/hour)\n', n_ep,  m_ep,  iqr_ep);
fprintf('NESD:          n=%d  median=%.3g  IQR=%.3g (spikes/hour)\n', n_nes, m_nes, iqr_nes);
fprintf('rank-sum p=%.3g\n', p_rankSum_B);
fprintf('\n=== Panel C: General vs Temporal vs Frontal — %s ===\n', segLabel);
disp(SubtypeStatsTable(:,{'EpiType4','Median','IQR','GroupCount'}));
for i = 1:3
    fprintf('Pair %s vs %s: p=%.3g  (Bonferroni: %.3g)\n', ...
        SubtypePairs(i,1), SubtypePairs(i,2), p_pair_raw(i), p_pair_bonf(i));
end

%% ======================= FIGURE 2: SPEARMAN (log–log, fixed axes) =======================
PatientSpikeSz_All   = Views.PatientSpikeSz_All;
PatientSpikeSz_Typed = Views.PatientSpikeSz_Typed;

% Overall
x_all = PatientSpikeSz_All.MeanSpikeRate_perMin;
y_all = PatientSpikeSz_All.MeanSzFreq;
mask_all = isfinite(x_all) & isfinite(y_all);
[rs_all, p_all] = corr(x_all(mask_all), y_all(mask_all), 'Type','Spearman','Rows','complete');
n_all = sum(mask_all);

% By group
rowsOut = {};
for g = canonical3
    % get rows for this group
    m = (PatientSpikeSz_Typed.EpiType3 == g);
    x = PatientSpikeSz_Typed.MeanSpikeRate_perMin(m);
    y = PatientSpikeSz_Typed.MeanSzFreq(m);
    mask = isfinite(x) & isfinite(y);
    n = sum(mask);
    if n >= 3
        [rs, p] = corr(x(mask), y(mask), 'Type','Spearman','Rows','complete');
    else
        rs = NaN; p = NaN;
    end
    rowsOut(end+1,:) = {char(g), n, rs, p}; %#ok<SAGROW>
end
SpearmanResults = cell2table(rowsOut, 'VariableNames', {'Group','N','Spearman_r','p_raw'});
k = height(SpearmanResults);
SpearmanResults.p_bonf = min(SpearmanResults.p_raw * k, 1);

disp('=== Spearman correlations: SpikeRate_perMin vs MeanSzFreq (per month) ===');
disp([table("Overall (all epilepsy)", n_all, rs_all, p_all, NaN, ...
      'VariableNames', {'Group','N','Spearman_r','p_raw','p_bonf'}); SpearmanResults])

% Log transforms with eps guards
minpos_rate = min(PatientSpikeSz_All.MeanSpikeRate_perMin(PatientSpikeSz_All.MeanSpikeRate_perMin>0));
if isempty(minpos_rate), minpos_rate=1e-6; end
minpos_sz   = min(PatientSpikeSz_All.MeanSzFreq(PatientSpikeSz_All.MeanSzFreq>0));
if isempty(minpos_sz),   minpos_sz=1e-6; end
eps_rate = 0.5*minpos_rate;
eps_sz   = 0.5*minpos_sz;

Tall = table;
Tall.logSpikeRate = log10(PatientSpikeSz_All.MeanSpikeRate_perMin + ...
                          (PatientSpikeSz_All.MeanSpikeRate_perMin<=0).*eps_rate);
Tall.logSzFreq    = log10(PatientSpikeSz_All.MeanSzFreq + ...
                          (PatientSpikeSz_All.MeanSzFreq<=0).*eps_sz);

T = table;
T.EpiType3      = categorical(PatientSpikeSz_Typed.EpiType3, canonical3);
T.logSpikeRate  = log10(PatientSpikeSz_Typed.MeanSpikeRate_perMin + ...
                        (PatientSpikeSz_Typed.MeanSpikeRate_perMin<=0).*eps_rate);
T.logSzFreq     = log10(PatientSpikeSz_Typed.MeanSzFreq + ...
                        (PatientSpikeSz_Typed.MeanSzFreq<=0).*eps_sz);
ok = isfinite(T.logSpikeRate) & isfinite(T.logSzFreq) & ~ismissing(T.EpiType3);
T = T(ok,:);
presentCats = categories(removecats(T.EpiType3));

isZeroSz_all   = (PatientSpikeSz_All.MeanSzFreq==0);
isZeroRate_all = (PatientSpikeSz_All.MeanSpikeRate_perMin==0);
onlySz_all     =  isZeroSz_all & ~isZeroRate_all;
onlyRate_all   = ~isZeroSz_all &  isZeroRate_all;
bothZero_all   =  isZeroSz_all &  isZeroRate_all;
nonZero_all    = ~(isZeroSz_all | isZeroRate_all);

xZero = log10(eps_sz);
yZero = log10(eps_rate);
xLims = spearman_xLims;
yLims = spearman_yLims;
xrange = diff(xLims);
yrange = diff(yLims);
dx = 0.008*xrange;
dy = 0.008*yrange;

% ---- Draw figure ----
f2 = figure('Color','w','Position',[60 60 1200 900]);
tiledlayout(f2,2,2,'Padding','compact','TileSpacing','compact');
fontL = 20;

% A. Overall
axA2 = nexttile(1); hold(axA2,'on'); grid(axA2,'on'); box(axA2,'off');
xline(axA2, xZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
yline(axA2, yZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
baseColor = [0.4 0.4 0.4];
scatter(axA2, Tall.logSzFreq(nonZero_all), Tall.logSpikeRate(nonZero_all), 14, ...
    baseColor, 'filled', 'MarkerFaceAlpha', 0.25);
plot(axA2, Tall.logSzFreq(onlySz_all)+dx,   Tall.logSpikeRate(onlySz_all),   '*', ...
    'Color',baseColor,'MarkerSize',7,'LineWidth',1);
plot(axA2, Tall.logSzFreq(onlyRate_all),    Tall.logSpikeRate(onlyRate_all)+dy,'*', ...
    'Color',baseColor,'MarkerSize',8,'LineWidth',1);
plot(axA2, Tall.logSzFreq(bothZero_all)+dx, Tall.logSpikeRate(bothZero_all)+dy,'*', ...
    'Color',baseColor,'MarkerSize',8,'LineWidth',1.2);

X = [ones(sum(nonZero_all),1), Tall.logSzFreq(nonZero_all)];
b = X \ Tall.logSpikeRate(nonZero_all);
xgrid = linspace(xLims(1), xLims(2), 300)';
plot(axA2, xgrid, b(1) + b(2)*xgrid, 'k-', 'LineWidth', 2);

xlim(axA2, xLims); ylim(axA2, yLims);
xlabel(axA2,'log_{10} Seizures per month','FontSize',fontL);
ylabel(axA2,'log_{10} Spikes per minute','FontSize',fontL);
title(axA2, sprintf('A. All epilepsy  (N=%d) — %s', height(Tall), segLabel), ...
    'FontSize',fontL, 'FontWeight','bold');
txtA = sprintf('Spearman r=%.3f, p=%.3g', rs_all, p_all);
text(axA2, 0.98, 0.95, txtA, 'Units','normalized', ...
     'HorizontalAlignment','right','VerticalAlignment','top', ...
     'FontSize',fontL-2,'FontWeight','bold');
set(axA2,'FontSize',fontL);

% B/C/D panels in order: Frontal, Temporal, General
panelOrder  = {'Frontal','Temporal','General'};
panelTitle  = {'B. Frontal','C. Temporal','D. General'};

for p = 1:3
    ax = nexttile(p+1); hold(ax,'on'); grid(ax,'on'); box(ax,'off');

    if ~ismember(panelOrder{p}, presentCats)
        axis(ax,'off'); continue;
    end
    gStr   = string(panelOrder{p});
    tgtCat = categorical(gStr, categories(PatientSpikeSz_Typed.EpiType3));
    gi     = find(strcmp(presentCats, panelOrder{p}), 1);
    col    = lines(3); col = col(min(gi, size(col,1)),:);

    isZeroSz_g   = (PatientSpikeSz_Typed.MeanSzFreq==0);
    isZeroRate_g = (PatientSpikeSz_Typed.MeanSpikeRate_perMin==0);
    ok_g = ~isnan(PatientSpikeSz_Typed.MeanSzFreq) & ...
           ~isnan(PatientSpikeSz_Typed.MeanSpikeRate_perMin) & ...
           (PatientSpikeSz_Typed.EpiType3 == tgtCat);

    logX = log10(PatientSpikeSz_Typed.MeanSzFreq(ok_g) + ...
                 (PatientSpikeSz_Typed.MeanSzFreq(ok_g)<=0).*eps_sz);
    logY = log10(PatientSpikeSz_Typed.MeanSpikeRate_perMin(ok_g) + ...
                 (PatientSpikeSz_Typed.MeanSpikeRate_perMin(ok_g)<=0).*eps_rate);

    onlySz = isZeroSz_g(ok_g) & ~isZeroRate_g(ok_g);
    onlyRt = ~isZeroSz_g(ok_g) & isZeroRate_g(ok_g);
    bothZ  = isZeroSz_g(ok_g) & isZeroRate_g(ok_g);
    nonZ   = ~(onlySz | onlyRt | bothZ);

    xline(ax, xZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
    yline(ax, yZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    scatter(ax, logX(nonZ), logY(nonZ), 18, col, 'filled', 'MarkerFaceAlpha', 0.35);
    if any(onlySz)
        plot(ax, logX(onlySz)+dx, logY(onlySz),   '*','Color',col,'MarkerSize',8,'LineWidth',1.1);
    end
    if any(onlyRt)
        plot(ax, logX(onlyRt),    logY(onlyRt)+dy,'*','Color',col,'MarkerSize',8,'LineWidth',1.1);
    end
    if any(bothZ)
        plot(ax, logX(bothZ)+dx,  logY(bothZ)+dy, '*','Color',col,'MarkerSize',9,'LineWidth',1.2);
    end

    if nnz(nonZ) >= 3
        Xg = [ones(nnz(nonZ),1), logX(nonZ)];
        bg = Xg \ logY(nonZ);
        xg = linspace(xLims(1), xLims(2), 250)';
        plot(ax, xg, bg(1)+bg(2)*xg, '-', 'Color', col, 'LineWidth', 2);
    end

    xlim(ax, xLims); ylim(ax, yLims);

    row = SpearmanResults(strcmp(SpearmanResults.Group, string(panelOrder{p})), :);
    if ~isempty(row) && row.N >= 3 && isfinite(row.Spearman_r)
        txt = sprintf('Spearman r=%.3f, p_{bonf}=%.3g', row.Spearman_r, row.p_bonf);
    else
        txt = 'Insufficient data';
    end
    nNow = sum(T.EpiType3==tgtCat);
    title(ax, sprintf('%s (N=%d)', panelTitle{p}, nNow), 'FontSize', fontL, 'FontWeight','bold');
    text(ax, 0.98, 0.95, txt, 'Units','normalized', ...
         'HorizontalAlignment','right','VerticalAlignment','top', ...
         'FontSize', fontL-3, 'FontWeight','bold');
    set(ax,'FontSize',fontL);
end

if ~exist(fileparts(fig2_out),'dir'), mkdir(fileparts(fig2_out)); end
exportgraphics(f2, fig2_out, 'Resolution', 300);
fprintf('Saved Fig 2 (Spearman): %s\n', fig2_out);

if ~isempty(resultsCsv)
    if ~exist(fileparts(resultsCsv),'dir'), mkdir(fileparts(resultsCsv)); end
    ResultsOut = [table("Overall (all epilepsy)", n_all, rs_all, p_all, NaN, ...
        'VariableNames', {'Group','N','Spearman_r','p_raw','p_bonf'}); SpearmanResults];
    writetable(ResultsOut, resultsCsv);
    fprintf('Saved results to: %s\n', resultsCsv);
end

%% ======================= FIGURE 3: 2x2 PRESENCE TABLE & CHI-SQUARE =======================
make_spike_sz_presence_chi2(PatientSpikeSz_All);

%% ======================= FIGURE 4: SPEARMAN — NON-ZERO ONLY =======================
fig4_out = '../figures/spearman_spikerate_szfreq_log10_nonzero.png';
make_spearman_nonzero_fig(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
    canonical3, spearman_xLims, spearman_yLims, fig4_out);


%% ======================= HELPER FUNCTIONS ===================================

function S = ensure_spikerates(S, countCol, durCol)
% Ensure S has SpikeRate_perHour and SpikeRate_perMin using the selected segment.

if ~isnumeric(S.(countCol)), S.(countCol) = double(S.(countCol)); end
if ~isnumeric(S.(durCol)),   S.(durCol)   = double(S.(durCol));   end
badDur = ~isfinite(S.(durCol)) | S.(durCol) <= 0;
if any(badDur)
    warning('Found %d rows with nonpositive/NaN %s; spike rates set to NaN there.', nnz(badDur), char(durCol));
end
S.SpikeRate_perHour = nan(height(S),1);
ok = ~badDur;
S.SpikeRate_perHour(ok) = (S.(countCol)(ok) ./ S.(durCol)(ok)) * 3600;
S.SpikeRate_perMin  = S.SpikeRate_perHour / 60;
S.SpikeRate_perMin(S.SpikeRate_perMin > -1e-12 & S.SpikeRate_perMin < 0) = 0;
S.SpikeRate_perMin(S.SpikeRate_perMin < -1e-12) = NaN;
end


function [PatientTypingAll, SzFreqPerPatient] = build_patient_metrics_from_report(R, canonical3,hassz_replace)

% Erin checked this function 11/20

% Build patient-level typing (Type + Specific → EpiType3) and MeanSzFreq once.
if ~isstring(R.sz_freqs),         R.sz_freqs         = string(R.sz_freqs);         end
if ~isstring(R.visit_dates_deid), R.visit_dates_deid = string(R.visit_dates_deid); end
if ~isstring(R.visit_hasSz)
    R.visit_hasSz = string(R.visit_hasSz);
end

% ----- seizure frequency per patient (MeanSzFreq) -----
PV = table('Size',[0 4], 'VariableTypes',{'double','datetime','double','double'}, ...
           'VariableNames',{'Patient','VisitDate','Freq','HasSz'});

% Loop over eegs in redcap report
for j = 1:height(R)

    % Get patient id
    pid = double(R.patient_id(j));

    % date formatting
    ds = strtrim(R.visit_dates_deid(j)); dd = string([]);
    if strlength(ds)>0 && ds~="[]" && ds~=""
        try
            dd = string(jsondecode(char(ds)));
        catch
            dd = string(regexp(ds,'\d{4}-\d{2}-\d{2}','match'));
        end
    end
    try
        d = datetime(dd,'InputFormat','yyyy-MM-dd');
    catch
        d = datetime(dd);
    end

    % seizure frequency formatting
    s = strtrim(R.sz_freqs(j)); 
    if strlength(s)>0 && s~="[]" && s~=""
        s = regexprep(s,'null','NaN','ignorecase');
        v = double(jsondecode(char(s)));
    else
        v = [];
    end
    v(~isfinite(v)) = NaN; v(v<0) = NaN; % -2 turns into NaN here

    % hasSz
    hs = strtrim(R.visit_hasSz(j)); h = nan(0,1);
    if strlength(hs)>0 && hs~="[]" && hs~=""
        h = double(jsondecode(char(hs)));
    else
        h = [];
    end
    h(h==2) = nan; % set 2 to nan

    % number of visits should equal number of sz freqs and number of has sz
    assert(length(d)==length(v))
    assert(length(d) == length(h))
    n = length(h);
    if n==0, continue; end
    PV = [PV; table(repmat(pid,n,1), d(1:n), v(1:n), h(1:n), ...
        'VariableNames', PV.Properties.VariableNames)]; %#ok<AGROW>
end

% This bit of logic handles the fact that every eeg row lists ALL visits
% for that patient, and so this returns just unique visits by grouping by
% patient id + visit date
[gv, pid_keys, date_keys] = findgroups(PV.Patient, PV.VisitDate); % group by patient and visit date

% ---- Sanity check: all finite frequency values and HasSz in each group must be identical ----
numGroups = max(gv);
for g = 1:numGroups
    xf = PV.Freq(gv == g);
    xf = xf(isfinite(xf));

    if numel(xf) > 1
        % tolerance for float noise
        if max(xf) - min(xf) > 1e-12
            error('Freq mismatch in group %d (Patient %g, Date %s): values = %s', ...
                g, pid_keys(g), string(date_keys(g)), mat2str(xf));
        end
    end

    xh = PV.HasSz(gv == g);
    xh = xh(isfinite(xh));   % ignore NaN (missing)
    if numel(xh) > 1
        if max(xh) ~= min(xh)
            error('HasSz mismatch in group %d (Patient %g, Date %s): values = %s', ...
                g, pid_keys(g), string(date_keys(g)), mat2str(xh));
        end
    end

end

Freq_agg = splitapply(@(x) mean(x(isfinite(x)),'omitnan'), PV.Freq, gv); % Get a single unique sz frequency for each visit
Has_agg  = splitapply(@(x) max_hasSz(x(isfinite(x))), PV.HasSz, gv); % get a single unique has sz for each visit

% build table for unique visits - I can check that this agrees with redcap
% (I did a few spot checks)
Vuniq = table(pid_keys, date_keys, Freq_agg, Has_agg, ...
    'VariableNames', {'Patient','VisitDate','Freq','HasSz'});
Vuniq.Freq(Vuniq.Freq<0) = NaN; % negative sz frequencies should be nans

if hassz_replace == 1
    % If visit frequency is NaN but Has Sz is 0, change visit frequency to zero!!
    Vuniq.Freq(~isfinite(Vuniq.Freq) & Vuniq.HasSz==0) = 0; 
end

% Now group by patient
[gpV,pidsV] = findgroups(Vuniq.Patient);

% Get the mean sz frequency across visits for each patient, removing nans
MeanSzFreq  = splitapply(@(x) mean(x,'omitnan'), Vuniq.Freq, gpV);

% Get mean sz frequency per patient
SzFreqPerPatient = table(pidsV, MeanSzFreq, 'VariableNames', {'Patient','MeanSzFreq'});

% ----- typing (Type + Specific → EpiType3) -----
RtType = sortrows(R(~ismissing(R.epilepsy_type) & ...
                    strlength(strtrim(R.epilepsy_type))>0, ...
                    {'patient_id','epilepsy_type'}), 'patient_id');
[uid_t, ia_t] = unique(RtType.patient_id,'stable');
PatientTypingAll = table(double(uid_t), string(RtType.epilepsy_type(ia_t)), ...
                         'VariableNames', {'Patient','EpilepsyType'});

RtSpec = sortrows(R(~ismissing(R.epilepsy_specific) & ...
                    strlength(strtrim(R.epilepsy_specific))>0, ...
                    {'patient_id','epilepsy_specific'}), 'patient_id');
[uid_s, ia_s] = unique(RtSpec.patient_id,'stable');
PatientTypingAll = outerjoin(PatientTypingAll, ...
    table(double(uid_s), string(RtSpec.epilepsy_specific(ia_s)), ...
          'VariableNames', {'Patient','EpilepsySpecific'}), ...
    'Keys','Patient','MergeKeys',true);

spec_norm = lower(strtrim(PatientTypingAll.EpilepsySpecific));
type_norm = lower(strtrim(PatientTypingAll.EpilepsyType));
isTemporal = contains(spec_norm,"temporal");
isFrontal  = contains(spec_norm,"frontal");
isGeneral  = contains(type_norm,"general");

SpecCanon = strings(height(PatientTypingAll),1);
SpecCanon(isTemporal) = "Temporal";
SpecCanon(isFrontal)  = "Frontal";
SpecCanon((strlength(SpecCanon)==0) & isGeneral) = "General";
PatientTypingAll.EpiType3 = categorical(SpecCanon, canonical3);
end

function out = max_hasSz(x)
    xf = x(isfinite(x));
    if isempty(xf)
        out = NaN;   % this visit has no usable hasSz info
    else
        out = max(xf);
    end
end

function Views = build_filtered_view(SessionsIn, ReportIn, PatientTypingAll, SzFreqPerPatient, ...
                                     only_amb, countCol, durCol, NESD_LABEL, badTypes, canonical3)
% Erin checked this function 11/21

% Returns all filtered, ready-to-plot pieces in one struct.
% Note: OUTPATIENT restriction is assumed to be already applied up front.

% --- apply only_amb to SessionsIn ---
SessionsFiltered = SessionsIn;
if only_amb == 1
    SessionsFiltered(SessionsFiltered.Duration_sec < 12*3600,:) = [];
elseif only_amb == 2
    SessionsFiltered(SessionsFiltered.Duration_sec > 12*3600,:) = [];
end
if ~isnumeric(SessionsFiltered.Patient)
    SessionsFiltered.Patient = double(str2double(string(SessionsFiltered.Patient)));
end
if ~isnumeric(SessionsFiltered.Session)
    SessionsFiltered.Session = double(str2double(string(SessionsFiltered.Session)));
end

% --- restrict report rows to sessions in SessionsFiltered ---
RSess_keep = unique(SessionsFiltered(:,{'Patient','Session'}));
ReportForKeptSessions = innerjoin(ReportIn, RSess_keep, ...
    'LeftKeys', {'patient_id','session_number'}, ...
    'RightKeys', {'Patient','Session'});
ReportForKeptSessions.Patient = ReportForKeptSessions.patient_id;
ReportForKeptSessions.Session = ReportForKeptSessions.session_number;

% --- patient typing for filtered set (subset of PatientTypingAll) ---
PatientsKept       = unique(SessionsFiltered(:,{'Patient'}));
TypingFiltered     = innerjoin(PatientTypingAll, PatientsKept, 'Keys','Patient');

% --- per patient mean spike rate (selected segment) ---
[gpS, pidsS] = findgroups(SessionsFiltered.Patient);
MeanSpikeRate_perHour = splitapply(@(x) mean(x,'omitnan'), SessionsFiltered.SpikeRate_perHour, gpS);
MeanSpikeRate_perMin  = splitapply(@(x) mean(x,'omitnan'), SessionsFiltered.SpikeRate_perMin,  gpS);
PatientLevelSpikeRates = table(pidsS, MeanSpikeRate_perHour, MeanSpikeRate_perMin, ...
                               'VariableNames', {'Patient','MeanSpikeRate_perHour','MeanSpikeRate_perMin'});

% attach typing to patient-level spike rates
PatientLevelSpikeRates = innerjoin(PatientLevelSpikeRates, TypingFiltered, 'Keys','Patient');

% --- per-session spikes/min for Panel C (use selected segment columns) ---
[grpPS, pKeys, sKeys] = findgroups(SessionsFiltered.Patient, SessionsFiltered.Session);
nanSum = @(x) sum(x(~isnan(x)));
Sum_Total = splitapply(nanSum, SessionsFiltered.(countCol), grpPS);
Sum_Dur   = splitapply(nanSum, SessionsFiltered.(durCol),   grpPS);
SessionLevelSpikeRates = table(pKeys, sKeys, Sum_Total, Sum_Dur, ...
    'VariableNames', {'Patient','Session','Sum_Total','Sum_Dur'});
SessionLevelSpikeRates.SpikesPerMin = nan(height(SessionLevelSpikeRates),1);
validDur = SessionLevelSpikeRates.Sum_Dur > 0;
SessionLevelSpikeRates.SpikesPerMin(validDur) = ...
    (SessionLevelSpikeRates.Sum_Total(validDur) ./ SessionLevelSpikeRates.Sum_Dur(validDur)) * 60;

% --- Build PatientSpikeSz (for Spearman) ---
SzFreqFiltered = innerjoin(SzFreqPerPatient, PatientsKept, 'Keys','Patient');
PatientSpikeSz_All = innerjoin( ...
    PatientLevelSpikeRates(:,{'Patient','MeanSpikeRate_perMin'}), ...
    SzFreqFiltered, 'Keys','Patient');
PatientSpikeSz_All = PatientSpikeSz_All(isfinite(PatientSpikeSz_All.MeanSzFreq) & ...
                                        isfinite(PatientSpikeSz_All.MeanSpikeRate_perMin), :);

keep3 = ~ismissing(TypingFiltered.EpiType3) & ismember(TypingFiltered.EpiType3, canonical3);
PatientSpikeSz_Typed = innerjoin( ...
    PatientLevelSpikeRates(:,{'Patient','MeanSpikeRate_perMin'}), ...
    innerjoin(SzFreqFiltered, TypingFiltered(keep3,{'Patient','EpiType3'}),'Keys','Patient'), ...
    'Keys','Patient');
PatientSpikeSz_Typed = PatientSpikeSz_Typed(isfinite(PatientSpikeSz_Typed.MeanSzFreq) & ...
                                            isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perMin) & ...
                                            ~ismissing(PatientSpikeSz_Typed.EpiType3), :);

% --- Precompute items needed for Fig 1 summaries ---
etype_norm = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));
IsEpilepsyMask = ~(etype_norm==lower(NESD_LABEL)) & ...
                 ~ismember(etype_norm, badTypes) & ...
                 strlength(PatientLevelSpikeRates.EpilepsyType)>0;
IsNESDMask     =  (etype_norm==lower(NESD_LABEL));

EpiType4 = strings(height(PatientLevelSpikeRates),1);
et_low = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));
es_low = lower(strtrim(string(PatientLevelSpikeRates.EpilepsySpecific)));
isGeneralType = contains(et_low,"general");
isTemporal    = contains(es_low,"temporal");
isFrontal     = contains(es_low,"frontal");
EpiType4(isGeneralType) = "General";
EpiType4(~isGeneralType & isTemporal) = "Temporal";
EpiType4(~isGeneralType & ~isTemporal & isFrontal) = "Frontal";
PatientLevelSpikeRates.EpiType4 = EpiType4;

InCanonical3Mask = IsEpilepsyMask & ismember(PatientLevelSpikeRates.EpiType4, canonical3);

Canonical3_SubsetTable = PatientLevelSpikeRates(InCanonical3Mask, {'EpiType4','MeanSpikeRate_perHour'});
Canonical3_SubsetTable.EpiType4 = categorical(string(Canonical3_SubsetTable.EpiType4), canonical3);
[g3, cats3]  = findgroups(Canonical3_SubsetTable.EpiType4);
medVals      = splitapply(@(x) median(x,'omitnan'), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
iqrVals      = splitapply(@(x) prctile(x,75) - prctile(x,25), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
nVals        = splitapply(@(x) sum(isfinite(x)), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
Canonical3_Stats = table(cats3, nVals, medVals, iqrVals, ...
    'VariableNames', {'EpiType4','GroupCount','Median','IQR'});

Canonical3_Pairs = ["General","Temporal"; "General","Frontal"; "Temporal","Frontal"];
p_pair = NaN(3,1);
for i = 1:3
    A = Canonical3_Pairs(i,1); B = Canonical3_Pairs(i,2);
    xa = Canonical3_SubsetTable.MeanSpikeRate_perHour(Canonical3_SubsetTable.EpiType4==A);
    xb = Canonical3_SubsetTable.MeanSpikeRate_perHour(Canonical3_SubsetTable.EpiType4==B);
    if nnz(isfinite(xa))>=3 && nnz(isfinite(xb))>=3
        p_pair(i) = ranksum(xa, xb, 'method','approx');
    end
end
p_pair_bonf = min(p_pair*3, 1);

% ---------- Pack ----------
Views.SessionsForFigures     = SessionsFiltered;
Views.ReportForKeptSessions  = ReportForKeptSessions;
Views.PatientTypingFiltered  = TypingFiltered;

Views.SessionLevelSpikeRates = SessionLevelSpikeRates;
Views.PatientLevelSpikeRates = PatientLevelSpikeRates;

Views.PatientSpikeSz_All     = PatientSpikeSz_All;
Views.PatientSpikeSz_Typed   = PatientSpikeSz_Typed;

Views.IsEpilepsyMask         = IsEpilepsyMask;
Views.IsNESDMask             = IsNESDMask;

Views.Canonical3_SubsetTable = Canonical3_SubsetTable;
Views.Canonical3_Stats       = Canonical3_Stats;
Views.Canonical3_Pairs       = Canonical3_Pairs;
Views.PvalsPairwise          = p_pair;
Views.PvalsPairwiseBonf      = p_pair_bonf;
end

function compare_spikerate_by_modality(Sessions, eps_per_min)
% Box+swarm plot of log10(spikes/min) for Ambulatory vs Routine EEGs (FULL view)
assert(ismember('Duration_sec',Sessions.Properties.VariableNames) && ...
       ismember('SpikeRate_perMin',Sessions.Properties.VariableNames), ...
    'Sessions must contain Duration_sec and SpikeRate_perMin.');
isAmb  = Sessions.Duration_sec >= 12*3600;
isRout = Sessions.Duration_sec <= 12*3600;
xAmb   = Sessions.SpikeRate_perMin(isAmb);
xRout  = Sessions.SpikeRate_perMin(isRout);
xAmbLog  = log10(xAmb  + (xAmb<=0).*eps_per_min);
xRoutLog = log10(xRout + (xRout<=0).*eps_per_min);
[p_rs,~] = ranksum(xAmbLog,xRoutLog,'method','approx');

figure('Color','w','Position',[200 200 600 500]);
hold on; grid on; box off;
cats = categorical({'Routine','Ambulatory'});
vals = [xRoutLog; xAmbLog];
grps = [repmat(cats(1),numel(xRoutLog),1); repmat(cats(2),numel(xAmbLog),1)];
boxchart(grps,vals,'BoxFaceAlpha',0.25);
swarmchart(grps,vals,18,'filled','MarkerFaceAlpha',0.45);
yline(log10(eps_per_min),':','Color',[0.4 0.4 0.4],'LineWidth',1.2);
ylabel('log_{10}(spikes/min)','FontSize',14);
title(sprintf('Spike rate by modality  (rank-sum p = %.3g)',p_rs),'FontSize',16);
set(gca,'FontSize',14);

fprintf('Spike-rate comparison:\n  Routine N=%d, Ambulatory N=%d, rank-sum p=%.3g\n', ...
        numel(xRoutLog), numel(xAmbLog), p_rs);
end
function compare_szfreq_by_modality(Sessions, SzFreqPerPatient)
% Box+swarm plot of log10(seizures/month) for Ambulatory vs Routine EEGs (FULL view)
Sess = unique(Sessions(:,{'Patient','Session','Duration_sec'}));
Sess.Modality = repmat("Routine",height(Sess),1);
Sess.Modality(Sess.Duration_sec >= 12*3600) = "Ambulatory";
JR = innerjoin(Sess,SzFreqPerPatient(:,{'Patient','MeanSzFreq'}),'Keys','Patient');
fAmb  = JR.MeanSzFreq(JR.Modality=="Ambulatory");
fRout = JR.MeanSzFreq(JR.Modality=="Routine");
fAmb(fAmb<=0)  = NaN;
fRout(fRout<=0) = NaN;
eps_sz = 0.5*min([min(fAmb,[],'omitnan'),min(fRout,[],'omitnan'),1e-6]);
xAmb  = log10(fAmb  + (fAmb<=0).*eps_sz);
xRout = log10(fRout + (fRout<=0).*eps_sz);
[p_rs,~] = ranksum(xAmb,xRout,'method','approx');

figure('Color','w','Position',[200 200 600 500]);
hold on; grid on; box off;
cats = categorical({'Routine','Ambulatory'});
vals = [xRout; xAmb];
grps = [repmat(cats(1),numel(xRout),1); repmat(cats(2),numel(xAmb),1)];
boxchart(grps,vals,'BoxFaceAlpha',0.25);
swarmchart(grps,vals,18,'filled','MarkerFaceAlpha',0.45);
ylabel('log_{10}(seizures/month)','FontSize',14);
title(sprintf('Seizure frequency by modality  (rank-sum p = %.3g)',p_rs),'FontSize',16);
set(gca,'FontSize',14);

fprintf('Seizure-frequency comparison:\n  Routine N=%d, Ambulatory N=%d, rank-sum p=%.3g\n', ...
        numel(xRout), numel(xAmb), p_rs);
end

%% ======================= HELPERS =======================
function ylog = to_log10_per_min(x_per_hour, eps_per_min)
pm = double(x_per_hour) / 60;
pm(~isfinite(pm) | pm <= 0) = eps_per_min;
ylog = log10(pm);
end

function ylog = to_log10_from_per_min(x_per_min, eps_per_min)
pm = double(x_per_min);
pm(~isfinite(pm) | pm <= 0) = eps_per_min;
ylog = log10(pm);
end

function add_sigbar(ax, x1, x2, y, ptext)
tick = 0.03 * diff(ax.YLim);
plot(ax, [x1 x1 x2 x2], [y- tick, y, y, y- tick], 'k-', 'LineWidth', 1.3);
text(ax, mean([x1 x2]), y + 0.02*diff(ax.YLim), ptext, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);
end

function make_spike_sz_presence_chi2(PatientSpikeSz_All)
% make_spike_sz_presence_chi2
%   Builds a 2x2 contingency table for:
%       spikes present = MeanSpikeRate_perMin > 0
%       seizures present = MeanSzFreq > 0
%   Runs chi-square test of independence and shows a figure + console stats.

rate = PatientSpikeSz_All.MeanSpikeRate_perMin;
freq = PatientSpikeSz_All.MeanSzFreq;

valid = isfinite(rate) & isfinite(freq);
rate = rate(valid);
freq = freq(valid);

hasSpikes = rate > 0;
hasSz     = freq > 0;

x = categorical(hasSpikes, [false true], {'No spikes','Spikes>0'});
y = categorical(hasSz,     [false true], {'No seizures','Seizures>0'});

[tbl, chi2, p, ~] = crosstab(x, y);

% Console output
fprintf('\n=== 2x2 presence table: spikes>0 vs seizures>0 (patient-level) ===\n');
rowNames = {'NoSpikes','SpikesGT0'};
colNames = {'NoSeizures','SeizuresGT0'};
disp(array2table(tbl, 'RowNames', rowNames, 'VariableNames', colNames));
df = (size(tbl,1)-1) * (size(tbl,2)-1);
fprintf('Chi-square(%d) = %.3f, p = %.3g\n', df, chi2, p);

% Figure
f = figure('Color','w','Position',[200 200 480 420]);
imagesc(tbl);
colorbar;
axis equal tight;
set(gca,'XTick',1:2,'XTickLabel',colNames, ...
        'YTick',1:2,'YTickLabel',rowNames, ...
        'FontSize',14);
xlabel('Seizures (MeanSzFreq>0)');
ylabel('Spikes (MeanSpikeRate>0)');
title(sprintf('Spikes>0 vs Seizures>0  (\\chi^2=%.2f, p=%.3g)', chi2, p), ...
      'FontSize',16,'FontWeight','bold');

% overlay counts
for i = 1:2
    for j = 1:2
        text(j, i, sprintf('%d', tbl(i,j)), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',16,'FontWeight','bold','Color','k');
    end
end

end

function make_spearman_nonzero_fig(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
                                   canonical3, xLims, yLims, outPng)
% make_spearman_nonzero_fig
%   Repeats the spike-vs-seizure Spearman correlations using ONLY
%   patients with non-zero spikes and non-zero seizures.
%
%   Uses log10(MeanSzFreq) vs log10(MeanSpikeRate_perMin)
%   Panels: A = all epilepsy, B/C/D = Frontal / Temporal / General.

fontL = 20;

% Ensure canonical3 is a cell array of char for categorical()
if isstring(canonical3)
    canonical3 = cellstr(canonical3);
end

%% --- Overall (all epilepsy) ---
rate_all = PatientSpikeSz_All.MeanSpikeRate_perMin;
freq_all = PatientSpikeSz_All.MeanSzFreq;

mask_all = isfinite(rate_all) & isfinite(freq_all) & ...
           (rate_all > 0) & (freq_all > 0);

x_all = log10(freq_all(mask_all));
y_all = log10(rate_all(mask_all));
n_all = numel(x_all);

if n_all >= 3
    [rs_all, p_all] = corr(x_all, y_all, 'Type','Spearman','Rows','complete');
else
    rs_all = NaN; p_all = NaN;
end

%% --- By group (EpiType3) with non-zero only ---
rowsOut = {};
for g = canonical3
    % g is a cellstr element; convert to char for comparison
    gStr = char(g);

    m = strcmp(string(PatientSpikeSz_Typed.EpiType3), gStr) & ...
        isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perMin) & ...
        isfinite(PatientSpikeSz_Typed.MeanSzFreq) & ...
        (PatientSpikeSz_Typed.MeanSpikeRate_perMin > 0) & ...
        (PatientSpikeSz_Typed.MeanSzFreq > 0);

    x = log10(PatientSpikeSz_Typed.MeanSzFreq(m));
    y = log10(PatientSpikeSz_Typed.MeanSpikeRate_perMin(m));
    n = numel(x);
    if n >= 3
        [rs, p] = corr(x, y, 'Type','Spearman','Rows','complete');
    else
        rs = NaN; p = NaN;
    end
    rowsOut(end+1,:) = {gStr, n, rs, p}; %#ok<SAGROW>
end
SpearmanResults_nz = cell2table(rowsOut, 'VariableNames', {'Group','N','Spearman_r','p_raw'});
k = height(SpearmanResults_nz);
SpearmanResults_nz.p_bonf = min(SpearmanResults_nz.p_raw * k, 1);

% Console summary
fprintf('\n=== Spearman (NON-ZERO ONLY): log10(spikes/min) vs log10(seizures/month) ===\n');
disp([table("Overall (all epilepsy, non-zero only)", n_all, rs_all, p_all, NaN, ...
      'VariableNames', {'Group','N','Spearman_r','p_raw','p_bonf'}); SpearmanResults_nz])

%% --- Build typed table for plotting ---
mask_typed = isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perMin) & ...
             isfinite(PatientSpikeSz_Typed.MeanSzFreq) & ...
             (PatientSpikeSz_Typed.MeanSpikeRate_perMin > 0) & ...
             (PatientSpikeSz_Typed.MeanSzFreq > 0) & ...
             ~ismissing(PatientSpikeSz_Typed.EpiType3);

T = table;
T.EpiType3     = categorical(string(PatientSpikeSz_Typed.EpiType3(mask_typed)), canonical3);
T.logSpikeRate = log10(PatientSpikeSz_Typed.MeanSpikeRate_perMin(mask_typed));
T.logSzFreq    = log10(PatientSpikeSz_Typed.MeanSzFreq(mask_typed));

presentCats = categories(removecats(T.EpiType3));

%% --- Draw figure ---
f = figure('Color','w','Position',[60 60 1200 900]);
tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');

% A. Overall
axA = nexttile(1); hold(axA,'on'); grid(axA,'on'); box(axA,'off');
scatter(axA, x_all, y_all, 18, [0.3 0.3 0.3], ...
    'filled','MarkerFaceAlpha',0.35);
if n_all >= 3
    X = [ones(n_all,1), x_all(:)];
    b = X \ y_all(:);
    xgrid = linspace(xLims(1), xLims(2), 300)';
    plot(axA, xgrid, b(1) + b(2)*xgrid, 'k-', 'LineWidth', 2);
end
xlim(axA, xLims); ylim(axA, yLims);
xlabel(axA,'log_{10} Seizures per month','FontSize',fontL);
ylabel(axA,'log_{10} Spikes per minute','FontSize',fontL);
title(axA, sprintf('A. All epilepsy (N=%d, non-zero only)', n_all), ...
    'FontSize',fontL,'FontWeight','bold');
txtA = sprintf('Spearman r=%.3f, p=%.3g', rs_all, p_all);
text(axA, 0.98, 0.95, txtA, 'Units','normalized', ...
     'HorizontalAlignment','right','VerticalAlignment','top', ...
     'FontSize',fontL-2,'FontWeight','bold');
set(axA,'FontSize',fontL);

% B/C/D panels: Frontal, Temporal, General
panelOrder = {'Frontal','Temporal','General'};
panelTitle = {'B. Frontal','C. Temporal','D. General'};
cols = lines(3);

for p = 1:3
    ax = nexttile(p+1); hold(ax,'on'); grid(ax,'on'); box(ax,'off');

    if ~ismember(panelOrder{p}, presentCats)
        axis(ax,'off');
        continue;
    end

    % One-category categorical target
    tgt = categorical(panelOrder(p), canonical3);
    m   = (T.EpiType3 == tgt);
    xg  = T.logSzFreq(m);
    yg  = T.logSpikeRate(m);
    nNow = numel(xg);

    if nNow == 0
        axis(ax,'off');
        continue;
    end

    col = cols(min(p,size(cols,1)),:);
    scatter(ax, xg, yg, 24, col, 'filled','MarkerFaceAlpha',0.4);

    if nNow >= 3
        Xg = [ones(nNow,1), xg(:)];
        bg = Xg \ yg(:);
        xgrid = linspace(xLims(1), xLims(2), 250)';
        plot(ax, xgrid, bg(1) + bg(2)*xgrid, '-', 'Color', col, 'LineWidth', 2);
    end

    xlim(ax, xLims); ylim(ax, yLims);
    row = SpearmanResults_nz(strcmp(SpearmanResults_nz.Group, panelOrder{p}), :);
    if ~isempty(row) && row.N >= 3 && isfinite(row.Spearman_r)
        txt = sprintf('Spearman r=%.3f, p_{bonf}=%.3g', row.Spearman_r, row.p_bonf);
    else
        txt = 'Insufficient data';
    end

    title(ax, sprintf('%s (N=%d, non-zero only)', panelTitle{p}, nNow), ...
        'FontSize',fontL,'FontWeight','bold');
    text(ax, 0.98, 0.95, txt, 'Units','normalized', ...
         'HorizontalAlignment','right','VerticalAlignment','top', ...
         'FontSize',fontL-3,'FontWeight','bold');
    set(ax,'FontSize',fontL);
end

if nargin >= 6 && ~isempty(outPng)
    if ~exist(fileparts(outPng),'dir'), mkdir(fileparts(outPng)); end
    exportgraphics(f, outPng, 'Resolution', 300);
    fprintf('Saved Fig (Spearman non-zero only): %s\n', outPng);
end

end
