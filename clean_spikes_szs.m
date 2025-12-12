
%% ======================= CONFIG =======================

% Paths to data
spikeSummaryMultiCsv = '../data/SN_counts/spike_counts_summary_multiThresh.csv';
reportCsv            = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-12-11_1316.csv';

% Max number of hours to call something a routine (as opposed to
% ambulatory)
MAX_ROUTINE_HOURS = 4;

% Label constants
NESD_LABEL = "Non-Epileptic Seizure Disorder";
badTypes   = lower(["Uncertain if Epilepsy","Unknown or MRN not found",""]); % excluded from 'Epilepsy'
canonical3 = ["General","Temporal","Frontal"];   % canonical 3 subtypes of epilepsy

% Plotting params (shared)
EPS_RATE = 30e-3;                 % → dotted "zero" line at y = -3 in log10(spikes/min)
Y_ZERO      = log10(EPS_RATE);
Y_LIMS      = [-2 4];             % fixed y-lims for box/swarm plots

% Spearman-figure axes (fixed across panels)
spearman_xLims = [-3.5, 4];         % log10(seizures/month)
spearman_yLims = [-1.5, 3];           % log10(spikes/min)

% Outputs
figHasSzFrac_out = '../figures/spikerate_by_HasSzFraction.png';
fig1_out   = '../figures/Fig1.png';
fig2_out   = '../figures/Fig2.png'; 
figS1_out   = '../figures/FigS1.png';
figS2_out = '../figures/FigS2.png';
figS3_out = '../figures/FigS3.png';
Table1Csv = fullfile('..','output','Table1.csv');
resultsHtml = '../output/results_summary.html';


%% ======================= LOAD REDCAP REPORT  =======================
% ---------- Load raw report (for outpatient rule + everything else) ----------
ReportTable = readtable(reportCsv,'TextType','string','VariableNamingRule','preserve');
acqStr_report = lower(strtrim(string(ReportTable.acquired_on)));

%% 1) Allowed visit types
allowable_visits = ["CONSULT VISIT","ESTABLISHED PATIENT VISIT",...
    "FOLLOW-UP PATIENT CLINIC",...
    "NEW PATIENT CLINIC","NEW PATIENT VISIT",...
    "NPV MANAGEMENT DURING COVID-19",...
    "NPV NEUROLOGY",...
    "RETURN ANNUAL VISIT","RETURN PATIENT EXTENDED",...
    "RETURN PATIENT VISIT","RPV MANAGEMENT DURING COVID-19","TELEHEALTH VIDEO VISIT RETURN",...
    ];

%% ======================= VISIT-TYPE FILTER =======================
% Drop any clinic visits whose visit_type is NOT in allowable_visits.
% This rewrites visit_dates_deid, visit_hasSz, sz_freqs, and visit_type
% so that all downstream code only "sees" allowed visit types.

totalVisits_before = 0;
totalVisits_after  = 0;

% ---------- Clean rows where visit_type is only [null] ----------
vt_raw_all = strtrim(string(ReportTable.visit_type));
mask_null_only = (vt_raw_all == "[null]") | (vt_raw_all == "null");
ReportTable.visit_type(mask_null_only)        = "[]";
ReportTable.visit_dates_deid(mask_null_only)  = "[]";
ReportTable.sz_freqs(mask_null_only)          = "[]";
ReportTable.visit_hasSz(mask_null_only)       = "[]";

% loop over redcap report
for i = 1:height(ReportTable)

    % Raw JSON-ish strings for this EEG row with visit-level data
    vt_raw    = strtrim(string(ReportTable.visit_type(i)));
    dates_raw = strtrim(string(ReportTable.visit_dates_deid(i)));
    sz_raw    = strtrim(string(ReportTable.sz_freqs(i)));
    hs_raw    = strtrim(string(ReportTable.visit_hasSz(i)));

    % If no usable visit_type content → wipe all visit arrays
    if strlength(vt_raw)==0 || vt_raw=="[]" || vt_raw=="<missing>"
        % Normalize to empty arrays
        ReportTable.visit_type(i)        = "[]";
        ReportTable.visit_dates_deid(i)  = "[]";
        ReportTable.sz_freqs(i)          = "[]";
        ReportTable.visit_hasSz(i)       = "[]";
        continue
    end


    % ---------- 1) Decode visit_type into a string array vt ----------
    vt_dec = jsondecode(char(vt_raw));   % can be cell, string, numeric, etc.
    

    % Convert vt_dec to a column string array vt(:)
    if iscell(vt_dec)
        % Typical case: cell array like {'EEG'; []; 'NEW PATIENT VISIT'}
        n = numel(vt_dec);
        vt = strings(n,1);
        for k = 1:n
            x = vt_dec{k};
            if ischar(x) || (isstring(x) && isscalar(x))
                vt(k) = string(x);
            else
                % null / empty / non-string → treat as empty label
                vt(k) = "";
            end
        end
    elseif ischar(vt_dec) || (isstring(vt_dec) && isscalar(vt_dec))
        % Single string, e.g. "EEG"
        vt = string(vt_dec);
        n  = numel(vt);
    elseif isnumeric(vt_dec)
        if isscalar(vt_dec) && isnan(vt_dec)
            % Entire visit_type unusable → wipe all visit arrays
            ReportTable.visit_type(i)        = "[]";
            ReportTable.visit_dates_deid(i)  = "[]";
            ReportTable.sz_freqs(i)          = "[]";
            ReportTable.visit_hasSz(i)       = "[]";
            continue
        else
            % numeric array (rare) → convert to strings
            vt = string(vt_dec(:));
            n  = numel(vt);
        end

    else
        error('weird format')
    end

    % ---------- 2) Decode the other arrays ----------
    % Dates
    if strlength(dates_raw)>0 && dates_raw~="[]"
        dates_cell = jsondecode(char(dates_raw));    % usually a cell array of date strings
    else
        dates_cell = {};
    end

    % sz_freqs (may contain nulls → NaN)
    if strlength(sz_raw)>0 && sz_raw~="[]"
        sz_cell = jsondecode(char(sz_raw));
    else
        sz_cell = [];
    end

    % visit_hasSz
    if strlength(hs_raw)>0 && hs_raw~="[]"
        hs_cell = jsondecode(char(hs_raw));
    else
        hs_cell = [];
    end

    % ---------- 3) Align lengths (take the minimum across arrays) ----------
    len_vt    = numel(vt);
    len_dates = numel(dates_cell);
    len_sz    = numel(sz_cell);
    len_hs    = numel(hs_cell);

    n_all = len_vt;

    if n_all == 0
        % No usable visits → wipe the row to consistent empty lists
        ReportTable.visit_type(i)        = "[]";
        ReportTable.visit_dates_deid(i)  = "[]";
        ReportTable.sz_freqs(i)          = "[]";
        ReportTable.visit_hasSz(i)       = "[]";
        continue
    end


    % If not all equal, throw an error
    if ~(len_vt==len_dates && len_vt==len_sz && len_vt==len_hs)
        error('mismatched lengths')
    end

    vt_use    = vt(1:n_all);
    dates_use = dates_cell(1:n_all);
    sz_use    = sz_cell(1:n_all);
    hs_use    = hs_cell(1:n_all);

    totalVisits_before = totalVisits_before + n_all;

    % ---------- 4) Keep only allowable visit types ----------
    keepMask = ismember(vt_use, allowable_visits);

    if ~any(keepMask)
        % No allowable visits for this EEG row → set arrays to empty JSON lists
        ReportTable.visit_dates_deid(i) = "[]";
        ReportTable.visit_hasSz(i)      = "[]";
        ReportTable.sz_freqs(i)         = "[]";
        ReportTable.visit_type(i)       = "[]";
        continue
    end

    % only keep those where keepMask == 1
    vt_filt    = cellstr(vt_use(keepMask));   % back to cellstr for jsonencode
    dates_filt = dates_use(keepMask);
    sz_filt    = sz_use(keepMask);
    hs_filt    = hs_use(keepMask);

    totalVisits_after = totalVisits_after + numel(vt_filt);

    % ---------- 5) Encode back to JSON strings ----------
    ReportTable.visit_dates_deid(i) = string(jsonencode(dates_filt));
    ReportTable.visit_hasSz(i)      = string(jsonencode(hs_filt));
    ReportTable.sz_freqs(i)         = string(jsonencode(sz_filt));
    ReportTable.visit_type(i)       = string(jsonencode(vt_filt));

end

fprintf('[Visit-type filter] Total clinic visits before filter: %d\n', totalVisits_before);
fprintf('[Visit-type filter] Total clinic visits after filter:  %d (kept %.1f%%)\n', ...
    totalVisits_after, 100*totalVisits_after/max(1,totalVisits_before));

%% ======================= Load spike summary csv =======================
SpikeSummaryTable = readtable(spikeSummaryMultiCsv,'TextType','string','VariableNamingRule','preserve');
% Get columns for spike rates
countCol = "count_0_46"; % this is the default spike prob threshold for SN2
durCol   = "Duration_sec";

% Convert spike count to spike rate in spikes/hour
SpikeSummaryTable.SpikeRate_perHour = nan(height(SpikeSummaryTable),1);
SpikeSummaryTable.SpikeRate_perHour = SpikeSummaryTable.(countCol)./SpikeSummaryTable.(durCol) * 3600;


%% ======================= OUTPATIENT FILTER =======================
fprintf('\n=== BUILDING OUTPATIENT KEYS (site rule OR report_patient_class OR jay_in_or_out) ===\n');

% ---------- 1) Outpatient by acquired_on (Spe/Radnor) ----------
OutptBySite = table('Size',[0 2], ...
                     'VariableTypes',{'double','double'}, ...
                     'VariableNames',{'Patient','Session'});

isOutpt_site = (contains(acqStr_report,"spe",'IgnoreCase',true) | contains(acqStr_report,"radnor",'IgnoreCase',true));

if any(isOutpt_site)
    OutptBySite = unique(ReportTable(isOutpt_site, {'patient_id','session_number'}));
    OutptBySite.Properties.VariableNames = {'Patient','Session'};
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

% ---------- 5) Apply filter to spike summary + report: OUTPATIENT + ROUTINE ONLY ----------
isRoutine = isfinite(SpikeSummaryTable.(durCol)) & ...
            SpikeSummaryTable.(durCol) <= MAX_ROUTINE_HOURS * 3600;

RoutineKeys = unique(SpikeSummaryTable(isRoutine, {'Patient','Session'}));

% Intersect outpatient keys with routine keys
OutptRoutineKeys = innerjoin(OutptKeys, RoutineKeys, 'Keys', {'Patient','Session'});

% Now restrict both SpikeSummaryTable and ReportTable to outpatient + routine
SpikeSummaryTable = innerjoin(SpikeSummaryTable, OutptRoutineKeys, ...
    'Keys', {'Patient','Session'});

ReportTable = innerjoin(ReportTable, OutptRoutineKeys, ...
    'LeftKeys', {'patient_id','session_number'}, ...
    'RightKeys', {'Patient','Session'});

fprintf(['\n[Outpatient filter] Kept ONLY outpatient routine EEGs (Duration_sec <= 4h)\n' ...
         '  Outpatient defined as ANY of:\n' ...
         '    - acquired_on containing "Spe" or "Radnor" (case-insensitive) OR\n' ...
         '    - report_PATIENT_CLASS == "Outpatient" OR\n' ...
         '    - jay_in_or_out == "out"\n']);
fprintf('[Outpatient+routine filter] Kept %d spike rows and %d report rows.\n', ...
    height(SpikeSummaryTable), height(ReportTable));


%% ======================= BUILD PATIENT-LEVEL METRICS ONCE =======================
[PatientTypingAll, SzFreqPerPatient] = build_patient_metrics_from_report(ReportTable, canonical3);

%% ======================= RESTRICT TO STUDY TYPE (ROUTINE-ONLY VIEW) =======================
Views = build_filtered_view(SpikeSummaryTable, ReportTable, PatientTypingAll, SzFreqPerPatient, ...
                            NESD_LABEL, badTypes, canonical3);

%% ======================= PATIENT-LEVEL: SEX (EPILEPSY ONLY) SPIKES + SEIZURES =======================
% Two-panel figure:
%   A) Patient mean spike rate across EEGs (log10 spikes/hour)
%   B) Patient mean seizure frequency across visits (Rule 1) (log10 seizures/month)
%{
SEX_VAR = 'nlp_gender';  % expects 'F' / 'M'
figSex_patient_epi_out = '../figures/sex_patientLevel_epilepsyOnly_spikes_and_seizures.png';

% --- Inputs already in your pipeline ---
PL = Views.PatientLevelSpikeRates;   % one row per patient, includes MeanSpikeRate_perHour
Rk = Views.ReportForKeptSessions;    % report rows for kept sessions (has Patient, nlp_gender)
isEpiMask = Views.IsEpilepsyMask;

% --- 1) Epilepsy-only patient list + spike-rate table ---
epiPatients = PL.Patient(isEpiMask);
PL_epi = PL(isEpiMask, {'Patient','MeanSpikeRate_perHour'});

% --- 2) Collapse sex to ONE value per patient (error on conflicts) ---
sex_raw = upper(strtrim(string(Rk.(SEX_VAR))));
[grpP, pid] = findgroups(double(Rk.Patient));

SexPerPatient = strings(max(grpP),1);
for g = 1:max(grpP)
    thisPid = pid(g);

    % only consider epilepsy patients
    if ~ismember(thisPid, epiPatients)
        SexPerPatient(g) = "";
        continue;
    end

    x = sex_raw(grpP==g);
    x = x(~ismissing(x) & strlength(x)>0);

    if isempty(x)
        SexPerPatient(g) = "";
    else
        if any(x ~= x(1))
            error('Conflicting nlp_gender within epilepsy patient %g: %s', ...
                thisPid, strjoin(unique(x), ","));
        end
        SexPerPatient(g) = x(1);
    end
end

SexTbl = table(pid, SexPerPatient, 'VariableNames', {'Patient','SexCode'});

% --- 3) Join sex onto epilepsy-only patient table (spikes) ---
Jspk = innerjoin(PL_epi, SexTbl, 'Keys','Patient');

isF_spk = (Jspk.SexCode == "F");
isM_spk = (Jspk.SexCode == "M");
isU_spk = ~(isF_spk | isM_spk);

x_f_spk = Jspk.MeanSpikeRate_perHour(isF_spk);
x_m_spk = Jspk.MeanSpikeRate_perHour(isM_spk);

x_f_spk = x_f_spk(isfinite(x_f_spk));
x_m_spk = x_m_spk(isfinite(x_m_spk));

n_f_spk = numel(x_f_spk);
n_m_spk = numel(x_m_spk);
n_u_spk = nnz(isU_spk);

% --- 4) Patient-level seizure frequency (Rule 1), epilepsy-only, join sex ---
% Build unique-visit table (Rule 1 replacement already applied inside)
Vuniq_R1 = build_visit_level_freq_R1(ReportTable);  % Patient, VisitDate, Freq_R1

% Restrict to epilepsy patients and compute per-patient mean seizure frequency
Vepi = innerjoin(Vuniq_R1, table(epiPatients,'VariableNames',{'Patient'}), 'Keys','Patient');

[gv, pkeys] = findgroups(Vepi.Patient);
MeanSzFreq_R1 = splitapply(@(x) mean(x,'omitnan'), Vepi.Freq_R1, gv);

SzP = table(pkeys, MeanSzFreq_R1, 'VariableNames', {'Patient','MeanSzFreq_R1'});

% Attach sex
Jsz = innerjoin(SzP, SexTbl, 'Keys','Patient');

isF_sz = (Jsz.SexCode == "F");
isM_sz = (Jsz.SexCode == "M");
isU_sz = ~(isF_sz | isM_sz);

x_f_sz = Jsz.MeanSzFreq_R1(isF_sz);
x_m_sz = Jsz.MeanSzFreq_R1(isM_sz);

x_f_sz = x_f_sz(isfinite(x_f_sz));
x_m_sz = x_m_sz(isfinite(x_m_sz));

n_f_sz = numel(x_f_sz);
n_m_sz = numel(x_m_sz);
n_u_sz = nnz(isU_sz);

% --- 5) Stats (rank-sum + Cliff’s delta) ---
% Spikes
p_spk = NaN; d_spk = NaN; med_f_spk=NaN; med_m_spk=NaN; iqr_f_spk=[NaN NaN]; iqr_m_spk=[NaN NaN];
if n_f_spk >= 3 && n_m_spk >= 3
    p_spk = ranksum(x_f_spk, x_m_spk, 'method','approx');
    d_spk = cliff_delta(x_f_spk, x_m_spk);   % positive => Women > Men
    med_f_spk = median(x_f_spk,'omitnan'); iqr_f_spk = prctile(x_f_spk,[25,75]);
    med_m_spk = median(x_m_spk,'omitnan'); iqr_m_spk = prctile(x_m_spk,[25,75]);
else
    warning('Not enough epilepsy patients with sex F/M for spikes (F=%d, M=%d). Unknown/other=%d', ...
        n_f_spk, n_m_spk, n_u_spk);
end

% Seizures
p_sz = NaN; d_sz = NaN; med_f_sz=NaN; med_m_sz=NaN; iqr_f_sz=[NaN NaN]; iqr_m_sz=[NaN NaN];
if n_f_sz >= 3 && n_m_sz >= 3
    p_sz = ranksum(x_f_sz, x_m_sz, 'method','approx');
    d_sz = cliff_delta(x_f_sz, x_m_sz);      % positive => Women > Men
    med_f_sz = median(x_f_sz,'omitnan'); iqr_f_sz = prctile(x_f_sz,[25,75]);
    med_m_sz = median(x_m_sz,'omitnan'); iqr_m_sz = prctile(x_m_sz,[25,75]);
else
    warning('Not enough epilepsy patients with sex F/M for seizures (F=%d, M=%d). Unknown/other=%d', ...
        n_f_sz, n_m_sz, n_u_sz);
end

fprintf('\n=== Epilepsy patients: patient-level by sex ===\n');
fprintf('SPIKES (mean across EEGs):\n');
fprintf('  Women: N=%d, median=%.2f (%.2f–%.2f) spikes/hour\n', n_f_spk, med_f_spk, iqr_f_spk(1), iqr_f_spk(2));
fprintf('  Men:   N=%d, median=%.2f (%.2f–%.2f) spikes/hour\n', n_m_spk, med_m_spk, iqr_m_spk(1), iqr_m_spk(2));
fprintf('  ranksum p=%.3g, Cliff''s delta=%.2f (Women vs Men)\n', p_spk, d_spk);
fprintf('SEIZURES (mean across visits, Rule 1):\n');
fprintf('  Women: N=%d, median=%.2f (%.2f–%.2f) seizures/month\n', n_f_sz, med_f_sz, iqr_f_sz(1), iqr_f_sz(2));
fprintf('  Men:   N=%d, median=%.2f (%.2f–%.2f) seizures/month\n', n_m_sz, med_m_sz, iqr_m_sz(1), iqr_m_sz(2));
fprintf('  ranksum p=%.3g, Cliff''s delta=%.2f (Women vs Men)\n', p_sz, d_sz);

% --- 6) Plot: two subplots (A spikes, B seizures) ---
EPS_FREQ = 1e-3;
Y_SZ_ZERO = log10(EPS_FREQ);
Y_SZ_LIMS = spearman_xLims;  % reuse [-3.5 4]

fSex = figure('Color','w','Position',[80 80 1200 520]);
tiledlayout(fSex,1,2,'TileSpacing','compact','Padding','compact');

% A) Spikes
axA = nexttile; hold(axA,'on'); box(axA,'off'); grid(axA,'on');

Y_A = to_log10_per_hour([x_f_spk(:); x_m_spk(:)], EPS_RATE);
Y_A = add_y_jitter_eps(Y_A, Y_ZERO, Y_LIMS, 0.02);
G_A = [repmat("Women", n_f_spk, 1); repmat("Men", n_m_spk, 1)];

boxchart(axA, categorical(G_A), Y_A, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(axA, categorical(G_A), Y_A, 18, 'filled','MarkerFaceAlpha',0.25);

yline(axA, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axA, Y_LIMS);
ylabel(axA, 'log_{10}(spikes/hour)');
title(axA, 'A. Mean spike rate by sex');

if isfinite(p_spk)
    ySig = Y_LIMS(2) - 0.08 * range(Y_LIMS);
    if p_spk < 0.001, add_sigbar(axA, 1, 2, ySig, 'p < 0.001');
    else,            add_sigbar(axA, 1, 2, ySig, sprintf('p = %.3g', p_spk));
    end
end

labsA = string(axA.XTickLabel);
labsA(labsA=="Women") = sprintf('Women (N=%d)', n_f_spk);
labsA(labsA=="Men")   = sprintf('Men (N=%d)',   n_m_spk);
axA.XTickLabel = labsA;
axA.XTickLabelRotation = 20;
set(axA,'FontSize',20);

% B) Seizures (Rule 1)
axB = nexttile; hold(axB,'on'); box(axB,'off'); grid(axB,'on');

Y_B = [to_log10_per_month(x_f_sz(:), EPS_FREQ); to_log10_per_month(x_m_sz(:), EPS_FREQ)];
Y_B = add_y_jitter_eps(Y_B, Y_SZ_ZERO, Y_SZ_LIMS, 0.02);
G_B = [repmat("Women", n_f_sz, 1); repmat("Men", n_m_sz, 1)];

boxchart(axB, categorical(G_B), Y_B, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(axB, categorical(G_B), Y_B, 18, 'filled','MarkerFaceAlpha',0.25);

yline(axB, Y_SZ_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axB, Y_SZ_LIMS);
ylabel(axB, 'log_{10}(seizures/month)');
title(axB, 'B. Mean seizure frequency by sex');

if isfinite(p_sz)
    ySig = Y_SZ_LIMS(2) - 0.08 * range(Y_SZ_LIMS);
    if p_sz < 0.001, add_sigbar(axB, 1, 2, ySig, 'p < 0.001');
    else,           add_sigbar(axB, 1, 2, ySig, sprintf('p = %.3g', p_sz));
    end
end

labsB = string(axB.XTickLabel);
labsB(labsB=="Women") = sprintf('Women (N=%d)', n_f_sz);
labsB(labsB=="Men")   = sprintf('Men (N=%d)',   n_m_sz);
axB.XTickLabel = labsB;
axB.XTickLabelRotation = 20;
set(axB,'FontSize',20);

% Save
if ~exist(fileparts(figSex_patient_epi_out),'dir')
    mkdir(fileparts(figSex_patient_epi_out));
end
exportgraphics(fSex, figSex_patient_epi_out, 'Resolution', 300);
fprintf('Saved epilepsy-only patient-level sex (spikes+seizures) figure: %s\n', figSex_patient_epi_out);
%}



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
% This first chunk of code is to figure out if the report listed spikes,
% reconciling two different types of reports (main and jay)
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

% Determine if spikes are present according to jay database
presentJay_any  = isF_P | isM_P | isG_P;           % any jay present
allJay_absent   = isF_A & isM_A & isG_A;          % all three explicitly "absent"
allJay_present  = isF_P & isM_P & isG_P;          % all three explicitly "present"

blankMain   = ~(isMainPresent | isMainAbsent);    % nothing interpretable in main col
blankJayAll = ~(isF_P | isF_A) & ~(isM_P | isM_A) & ~(isG_P | isG_A);  % all jay blank/uninterpretable

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



JoinA = innerjoin(SessionLevelSpikeRates(:,{'Patient','Session','SpikesPerHour'}), ...
                  ReportSlim, 'Keys', {'Patient','Session'});

x_abs = JoinA.SpikesPerHour(JoinA.ReportStatus=="absent");
x_pre = JoinA.SpikesPerHour(JoinA.ReportStatus=="present");
p_rankSum_A = ranksum(x_abs, x_pre, 'method','approx');
m_abs  = median(x_abs,'omitnan');   iqr_abs  = prctile(x_abs,[25,75]);
m_pre  = median(x_pre,'omitnan');   iqr_pre  = prctile(x_pre,[25,75]);
effectA_cliff = cliff_delta(x_pre, x_abs);


Y_A = to_log10_per_hour([x_abs(:); x_pre(:)], EPS_RATE);
Y_A = add_y_jitter_eps(Y_A, Y_ZERO, Y_LIMS, 0.02);  % ~±1% of y-range
G_A = [repmat("Absent", numel(x_abs), 1); repmat("Present", numel(x_pre), 1)];


% ---- Panel B calculations: Epilepsy (any) vs NESD (patient-level means) ----
x_ep  = PatientLevelSpikeRates.MeanSpikeRate_perHour(isEpilepsyMask); % combined gen + focal, focal, gen, unclassified (leaves out uncertain, unknown, PNEE)
x_nes = PatientLevelSpikeRates.MeanSpikeRate_perHour(isNESDMask); % just pnee
n_ep  = nnz(isfinite(x_ep)); % how many finite numerical values
n_nes = nnz(isfinite(x_nes));
p_rankSum_B = ranksum(x_ep, x_nes, 'method','approx'); % rank sum test
m_ep  = median(x_ep,'omitnan');   iqr_ep  = prctile(x_ep,[25,75]); % median and iqr
m_nes = median(x_nes,'omitnan');  iqr_nes = prctile(x_nes,[25,75]);
effectB_cliff = cliff_delta(x_ep, x_nes);
Y_B = [to_log10_per_hour(x_ep(:), EPS_RATE); to_log10_per_hour(x_nes(:), EPS_RATE)];
Y_B = add_y_jitter_eps(Y_B, Y_ZERO, Y_LIMS, 0.02);
G_B = [repmat("Epilepsy", n_ep, 1); repmat("NESD", n_nes, 1)];


% ---- Panel C calculations: General vs Temporal vs Frontal (patient-level) ----
Y_C = to_log10_per_hour(SubtypeSubsetTable.MeanSpikeRate_perHour, EPS_RATE);
Y_C = add_y_jitter_eps(Y_C, Y_ZERO, Y_LIMS, 0.02);

% ---- Panel C omnibus test: Kruskal–Wallis ----
groupsC = SubtypeSubsetTable.EpiType4;   % categorical with levels: General, Temporal, Frontal
[p_kw_C, tbl_kw_C, stats_kw_C] = kruskalwallis( ...
        SubtypeSubsetTable.MeanSpikeRate_perHour, groupsC, 'off');

SS_total = tbl_kw_C{end, 2};
SS_group = tbl_kw_C{2, 2};
eta2_kw_C = SS_group / SS_total;

% ---- Draw figure ----
f1 = figure('Color','w','Position',[60 60 1500 520]);
tiledlayout(f1,1,3,'TileSpacing','compact','Padding','loose');

% A) Report Present vs Absent
axA = nexttile; hold(axA,'on'); box(axA,'off'); grid(axA,'on');
boxchart(axA, categorical(G_A), Y_A, 'BoxFaceAlpha',0.25,'MarkerStyle', 'none');
swarmchart(axA, categorical(G_A), Y_A, 18, 'filled','MarkerFaceAlpha',0.25);    
yline(axA, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axA, Y_LIMS);
ylabel(axA, 'log_{10}(spikes/hour)');
title(axA, 'A. Reported presence or absence of spikes');
if p_rankSum_A < 0.001
    add_sigbar(axA, 1, 2, Y_LIMS(2) - 0.08*range(Y_LIMS), sprintf('p < 0.001'));
else
    add_sigbar(axA, 1, 2, Y_LIMS(2) - 0.08*range(Y_LIMS), sprintf('p = %.3g', p_rankSum_A));
end

set(axA,'FontSize',20);
labelsA = string(axA.XTickLabel);  
labelsA(labelsA=="Absent") = sprintf('Absent (N=%d)', nnz(isfinite(x_abs)));
labelsA(labelsA=="Present")     = sprintf('Present (N=%d)', nnz(isfinite(x_pre)));  
axA.XTickLabel = labelsA;
axA.XTickLabelRotation = 20;   

% B) Epilepsy vs NESD
axB = nexttile; hold(axB,'on'); box(axB,'off'); grid(axB,'on');
boxchart(axB, categorical(G_B), Y_B, 'BoxFaceAlpha',0.25,'MarkerStyle', 'none');
swarmchart(axB, categorical(G_B), Y_B, 18, 'filled','MarkerFaceAlpha',0.25);    
yline(axB, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axB, Y_LIMS);
ylabel(axB, 'log_{10}(spikes/hour)');
title(axB, 'B. Epilepsy versus NESD');
if p_rankSum_B < 0.001
    add_sigbar(axB, 1, 2, Y_LIMS(2) - 0.08*range(Y_LIMS), sprintf('p < 0.001'));
else
    add_sigbar(axB, 1, 2, Y_LIMS(2) - 0.08*range(Y_LIMS), sprintf('p = %.3g', p_rankSum_B));
end
set(axB,'FontSize',20);


labelsB = string(axB.XTickLabel); 
labelsB(labelsB=="Epilepsy") = sprintf('Epilepsy (N=%d)', n_ep);
labelsB(labelsB=="NESD")     = sprintf('NESD (N=%d)', n_nes);  
axB.XTickLabel = labelsB;
axB.XTickLabelRotation = 20;  


% C) General vs Temporal vs Frontal
axC = nexttile; hold(axC,'on'); box(axC,'off'); grid(axC,'on');
boxchart(axC, SubtypeSubsetTable.EpiType4, Y_C, 'BoxFaceAlpha',0.25,'MarkerStyle', 'none');
swarmchart(axC, SubtypeSubsetTable.EpiType4, Y_C, 18, 'filled','MarkerFaceAlpha',0.25);   
yline(axC, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axC, Y_LIMS);
ylabel(axC, 'log_{10}(spikes/hour)');
title(axC, 'C. Epilepsy subtype');
set(axC,'FontSize',20);


set(axC,'FontSize',20);
labelsC = string(axC.XTickLabel);  

for i = 1:height(SubtypeStatsTable)
    lab   = string(SubtypeStatsTable.EpiType4(i));
    nHere = SubtypeStatsTable.GroupCount(i);
    labelsC(labelsC==lab) = sprintf('%s (N=%d)', lab, nHere);
end

axC.XTickLabel = labelsC;
axC.XTickLabelRotation = 20;

% Pairwise bars with Bonferroni-coded labels
yTop  = Y_LIMS(2);
yStep = 0.09 * range(Y_LIMS);
y0    = yTop - 0.09 * range(Y_LIMS);
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



%% ======================= FIGURE 2: SPEARMAN (log–log, fixed axes) =======================
PatientSpikeSz_All   = Views.PatientSpikeSz_All;
PatientSpikeSz_Typed = Views.PatientSpikeSz_Typed;

% Figure 2
[SpearmanResults_main, rs_all_main, p_all_main, n_all_main] = ...
    spearman_plotting_function(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
        canonical3, spearman_xLims, spearman_yLims, ...
        fig2_out, ...
        'MeanSzFreq', '', false);

% Supplemental figure 1
[SpearmanResults_nz, rs_all_nz, p_all_nz, n_all_nz] = ...
    spearman_plotting_function(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
        canonical3, spearman_xLims, spearman_yLims, ...
        figS1_out, ...
        'MeanSzFreq', ' (non-zero spikes and seizures)', true);


%% ======================= FIGURE S2: SZ FREQ BY REPORTED SPIKES ACROSS EEGs =======================
% Goal:
%   Compare mean seizure frequency across visits (Rule 1 only, MeanSzFreq)
%   between:
%     - Group 1: patients for whom ALL EEGs with reports are reported to have NO spikes
%     - Group 2: patients for whom AT LEAST ONE EEG with a report is reported to have spikes.
%
%   Analysis: Wilcoxon rank-sum test
%   Plot: log10(seizures per month) for the two groups (Fig. S3)


% Use the same ReportSlim logic as Control Panel A (only rows with resolved present/absent)
% ReportSlim has columns: Patient, Session, ReportStatus (categorical "absent"/"present")
% If you moved the creation of ReportSlim, ensure it exists before this block.

RS = ReportSlim;  % shorthand

% Group by patient: determine which patients ever have "present" vs only "absent"
[gpRS, pidRS] = findgroups(RS.Patient);

hasPresent = splitapply(@(x) any(x == "present"), RS.ReportStatus, gpRS);
hasAbsent  = splitapply(@(x) any(x == "absent"),  RS.ReportStatus, gpRS);

% Build a per-patient table for grouping
ReportPerPatient = table(pidRS, hasPresent, hasAbsent, ...
    'VariableNames', {'Patient','HasPresent','HasAbsent'});

% Find epilepsy patients
EpPatients = Views.PatientLevelSpikeRates.Patient(Views.IsEpilepsyMask);
EpTable    = table(EpPatients, 'VariableNames', {'Patient'});

% Join with seizure frequencies (Rule 1 only: MeanSzFreq)
% SzFreqPerPatient has columns: Patient, MeanSzFreq, MeanSzFreq_rule12, MeanSzFreq_raw
S2 = innerjoin(SzFreqPerPatient, ReportPerPatient, 'Keys','Patient');

% Restrict to epilepsy-only
S2 = innerjoin(S2, EpTable, 'Keys','Patient');

% Define groups:
%   Group 1: ALL EEGs with reports have "absent" AND none have "present"
%   Group 2: AT LEAST ONE EEG with "present" (regardless of whether some are "absent")
mask_allAbsent  =  S2.HasAbsent & ~S2.HasPresent;
mask_anyPresent =  S2.HasPresent;

freq_allAbsent  = S2.MeanSzFreq(mask_allAbsent);
freq_anyPresent = S2.MeanSzFreq(mask_anyPresent);

% Drop NaNs (patients with no usable seizure frequency under Rule 1)
freq_allAbsent  = freq_allAbsent(isfinite(freq_allAbsent));
freq_anyPresent = freq_anyPresent(isfinite(freq_anyPresent));

n_allAbsent  = numel(freq_allAbsent);
n_anyPresent = numel(freq_anyPresent);

% Wilcoxon rank-sum test
p_rankSum_S2 = NaN;
if n_allAbsent >= 1 && n_anyPresent >= 1
    p_rankSum_S2 = ranksum(freq_allAbsent, freq_anyPresent, 'method','approx');
end

% Medians, IQRs, and effect size (Cliff's delta)
m_allAbsent   = median(freq_allAbsent,'omitnan');
iqr_allAbsent = prctile(freq_allAbsent,[25,75]);

m_anyPresent   = median(freq_anyPresent,'omitnan');
iqr_anyPresent = prctile(freq_anyPresent,[25,75]);

effectS2_cliff = cliff_delta(freq_anyPresent, freq_allAbsent);


% ----- Plot log10(seizures/month) with eps floor and jitter -----
% Use a small epsilon for seizure frequency, similar spirit to Spearman plotting
EPS_FREQ = 1e-3;   % 0.01 seizures/month floor
Y_S2_ZERO = log10(EPS_FREQ);
Y_S2_LIMS = spearman_xLims;   % reuse seizure-frequency axis limits

Y_S2 = [to_log10_per_month(freq_allAbsent,  EPS_FREQ); ...
        to_log10_per_month(freq_anyPresent, EPS_FREQ)];
Y_S2 = add_y_jitter_eps(Y_S2, Y_S2_ZERO, Y_S2_LIMS, 0.02);  % ~±1% of y-range

G_S2 = [repmat("All EEGs: no spikes",   n_allAbsent, 1); ...
        repmat("≥1 EEG: spikes present", n_anyPresent, 1)];

% Draw figure
fS2 = figure('Color','w','Position',[100 100 700 500]);
axS2 = axes(fS2); hold(axS2,'on'); box(axS2,'off'); grid(axS2,'on');

boxchart(axS2, categorical(G_S2), Y_S2, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(axS2, categorical(G_S2), Y_S2, 18, 'filled','MarkerFaceAlpha',0.25);

yline(axS2, Y_S2_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axS2, Y_S2_LIMS);
ylabel(axS2, 'log_{10}(seizures/month)');
title(axS2, 'Fig. S2. Mean seizure frequency by reported spikes across EEGs');

% Significance bar
if ~isnan(p_rankSum_S2)
    yTop  = Y_S2_LIMS(2);
    ySig  = yTop - 0.08 * range(Y_S2_LIMS);
    if p_rankSum_S2 < 0.001
        pLabel = 'p < 0.001';
    else
        pLabel = sprintf('p = %.3g', p_rankSum_S2);
    end
    add_sigbar(axS2, 1, 2, ySig, pLabel);
end

set(axS2,'FontSize',20);

% Custom x tick labels with N
labelsS2 = string(axS2.XTickLabel);
labelsS2(labelsS2=="All EEGs: no spikes") = ...
    sprintf('All EEGs: no spikes (N=%d)', n_allAbsent);
labelsS2(labelsS2=="≥1 EEG: spikes present") = ...
    sprintf('\x2265 1 EEG: spikes present (N=%d)', n_anyPresent);  % ≥ if it renders; otherwise ">= 1"
axS2.XTickLabel = labelsS2;
axS2.XTickLabelRotation = 15;

% Save figure
if ~exist(fileparts(figS2_out),'dir'), mkdir(fileparts(figS2_out)); end
exportgraphics(fS2, figS2_out, 'Resolution', 300);
fprintf('Saved Fig S2: %s\n', figS2_out);

%% ======================= TABLE 1: COHORT CHARACTERISTICS =======================
% Build Table 1 for all patients with >=1 outpatient routine EEG

BIRTH_VAR = 'deid_birth_date';   % de-identified birth date (M/DD/YY)
SEX_VAR   = 'nlp_gender';        % 'F' / 'M' from NLP

Rk = Views.ReportForKeptSessions;        % outpatient + routine + kept sessions
PL = Views.PatientLevelSpikeRates;       % per-patient mean spike rate + typing
SzP = SzFreqPerPatient;                  % per-patient mean seizure freq (Rule 1)


% ---------- 1. Base cohort: one row per patient ----------
AllPatients = PL.Patient;           % these are all patients with >=1 outpatient routine EEG
N_total     = numel(AllPatients);

% ---------- 2. Age at first visit (from de-identified birth date) ----------
% Age at first visit is defined as the number of years between
% deid_birth_date and 01-Jan-2000 (the de-identified reference date).

birth_str = strtrim(string(Rk.(BIRTH_VAR)));

% Treat empty/null-like entries as missing
isMissingBirth = birth_str == "" | birth_str == "null" | birth_str == "[null]";

birth_dt = NaT(size(birth_str));
birth_dt(~isMissingBirth) = datetime( ...
    birth_str(~isMissingBirth), ...
    'InputFormat','yyyy-MM-dd');  

refDate = datetime(2000,1,1);   % first-visit reference date in de-identified time

age_first_vec = NaN(size(birth_dt));
validBirth = ~isnat(birth_dt);
age_first_vec(validBirth) = days(refDate - birth_dt(validBirth)) / 365.25;

% Collapse to one age per patient (should all be the same within a patient)
[grpAge, pidAge] = findgroups(Rk.Patient);
AgeFirst = splitapply(@local_min_omitnan, age_first_vec, grpAge);

AgeTable = table(pidAge, AgeFirst, ...
    'VariableNames', {'Patient','AgeFirst'});


% Join to main patient list, preserving PL.Patient order
AgeTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), AgeTable, ...
                     'Keys','Patient', 'RightVariables','AgeFirst');

age_vec = AgeTable.AgeFirst;
age_med = median(age_vec,'omitnan');
age_q   = prctile(age_vec,[25,75]);

% ---------- 3. Sex (Women / Men / Unknown/Other) ----------
% Use nlp_gender: "F" / "M" per report row, collapsed to 1 value per patient

sex_raw = upper(strtrim(string(Rk.(SEX_VAR))));  % expect "F", "M", or missing

[grpSex, pidSex] = findgroups(Rk.Patient);
sex_per_patient = splitapply(@(s) local_first_nonmissing(s), sex_raw, grpSex);

SexTable = table(pidSex, sex_per_patient, ...
    'VariableNames', {'Patient','SexCode'});
SexTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), SexTable, ...
                     'Keys','Patient','RightVariables','SexCode');

isFemale    = SexTable.SexCode == "F";
isMale      = SexTable.SexCode == "M";
isUnknownSex = ~(isFemale | isMale);

n_female  = nnz(isFemale);
n_male    = nnz(isMale);
n_sex_unk = nnz(isUnknownSex);


% ---------- 4. Epilepsy diagnosis (Epilepsy / PNES / Unknown) ----------
IsEpilepsyMask = Views.IsEpilepsyMask;   % same length as PL
IsNESDMask     = Views.IsNESDMask;

n_epi   = nnz(IsEpilepsyMask);
n_pnes  = nnz(IsNESDMask);
n_diag_unknown = N_total - n_epi - n_pnes;

% ---------- 5. Epilepsy subtype ----------
% Temporal / Frontal / Generalized / Other / Unknown
E4 = string(PL.EpiType4);       % "General","Temporal","Frontal","" etc.
isTemp   = (E4 == "Temporal")   & IsEpilepsyMask;
isFront  = (E4 == "Frontal")    & IsEpilepsyMask;
isGen    = (E4 == "General")    & IsEpilepsyMask;

% "Unknown" subtype 
typing = PL.EpilepsySpecific;
isSubtypeUnknown = (typing == "Unknown or MRN not found" | ...
    typing == "Unclassified or Unspecified") & ~isTemp & ~isFront & ~isGen & isEpilepsyMask;

% "Other" = epilepsy, but not Temporal/Frontal/Generalized, and has some epilepsy label
isOtherSubtype = IsEpilepsyMask & ~isSubtypeUnknown & ~(isTemp | isFront | isGen);

n_temp    = nnz(isTemp);
n_front   = nnz(isFront);
n_gen     = nnz(isGen);
n_other   = nnz(isOtherSubtype);
n_sub_unk = nnz(isSubtypeUnknown);

assert(n_temp+n_front+n_gen+n_other+n_sub_unk == sum(IsEpilepsyMask==1))

% ---------- 6. Number of clinic visits per patient (median, IQR) ----------
% Use visit-level table from build_visit_level_freq_R1 (Rule 1 only)
Vuniq_all = build_visit_level_freq_R1(ReportTable);  % all patients with outpatient routine EEGs

[grpV, pidV] = findgroups(Vuniq_all.Patient);
nVisits = splitapply(@(x) numel(unique(x)), Vuniq_all.VisitDate, grpV);

VisitsTable = table(pidV, nVisits, 'VariableNames', {'Patient','NumVisits'});
VisitsTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), VisitsTable, ...
                        'Keys','Patient','RightVariables','NumVisits');

vis_vec = VisitsTable.NumVisits;
vis_med = median(vis_vec,'omitnan');
vis_q   = prctile(vis_vec,[25,75]);

% ---------- 7. Number of EEGs per patient (median, IQR) ----------
Sess = Views.SessionsForFigures;   % outpatient + routine sessions
[grpS, pidS] = findgroups(Sess.Patient);
nEEG = splitapply(@(x) numel(unique(x)), Sess.Session, grpS);

EEGTable = table(pidS, nEEG, 'VariableNames', {'Patient','NumEEG'});
EEGTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), EEGTable, ...
                     'Keys','Patient','RightVariables','NumEEG');

eeg_vec = EEGTable.NumEEG;
eeg_med = median(eeg_vec,'omitnan');
eeg_q   = prctile(eeg_vec,[25,75]);

% ---------- 8. Mean seizure frequency per patient (median, IQR) ----------
SzP_join = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), SzP, 'Keys','Patient');
sf_vec   = SzP_join.MeanSzFreq;
sf_vec   = sf_vec(isfinite(sf_vec));
sf_med   = median(sf_vec,'omitnan');
sf_q     = prctile(sf_vec,[25,75]);

% ---------- 9. Mean spike rate per patient (median, IQR) ----------
% Use MeanSpikeRate_perHour (spikes/hour)
sr_vec = PL.MeanSpikeRate_perHour;
sr_vec = sr_vec(isfinite(sr_vec));
sr_med = median(sr_vec,'omitnan');
sr_q   = prctile(sr_vec,[25,75]);

% ---------- 10. Reported spikes per EEG (Present / Absent / Unknown) ----------
% Use ReportSlim and same logic as Fig S3 but per patient
% Start from all outpatient routine sessions with spike info logic applied
RS_all = ReportForKeptSessions(:, {'Patient','Session','ReportStatus'});

% Turn missing into an explicit "unknown" category
status_str = string(RS_all.ReportStatus);
status_str(ismissing(RS_all.ReportStatus)) = "unknown";

RS_all.SpikeStatus = categorical(status_str, ["absent","present","unknown"]);

n_rep_present = nnz(RS_all.SpikeStatus == "present");
n_rep_absent  = nnz(RS_all.SpikeStatus == "absent");
n_rep_unknown = nnz(RS_all.SpikeStatus == "unknown");
n_eegs_all = height(RS_all);

assert(n_rep_present + n_rep_absent + n_rep_unknown == height(RS_all));

% ---------- 11. Build Table 1 ----------
Var   = {};
Cat   = {};
Ncol  = [];
Pcol  = [];
Summ  = {};

% 1) Total N
Var{end+1}  = 'Total N patients with ≥1 outpatient routine EEG';
Cat{end+1}  = '';
Ncol(end+1) = N_total;
Pcol(end+1) = 100;
Summ{end+1} = "";

% 2) Age at first visit
ageStr = sprintf('%.1f (%.1f-%.1f)', age_med, age_q(1), age_q(2));
Var{end+1}  = 'Age at first visit (years)';
Cat{end+1}  = '';
Ncol(end+1) = NaN;
Pcol(end+1) = NaN;
Summ{end+1} = ageStr;

% 3) Sex
Var{end+1}  = 'Sex';
Cat{end+1}  = 'Women';
Ncol(end+1) = n_female;
Pcol(end+1) = 100 * n_female / N_total;
Summ{end+1} = "";

Var{end+1}  = 'Sex';
Cat{end+1}  = 'Men';
Ncol(end+1) = n_male;
Pcol(end+1) = 100 * n_male / N_total;
Summ{end+1} = "";

Var{end+1}  = 'Sex';
Cat{end+1}  = 'Unknown/Other';
Ncol(end+1) = n_sex_unk;
Pcol(end+1) = 100 * n_sex_unk / N_total;
Summ{end+1} = "";

% 4) Epilepsy diagnosis
Var{end+1}  = 'Epilepsy diagnosis';
Cat{end+1}  = 'Epilepsy';
Ncol(end+1) = n_epi;
Pcol(end+1) = 100 * n_epi / N_total;
Summ{end+1} = "";

Var{end+1}  = 'Epilepsy diagnosis';
Cat{end+1}  = 'PNES';
Ncol(end+1) = n_pnes;
Pcol(end+1) = 100 * n_pnes / N_total;
Summ{end+1} = "";

Var{end+1}  = 'Epilepsy diagnosis';
Cat{end+1}  = 'Unknown';
Ncol(end+1) = n_diag_unknown;
Pcol(end+1) = 100 * n_diag_unknown / N_total;
Summ{end+1} = "";

% 5) Epilepsy subtype
Var{end+1}  = 'Epilepsy subtype';
Cat{end+1}  = 'Temporal lobe';
Ncol(end+1) = n_temp;
Pcol(end+1) = 100 * n_temp / n_epi;
Summ{end+1} = "";

Var{end+1}  = 'Epilepsy subtype';
Cat{end+1}  = 'Frontal lobe';
Ncol(end+1) = n_front;
Pcol(end+1) = 100 * n_front / n_epi;
Summ{end+1} = "";

Var{end+1}  = 'Epilepsy subtype';
Cat{end+1}  = 'Generalized';
Ncol(end+1) = n_gen;
Pcol(end+1) = 100 * n_gen / n_epi;
Summ{end+1} = "";

Var{end+1}  = 'Epilepsy subtype';
Cat{end+1}  = 'Other';
Ncol(end+1) = n_other;
Pcol(end+1) = 100 * n_other / n_epi;
Summ{end+1} = "";

Var{end+1}  = 'Epilepsy subtype';
Cat{end+1}  = 'Unknown';
Ncol(end+1) = n_sub_unk;
Pcol(end+1) = 100 * n_sub_unk / n_epi;
Summ{end+1} = "";

% 6) Number of clinic visits
visStr = sprintf('%.1f (%.1f-%.1f)', vis_med, vis_q(1), vis_q(2));
Var{end+1}  = 'Number of clinic visits';
Cat{end+1}  = '';
Ncol(end+1) = NaN;
Pcol(end+1) = NaN;
Summ{end+1} = visStr;

% 7) Number of EEGs
eegStr = sprintf('%.1f (%1.1f-%1.1f)', eeg_med, eeg_q(1), eeg_q(2));
Var{end+1}  = 'Number of EEGs';
Cat{end+1}  = '';
Ncol(end+1) = NaN;
Pcol(end+1) = NaN;
Summ{end+1} = eegStr;

% 8) Mean seizure frequency
sfStr = sprintf('%.2f (%.2f-%.2f)', sf_med, sf_q(1), sf_q(2));
Var{end+1}  = 'Mean seizure frequency (seizures/month)';
Cat{end+1}  = '';
Ncol(end+1) = NaN;
Pcol(end+1) = NaN;
Summ{end+1} = sfStr;

% 9) Mean spike rate
srStr = sprintf('%.2f (%.2f-%.2f)', sr_med, sr_q(1), sr_q(2));
Var{end+1}  = 'Mean spike rate (spikes/hour)';
Cat{end+1}  = '';
Ncol(end+1) = NaN;
Pcol(end+1) = NaN;
Summ{end+1} = srStr;

% 10) Reported spikes
Var{end+1}  = 'Reported spikes';
Cat{end+1}  = 'Present';
Ncol(end+1) = n_rep_present;
Pcol(end+1) = 100 * n_rep_present / n_eegs_all;
Summ{end+1} = "";

Var{end+1}  = 'Reported spikes';
Cat{end+1}  = 'Absent';
Ncol(end+1) = n_rep_absent;
Pcol(end+1) = 100 * n_rep_absent / n_eegs_all;
Summ{end+1} = "";

Var{end+1}  = 'Reported spikes';
Cat{end+1}  = 'Unknown';
Ncol(end+1) = n_rep_unknown;
Pcol(end+1) = 100 * n_rep_unknown / n_eegs_all;
Summ{end+1} = "";

Table1 = table(string(Var(:)), string(Cat(:)), Ncol(:), Pcol(:), string(Summ(:)), ...
    'VariableNames', {'Variable','Category','N','Percent','Summary'});

Table1.Percent = round(Table1.Percent, 1);

% ====== FLATTEN Table1 TO TWO-COLUMN DISPLAY TABLE (WITH HEADINGS) ======

% Categorical variables that should get a heading row with blank stats
catVars = [
    "Sex"
    "Epilepsy diagnosis"
    "Epilepsy subtype"
    "Reported spikes"
];

% We will iterate by variable name in the order they appear
uVars = unique(string(Table1.Variable), 'stable');

OutVar  = {};
OutStat = {};

for v = 1:numel(uVars)
    varName = uVars(v);
    maskVar = string(Table1.Variable) == varName;

    % Extract all rows for this variable
    Tsub = Table1(maskVar, :);

    % Is this variable treated as categorical?
    if ismember(varName, catVars)

        % ---- Heading row ----
        OutVar{end+1,1}  = varName;
        OutStat{end+1,1} = "";

        % ---- Subcategory rows ----
        % Keep only rows with a non-empty Category
        catMask = strlength(strtrim(string(Tsub.Category))) > 0;
        Tcats   = Tsub(catMask, :);

        for j = 1:height(Tcats)
            catName = string(Tcats.Category(j));
            Ni      = Tcats.N(j);
            Pi      = Tcats.Percent(j);

            OutVar{end+1,1} = "    " + catName;

            if isfinite(Ni) && isfinite(Pi)
                OutStat{end+1,1} = sprintf('%d (%.1f%%)', Ni, Pi);
            elseif isfinite(Ni)
                OutStat{end+1,1} = sprintf('%d', Ni);
            else
                OutStat{end+1,1} = "";
            end
        end

    else
        % ---- Continuous (or single-row) variable ----
        % Take the first row for this variable
        row = Tsub(1,:);

        varLabel = varName;
        summ     = string(row.Summary);
        Ni       = row.N;
        Pi       = row.Percent;

        OutVar{end+1,1} = varLabel;

        if strlength(strtrim(summ)) > 0
            % e.g., "45.0 (32.0–58.0)"
            OutStat{end+1,1} = summ;
        elseif isfinite(Ni) && ~isfinite(Pi)
            % e.g., total N patients
            OutStat{end+1,1} = sprintf('%d', Ni);
        elseif isfinite(Ni) && isfinite(Pi)
            % fallback: N (%)
            OutStat{end+1,1} = sprintf('%d (%.1f%%)', Ni, Pi);
        else
            OutStat{end+1,1} = "";
        end
    end
end

Table1_flat = table(string(OutVar), string(OutStat), ...
    'VariableNames', {'Variable','Statistic'});
writetable(Table1_flat,Table1Csv);


%% ======================= PAIRED ANALYSIS =======================

%% -------- PARAMETERS FOR VISIT–EEG MATCHING --------
maxDaysBefore = 30;    % visit cannot be more than 30 days BEFORE EEG
maxDaysAfter  = 30;    % visit cannot be more than 30 days AFTER EEG
min_abs_diff_spikes = 0; % 0 means that high spike eeg must have more spikes than low spike eeg

%% -------- 1) Build visit-level frequency table (Rule 1 only) --------
% Helper function should return:
%   Vuniq_R1: Patient (double), VisitDate (datetime), Freq_R1 (double)
Vuniq_R1 = build_visit_level_freq_R1(ReportTable);

% Restrict to epilepsy patients based on Views
EpPatients = Views.PatientLevelSpikeRates.Patient(Views.IsEpilepsyMask);
EpPatients = unique(EpPatients);
EpTable    = table(EpPatients, 'VariableNames', {'Patient'});

Vuniq_R1 = innerjoin(Vuniq_R1, EpTable, 'Keys', 'Patient');

fprintf('[Spike-first (closest-visit) analysis] Visit-level Rule 1 table: %d rows for epilepsy patients.\n', ...
    height(Vuniq_R1));

%% -------- 2) Build session-level spike rates + EEG dates --------
% Session-level spike rates (outpatient + routine only) from Views:
SessRates = Views.SessionLevelSpikeRates;   % Patient, Session, SpikesPerHour, SpikesPerMin

% EEG start times live in the report table for the kept sessions:
Rk      = Views.ReportForKeptSessions;      % has Patient, Session, start_time_deid, etc.
EEG_raw = Rk.start_time_deid;
if isdatetime(EEG_raw)
    EEG_dt = EEG_raw;
else
    EEG_dt = datetime(strtrim(string(EEG_raw)), ...
                      'InputFormat','yyyy-MM-dd''T''HH:mm:ss');
end

SessionDates = table(Rk.Patient, Rk.Session, EEG_dt, ...
    'VariableNames', {'Patient','Session','EEG_Date'});

% Join spike rates + EEG dates (still outpatient + routine)
SessWithDate = innerjoin(SessRates, SessionDates, 'Keys', {'Patient','Session'});

% Restrict to epilepsy patients only
SessWithDate = innerjoin(SessWithDate, EpTable, 'Keys', 'Patient');

fprintf('[Spike-first (closest-visit) analysis] Session-level EEG table: %d epilepsy sessions.\n', ...
    height(SessWithDate));

%% -------- 3) Build all EEG–visit matches within the allowed window --------
Pairs = table('Size',[0 7], ...
    'VariableTypes', {'double','double','datetime','datetime','double','double','double'}, ...
    'VariableNames', {'Patient','Session','EEG_Date','VisitDate','GapDays', ...
                      'SpikeRate_perHour','Freq_R1'});

for pid = EpPatients'
    % All outpatient routine EEG sessions for this patient
    s_mask = (SessWithDate.Patient == pid);
    if ~any(s_mask), continue; end
    S_sub = SessWithDate(s_mask, :);

    % All allowable clinic visits for this patient (already filtered upstream)
    v_mask = (Vuniq_R1.Patient == pid);
    if ~any(v_mask), continue; end
    V_sub = Vuniq_R1(v_mask, :);

    for j = 1:height(S_sub)
        EEGd = S_sub.EEG_Date(j);

        % Positive gapDays = EEG after visit; negative = EEG before visit
        gapDays = days(EEGd - V_sub.VisitDate);
        keep    = (gapDays >= -maxDaysBefore) & (gapDays <= maxDaysAfter);

        if ~any(keep), continue; end

        tmp = table( ...
            repmat(pid, nnz(keep), 1), ...
            repmat(S_sub.Session(j), nnz(keep), 1), ...
            repmat(EEGd, nnz(keep), 1), ...
            V_sub.VisitDate(keep), ...
            gapDays(keep), ...
            repmat(S_sub.SpikesPerHour(j), nnz(keep), 1), ...
            V_sub.Freq_R1(keep), ...
            'VariableNames', Pairs.Properties.VariableNames);

        Pairs = [Pairs; tmp]; %#ok<AGROW>
    end
end

fprintf('\n[Spike-first (closest-visit) analysis] Built %d EEG–visit matches across %d epilepsy patients.\n', ...
    height(Pairs), numel(unique(Pairs.Patient)));

% Drop rows with non-finite frequencies or spike rates
Pairs = Pairs(isfinite(Pairs.Freq_R1) & isfinite(Pairs.SpikeRate_perHour), :);

if isempty(Pairs)
    warning('No finite EEG–visit matches for spike-first (closest-visit) analysis.');
end

%% -------- 4) For each (Patient, Session), keep ONLY the closest visit WITH finite seizure freq --------
[grpPS, ~, ~] = findgroups(Pairs.Patient, Pairs.Session);

keepPerEEG = false(height(Pairs),1);

for g = 1:max(grpPS)
    idx = find(grpPS == g);   % rows for one (Patient, Session)

    % Only consider visits with non-NaN seizure frequency
    idx_ok = idx(isfinite(Pairs.Freq_R1(idx)));

    if isempty(idx_ok)
        % No usable visit frequencies for this EEG → keep none (EEG will be dropped)
        continue;
    end

    if numel(idx_ok) == 1
        keepPerEEG(idx_ok) = true;
        continue;
    end

    % Choose closest among the usable visits
    [~, bestLocal] = min(abs(Pairs.GapDays(idx_ok)));
    keepPerEEG(idx_ok(bestLocal)) = true;
end

Pairs = Pairs(keepPerEEG, :);

fprintf('[Spike-first (closest finite-visit) analysis] After closest-visit filter: %d EEG–visit pairs remain.\n', ...
    height(Pairs));
SessAgg = Pairs;


%% -------- 5) For each patient, pick LOWEST- and HIGHEST-spike EEG --------
[grpP, uPatients] = findgroups(SessAgg.Patient);

nP           = numel(uPatients);
Spike_low    = nan(nP,1);
Spike_high   = nan(nP,1);
Freq_low     = nan(nP,1);
Freq_high    = nan(nP,1);
EEGDate_low  = NaT(nP,1);
EEGDate_high = NaT(nP,1);
Visit_low    = NaT(nP,1);
Visit_high   = NaT(nP,1);
Sess_low_id  = nan(nP,1);
Sess_high_id = nan(nP,1);

for k = 1:nP
    idx = find(grpP == k);
    if numel(idx) < 2
        % Need at least 2 EEG sessions with matched visits for this patient
        continue;
    end

    r = SessAgg.SpikeRate_perHour(idx);   % spike rate per EEG
    f = SessAgg.Freq_R1(idx);         % seizure freq per EEG (from its single visit)

    % Require at least two sessions with finite seizure frequency & spike rate
    finiteMask = isfinite(r) & isfinite(f);
    if nnz(finiteMask) < 2
        continue;
    end

    idx_use = idx(finiteMask);
    r_use   = r(finiteMask);
    f_use   = f(finiteMask);

    % Identify lowest-spike and highest-spike EEG within this patient
    [~, iLowLocal]  = min(r_use);
    [~, iHighLocal] = max(r_use);

    if iLowLocal == iHighLocal
        % All spike rates identical for this patient's usable sessions
        continue;
    end

    if r_use(iHighLocal) - r_use(iLowLocal) <= min_abs_diff_spikes
        % Spike rate between high and low spike freq EEG too similar
        continue;
    end

    idx_low  = idx_use(iLowLocal);
    idx_high = idx_use(iHighLocal);

    % ---- NEW: Exclude patients where the SAME visit is paired to both EEGs ----
    if SessAgg.VisitDate(idx_low) == SessAgg.VisitDate(idx_high)
        % This visit would contribute to both low- and high-spike EEG → skip patient
        continue;
    end

    % Store values
    Spike_low(k)    = SessAgg.SpikeRate_perHour(idx_low);
    Spike_high(k)   = SessAgg.SpikeRate_perHour(idx_high);
    Freq_low(k)     = SessAgg.Freq_R1(idx_low);
    Freq_high(k)    = SessAgg.Freq_R1(idx_high);
    EEGDate_low(k)  = SessAgg.EEG_Date(idx_low);
    EEGDate_high(k) = SessAgg.EEG_Date(idx_high);
    Visit_low(k)    = SessAgg.VisitDate(idx_low);
    Visit_high(k)   = SessAgg.VisitDate(idx_high);
    Sess_low_id(k)  = SessAgg.Session(idx_low);
    Sess_high_id(k) = SessAgg.Session(idx_high);
end

%% -------- 6) Keep only valid pairs and run paired Wilcoxon --------
valid = isfinite(Spike_low) & isfinite(Spike_high) & ...
        isfinite(Freq_low)  & isfinite(Freq_high);

Spike_low    = Spike_low(valid);
Spike_high   = Spike_high(valid);
Freq_low     = Freq_low(valid);
Freq_high    = Freq_high(valid);
EEGDate_low  = EEGDate_low(valid);
EEGDate_high = EEGDate_high(valid);
Visit_low    = Visit_low(valid);
Visit_high   = Visit_high(valid);
Sess_low_id  = Sess_low_id(valid);
Sess_high_id = Sess_high_id(valid);
uPatients    = uPatients(valid);

nPairs = numel(Freq_low);
fprintf('[Spike-first (closest-visit) analysis] %d patients have valid low/high-spike EEG pairs with distinct visits.\n', ...
    nPairs);

% Build a table so you can inspect per-patient data used in the analysis
PairedEEGTable = table(uPatients, Sess_low_id, Sess_high_id, ...
    EEGDate_low, EEGDate_high, Visit_low, Visit_high, ...
    Spike_low, Spike_high, Freq_low, Freq_high, ...
    'VariableNames', {'Patient','Session_low','Session_high', ...
                      'EEGDate_low','EEGDate_high', ...
                      'VisitDate_low','VisitDate_high', ...
                      'SpikeRate_low_perHour','SpikeRate_high_perHour', ...
                      'Freq_low_perMonth','Freq_high_perMonth'});

fprintf('[Spike-first (closest-visit) analysis] Created PairedEEGTable with %d rows.\n', ...
    height(PairedEEGTable));


% Paired Wilcoxon (signed-rank) on seizure frequency
[p_signed, ~, stats_signed] = signrank(Freq_high, Freq_low, 'method','approx');
medDiffFreq = median(Freq_high - Freq_low, 'omitnan');

fprintf('\n=== Spike-first (closest-visit) analysis: seizure frequency at HIGH- vs LOW-spike EEGs ===\n');
fprintf('Number of patients with valid pairs: %d\n', nPairs);
fprintf('Median seizure freq at LOW-spike EEG:  %.3f seizures/month\n', ...
    median(Freq_low,'omitnan'));
fprintf('Median seizure freq at HIGH-spike EEG: %.3f seizures/month\n', ...
    median(Freq_high,'omitnan'));
fprintf('Median paired difference (HIGH - LOW):  %.3f seizures/month\n', medDiffFreq);

if isfield(stats_signed,'zval') && ~isempty(stats_signed.zval)
    z        = stats_signed.zval;
    r_effect = z / sqrt(nPairs);  % simple effect size
    fprintf('Wilcoxon signed-rank p = %.3g, z = %.3f, r ≈ %.3f\n', ...
        p_signed, z, r_effect);
else
    fprintf('Wilcoxon signed-rank p = %.3g\n', p_signed);
end

% Checking correlation
if 0
    delta_spikes = Spike_high-Spike_low;
    delta_freq = Freq_high-Freq_low;
    figure
    plot(delta_spikes,delta_freq,'o');
    [r,p]=corr(delta_spikes,delta_freq,"Type","Spearman");
end

%% -------- 8) Paired plot of seizure frequency --------
figP = figure('Color','w','Position',[100 100 650 500]);
axP  = axes(figP); hold(axP,'on'); box(axP,'off'); grid(axP,'on');

low_freq_log = to_log10_per_month(Freq_low,  EPS_FREQ);
high_freq_log = to_log10_per_month(Freq_high,  EPS_FREQ);
Y_S3_ZERO = log10(EPS_FREQ);

% Light grey paired lines
for i = 1:nPairs
    plot(axP, [1 2], [low_freq_log(i), high_freq_log(i)], '-', ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
end

hold on
yline(axP, Y_S3_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

% Overlay points (jitter a bit horizontally so they don't sit exactly on 1/2)
x_low  = 1 + (rand(nPairs,1)-0.5)*0.04;
x_high = 2 + (rand(nPairs,1)-0.5)*0.04;
scatter(axP, x_low,  low_freq_log,  30, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(axP, x_high, high_freq_log, 30, 'filled', 'MarkerFaceAlpha', 0.6);

xlim(axP, [0.5 2.5]);
xticks(axP, [1 2]);
xticklabels(axP, {'Low-spike EEG','High-spike EEG'});
ylabel(axP, 'log_{10}(seizures/month)', 'FontSize', 14);


% Add a significance bar using your add_sigbar helper if available

yMax = max([low_freq_log; high_freq_log]) * 1.4;
yMin = min([low_freq_log; high_freq_log]) * 1.2;
ylim(axP, [yMin, yMax]);
ySig = yMax - 0.2*(yMax);
if p_signed < 0.001
    pLabel = 'p < 0.001';
else
    pLabel = sprintf('p = %.2g', p_signed);
end
add_sigbar(axP, 1, 2, ySig, pLabel);
title(axP, sprintf('Paired seizure frequencies (N = %d patients)', nPairs), ...
    'FontSize', 20);

set(axP,'FontSize',20);
exportgraphics(figP, figS3_out, 'Resolution', 300);


delta_spikes = PairedEEGTable.SpikeRate_high_perHour-PairedEEGTable.SpikeRate_low_perHour;
delta_szs = PairedEEGTable.Freq_high_perMonth-PairedEEGTable.Freq_low_perMonth;

%% ======================= WRITE HTML RESULTS SUMMARY =======================
if ~exist(fileparts(resultsHtml),'dir')
    mkdir(fileparts(resultsHtml));
end

fid = fopen(resultsHtml,'w');
if fid == -1
    error('Could not open HTML results file: %s', resultsHtml);
end

fprintf(fid, '<html><head><meta charset="UTF-8">\n');
fprintf(fid, '<title>Spike Rate & Seizure Frequency Results</title></head><body>\n');
fprintf(fid, '<h1>Spike Rate and Seizure Frequency Results Summary</h1>\n');

%% Table 1 general overview

fprintf(fid, '<h2>Cohort summary</h2>\n');
fprintf(fid,['<p>We included %d patients with at least one outpatient routine EEG '...
    '(%d EEGs). %d patients (%1.1f%%) carried a '...
    'diagnosis of epilepsy. Median (IQR) monthly seizure frequency was '...
    '%1.2f (%1.2f-%1.2f), and median (IQR) spikes/hour across EEGs '...
    'was %1.2f (%1.2f-%1.2f) (Table 1).</p>'],...
    N_total,n_eegs_all,n_epi,n_epi/N_total*100,...
    sf_med,sf_q(1), sf_q(2),...
    sr_med, sr_q(1), sr_q(2));


%% ---- Figure 1: Control Panels ----
fprintf(fid, '<h2>Spike rates by patient groups</h2>\n');


% Panel A
pA_str = format_p_html(p_rankSum_A);
fprintf(fid,['<p>Automatically detected spike rates were higher in EEGs with clinically-reported spikes '...
            '(median spike rate %.2f (%.2f-%.2f) spikes/hour) than '...
            'those without reported spikes (%.2f (%.2f-%.2f) '...
            'spikes/hour) (%s, Cliff''s &delta; = %.2f; Fig. 1A). '],...
            m_pre,iqr_pre(1),iqr_pre(2),m_abs,...
            iqr_abs(1),iqr_abs(2),pA_str,effectA_cliff);


% Panel B
pB_str = format_p_html(p_rankSum_B);
fprintf(fid,['Patients with epilepsy also had higher spike rates '...
            '(%.2f (%.2f-%.2f) spikes/hour) than those '...
            'with NESD (%.2f (%.2f-%.2f) spikes/hour) '...
            '(%s, &delta; = %.2f; Fig. 1B). '],...
            m_ep,iqr_ep(1),iqr_ep(2),m_nes,...
            iqr_nes(1),iqr_nes(2),pB_str,effectB_cliff);

% Panel C
fprintf(fid,['Spike rates differed across epilepsy subtypes (Kruskal–Wallis '...
            '%s, η² ≈ %.3f). '],format_p_html(p_kw_C),eta2_kw_C);

fprintf(fid,['Generalized epilepsy demonstrated higher spike rates '...
            '(%.2f (%.2f-%.2f)) than temporal '...
            '(%.2f (%.2f-%.2f); Bonferroni-adjusted '...
            '%s) and frontal '...
            '(%.2f (%.2f-%.2f); '...
            '%s). The temporal versus frontal comparison '...
            'was not significant (%s; Fig. 1C).</p>'],...
            SubtypeStatsTable.Median(1),SubtypeStatsTable.P25(1),SubtypeStatsTable.P75(1),...
            SubtypeStatsTable.Median(2),SubtypeStatsTable.P25(2),SubtypeStatsTable.P75(2),...
            format_p_html(p_pair_bonf(1)),...
            SubtypeStatsTable.Median(3),SubtypeStatsTable.P25(3),SubtypeStatsTable.P75(3),...
            format_p_html(p_pair_bonf(2)),...
            format_p_html(p_pair_bonf(3)));


% Panel C



%% ---- Figure 2: Spearman correlations ----
fprintf(fid, '<h2>Relationship between spike rate and seizure frequency</h2>\n');

fprintf(fid,['<p>Across all epilepsy patients (N = %d), spike rate '...
    'and seizure frequency were positively correlated '...
    '(&rho; = %.2f, %s).'],n_all_main,rs_all_main,format_p_html(p_all_main));

fprintf(fid,[' Subtype-specific correlations were significant for generalized epilepsy '...
    '(N = %d, &rho; = %.2f, Bonferroni-adjusted %s) and '...
    'temporal lobe epilepsy (N = %d, &rho; = %.2f, %s), but not frontal '...
    'lobe epilepsy (N = %d, &rho; = %.2f, %s; Fig 2A-D). '],...
    SpearmanResults_main.N(1),SpearmanResults_main.Spearman_r(1),format_p_html(SpearmanResults_main.p_bonf(1)),...
    SpearmanResults_main.N(2),SpearmanResults_main.Spearman_r(2),format_p_html(SpearmanResults_main.p_bonf(2)),...
    SpearmanResults_main.N(3),SpearmanResults_main.Spearman_r(3),format_p_html(SpearmanResults_main.p_bonf(3)));

fprintf(fid,['Results were consistent when restricting analyses to '...
    'patients with non-zero spike rates and seizure frequencies (Fig. S1). ']);

fprintf(fid,['Patients who ever had spikes reported on at least one EEG also had '...
    'higher mean seizure frequencies than patients whose EEGs consistently '...
    'lacked spikes (Fig. S2). ']);

fprintf(fid,['Together, these findings indicate that spike burden reflects seizure severity at the population level. ']);


%% Paired analysis

fprintf(fid,['Among %d patients with multiple EEGs-clinic visit pairs, '...
    'seizure frequency did not differ significantly between periods '...
    'near low- versus high-spike-rate EEGs '...
    '(%1.1f (%1.1f-%1.1f) vs %1.1f (%1.1f-%1.1f) seizures/month; '...
    'W = %1.1f, %s; Fig. S3).</p>'],...
    nPairs,median(Freq_low,'omitnan'),prctile(Freq_low,25),prctile(Freq_low,75),...
    median(Freq_high,'omitnan'),prctile(Freq_high,25),prctile(Freq_high,75),...
    stats_signed.signedrank,format_p_html(p_signed));

fprintf(fid, '</body></html>\n');
fclose(fid);


%% ======================= HELPER FUNCTIONS ===================================




function [PatientTypingAll, SzFreqPerPatient] = build_patient_metrics_from_report(R, canonical3)

% ----- seizure frequency per patient (MeanSzFreq) -----
PV = table('Size',[0 4], 'VariableTypes',{'double','datetime','double','double'}, ...
           'VariableNames',{'Patient','VisitDate','Freq','HasSz'});

% Loop over eegs in redcap report to get visit level info
for j = 1:height(R)

    % Get patient id
    pid = double(R.patient_id(j));

    % date formatting
    ds = strtrim(R.visit_dates_deid(j)); 
    dd = string(jsondecode(char(ds)));
    d = datetime(dd,'InputFormat','yyyy-MM-dd');

    % seizure frequency formatting
    s = strtrim(R.sz_freqs(j)); 
    s = regexprep(s,'null','NaN','ignorecase');
    v = double(jsondecode(char(s)));
    v(~isfinite(v)) = NaN; 
    v(v<0) = NaN; % -2 turns into NaN here

    % hasSz
    hs = strtrim(R.visit_hasSz(j)); 
    h = double(jsondecode(char(hs)));
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
Has_agg  = splitapply(@(x) max_hasSz(x(isfinite(x))), PV.HasSz, gv);    % get a single unique has sz for each visit

% build table for unique visits - can be checked against redcap
Vuniq = table(pid_keys, date_keys, Freq_agg, Has_agg, ...
    'VariableNames', {'Patient','VisitDate','Freq','HasSz'});
Vuniq.Freq(Vuniq.Freq<0) = NaN; % negative sz frequencies should be nans
% -------- keep three versions of Freq ----------
Vuniq.FreqNoReplace = Vuniq.Freq;          % raw frequency (for QC only)
Vuniq.Freq_R1       = Vuniq.FreqNoReplace; % Rule 1 only
Vuniq.Freq_R12      = Vuniq.FreqNoReplace; % Rule 1 + Rule 2

% =========== If HasSz==0 & missing Freq -> seizure frequency = 0 ====================
mask_rule1 = ~isfinite(Vuniq.Freq_R1) & Vuniq.HasSz==0;
Vuniq.Freq_R1(mask_rule1)  = 0;

% Now group by patient
[gpV,pidsV] = findgroups(Vuniq.Patient);

% Patient-level mean sz frequency:
MeanSzFreq        = splitapply(@(x) mean(x,'omitnan'), Vuniq.Freq_R1,  gpV);

SzFreqPerPatient = table(pidsV, MeanSzFreq, ...
    'VariableNames', {'Patient','MeanSzFreq'});

% 4) Fraction of visits with HasSz==1 (among HasSz==0 or 1)
FracVisits_HasSz1 = splitapply(@local_frac, Vuniq.HasSz, gpV);

SzFreqPerPatient.FracVisits_HasSz1 = FracVisits_HasSz1;


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



function Views = build_filtered_view(SessionsFiltered, ReportIn, PatientTypingAll, SzFreqPerPatient, ...
                                     NESD_LABEL, badTypes, canonical3)
% build_filtered_view
%   Bundle all "view building" into one place:
%     - Filter report rows to sessions in SessionsIn
%     - Compute patient-level spike rates + attach typing
%     - Mark epilepsy vs NESD
%     - Compute per-session spike rates (spikes/min)
%     - Build spike+seizure tables for Spearman plots (All + Typed)
%     - Build canonical3 subtype tables + stats + pairwise p-values
%
%   Inputs:
%     SessionsIn        : session-level spike table (already outpatient-filtered)
%     ReportIn          : REDCap report table (already outpatient-filtered)
%     PatientTypingAll  : patient-level typing (EpilepsyType / EpilepsySpecific / EpiType3)
%     SzFreqPerPatient  : patient-level seizure frequency (MeanSzFreq / MeanSzFreq_raw)
%     countCol, durCol  : names used to recompute per-session rates
%     NESD_LABEL        : string, e.g. "Non-Epileptic Seizure Disorder"
%     badTypes          : lower-case strings to exclude from "epilepsy"
%     canonical3        : ["General","Temporal","Frontal"]
%
%   Output:
%     Views (struct) with fields:
%       - SessionsForFigures, ReportForKeptSessions, PatientTypingFiltered
%       - SessionLevelSpikeRates, PatientLevelSpikeRates
%       - PatientSpikeSz_All, PatientSpikeSz_Typed
%       - IsEpilepsyMask, IsNESDMask
%       - Canonical3_SubsetTable, Canonical3_Stats, Canonical3_Pairs
%       - PvalsPairwise, PvalsPairwiseBonf

    %% ---------- 1) Basic filtering: sessions + reports ----------
    % Keep only report rows that correspond to kept sessions
    SessKeys = unique(SessionsFiltered(:,{'Patient','Session'}));
    ReportForKeptSessions = innerjoin(ReportIn, SessKeys, ...
        'LeftKeys', {'patient_id','session_number'}, ...
        'RightKeys', {'Patient','Session'});

    % Promote Patient/Session columns on report for convenience
    ReportForKeptSessions.Patient = ReportForKeptSessions.patient_id;
    ReportForKeptSessions.Session = ReportForKeptSessions.session_number;

    % Patients actually present in this view
    PatientsKept   = unique(SessionsFiltered(:,{'Patient'}));
    TypingFiltered = innerjoin(PatientTypingAll, PatientsKept, 'Keys','Patient');


    %% ---------- 2) Patient-level mean spike rates + typing ----------
    [gpS, pidsS] = findgroups(SessionsFiltered.Patient);

    MeanSpikeRate_perHour = splitapply(@(x) mean(x,'omitnan'), ...
        SessionsFiltered.SpikeRate_perHour, gpS);

    PatientLevelSpikeRates = table(pidsS, MeanSpikeRate_perHour,  ...
        'VariableNames', {'Patient','MeanSpikeRate_perHour'});

    % Attach epilepsy typing
    PatientLevelSpikeRates = innerjoin(PatientLevelSpikeRates, TypingFiltered, 'Keys','Patient');


    %% ---------- 3) Epilepsy vs NESD masks ----------
    etype_norm = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));

    IsNESDMask = (etype_norm == lower(strtrim(NESD_LABEL)));
    IsBadType  = ismember(etype_norm, badTypes) | ...
                 ismissing(etype_norm) | (strlength(etype_norm)==0);

    % "Epilepsy" here = not NESD and not any of the excluded/bad labels
    IsEpilepsyMask = ~IsNESDMask & ~IsBadType;

    %% ---------- 4) Per-session spike rates (already 1 row per EEG) ----------
    SessionLevelSpikeRates = SessionsFiltered(:, {'Patient','Session','SpikeRate_perHour'});
    SessionLevelSpikeRates.Properties.VariableNames{'SpikeRate_perHour'} = 'SpikesPerHour';

    %% ---------- 5) Spike + seizure tables for Spearman ----------
    % Start from seizure frequencies for patients in this view
    SzFreqFiltered = innerjoin(SzFreqPerPatient, PatientsKept, 'Keys','Patient');

    % Further restrict to epilepsy patients (for "all epilepsy" panels)
    EpilepsyPatients = table( ...
        PatientLevelSpikeRates.Patient(IsEpilepsyMask), ...
        'VariableNames', {'Patient'});

    SzFreqEpilepsy = innerjoin(SzFreqFiltered, EpilepsyPatients, 'Keys','Patient');

    % ---- 5A. All epilepsy: spike rate + sz frequency ----
    PatientSpikeSz_All = innerjoin( ...
        PatientLevelSpikeRates(IsEpilepsyMask, {'Patient','MeanSpikeRate_perHour'}), ...
        SzFreqEpilepsy, 'Keys','Patient');

    % Require finite spike rate and at least one non-missing sz frequency
    keepAll = isfinite(PatientSpikeSz_All.MeanSpikeRate_perHour) & ...
          (isfinite(PatientSpikeSz_All.MeanSzFreq));

    PatientSpikeSz_All = PatientSpikeSz_All(keepAll, :);


    % ---- 5B. Typed (Frontal/Temporal/General) subset ----
    keepCanon3 = ~ismissing(TypingFiltered.EpiType3) & ...
                 ismember(TypingFiltered.EpiType3, canonical3);

    SzFreqCanon = innerjoin(SzFreqEpilepsy, ...
        TypingFiltered(keepCanon3, {'Patient','EpiType3'}), 'Keys','Patient');

    PatientSpikeSz_Typed = innerjoin( ...
        PatientLevelSpikeRates(IsEpilepsyMask, {'Patient','MeanSpikeRate_perHour'}), ...
        SzFreqCanon, 'Keys','Patient');

    keepTyped = isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perHour) & ...
                (isfinite(PatientSpikeSz_Typed.MeanSzFreq)) & ...
                ~ismissing(PatientSpikeSz_Typed.EpiType3);
    PatientSpikeSz_Typed = PatientSpikeSz_Typed(keepTyped, :);


    %% ---------- 6) Canonical3 subtype stats (General/Temporal/Frontal) ----------
    % Build EpiType4 = {General, Temporal, Frontal, (other/empty)}
    % Do the stats comparing spike rates across epilepsy groups
    EpiType4 = strings(height(PatientLevelSpikeRates),1);
    et_low = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));
    es_low = lower(strtrim(string(PatientLevelSpikeRates.EpilepsySpecific)));

    isGeneralType = contains(et_low, "general");
    isTemporal    = contains(es_low, "temporal");
    isFrontal     = contains(es_low, "frontal");

    EpiType4(isGeneralType)                            = "General";
    EpiType4(~isGeneralType & isTemporal)              = "Temporal";
    EpiType4(~isGeneralType & ~isTemporal & isFrontal) = "Frontal";

    PatientLevelSpikeRates.EpiType4 = EpiType4;

    InCanonical3Mask = IsEpilepsyMask & ismember(PatientLevelSpikeRates.EpiType4, canonical3);

    Canonical3_SubsetTable = PatientLevelSpikeRates(InCanonical3Mask, ...
        {'EpiType4','MeanSpikeRate_perHour'});
    Canonical3_SubsetTable.EpiType4 = categorical(string(Canonical3_SubsetTable.EpiType4), canonical3);

    [g3, cats3]  = findgroups(Canonical3_SubsetTable.EpiType4);
    medVals      = splitapply(@(x) median(x,'omitnan'), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
    iqrVals      = splitapply(@(x) [prctile(x,25), prctile(x,75)], Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
    nVals        = splitapply(@(x) sum(isfinite(x)), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);

    Canonical3_Stats = table(cats3, nVals, medVals, ...
    iqrVals(:,1), iqrVals(:,2), ...
    'VariableNames', {'EpiType4','GroupCount','Median','P25','P75'});


    Canonical3_Pairs = ["General","Temporal";
                        "General","Frontal";
                        "Temporal","Frontal"];

    p_pair = NaN(3,1);
    for i = 1:3
        A = Canonical3_Pairs(i,1);
        B = Canonical3_Pairs(i,2);
        xa = Canonical3_SubsetTable.MeanSpikeRate_perHour(Canonical3_SubsetTable.EpiType4==A);
        xb = Canonical3_SubsetTable.MeanSpikeRate_perHour(Canonical3_SubsetTable.EpiType4==B);
        if nnz(isfinite(xa))>=3 && nnz(isfinite(xb))>=3
            p_pair(i) = ranksum(xa, xb, 'method','approx');
        end
    end
    p_pair_bonf = min(p_pair * 3, 1); % just multiplying by 3 for bonferroni correction


    %% ---------- 7) Pack outputs ----------
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



%% ======================= HELPERS =======================

function ylog = to_log10_per_hour(x_per_hour, eps_rate)
    rate_hr = double(x_per_hour);
    rate_hr(~isfinite(rate_hr) | rate_hr <= 0) = eps_rate;
    ylog = log10(rate_hr);
end

%

function add_sigbar(ax, x1, x2, y, ptext)
tick = 0.03 * diff(ax.YLim);
plot(ax, [x1 x1 x2 x2], [y- tick, y, y, y- tick], 'k-', 'LineWidth', 1.3);
text(ax, mean([x1 x2]), y + 0.02*diff(ax.YLim), ptext, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);
end




function [SpearmanResults, rs_all, p_all, n_all] = ...
    spearman_plotting_function(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
                              canonical3, spearman_xLims, spearman_yLims, ...
                              fig_out, ...
                              freqFieldName, labelSuffix, nonZeroOnly)
% freqFieldName: 'MeanSzFreq' or 'MeanSzFreq_raw'
% nonZeroOnly: true  -> require spikes>0 AND seizures>0
%               false -> allow zeros, use eps guards

%% ---- Overall: define mask once and reuse it everywhere ----
x_all = PatientSpikeSz_All.MeanSpikeRate_perHour;
y_all = PatientSpikeSz_All.(freqFieldName); % mean sz frequency

mask_all = isfinite(x_all) & isfinite(y_all);

if nonZeroOnly
    mask_all = mask_all & (x_all > 0) & (y_all > 0);
end

n_all = sum(mask_all);


% Correlation uses ONLY rows with finite x,y
[rs_all, p_all] = corr(x_all(mask_all), y_all(mask_all), ...
    'Type','Spearman','Rows','complete');


%% ---- By group (Frontal / Temporal / General) ----
rowsOut = {};
for g = canonical3
    mBase = (PatientSpikeSz_Typed.EpiType3 == g);
    x = PatientSpikeSz_Typed.MeanSpikeRate_perHour(mBase);
    y = PatientSpikeSz_Typed.(freqFieldName)(mBase);

    mask = isfinite(x) & isfinite(y);
    if nonZeroOnly
        mask = mask & (x > 0) & (y > 0);
    end

    n = sum(mask);
    [rs, p] = corr(x(mask), y(mask), 'Type','Spearman','Rows','complete');
 
    rowsOut(end+1,:) = {char(g), n, rs, p}; %#ok<SAGROW>
end

SpearmanResults = cell2table(rowsOut, ...
    'VariableNames', {'Group','N','Spearman_r','p_raw'});
k = height(SpearmanResults);
SpearmanResults.p_bonf = min(SpearmanResults.p_raw * k, 1);

%% ---- Log transforms with eps guards: use only mask_all rows ----
x_used = x_all(mask_all);
y_used = y_all(mask_all);

% Small positive epsilons based on the used data
minpos_rate = min(x_used(x_used > 0));
if isempty(minpos_rate), minpos_rate = 1e-6; end

minpos_sz = min(y_used(y_used > 0));
if isempty(minpos_sz),   minpos_sz   = 1e-6; end

eps_rate = 0.5 * minpos_rate;
eps_sz   = 0.5 * minpos_sz;

% Overall table (ONLY rows used in corr)
Tall = table;
Tall.SpikeRate_perHour = x_used;
Tall.SzFreq           = y_used;

Tall.logSpikeRate = log10(Tall.SpikeRate_perHour + ...
                      (Tall.SpikeRate_perHour<=0).*eps_rate);
Tall.logSzFreq    = log10(Tall.SzFreq + ...
                      (Tall.SzFreq<=0).*eps_sz);

% Zero/non-zero classification on the same used rows
isZeroSz_all   = (Tall.SzFreq == 0);
isZeroRate_all = (Tall.SpikeRate_perHour == 0);
onlySz_all     =  isZeroSz_all & ~isZeroRate_all;
onlyRate_all   = ~isZeroSz_all &  isZeroRate_all;
bothZero_all   =  isZeroSz_all &  isZeroRate_all;
nonZero_all    = ~(isZeroSz_all | isZeroRate_all);

% ---- Typed table: unchanged, but already uses finite mask ----
T = table;
T.EpiType3      = categorical(PatientSpikeSz_Typed.EpiType3, canonical3);
T.logSpikeRate  = log10(PatientSpikeSz_Typed.MeanSpikeRate_perHour + ...
                    (PatientSpikeSz_Typed.MeanSpikeRate_perHour<=0).*eps_rate);
T.logSzFreq     = log10(PatientSpikeSz_Typed.(freqFieldName) + ...
                    (PatientSpikeSz_Typed.(freqFieldName)<=0).*eps_sz);
ok = isfinite(T.logSpikeRate) & isfinite(T.logSzFreq) & ~ismissing(T.EpiType3);
if nonZeroOnly
    ok = ok & (PatientSpikeSz_Typed.MeanSpikeRate_perHour>0) & ...
             (PatientSpikeSz_Typed.(freqFieldName)>0);
end
T = T(ok,:);
presentCats = categories(removecats(T.EpiType3));

% ---- Axes helpers ----
xZero = log10(eps_sz);
yZero = log10(eps_rate);
xLims = spearman_xLims;
yLims = spearman_yLims;


%% ---- Draw figure ----
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
plot(axA2, Tall.logSzFreq(onlySz_all),   Tall.logSpikeRate(onlySz_all),   '*', ...
    'Color',baseColor,'MarkerSize',7,'LineWidth',1);
plot(axA2, Tall.logSzFreq(onlyRate_all),    Tall.logSpikeRate(onlyRate_all),'*', ...
    'Color',baseColor,'MarkerSize',8,'LineWidth',1);
plot(axA2, Tall.logSzFreq(bothZero_all), Tall.logSpikeRate(bothZero_all),'*', ...
    'Color',baseColor,'MarkerSize',8,'LineWidth',1.2);

% Regression line on the same used rows
if n_all >= 3
    X = [ones(n_all,1), Tall.logSzFreq];
    b = X \ Tall.logSpikeRate;
    xgrid = linspace(xLims(1), xLims(2), 300)';
    plot(axA2, xgrid, b(1) + b(2)*xgrid, 'k-', 'LineWidth', 2);
end

xlim(axA2, xLims); ylim(axA2, yLims);
xlabel(axA2,'log_{10} Seizures per month','FontSize',fontL);
ylabel(axA2,'log_{10} Spikes per hour','FontSize',fontL);

title(axA2, sprintf('A. All epilepsy%s (N=%d)', labelSuffix, n_all), ...
    'FontSize',fontL, 'FontWeight','bold');

if p_all < 0.001
    txtA = sprintf('ρ=%.2f, p<0.001', rs_all);
else
    txtA = sprintf('ρ=%.2f, p=%.2g', rs_all, p_all);
end
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

    % mask for this group using Typed table
    idx = (PatientSpikeSz_Typed.EpiType3 == tgtCat) & ...
      isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perHour) & ...
      isfinite(PatientSpikeSz_Typed.(freqFieldName));

    if nonZeroOnly
        idx = idx & (PatientSpikeSz_Typed.MeanSpikeRate_perHour > 0) & ...
                   (PatientSpikeSz_Typed.(freqFieldName) > 0);
    end


    logX = log10(PatientSpikeSz_Typed.(freqFieldName)(idx) + ...
                 (PatientSpikeSz_Typed.(freqFieldName)(idx)<=0).*eps_sz);
    logY = log10(PatientSpikeSz_Typed.MeanSpikeRate_perHour(idx) + ...
                 (PatientSpikeSz_Typed.MeanSpikeRate_perHour(idx)<=0).*eps_rate);

    isZeroSz_g   = (PatientSpikeSz_Typed.(freqFieldName)==0);
    isZeroRate_g = (PatientSpikeSz_Typed.MeanSpikeRate_perHour==0);
    ok_g = idx;
    onlySz = isZeroSz_g(ok_g) & ~isZeroRate_g(ok_g);
    onlyRt = ~isZeroSz_g(ok_g) & isZeroRate_g(ok_g);
    bothZ  = isZeroSz_g(ok_g) & isZeroRate_g(ok_g);
    nonZ   = ~(onlySz | onlyRt | bothZ);

    xline(ax, xZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
    yline(ax, yZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    scatter(ax, logX(nonZ), logY(nonZ), 18, col, 'filled', 'MarkerFaceAlpha', 0.35);
    if any(onlySz)
        plot(ax, logX(onlySz), logY(onlySz),   '*','Color',col,'MarkerSize',8,'LineWidth',1.1);
    end
    if any(onlyRt)
        plot(ax, logX(onlyRt),    logY(onlyRt),'*','Color',col,'MarkerSize',8,'LineWidth',1.1);
    end
    if any(bothZ)
        plot(ax, logX(bothZ),  logY(bothZ), '*','Color',col,'MarkerSize',9,'LineWidth',1.2);
    end

    if nnz(nonZ) >= 3
        Xg = [ones(nnz(nonZ),1), logX(nonZ)];
        bg = Xg \ logY(nonZ);
        xg = linspace(xLims(1), xLims(2), 250)';
        plot(ax, xg, bg(1)+bg(2)*xg, '-', 'Color', col, 'LineWidth', 2);
    end

    xlim(ax, xLims); ylim(ax, yLims);

    row = SpearmanResults(strcmp(SpearmanResults.Group, string(panelOrder{p})), :);
    if row.p_bonf < 0.001
        txt = sprintf('ρ=%.2f, p_{bonf}<0.001', row.Spearman_r);
    else
        txt = sprintf('ρ=%.2f, p_{bonf}=%.2g', row.Spearman_r, row.p_bonf);
    end

    xlabel(ax,'log_{10} Seizures per month','FontSize',fontL);
    ylabel(ax,'log_{10} Spikes per hour','FontSize',fontL);

    
    nNow = sum(T.EpiType3==tgtCat);
    title(ax, sprintf('%s%s (N=%d)', panelTitle{p}, labelSuffix, nNow), ...
        'FontSize', fontL, 'FontWeight','bold');
    text(ax, 0.98, 0.95, txt, 'Units','normalized', ...
         'HorizontalAlignment','right','VerticalAlignment','top', ...
         'FontSize', fontL-3, 'FontWeight','bold');
    set(ax,'FontSize',fontL);
end

if ~exist(fileparts(fig_out),'dir'), mkdir(fileparts(fig_out)); end
exportgraphics(f2, fig_out, 'Resolution', 300);
fprintf('Saved Spearman figure: %s\n', fig_out);


end

function Yj = add_y_jitter_eps(Y, Y_ZERO, Y_LIMS, frac)
% add_y_jitter_eps
%   Adds tiny vertical jitter ONLY to points that sit exactly at Y_ZERO
%   (i.e., the eps floor), so overlapping zero-spike points become visible.
%
%   frac controls how big the jitter is as a fraction of the y-range
%   (e.g., 0.01 → ±0.5% of range, 0.02 → ±1% of range).

    if nargin < 4
        frac = 0.01;  % default: small jitter
    end

    Yj = Y;
    mask = abs(Y - Y_ZERO) < 1e-9;  % points at eps floor

    if any(mask)
        amp = frac * diff(Y_LIMS);  % total jitter range
        Yj(mask) = Yj(mask) + (rand(sum(mask),1) - 0.5) * amp;
    end
end

function pStr = format_p_html(p)

    if isnan(p)
        pStr = 'p = NaN';
        return;
    end

    % Rule 1: very small p-values
    if p < 0.001
        pStr = 'p &lt; 0.001';
        return;
    end

    % Rule 2: for p < 0.01 use 2 significant digits, fixed formatting
    if p < 0.01
        % 2 significant digits
        s = sprintf('%.2g', p);   % e.g., "0.0073"
        % Ensure a leading zero (sometimes %.2g gives ".0073")
        if startsWith(s, '.')
            s = ['0' s];
        end
        pStr = ['p = ' s];
        return;
    end

    % Rule 3: p >= 0.01 → always show 2 decimals (e.g., 0.60)
    pStr = sprintf('p = %.2f', p);
end


function d = cliff_delta(x1, x2)
% cliff_delta
%   Computes Cliff's delta using ranksum-based U statistic.
%   Positive d means x1 tends to be larger than x2.

    x1 = x1(:); x2 = x2(:);
    x1 = x1(isfinite(x1));
    x2 = x2(isfinite(x2));

    n1 = numel(x1);
    n2 = numel(x2);

    if n1 == 0 || n2 == 0
        d = NaN;
        return;
    end

    % ranksum returns sum of ranks for first sample (x1)
    [~,~,stats] = ranksum(x1, x2, 'method','approx');
    R1 = stats.ranksum;
    U1 = R1 - n1*(n1+1)/2;

    % Cliff's delta in terms of U
    d = (2 * U1 / (n1 * n2)) - 1;
end

function ylog = to_log10_per_month(freq_per_month, eps_freq)
% to_log10_per_month
%   Convert seizure frequency (per month) to log10 scale with a small
%   epsilon floor for zeros and nonpositive values.

    f = double(freq_per_month);
    f(~isfinite(f) | f <= 0) = eps_freq;
    ylog = log10(f);
end


function Vuniq_R1 = build_visit_level_freq_R1(R)
% build_visit_level_freq_R1
%
%   Output:
%     Vuniq_R1 (table with columns Patient, VisitDate, Freq_R1)

PV = table('Size',[0 4], ...
    'VariableTypes',{'double','datetime','double','double'}, ...
    'VariableNames',{'Patient','VisitDate','Freq','HasSz'});

for j = 1:height(R)
    pid = double(R.patient_id(j));

    % visit dates
    ds = strtrim(R.visit_dates_deid(j));
    if ds == "[]" || strlength(ds)==0
        continue;
    end
    dd = string(jsondecode(char(ds)));
    d  = datetime(dd,'InputFormat','yyyy-MM-dd');

    % seizure frequency
    s = strtrim(R.sz_freqs(j));
    s = regexprep(s,'null','NaN','ignorecase');
    v = double(jsondecode(char(s)));
    v(~isfinite(v)) = NaN;
    v(v<0)          = NaN;

    % hasSz
    hs = strtrim(R.visit_hasSz(j));
    h  = double(jsondecode(char(hs)));
    h(h==2) = NaN;  % set "2" to NaN

    assert(numel(d) == numel(v) && numel(d) == numel(h), ...
        'visit_dates_deid, sz_freqs, and visit_hasSz must have same length per row.');

    n = numel(h);
    if n==0, continue; end

    PV = [PV; table(repmat(pid,n,1), d(:), v(:), h(:), ...
        'VariableNames', PV.Properties.VariableNames)]; %#ok<AGROW>
end

[gv, pid_keys, date_keys] = findgroups(PV.Patient, PV.VisitDate);

Freq_agg = splitapply(@(x) mean(x(isfinite(x)),'omitnan'), PV.Freq, gv);
Has_agg  = splitapply(@(x) max_hasSz(x(isfinite(x))),    PV.HasSz, gv);

Vuniq = table(pid_keys, date_keys, Freq_agg, Has_agg, ...
    'VariableNames', {'Patient','VisitDate','Freq','HasSz'});

Vuniq.Freq(Vuniq.Freq < 0) = NaN;

% --- Rule 1: HasSz==0 & missing Freq -> 0 ---
Vuniq.Freq_R1 = Vuniq.Freq;
mask_rule1 = ~isfinite(Vuniq.Freq_R1) & (Vuniq.HasSz == 0);
Vuniq.Freq_R1(mask_rule1) = 0;

Vuniq_R1 = Vuniq(:, {'Patient','VisitDate','Freq_R1'});
end

function out = local_first_nonmissing(s)
    s = s(:);
    s = s(~ismissing(s) & strlength(s)>0);
    if isempty(s)
        out = "";
    else
        out = s(1);
    end
end

function out = local_frac(hsVec)
% local_frac
%   Given a vector of HasSz values for one patient across visits,
%   compute fraction of visits with HasSz==1 among visits with HasSz==0 or 1.

    hs = hsVec(isfinite(hsVec));   % keep only finite (should be 0/1)
    if isempty(hs)
        out = NaN;
        return;
    end

    validMask = (hs==0 | hs==1);
    denom = nnz(validMask);
    if denom == 0
        out = NaN;
    else
        out = nnz(hs(validMask)==1) / denom;
    end
end

function out = local_min_omitnan(a)
% local_min_omitnan
%   Returns the minimum of a vector, ignoring NaNs.
%   If all entries are NaN (or non-finite), returns NaN (scalar),
%   which keeps splitapply happy.

    a = a(isfinite(a));
    if isempty(a)
        out = NaN;
    else
        out = min(a);
    end
end
