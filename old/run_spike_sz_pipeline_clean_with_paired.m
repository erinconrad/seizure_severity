function run_spike_sz_pipeline_clean_with_paired()

%% ======================= RNG =======================
rng(0);

%% ======================= PATHS =======================
% Inputs
spikeSummaryMultiCsv = '../data/SN_counts/spike_counts_summary_multiThresh.csv';
reportCsv            = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-12-11_1316.csv';

% Outputs
fig1_out   = '../figures/Fig1.png';
fig2_out   = '../figures/Fig2.png';
figS1_out  = '../figures/FigS1.png';
figS2_out  = '../figures/FigS2.png';
figS3_out  = '../figures/FigS3.png';
Table1Csv  = fullfile('..','output','Table1.csv');
resultsHtml = '../output/results_summary.html';

%% ======================= PARAMETERS ===================
MAX_ROUTINE_HOURS = 4;

NESD_LABEL = "Non-Epileptic Seizure Disorder";
badTypes   = lower(["Uncertain if Epilepsy","Unknown or MRN not found",""]); % excluded from epilepsy
canonical3 = ["General","Temporal","Frontal"];

EPS_RATE = 30e-3;               % spikes/hour floor for log plotting
Y_ZERO  = log10(EPS_RATE);
Y_LIMS  = [-2 4];
TITLE_Y_OFFSET = 0.02;

spearman_xLims = [-3.5, 4];     % log10(seizures/month)
spearman_yLims = [-1.5, 3];     % log10(spikes/hour)

nBoot = 5000;
alpha = 0.05;

% SpikeNet threshold columns
countCol = "count_0_46";
durCol   = "Duration_sec";

% Allowed outpatient clinic visit types (other things are telephone calls,
% care management encounters, etc.)
allowable_visits = [
    "CONSULT VISIT","ESTABLISHED PATIENT VISIT",...
    "FOLLOW-UP PATIENT CLINIC","NEW PATIENT CLINIC","NEW PATIENT VISIT",...
    "NPV MANAGEMENT DURING COVID-19","NPV NEUROLOGY",...
    "RETURN ANNUAL VISIT","RETURN PATIENT EXTENDED","RETURN PATIENT VISIT",...
    "RPV MANAGEMENT DURING COVID-19","TELEHEALTH VIDEO VISIT RETURN"
];

%% ======================= LOAD SPIKE SUMMARY =======================
SpikeSummaryTable = readtable(spikeSummaryMultiCsv,'TextType','string','VariableNamingRule','preserve');

require_cols(SpikeSummaryTable, ["Patient","Session",countCol,durCol], "SpikeSummaryTable");

SpikeSummaryTable.SpikeRate_perHour = SpikeSummaryTable.(countCol) ./ SpikeSummaryTable.(durCol) * 3600; % get spikes/hour

%% ======================= LOAD REPORT =======================
ReportTable = readtable(reportCsv,'TextType','string','VariableNamingRule','preserve');

require_cols(ReportTable, ...
    ["patient_id","session_number","acquired_on", ...
     "report_PATIENT_CLASS","jay_in_or_out", ...
     "visit_type","visit_dates_deid","sz_freqs","visit_hasSz", ...
     "epilepsy_type","epilepsy_specific","nlp_gender","deid_birth_date","start_time_deid", ...
     "report_SPORADIC_EPILEPTIFORM_DISCHARGES","jay_focal_epi","jay_multifocal_epi","jay_gen_epi"], ...
    "ReportTable");

%% ======================= REMOVE NON-OUTPATIENT VISITS =======================
ReportTable = filter_visit_arrays_by_type(ReportTable, allowable_visits);

%% ======================= OUTPATIENT + ROUTINE EEG FILTER =======================
[SpikeSummaryTable, ReportTable] = filter_outpatient_routine( ...
    SpikeSummaryTable, ReportTable, durCol, MAX_ROUTINE_HOURS); % this removes EEGs that are inpatient or too long

% Sanity: keys match 1:1 (both the spike summary and report table should
% only have one row for each unique patient-session combo)
assert_unique_keys(SpikeSummaryTable, "Patient","Session","SpikeSummaryTable");
assert_unique_keys(ReportTable, "patient_id","session_number","ReportTable");

%% ======================= BUILD VISIT-LEVEL + PATIENT-LEVEL METRICS =======================
% Visit-level: Patient, VisitDate, Freq, HasSz, Freq_R1
Vuniq = build_visit_level_table_R1(ReportTable);

% Patient typing from report 
PatientTypingAll = build_patient_typing_from_report(ReportTable, canonical3);

% Patient seizure metrics: MeanSzFreq + fraction hasSz==1 among valid visits
SzFreqPerPatient = build_patient_seizure_metrics(Vuniq);

%% ======================= BUILD VIEW BUNDLE =======================
Views = build_filtered_view( ...
    SpikeSummaryTable, ReportTable, PatientTypingAll, SzFreqPerPatient, ...
    NESD_LABEL, badTypes, canonical3);

%% ======================= FIGURE 1 (CONTROL PANELS) =======================
[Fig1, Fig1Stats] = make_fig1_controls(Views, EPS_RATE, Y_ZERO, Y_LIMS, TITLE_Y_OFFSET, nBoot, alpha);
save_fig(Fig1, fig1_out);
fprintf('Saved Fig 1: %s\n', fig1_out);

%% ======================= FIGURE 2 + FIG S1 (SPEARMAN) =======================
PatientSpikeSz_All   = Views.PatientSpikeSz_All;
PatientSpikeSz_Typed = Views.PatientSpikeSz_Typed;

[SpearmanResults_main, rs_all_main, p_all_main, n_all_main, rho_lo_main, rho_hi_main, subtype_ci_main] = ...
    spearman_plotting_function(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
        canonical3, spearman_xLims, spearman_yLims, ...
        fig2_out, ...
        'MeanSzFreq', '', false);

[~, rs_all_nz, p_all_nz, n_all_nz, rho_lo_nz, rho_hi_nz, subtype_ci_nz] = ...
    spearman_plotting_function(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
        canonical3, spearman_xLims, spearman_yLims, ...
        figS1_out, ...
        'MeanSzFreq', ' (non-zero only)', true);

fprintf('Saved Fig 2:  %s\n', fig2_out);
fprintf('Saved Fig S1: %s\n', figS1_out);

%% ======================= FIG S2 (SZ FREQ BY REPORTED SPIKES ACROSS EEGs) =======================
FigS2 = make_figS2_sz_by_reported_spikes(Views, SzFreqPerPatient, nBoot, alpha, spearman_xLims);
save_fig(FigS2, figS2_out);
fprintf('Saved Fig S2: %s\n', figS2_out);

%% ======================= TABLE 1 =======================
Table1_flat = build_table1_flat(Views, SzFreqPerPatient, Vuniq, EPS_RATE, nBoot, alpha);
if ~exist(fileparts(Table1Csv),'dir'), mkdir(fileparts(Table1Csv)); end
writetable(Table1_flat, Table1Csv);
fprintf('Wrote Table 1 CSV: %s\n', Table1Csv);

%% ======================= PAIRED ANALYSIS + SENSITIVITY (FIG S3) =======================
[FigS3, PairedStats] = make_figS3_paired_and_sensitivity(Views, Vuniq, spearman_xLims, nBoot, alpha);
save_fig(FigS3, figS3_out);
fprintf('Saved Fig S3: %s\n', figS3_out);

%% ======================= WRITE HTML SUMMARY =======================
write_results_html(resultsHtml, Views, SzFreqPerPatient, ...
    Fig1Stats, ...
    SpearmanResults_main, rs_all_main, p_all_main, n_all_main, rho_lo_main, rho_hi_main, subtype_ci_main, ...
    Views.ReportForKeptSessions, ...
    PairedStats);

fprintf('Wrote HTML summary: %s\n', resultsHtml);

end

%% =====================================================================
%% ======================= CORE PIPELINE HELPERS =======================
%% =====================================================================

function require_cols(T, cols, name)
missing = setdiff(string(cols), string(T.Properties.VariableNames));
if ~isempty(missing)
    error('%s is missing required columns: %s', name, strjoin(missing, ", "));
end
end

function assert_unique_keys(T, pidCol, sesCol, name)
key = string(T.(pidCol)) + "|" + string(T.(sesCol));
u = unique(key);
if numel(u) ~= numel(key)
    error('%s has duplicated (Patient,Session) keys (%d duplicates).', name, numel(key)-numel(u));
end
end

function save_fig(figH, outPath)
if ~exist(fileparts(outPath),'dir'), mkdir(fileparts(outPath)); end
exportgraphics(figH, outPath, 'Resolution', 300);
end

function R = filter_visit_arrays_by_type(R, allowable_visits)
% Rewrites visit_type / visit_dates_deid / sz_freqs / visit_hasSz in-place
% keeping only entries whose visit_type is in allowable_visits.

vt_raw_all = strtrim(string(R.visit_type));
mask_null_only = (vt_raw_all == "[null]") | (vt_raw_all == "null");
R.visit_type(mask_null_only)        = "[]";
R.visit_dates_deid(mask_null_only)  = "[]";
R.sz_freqs(mask_null_only)          = "[]";
R.visit_hasSz(mask_null_only)       = "[]";

total_before = 0;
total_after  = 0;

for i = 1:height(R)
    vt_raw    = strtrim(string(R.visit_type(i)));
    dates_raw = strtrim(string(R.visit_dates_deid(i)));
    sz_raw    = strtrim(string(R.sz_freqs(i)));
    hs_raw    = strtrim(string(R.visit_hasSz(i)));

    if vt_raw=="" || vt_raw=="[]" || vt_raw=="<missing>"
        R.visit_type(i)        = "[]";
        R.visit_dates_deid(i)  = "[]";
        R.sz_freqs(i)          = "[]";
        R.visit_hasSz(i)       = "[]";
        continue
    end

    vt = json_to_string_array(vt_raw);

    dates = json_to_string_array(dates_raw);  % date strings
    sz    = json_to_double_array(sz_raw);     % doubles (NaN ok)
    hs    = json_to_double_array(hs_raw);     % doubles (0/1/2/NaN)

    % Strict length check
    if ~(numel(vt)==numel(dates) && numel(vt)==numel(sz) && numel(vt)==numel(hs))
        error('Row %d: visit arrays have mismatched lengths (vt=%d, dates=%d, sz=%d, hs=%d).', ...
            i, numel(vt), numel(dates), numel(sz), numel(hs));
    end

    n_all = numel(vt);
    total_before = total_before + n_all;

    keepMask = ismember(string(vt), allowable_visits);

    if ~any(keepMask)
        R.visit_type(i)        = "[]";
        R.visit_dates_deid(i)  = "[]";
        R.sz_freqs(i)          = "[]";
        R.visit_hasSz(i)       = "[]";
        continue
    end

    vt_f    = cellstr(string(vt(keepMask)));
    dates_f = cellstr(string(dates(keepMask)));
    sz_f    = sz(keepMask);
    hs_f    = hs(keepMask);

    total_after = total_after + numel(vt_f);

    R.visit_type(i)        = string(jsonencode(vt_f));
    R.visit_dates_deid(i)  = string(jsonencode(dates_f));
    R.sz_freqs(i)          = string(jsonencode(sz_f));
    R.visit_hasSz(i)       = string(jsonencode(hs_f));
end

fprintf('[Visit-type filter] Total clinic visits before filter: %d\n', total_before);
fprintf('[Visit-type filter] Total clinic visits after filter:  %d (kept %.1f%%)\n', ...
    total_after, 100*total_after/max(1,total_before)); % expect to have kept 92%
end

function [S, R] = filter_outpatient_routine(S, R, durCol, MAX_ROUTINE_HOURS)
nR0 = height(R);
nS0 = height(S);

% outpatient if acquired on an spe or radnor machine
acqStr = lower(strtrim(string(R.acquired_on)));
isOutpt_site  = contains(acqStr,"spe") | contains(acqStr,"radnor");

% outpatient if epic report says it's outpatient
classStr = lower(strtrim(string(R.report_PATIENT_CLASS)));
isOutpt_class = (classStr == "outpatient");

% outpatient if Jay report say it's outpatient
jayStr = lower(strtrim(string(R.jay_in_or_out)));
isOutpt_jay = (jayStr == "out");

OutptKeys = unique(R(isOutpt_site | isOutpt_class | isOutpt_jay, {'patient_id','session_number'}));
OutptKeys.Properties.VariableNames = {'Patient','Session'};
if isempty(OutptKeys)
    error('No outpatient sessions identified by site/class/jay flags.');
end

% require it's not too long (remove ambulatories)
isRoutine = isfinite(S.(durCol)) & S.(durCol) <= MAX_ROUTINE_HOURS*3600;
RoutineKeys = unique(S(isRoutine, {'Patient','Session'}));

OutptRoutineKeys = innerjoin(OutptKeys, RoutineKeys, 'Keys', {'Patient','Session'});

S = innerjoin(S, OutptRoutineKeys, 'Keys', {'Patient','Session'});
R = innerjoin(R, OutptRoutineKeys, ...
    'LeftKeys', {'patient_id','session_number'}, ...
    'RightKeys', {'Patient','Session'});

fprintf('[Outpatient+routine] Kept %d/%d spike rows (%.1f%%), %d/%d report rows (%.1f%%)\n', ...
    height(S), nS0, 100*height(S)/max(1,nS0), ...
    height(R), nR0, 100*height(R)/max(1,nR0));
end

function Vuniq = build_visit_level_table_R1(R)
% Builds PV by exploding JSON arrays for each row, then collapses to unique (Patient,VisitDate)
PV = table('Size',[0 4], ...
    'VariableTypes',{'double','datetime','double','double'}, ...
    'VariableNames',{'Patient','VisitDate','Freq','HasSz'});

for j = 1:height(R)
    pid = double(R.patient_id(j));

    ds = strtrim(string(R.visit_dates_deid(j)));
    if ds=="[]" || ds=="", continue; end
    dd = string(jsondecode(char(ds)));
    d  = datetime(dd,'InputFormat','yyyy-MM-dd');

    s = strtrim(string(R.sz_freqs(j)));
    s = regexprep(s,'null','NaN','ignorecase');
    v = double(jsondecode(char(s)));
    v(~isfinite(v)) = NaN;
    v(v<0) = NaN;

    hs = strtrim(string(R.visit_hasSz(j)));
    h  = double(jsondecode(char(hs)));
    h(h==2) = NaN;

    if ~(numel(d)==numel(v) && numel(d)==numel(h))
        error('Row %d (patient %g): visit arrays mismatched after filtering.', j, pid);
    end

    n = numel(d);
    if n==0, continue; end

    PV = [PV; table(repmat(pid,n,1), d(:), v(:), h(:), ...
        'VariableNames', PV.Properties.VariableNames)]; %#ok<AGROW>
end

% Collapse duplicates
[gv, pid_keys, date_keys] = findgroups(PV.Patient, PV.VisitDate);

Freq_agg = splitapply(@(x) mean(x(isfinite(x)),'omitnan'), PV.Freq, gv);
Has_agg  = splitapply(@max_hasSz, PV.HasSz, gv);

Vuniq = table(pid_keys, date_keys, Freq_agg, Has_agg, ...
    'VariableNames', {'Patient','VisitDate','Freq','HasSz'});

% Rule 1: if HasSz==0 and Freq missing -> set Freq_R1=0
Vuniq.Freq_R1 = Vuniq.Freq;
mask_rule1 = ~isfinite(Vuniq.Freq_R1) & (Vuniq.HasSz==0);
Vuniq.Freq_R1(mask_rule1) = 0;
end

function SzP = build_patient_seizure_metrics(Vuniq)
[gp, pids] = findgroups(Vuniq.Patient);
MeanSzFreq = splitapply(@(x) mean(x,'omitnan'), Vuniq.Freq_R1, gp);
FracVisits_HasSz1 = splitapply(@local_frac_hasSz1, Vuniq.HasSz, gp);
SzP = table(pids, MeanSzFreq, FracVisits_HasSz1, ...
    'VariableNames', {'Patient','MeanSzFreq','FracVisits_HasSz1'});
end

function T = build_patient_typing_from_report(R, canonical3)
% Strict per-patient typing (errors on conflicts)
pid = double(R.patient_id);

etype = strtrim(string(R.epilepsy_type));
espec = strtrim(string(R.epilepsy_specific));

% Keep non-missing rows
okT = ~ismissing(etype) & strlength(etype)>0;
okS = ~ismissing(espec) & strlength(espec)>0;

% Collapse epilepsy_type
Ttype = table(pid(okT), etype(okT), 'VariableNames', {'Patient','EpilepsyType_raw'});
Ttype = sortrows(Ttype, 'Patient');
[uid, ~, g] = unique(Ttype.Patient, 'stable');
etype_one = strings(numel(uid),1);
for k=1:numel(uid)
    vals = unique(Ttype.EpilepsyType_raw(g==k));
    vals = vals(strlength(vals)>0);
    if isempty(vals), etype_one(k)=""; continue; end
    if numel(vals) > 1
        error('Conflicting epilepsy_type for Patient %g: %s', uid(k), strjoin(vals,", "));
    end
    etype_one(k) = vals(1);
end

% Collapse epilepsy_specific
Tspec = table(pid(okS), espec(okS), 'VariableNames', {'Patient','EpilepsySpecific_raw'});
Tspec = sortrows(Tspec, 'Patient');
[uidS, ~, gS] = unique(Tspec.Patient, 'stable');
espec_one = strings(numel(uidS),1);
for k=1:numel(uidS)
    vals = unique(Tspec.EpilepsySpecific_raw(gS==k));
    vals = vals(strlength(vals)>0);
    if isempty(vals), espec_one(k)=""; continue; end
    if numel(vals) > 1
        % You can relax this if you want, but conflicts here are usually real.
        error('Conflicting epilepsy_specific for Patient %g: %s', uidS(k), strjoin(vals,", "));
    end
    espec_one(k) = vals(1);
end

T = table(uid, etype_one, 'VariableNames', {'Patient','EpilepsyType'});
T = outerjoin(T, table(uidS, espec_one, 'VariableNames', {'Patient','EpilepsySpecific'}), ...
    'Keys','Patient','MergeKeys',true);

% Map to canonical3 using your rules
spec_norm = lower(strtrim(string(T.EpilepsySpecific)));
type_norm = lower(strtrim(string(T.EpilepsyType)));

E3 = strings(height(T),1);
E3(contains(spec_norm,"temporal")) = "Temporal";
E3(contains(spec_norm,"frontal"))  = "Frontal";
E3((strlength(E3)==0) & contains(type_norm,"general")) = "General";

T.EpiType3 = categorical(E3, canonical3);
end

%% =====================================================================
%% ======================= VIEW BUNDLE (CLEAN) =========================
%% =====================================================================

function Views = build_filtered_view(SessionsFiltered, ReportIn, PatientTypingAll, SzFreqPerPatient, ...
                                     NESD_LABEL, badTypes, canonical3)

SessKeys = unique(SessionsFiltered(:,{'Patient','Session'}));
ReportForKeptSessions = innerjoin(ReportIn, SessKeys, ...
    'LeftKeys', {'patient_id','session_number'}, ...
    'RightKeys', {'Patient','Session'});

ReportForKeptSessions.Patient = double(ReportForKeptSessions.patient_id);
ReportForKeptSessions.Session = double(ReportForKeptSessions.session_number);

PatientsKept = unique(SessionsFiltered.Patient);
TypingFiltered = innerjoin(PatientTypingAll, table(PatientsKept,'VariableNames',{'Patient'}), 'Keys','Patient');

% Patient-level mean spike rate
[gpS, pidsS] = findgroups(SessionsFiltered.Patient);
MeanSpikeRate_perHour = splitapply(@(x) mean(x,'omitnan'), SessionsFiltered.SpikeRate_perHour, gpS);
PatientLevelSpikeRates = table(pidsS, MeanSpikeRate_perHour, 'VariableNames', {'Patient','MeanSpikeRate_perHour'});
PatientLevelSpikeRates = innerjoin(PatientLevelSpikeRates, TypingFiltered, 'Keys','Patient');

% Epilepsy vs NESD masks
etype_norm = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));
IsNESDMask = (etype_norm == lower(strtrim(NESD_LABEL)));
IsBadType  = ismember(etype_norm, badTypes) | ismissing(etype_norm) | (strlength(etype_norm)==0);
IsEpilepsyMask = ~IsNESDMask & ~IsBadType;

% Session-level
SessionLevelSpikeRates = SessionsFiltered(:, {'Patient','Session','SpikeRate_perHour'});
SessionLevelSpikeRates.Properties.VariableNames{'SpikeRate_perHour'} = 'SpikesPerHour';

% Spearman tables
SzFreqFiltered = innerjoin(SzFreqPerPatient, table(PatientsKept,'VariableNames',{'Patient'}), 'Keys','Patient');
EpPatients = table(PatientLevelSpikeRates.Patient(IsEpilepsyMask), 'VariableNames', {'Patient'});
SzFreqEpilepsy = innerjoin(SzFreqFiltered, EpPatients, 'Keys','Patient');

PatientSpikeSz_All = innerjoin( ...
    PatientLevelSpikeRates(IsEpilepsyMask, {'Patient','MeanSpikeRate_perHour'}), ...
    SzFreqEpilepsy, 'Keys','Patient');

keepAll = isfinite(PatientSpikeSz_All.MeanSpikeRate_perHour) & isfinite(PatientSpikeSz_All.MeanSzFreq);
PatientSpikeSz_All = PatientSpikeSz_All(keepAll,:);

keepCanon3 = ~ismissing(TypingFiltered.EpiType3) & ismember(TypingFiltered.EpiType3, canonical3);
SzFreqCanon = innerjoin(SzFreqEpilepsy, TypingFiltered(keepCanon3, {'Patient','EpiType3'}), 'Keys','Patient');

PatientSpikeSz_Typed = innerjoin( ...
    PatientLevelSpikeRates(IsEpilepsyMask, {'Patient','MeanSpikeRate_perHour'}), ...
    SzFreqCanon, 'Keys','Patient');

keepTyped = isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perHour) & isfinite(PatientSpikeSz_Typed.MeanSzFreq) & ~ismissing(PatientSpikeSz_Typed.EpiType3);
PatientSpikeSz_Typed = PatientSpikeSz_Typed(keepTyped,:);

% Canonical3 group stats for Fig1C
EpiType4 = strings(height(PatientLevelSpikeRates),1);
et_low = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));
es_low = lower(strtrim(string(PatientLevelSpikeRates.EpilepsySpecific)));

isGeneralType = contains(et_low, "general");
isTemporal    = contains(es_low, "temporal");
isFrontal     = contains(es_low, "frontal");

EpiType4(isGeneralType)                            = "General";
EpiType4(~isGeneralType & isTemporal)              = "Temporal";
EpiType4(~isGeneralType & ~isTemporal & isFrontal) = "Frontal";

PatientLevelSpikeRates.EpiType4 = categorical(EpiType4, canonical3);

InCanonical3Mask = IsEpilepsyMask & ismember(PatientLevelSpikeRates.EpiType4, canonical3);
Canonical3_SubsetTable = PatientLevelSpikeRates(InCanonical3Mask, {'EpiType4','MeanSpikeRate_perHour'});

[g3, cats3] = findgroups(Canonical3_SubsetTable.EpiType4);
medVals = splitapply(@(x) median(x,'omitnan'), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
p25Vals = splitapply(@(x) prctile(x,25), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
p75Vals = splitapply(@(x) prctile(x,75), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
nVals   = splitapply(@(x) sum(isfinite(x)), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);

Canonical3_Stats = table(cats3, nVals, medVals, p25Vals, p75Vals, ...
    'VariableNames', {'EpiType4','GroupCount','Median','P25','P75'});

Canonical3_Pairs = ["General","Temporal";
                    "General","Frontal";
                    "Temporal","Frontal"];

p_pair = NaN(3,1);
for i = 1:3
    A = Canonical3_Pairs(i,1);
    B = Canonical3_Pairs(i,2);
    xa = Canonical3_SubsetTable.MeanSpikeRate_perHour(Canonical3_SubsetTable.EpiType4==categorical(A,canonical3));
    xb = Canonical3_SubsetTable.MeanSpikeRate_perHour(Canonical3_SubsetTable.EpiType4==categorical(B,canonical3));
    if nnz(isfinite(xa))>=3 && nnz(isfinite(xb))>=3
        p_pair(i) = ranksum(xa, xb, 'method','approx');
    end
end
p_pair_bonf = min(p_pair*3, 1);

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

%% =====================================================================
%% ======================= FIGURES + TABLES ============================
%% =====================================================================

function [f1, Fig1Stats] = make_fig1_controls(Views, EPS_RATE, Y_ZERO, Y_LIMS, TITLE_Y_OFFSET, nBoot, alpha)
SessionLevel = Views.SessionLevelSpikeRates;
Report = Views.ReportForKeptSessions;
PL = Views.PatientLevelSpikeRates;

n_ep  = nnz(Views.IsEpilepsyMask);
n_nes = nnz(Views.IsNESDMask);

% Panel A: reported spikes present/absent (resolved)
ReportSlim = resolve_reported_spike_status(Report);
JoinA = innerjoin(SessionLevel(:,{'Patient','Session','SpikesPerHour'}), ...
                  ReportSlim, 'Keys', {'Patient','Session'});

x_abs = JoinA.SpikesPerHour(JoinA.ReportStatus=="absent");
x_pre = JoinA.SpikesPerHour(JoinA.ReportStatus=="present");

pA = ranksum(x_abs, x_pre, 'method','approx');
effectA = cliff_delta(x_pre, x_abs);

[med_abs, lo_abs, hi_abs] = bootstrap_median_ci(x_abs, nBoot, alpha);
[med_pre, lo_pre, hi_pre] = bootstrap_median_ci(x_pre, nBoot, alpha);

Y_A = add_y_jitter_eps(to_log10_per_hour([x_abs; x_pre], EPS_RATE), log10(EPS_RATE), Y_LIMS, 0.02);
G_A = [repmat("Absent", numel(x_abs), 1); repmat("Present", numel(x_pre), 1)];

% Panel B: epilepsy vs NESD (patient means)
x_ep  = PL.MeanSpikeRate_perHour(Views.IsEpilepsyMask);
x_nes = PL.MeanSpikeRate_perHour(Views.IsNESDMask);

pB = ranksum(x_ep, x_nes, 'method','approx');
effectB = cliff_delta(x_ep, x_nes);

[med_ep, lo_ep, hi_ep] = bootstrap_median_ci(x_ep, nBoot, alpha);
[med_nes, lo_nes, hi_nes] = bootstrap_median_ci(x_nes, nBoot, alpha);

Y_B = add_y_jitter_eps([to_log10_per_hour(x_ep, EPS_RATE); to_log10_per_hour(x_nes, EPS_RATE)], ...
    Y_ZERO, Y_LIMS, 0.02);
G_B = [repmat("Epilepsy", nnz(isfinite(x_ep)), 1); repmat("NESD", nnz(isfinite(x_nes)), 1)];

% Panel C: subtype
Sub = Views.Canonical3_SubsetTable;
Y_C = add_y_jitter_eps(to_log10_per_hour(Sub.MeanSpikeRate_perHour, EPS_RATE), Y_ZERO, Y_LIMS, 0.02);

[p_kw, tbl_kw] = kruskalwallis(Sub.MeanSpikeRate_perHour, Sub.EpiType4, 'off');
SS_total = tbl_kw{end,2}; SS_group = tbl_kw{2,2};
eta2_kw = SS_group/SS_total;

% Draw
f1 = figure('Color','w','Position',[60 60 1500 520]);
tiledlayout(f1,1,3,'TileSpacing','compact','Padding','loose');

% A
axA = nexttile; hold(axA,'on'); box(axA,'off'); grid(axA,'on');
boxchart(axA, categorical(G_A), Y_A, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(axA, categorical(G_A), Y_A, 18, 'filled','MarkerFaceAlpha',0.18);
yline(axA, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axA, Y_LIMS);
ylabel(axA,'Spikes/hour (log scale)');
set_log10_ticks(axA,'y',EPS_RATE,Y_LIMS);
add_median_ci_overlay(axA, 1, med_abs, lo_abs, hi_abs, EPS_RATE);
add_median_ci_overlay(axA, 2, med_pre, lo_pre, hi_pre, EPS_RATE);
t = title(axA,'A. Reported presence or absence of spikes');
add_sigbar(axA,1,2, Y_LIMS(2)-0.08*range(Y_LIMS), p_label(pA));
set(axA,'FontSize',20);
t.Units='normalized'; t.Position(2)=t.Position(2)+TITLE_Y_OFFSET;

labelsA = string(axA.XTickLabel);  
labelsA(labelsA=="Absent") = sprintf('Absent (N=%d)', nnz(isfinite(x_abs)));
labelsA(labelsA=="Present")     = sprintf('Present (N=%d)', nnz(isfinite(x_pre)));  
axA.XTickLabel = labelsA;
axA.XTickLabelRotation = 20;   

% B
axB = nexttile; hold(axB,'on'); box(axB,'off'); grid(axB,'on');
boxchart(axB, categorical(G_B), Y_B, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(axB, categorical(G_B), Y_B, 18, 'filled','MarkerFaceAlpha',0.18);
yline(axB, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axB, Y_LIMS);
ylabel(axB,'Spikes/hour (log scale)');
set_log10_ticks(axB,'y',EPS_RATE,Y_LIMS);
add_median_ci_overlay(axB, 1, med_ep, lo_ep, hi_ep, EPS_RATE);
add_median_ci_overlay(axB, 2, med_nes, lo_nes, hi_nes, EPS_RATE);
t = title(axB,'B. Epilepsy versus NESD');
add_sigbar(axB,1,2, Y_LIMS(2)-0.08*range(Y_LIMS), p_label(pB));
set(axB,'FontSize',20);
t.Units='normalized'; t.Position(2)=t.Position(2)+TITLE_Y_OFFSET;
labelsB = string(axB.XTickLabel); 
labelsB(labelsB=="Epilepsy") = sprintf('Epilepsy (N=%d)', n_ep);
labelsB(labelsB=="NESD")     = sprintf('NESD (N=%d)', n_nes);  
axB.XTickLabel = labelsB;
axB.XTickLabelRotation = 20;  

% C
axC = nexttile; hold(axC,'on'); box(axC,'off'); grid(axC,'on');
boxchart(axC, Sub.EpiType4, Y_C, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(axC, Sub.EpiType4, Y_C, 18, 'filled','MarkerFaceAlpha',0.18);
yline(axC, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axC, Y_LIMS);
ylabel(axC,'Spikes/hour (log scale)');
set_log10_ticks(axC,'y',EPS_RATE,Y_LIMS);

% --- Add median + bootstrap CI overlays for Fig 1C (per subtype) ---
cats = categories(Sub.EpiType4);   % e.g., {'General','Temporal','Frontal'} in axis order
for k = 1:numel(cats)
    xg = Sub.MeanSpikeRate_perHour(Sub.EpiType4 == cats{k});
    [mg, log, hig] = bootstrap_median_ci(xg, nBoot, alpha);
    add_median_ci_overlay(axC, k, mg, log, hig, EPS_RATE);  % converts to log10 internally
end


t = title(axC, sprintf('C. Epilepsy subtype'));
set(axC,'FontSize',20);
t.Units='normalized'; t.Position(2)=t.Position(2)+TITLE_Y_OFFSET;

SubStats = Views.Canonical3_Stats;

labelsC = string(axC.XTickLabel);

for i = 1:height(SubStats)
    lab   = string(SubStats.EpiType4(i));
    nHere = SubStats.GroupCount(i);
    labelsC(labelsC == lab) = sprintf('%s (N=%d)', lab, nHere);
end

axC.XTickLabel = labelsC;
axC.XTickLabelRotation = 20;


% --- Pairwise significance bars for Fig 1C (Bonferroni-corrected) ---

SubtypePairs = Views.Canonical3_Pairs;      % e.g. ["General","Temporal"; ...]
p_pair_bonf  = Views.PvalsPairwiseBonf;     % length = nPairs
cats3 = categories(Sub.EpiType4);          % returns cellstr like {'General','Temporal','Frontal'}
cats3 = categorical(string(cats3));        % make it categorical for comparison

yTop  = Y_LIMS(2);
yStep = 0.08 * range(Y_LIMS);   % vertical spacing between bars
y0    = yTop - 0.05 * range(Y_LIMS);  % starting height

for i = 1:size(SubtypePairs,1)

    A = SubtypePairs(i,1);
    B = SubtypePairs(i,2);

    % x positions based on category order
    x1 = find(cats3 == categorical(string(A)));
    x2 = find(cats3 == categorical(string(B)));


    % Bonferroni-corrected p-value → label
    pval = p_pair_bonf(i);
    if isnan(pval)
        continue
    elseif pval < 1e-3
        lab = "***";
    elseif pval < 1e-2
        lab = "**";
    elseif pval < 5e-2
        lab = "*";
    else
        lab = "ns";  
    end

    % draw bar
    add_sigbar(axC, x1, x2, y0 - (i-1)*yStep, lab);
end

% ---------- stats bundle for HTML ----------
Fig1Stats = struct();

% Panel A
Fig1Stats.p_rankSum_A = pA;
Fig1Stats.effectA_cliff = effectA;
Fig1Stats.m_pre = med_pre; Fig1Stats.lo_pre = lo_pre; Fig1Stats.hi_pre = hi_pre;
Fig1Stats.m_abs = med_abs; Fig1Stats.lo_abs = lo_abs; Fig1Stats.hi_abs = hi_abs;

% Panel B
Fig1Stats.p_rankSum_B = pB;
Fig1Stats.effectB_cliff = effectB;
Fig1Stats.m_ep  = med_ep;  Fig1Stats.lo_ep  = lo_ep;  Fig1Stats.hi_ep  = hi_ep;
Fig1Stats.m_nes = med_nes; Fig1Stats.lo_nes = lo_nes; Fig1Stats.hi_nes = hi_nes;

% Panel C omnibus
Fig1Stats.p_kw_C = p_kw;
Fig1Stats.eta2_kw_C = eta2_kw;

% Panel C pairwise (already in Views)
Fig1Stats.p_pair_bonf = Views.PvalsPairwiseBonf;

% Panel C subtype medians + CI (bootstrap medians per group)
Sub = Views.Canonical3_SubsetTable;
subNames = string(categories(Sub.EpiType4));
subMed = nan(numel(subNames),1);
subCI  = nan(numel(subNames),3); % [median lo hi]
for k=1:numel(subNames)
    x = Sub.MeanSpikeRate_perHour(Sub.EpiType4 == subNames(k));
    [m, lo, hi] = bootstrap_median_ci(x, nBoot, alpha);
    subMed(k) = m;
    subCI(k,:) = [m lo hi];
end

% store in same “shape” you used before (median + CI)
Fig1Stats.SubtypeStatsTable = table(subNames, subMed, subCI(:,2), subCI(:,3), ...
    'VariableNames', {'Group','Median','CI_lo','CI_hi'});


end

function fS2 = make_figS2_sz_by_reported_spikes(Views, SzFreqPerPatient, nBoot, alpha, xLims_log10)
RS = resolve_reported_spike_status(Views.ReportForKeptSessions);

[gp, pid] = findgroups(RS.Patient);
hasPresent = splitapply(@(x) any(x=="present"), RS.ReportStatus, gp);
hasAbsent  = splitapply(@(x) any(x=="absent"),  RS.ReportStatus, gp);
RptP = table(pid, hasPresent, hasAbsent, 'VariableNames', {'Patient','HasPresent','HasAbsent'});

S2 = innerjoin(SzFreqPerPatient, RptP, 'Keys','Patient');

% epilepsy only
EpPatients = Views.PatientLevelSpikeRates.Patient(Views.IsEpilepsyMask);
S2 = innerjoin(S2, table(EpPatients,'VariableNames',{'Patient'}), 'Keys','Patient');

mask_allAbsent  = S2.HasAbsent & ~S2.HasPresent;
mask_anyPresent = S2.HasPresent;

freq_allAbsent  = S2.MeanSzFreq(mask_allAbsent);
freq_anyPresent = S2.MeanSzFreq(mask_anyPresent);

freq_allAbsent  = freq_allAbsent(isfinite(freq_allAbsent));
freq_anyPresent = freq_anyPresent(isfinite(freq_anyPresent));
n_allAbsent  = numel(freq_allAbsent);
n_anyPresent = numel(freq_anyPresent);


p = ranksum(freq_allAbsent, freq_anyPresent, 'method','approx');

EPS_FREQ = 1e-3;
Y_ZERO = log10(EPS_FREQ);

Y = [to_log10_per_month(freq_allAbsent, EPS_FREQ); to_log10_per_month(freq_anyPresent, EPS_FREQ)];
Y = add_y_jitter_eps(Y, Y_ZERO, xLims_log10, 0.02);
G = [repmat("All EEGs: no spikes", numel(freq_allAbsent), 1);
     repmat("≥1 EEG: spikes present", numel(freq_anyPresent), 1)];

[med1, lo1, hi1] = bootstrap_median_ci(freq_allAbsent, nBoot, alpha);
[med2, lo2, hi2] = bootstrap_median_ci(freq_anyPresent, nBoot, alpha);

fS2 = figure('Color','w','Position',[100 100 800 520]);
ax = axes(fS2); hold(ax,'on'); box(ax,'off'); grid(ax,'on');

boxchart(ax, categorical(G), Y, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(ax, categorical(G), Y, 18, 'filled','MarkerFaceAlpha',0.18);
yline(ax, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(ax, xLims_log10);
ylabel(ax,'Seizures/month (log scale)');
set_log10_ticks(ax,'y',EPS_FREQ,xLims_log10);

add_median_ci_overlay_month(ax,1,med1,lo1,hi1,EPS_FREQ);
add_median_ci_overlay_month(ax,2,med2,lo2,hi2,EPS_FREQ);

title(ax,'Mean seizure frequency by reported spikes across EEGs');
add_sigbar(ax,1,2, xLims_log10(2)-0.08*range(xLims_log10), p_label(p));
set(ax,'FontSize',20);

% ---- Add Ns to x tick labels ----
labels = string(ax.XTickLabel);

labels(labels=="All EEGs: no spikes") = ...
    sprintf('All EEGs: no spikes (N=%d)', n_allAbsent);

labels(labels=="≥1 EEG: spikes present") = ...
    sprintf('≥1 EEG: spikes present (N=%d)', n_anyPresent);

ax.XTickLabel = labels;
ax.XTickLabelRotation = 20;

end

function Table1_flat = build_table1_flat(Views, SzFreqPerPatient, Vuniq, EPS_RATE, nBoot, alpha)
Rk = Views.ReportForKeptSessions;
PL = Views.PatientLevelSpikeRates;

AllPatients = PL.Patient;
N_total = numel(AllPatients);

% Age
birth_str = strtrim(string(Rk.deid_birth_date));
isMiss = birth_str=="" | birth_str=="null" | birth_str=="[null]";
birth_dt = NaT(size(birth_str));
birth_dt(~isMiss) = datetime(birth_str(~isMiss),'InputFormat','yyyy-MM-dd');
refDate = datetime(2000,1,1);
age_first = NaN(size(birth_dt));
ok = ~isnat(birth_dt);
age_first(ok) = days(refDate - birth_dt(ok))/365.25;

[ga, pidA] = findgroups(Rk.Patient);
AgeFirst = splitapply(@local_min_omitnan, age_first, ga);
AgeTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), ...
                     table(pidA,AgeFirst,'VariableNames',{'Patient','AgeFirst'}), 'Keys','Patient');
age_vec = AgeTable.AgeFirst;
age_med = median(age_vec,'omitnan');
age_q = prctile(age_vec,[25,75]);

% Sex (first nonmissing)
sex_raw = upper(strtrim(string(Rk.nlp_gender)));
[gs, pidS] = findgroups(Rk.Patient);
sex_per = splitapply(@local_first_nonmissing, sex_raw, gs);
SexTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), ...
                     table(pidS,sex_per,'VariableNames',{'Patient','SexCode'}), 'Keys','Patient');
n_f = nnz(SexTable.SexCode=="F");
n_m = nnz(SexTable.SexCode=="M");
n_u = N_total - n_f - n_m;

% Dx
n_epi = nnz(Views.IsEpilepsyMask);
n_pnes = nnz(Views.IsNESDMask);
n_dx_unk = N_total - n_epi - n_pnes;

% Subtype (among epilepsy)
E4 = string(PL.EpiType4);
isTemp  = (E4=="Temporal") & Views.IsEpilepsyMask;
isFront = (E4=="Frontal")  & Views.IsEpilepsyMask;
isGen   = (E4=="General")  & Views.IsEpilepsyMask;
typing = strtrim(string(PL.EpilepsySpecific));
isSubUnk = Views.IsEpilepsyMask & ~(isTemp|isFront|isGen) & ...
    (typing=="Unknown or MRN not found" | typing=="Unclassified or Unspecified" | typing=="");
isOther = Views.IsEpilepsyMask & ~(isTemp|isFront|isGen) & ~isSubUnk;

n_temp = nnz(isTemp); n_front = nnz(isFront); n_gen = nnz(isGen);
n_other = nnz(isOther); n_subunk = nnz(isSubUnk);

% Visits per patient
V_all = Vuniq(:,{'Patient','VisitDate'});
[gv, pidV] = findgroups(V_all.Patient);
nVisits = splitapply(@(x) numel(unique(x)), V_all.VisitDate, gv);
VisitsTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), ...
    table(pidV,nVisits,'VariableNames',{'Patient','NumVisits'}), 'Keys','Patient');
vis_med = median(VisitsTable.NumVisits,'omitnan');
vis_q = prctile(VisitsTable.NumVisits,[25,75]);

% EEGs per patient
Sess = Views.SessionsForFigures;
[ge, pidE] = findgroups(Sess.Patient);
nEEG = splitapply(@(x) numel(unique(x)), Sess.Session, ge);
EEGTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), ...
    table(pidE,nEEG,'VariableNames',{'Patient','NumEEG'}), 'Keys','Patient');
eeg_med = median(EEGTable.NumEEG,'omitnan');
eeg_q = prctile(EEGTable.NumEEG,[25,75]);

% Mean seizure freq
SzJ = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), SzFreqPerPatient, 'Keys','Patient');
sf_vec = SzJ.MeanSzFreq; sf_vec = sf_vec(isfinite(sf_vec));
sf_med = median(sf_vec,'omitnan'); sf_q = prctile(sf_vec,[25,75]);
[~, sf_lo, sf_hi] = bootstrap_median_ci(sf_vec, nBoot, alpha);

% Mean spike rate
sr_vec = PL.MeanSpikeRate_perHour; sr_vec = sr_vec(isfinite(sr_vec));
sr_med = median(sr_vec,'omitnan'); sr_q = prctile(sr_vec,[25,75]);
[~, sr_lo, sr_hi] = bootstrap_median_ci(sr_vec, nBoot, alpha);

% Reported spikes per EEG
RS_all = resolve_reported_spike_status(Views.ReportForKeptSessions);
status = string(RS_all.ReportStatus);
status(ismissing(RS_all.ReportStatus)) = "unknown";
RS_all.SpikeStatus = categorical(status, ["absent","present","unknown"]);
n_rep_pre = nnz(RS_all.SpikeStatus=="present");
n_rep_abs = nnz(RS_all.SpikeStatus=="absent");
n_rep_unk = nnz(RS_all.SpikeStatus=="unknown");
n_eegs_all = height(RS_all);

% Build flattened 2-col table
OutVar = {};
OutStat = {};

push = @(v,s) assignin('caller','OutVar',[OutVar; {v}]); %#ok<NASGU>
% (Use explicit to avoid workspace tricks)
OutVar{end+1,1} = "Total N patients with ≥1 outpatient routine EEG";
OutStat{end+1,1} = sprintf('%d', N_total);

OutVar{end+1,1} = "Age at first visit (years)";
OutStat{end+1,1} = sprintf('%.1f (%.1f-%.1f)', age_med, age_q(1), age_q(2));

OutVar{end+1,1} = "Sex"; OutStat{end+1,1} = "";
OutVar{end+1,1} = "    Women"; OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_f, 100*n_f/N_total);
OutVar{end+1,1} = "    Men";   OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_m, 100*n_m/N_total);
OutVar{end+1,1} = "    Unknown/Other"; OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_u, 100*n_u/N_total);

OutVar{end+1,1} = "Epilepsy diagnosis"; OutStat{end+1,1} = "";
OutVar{end+1,1} = "    Epilepsy"; OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_epi, 100*n_epi/N_total);
OutVar{end+1,1} = "    PNES";     OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_pnes, 100*n_pnes/N_total);
OutVar{end+1,1} = "    Unknown";  OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_dx_unk, 100*n_dx_unk/N_total);

OutVar{end+1,1} = "Epilepsy subtype"; OutStat{end+1,1} = "";
OutVar{end+1,1} = "    Temporal lobe"; OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_temp, 100*n_temp/max(1,n_epi));
OutVar{end+1,1} = "    Frontal lobe";  OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_front, 100*n_front/max(1,n_epi));
OutVar{end+1,1} = "    Generalized";   OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_gen, 100*n_gen/max(1,n_epi));
OutVar{end+1,1} = "    Other";         OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_other, 100*n_other/max(1,n_epi));
OutVar{end+1,1} = "    Unknown";       OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_subunk, 100*n_subunk/max(1,n_epi));

OutVar{end+1,1} = "Number of clinic visits";
OutStat{end+1,1} = sprintf('%.1f (%.1f-%.1f)', vis_med, vis_q(1), vis_q(2));

OutVar{end+1,1} = "Number of EEGs";
OutStat{end+1,1} = sprintf('%.1f (%.1f-%.1f)', eeg_med, eeg_q(1), eeg_q(2));

OutVar{end+1,1} = "Mean seizure frequency (seizures/month)";
OutStat{end+1,1} = sprintf('%.2f (%.2f-%.2f); median CI [%.2f-%.2f]', sf_med, sf_q(1), sf_q(2), sf_lo, sf_hi);

OutVar{end+1,1} = "Mean spike rate (spikes/hour)";
OutStat{end+1,1} = sprintf('%.2f (%.2f-%.2f); median CI [%.2f-%.2f]', sr_med, sr_q(1), sr_q(2), sr_lo, sr_hi);

OutVar{end+1,1} = "Reported spikes"; OutStat{end+1,1} = "";
OutVar{end+1,1} = "    Present"; OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_rep_pre, 100*n_rep_pre/max(1,n_eegs_all));
OutVar{end+1,1} = "    Absent";  OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_rep_abs, 100*n_rep_abs/max(1,n_eegs_all));
OutVar{end+1,1} = "    Unknown"; OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_rep_unk, 100*n_rep_unk/max(1,n_eegs_all));

Table1_flat = table(string(OutVar), string(OutStat), 'VariableNames', {'Variable','Statistic'});
end

function [fS3, PairedStats] = make_figS3_paired_and_sensitivity(Views, Vuniq, yLims_log10, nBoot, alpha)
% Uses Vuniq (has VisitDate + Freq_R1) and session-level spikes+dates from report.
gapVec = [15 30 60];
nG = numel(gapVec);

SessRates = Views.SessionLevelSpikeRates;   % Patient, Session, SpikesPerHour
Rk = Views.ReportForKeptSessions;

% Parse EEG dates
EEG_raw = Rk.start_time_deid;
if isdatetime(EEG_raw)
    EEG_dt = EEG_raw;
else
    EEG_dt = datetime(strtrim(string(EEG_raw)), 'InputFormat',"yyyy-MM-dd'T'HH:mm:ss");
end
SessionDates = table(double(Rk.Patient), double(Rk.Session), EEG_dt, ...
    'VariableNames', {'Patient','Session','EEG_Date'});
SessWithDate = innerjoin(SessRates, SessionDates, 'Keys', {'Patient','Session'});

Vuniq_R1 = Vuniq(:,{'Patient','VisitDate','Freq_R1'});

EpPatients = Views.PatientLevelSpikeRates.Patient(Views.IsEpilepsyMask);

min_abs_diff_spikes = 0;

% Diagnostic where we are losing patients
Diag = paired_closest_visit_diag(SessWithDate, Vuniq_R1, EpPatients, 180, 180);


effectMat = nan(nG,nG);
pMat = nan(nG,nG);
nMat = nan(nG,nG);

for i=1:nG
    for j=1:nG
        [r_eff, p_signed, ~,nPairs] = run_paired_analysis_core( ...
            SessWithDate, Vuniq_R1, EpPatients, gapVec(i), gapVec(j), min_abs_diff_spikes);
        effectMat(i,j) = r_eff;
        pMat(i,j) = p_signed;
        nMat(i,j) = nPairs;
    end
end

% FDR
pVec = pMat(:);
valid = isfinite(pVec);
qVec = nan(size(pVec));
qVec(valid) = mafdr(pVec(valid), 'BHFDR', true);
qMat = reshape(qVec, size(pMat));

% Primary window 30/30
[r_eff_primary, p_primary, stats_signed, nPairs_primary, Freq_low, Freq_high] = run_paired_analysis_core( ...
    SessWithDate, Vuniq_R1, EpPatients, 30, 30, min_abs_diff_spikes);

% ---- Bundle paired stats for HTML ----
[medLow,  loLow,  hiLow]  = bootstrap_median_ci(Freq_low(isfinite(Freq_low)),  nBoot, alpha);
[medHigh, loHigh, hiHigh] = bootstrap_median_ci(Freq_high(isfinite(Freq_high)), nBoot, alpha);

PairedStats = struct();
PairedStats.nPairs = nPairs_primary;
PairedStats.p_signed = p_primary;
PairedStats.stats_signed = stats_signed;
PairedStats.r_effect = r_eff_primary;

PairedStats.medLow  = medLow;  PairedStats.loLow  = loLow;  PairedStats.hiLow  = hiLow;
PairedStats.medHigh = medHigh; PairedStats.loHigh = loHigh; PairedStats.hiHigh = hiHigh;



EPS_FREQ = 1e-3;
low_log  = to_log10_per_month(Freq_low,  EPS_FREQ);
high_log = to_log10_per_month(Freq_high, EPS_FREQ);
Y_ZERO = log10(EPS_FREQ);

fS3 = figure('Color','w','Position',[80 80 1400 550]);
tl = tiledlayout(fS3,1,2,'TileSpacing','compact','Padding','compact');

% Left: paired plot
axP = nexttile(tl,1); hold(axP,'on'); box(axP,'off'); grid(axP,'on');
for k=1:nPairs_primary
    plot(axP, [1 2], [low_log(k), high_log(k)], '-', 'Color',[0.7 0.7 0.7], 'LineWidth',1.2);
end
yline(axP, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

x_low  = 1 + (rand(nPairs_primary,1)-0.5)*0.04;
x_high = 2 + (rand(nPairs_primary,1)-0.5)*0.04;
scatter(axP, x_low,  low_log,  30, 'filled', 'MarkerFaceAlpha',0.6);
scatter(axP, x_high, high_log, 30, 'filled', 'MarkerFaceAlpha',0.6);

xlim(axP,[0.5 2.5]);
xticks(axP,[1 2]);
xticklabels(axP,{'Low-spike EEG','High-spike EEG'});
ylabel(axP,'Seizures/month (log scale)');
ylim(axP, [min([low_log;high_log])-0.3, max([low_log;high_log])+0.3]);
set_log10_ticks(axP,'y',EPS_FREQ, axP.YLim);

yDataMax = max([low_log(:); high_log(:)], [], 'omitnan');
yl = ylim(axP);
pad = 0.12 * range(yl);

if yDataMax + pad > yl(2)
    ylim(axP, [yl(1), yDataMax + pad]);
end

yBar = yDataMax + 0.06*range(ylim(axP));

add_sigbar(axP, 1, 2, yBar, sprintf('%s, r=%.2f', p_label(p_primary), r_eff_primary));


title(axP, sprintf('Paired seizure frequencies (N=%d)', nPairs_primary), 'FontSize',20);
set(axP,'FontSize',20);

% Right: sensitivity heatmap
axH = nexttile(tl,2);
imagesc(axH, effectMat);
axis(axH,'equal'); axis(axH,'tight');
colormap(axH, parula);
caxis(axH, [-0.5 0.5]);
cb = colorbar(axH);
ylabel(cb,'Effect size (Wilcoxon r)','FontSize',20);

xlabel(axH,'Max days after EEG','FontSize',20);
ylabel(axH,'Max days before EEG','FontSize',20);
axH.XTick = 1:nG; axH.YTick = 1:nG;
axH.XTickLabel = gapVec; axH.YTickLabel = gapVec;
title(axH,'Sensitivity to visit–EEG matching window','FontSize',20);

for ii=1:nG
    for jj=1:nG
        if ~isfinite(effectMat(ii,jj)) || ~isfinite(qMat(ii,jj)) || nMat(ii,jj)==0, continue; end
        text(axH, jj, ii, sprintf('N=%d\nr=%.2f', nMat(ii,jj), effectMat(ii,jj)), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',16,'Color','k','FontWeight','bold');
    end
end
set(axH,'FontSize',20);
end

function write_results_html(outPath, Views, SzFreqPerPatient, ...
    Fig1Stats, ...
    SpearmanResults_main, rs_all_main, p_all_main, n_all_main, rho_lo_main, rho_hi_main, subtype_ci_main, ...
    ReportForKeptSessions, ...
    PairedStats)

if ~exist(fileparts(outPath),'dir'), mkdir(fileparts(outPath)); end
fid = fopen(outPath,'w');
if fid==-1, error('Could not open %s', outPath); end

PL = Views.PatientLevelSpikeRates;
AllPatients = PL.Patient;
N_total = numel(AllPatients);

% EEG count (after outpatient+routine filtering)
n_eegs_all = height(ReportForKeptSessions);

% Epilepsy count
n_epi = nnz(Views.IsEpilepsyMask);

% cohort medians + CI
sf_vec = SzFreqPerPatient.MeanSzFreq;
sf_vec = sf_vec(isfinite(sf_vec));
sf_med = median(sf_vec,'omitnan');
[~, sf_ci_lo, sf_ci_hi] = bootstrap_median_ci(sf_vec, 5000, 0.05);

sr_vec = PL.MeanSpikeRate_perHour;
sr_vec = sr_vec(isfinite(sr_vec));
sr_med = median(sr_vec,'omitnan');
[~, sr_ci_lo, sr_ci_hi] = bootstrap_median_ci(sr_vec, 5000, 0.05);

fprintf(fid,'<html><head><meta charset="UTF-8"><title>Results</title></head><body>\n');

%% Table 1 general overview
fprintf(fid, '<h2>Cohort summary</h2>\n');
fprintf(fid,['<p>We included %d patients with at least one outpatient routine EEG '...
    '(%d EEGs). %d patients (%1.1f%%) carried a '...
    'diagnosis of epilepsy. Median [95%% CI] monthly seizure frequency was '...
    '%1.2f [%1.2f-%1.2f], and median [95%% CI] spikes/hour across EEGs '...
    'was %1.2f [%1.2f-%1.2f] (Table 1).</p>'],...
    N_total, n_eegs_all, n_epi, 100*n_epi/max(1,N_total), ...
    sf_med, sf_ci_lo, sf_ci_hi, ...
    sr_med, sr_ci_lo, sr_ci_hi);

%% ---- Figure 1: Control Panels ----
fprintf(fid, '<h2>Spike rates by patient groups</h2>\n');

% Panel A
pA_str = format_p_html(Fig1Stats.p_rankSum_A);
fprintf(fid,['<p>Automatically detected spike rates were higher in EEGs with clinically-reported spikes '...
            '(median spike rate %.2f [95%% CI %.2f-%.2f] spikes/hour) than '...
            'those without reported spikes (%.2f [%.2f-%.2f] '...
            'spikes/hour) (%s, Cliff''s &delta; = %.2f; Fig. 1A). '],...
            Fig1Stats.m_pre, Fig1Stats.lo_pre, Fig1Stats.hi_pre, ...
            Fig1Stats.m_abs, Fig1Stats.lo_abs, Fig1Stats.hi_abs, ...
            pA_str, Fig1Stats.effectA_cliff);

% Panel B
pB_str = format_p_html(Fig1Stats.p_rankSum_B);
fprintf(fid,['Patients with epilepsy also had higher spike rates '...
            '(%.2f [%.2f-%.2f] spikes/hour) than those '...
            'with NESD (%.2f [%.2f-%.2f] spikes/hour) '...
            '(%s, &delta; = %.2f; Fig. 1B). '],...
            Fig1Stats.m_ep, Fig1Stats.lo_ep, Fig1Stats.hi_ep, ...
            Fig1Stats.m_nes, Fig1Stats.lo_nes, Fig1Stats.hi_nes, ...
            pB_str, Fig1Stats.effectB_cliff);

% Panel C omnibus
fprintf(fid,['Spike rates differed across epilepsy subtypes (Kruskal–Wallis '...
            '%s, η² ≈ %.3f). '], format_p_html(Fig1Stats.p_kw_C), Fig1Stats.eta2_kw_C);

% Panel C pairwise narrative (assumes General, Temporal, Frontal exist)
Tsub = Fig1Stats.SubtypeStatsTable;
getRow = @(nm) Tsub(strcmp(string(Tsub.Group), string(nm)), :);

rG = getRow("General");
rT = getRow("Temporal");
rF = getRow("Frontal");

p_pair_bonf = Fig1Stats.p_pair_bonf; % [G vs T; G vs F; T vs F]

fprintf(fid,['Generalized epilepsy demonstrated higher spike rates '...
            '(%.2f [%.2f-%.2f]) than temporal '...
            '(%.2f [%.2f-%.2f]; Bonferroni-adjusted '...
            '%s) and frontal '...
            '(%.2f [%.2f-%.2f]; '...
            '%s). The temporal versus frontal comparison '...
            'was not significant (%s; Fig. 1C).</p>'],...
            rG.Median, rG.CI_lo, rG.CI_hi, ...
            rT.Median, rT.CI_lo, rT.CI_hi, format_p_html(p_pair_bonf(1)), ...
            rF.Median, rF.CI_lo, rF.CI_hi, format_p_html(p_pair_bonf(2)), ...
            format_p_html(p_pair_bonf(3)));

%% ---- Figure 2: Spearman correlations ----
fprintf(fid, '<h2>Relationship between spike rate and seizure frequency</h2>\n');

fprintf(fid,['<p>Across all epilepsy patients with documented seizure frequency (N = %d), spike rate '...
    'and seizure frequency were positively correlated '...
    '(&rho; = %.2f [95%% CI %.2f-%.2f], %s).'], ...
    n_all_main, rs_all_main, rho_lo_main, rho_hi_main, format_p_html(p_all_main));

fprintf(fid,[' Subtype-specific correlations were significant for generalized epilepsy '...
    '(N = %d, &rho; = %.2f [%.2f-%.2f], Bonferroni-adjusted %s) and '...
    'temporal lobe epilepsy (N = %d, &rho; = %.2f [%.2f-%.2f], %s), but not frontal '...
    'lobe epilepsy (N = %d, &rho; = %.2f [%.2f-%.2f], %s; Fig 2A-D). '],...
    SpearmanResults_main.N(1), SpearmanResults_main.Spearman_r(1), ...
    subtype_ci_main.ci_lo(1), subtype_ci_main.ci_hi(1), format_p_html(SpearmanResults_main.p_bonf(1)), ...
    SpearmanResults_main.N(2), SpearmanResults_main.Spearman_r(2), ...
    subtype_ci_main.ci_lo(2), subtype_ci_main.ci_hi(2), format_p_html(SpearmanResults_main.p_bonf(2)), ...
    SpearmanResults_main.N(3), SpearmanResults_main.Spearman_r(3), ...
    subtype_ci_main.ci_lo(3), subtype_ci_main.ci_hi(3), format_p_html(SpearmanResults_main.p_bonf(3)));

fprintf(fid,['Results were consistent when restricting analyses to non-zero data (Fig. S1). ']);
fprintf(fid,['Patients who ever had spikes reported on at least one EEG also had '...
    'higher mean seizure frequencies than patients whose EEGs consistently '...
    'lacked spikes (Fig. S2). ']);
fprintf(fid,['Together, these findings indicate that spike burden reflects seizure severity at the population level. ']);

%% ---- Paired analysis ----
fprintf(fid,['Among %d patients with multiple EEGs-clinic visit pairs, '...
    'seizure frequency did not differ significantly between periods '...
    'near low- versus high-spike-rate EEGs '...
    '(%1.1f [%1.1f-%1.1f)] vs %1.1f [%1.1f-%1.1f] seizures/month; '...
    'W = %1.1f, %s; Fig. S3).</p>'],...
    PairedStats.nPairs, ...
    PairedStats.medLow,  PairedStats.loLow,  PairedStats.hiLow, ...
    PairedStats.medHigh, PairedStats.loHigh, PairedStats.hiHigh, ...
    PairedStats.stats_signed.signedrank, format_p_html(PairedStats.p_signed));

fprintf(fid, '</body></html>\n');
fclose(fid);
end


%% =====================================================================
%% ======================= REPORT STATUS RESOLUTION =====================
%% =====================================================================

function ReportSlim = resolve_reported_spike_status(ReportForKeptSessions)
any_spikes = string(ReportForKeptSessions.report_SPORADIC_EPILEPTIFORM_DISCHARGES);
isMainPresent = (any_spikes=="present");
isMainAbsent  = (any_spikes=="absent");

rawF = lower(strtrim(string(ReportForKeptSessions.jay_focal_epi)));
rawM = lower(strtrim(string(ReportForKeptSessions.jay_multifocal_epi)));
rawG = lower(strtrim(string(ReportForKeptSessions.jay_gen_epi)));

isF_P = rawF=="present"; isF_A = rawF=="absent";
isM_P = rawM=="present"; isM_A = rawM=="absent";
isG_P = rawG=="present"; isG_A = rawG=="absent";

presentJay_any = isF_P | isM_P | isG_P;
allJay_absent  = isF_A & isM_A & isG_A;
allJay_present = isF_P & isM_P & isG_P;

blankMain   = ~(isMainPresent | isMainAbsent);
blankJayAll = ~(isF_P|isF_A) & ~(isM_P|isM_A) & ~(isG_P|isG_A);

disc1 = allJay_absent & isMainPresent;
disc2 = isMainAbsent & allJay_present;
if any(disc1|disc2)
    error('Discordant spike presence between main and jay_* columns.');
end

repCombined = strings(height(ReportForKeptSessions),1);
repCombined(isMainPresent | presentJay_any) = "present";
repCombined(allJay_absent & blankMain) = "absent";
repCombined(isMainAbsent & blankJayAll) = "absent";
repCombined(repCombined=="") = "unknown";

ReportSlim = ReportForKeptSessions(:, {'Patient','Session'});
ReportSlim.ReportStatus = categorical(repCombined, ["absent","present","unknown"]);

end

%% =====================================================================
%% ======================= PAIRED ANALYSIS CORE =========================
%% =====================================================================

function [r_effect, p_signed, stats_signed, nPairs, Freq_low, Freq_high] = run_paired_analysis_core(...
    SessWithDate, Vuniq_R1, EpPatients, maxDaysBefore, maxDaysAfter, min_abs_diff_spikes)

EpPatients = unique(double(EpPatients(:)));
SessWithDate = innerjoin(SessWithDate, table(EpPatients,'VariableNames',{'Patient'}), 'Keys','Patient');
Vuniq_R1     = innerjoin(Vuniq_R1,     table(EpPatients,'VariableNames',{'Patient'}), 'Keys','Patient');

Pairs = table('Size',[0 7], ...
    'VariableTypes', {'double','double','datetime','datetime','double','double','double'}, ...
    'VariableNames', {'Patient','Session','EEG_Date','VisitDate','GapDays','SpikeRate_perHour','Freq_R1'});

for pid = EpPatients'
    S_sub = SessWithDate(SessWithDate.Patient==pid,:);
    V_sub = Vuniq_R1(Vuniq_R1.Patient==pid,:);
    if isempty(S_sub) || isempty(V_sub), continue; end

    for j=1:height(S_sub)
        EEGd = S_sub.EEG_Date(j);
        gapDays = days(EEGd - V_sub.VisitDate);
        keep = (gapDays >= -maxDaysBefore) & (gapDays <= maxDaysAfter);
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

Pairs = Pairs(isfinite(Pairs.Freq_R1) & isfinite(Pairs.SpikeRate_perHour), :);
if isempty(Pairs)
    r_effect = NaN; p_signed = NaN; nPairs = 0; Freq_low=[]; Freq_high=[];
    return;
end

% closest visit per (Patient,Session)
[grpPS, ~, ~] = findgroups(Pairs.Patient, Pairs.Session);
keep = false(height(Pairs),1);
for g=1:max(grpPS)
    idx = find(grpPS==g);
    [~,best] = min(abs(Pairs.GapDays(idx)));
    keep(idx(best)) = true;
end
Pairs = Pairs(keep,:);

% choose low/high spike EEG per patient, require distinct visits
[gp, uP] = findgroups(Pairs.Patient);
Spike_low=nan(numel(uP),1); Spike_high=nan(numel(uP),1);
Freq_low=nan(numel(uP),1);  Freq_high=nan(numel(uP),1);
Visit_low=NaT(numel(uP),1); Visit_high=NaT(numel(uP),1);

for k=1:numel(uP)
    idx = find(gp==k);
    if numel(idx)<2, continue; end
    r = Pairs.SpikeRate_perHour(idx);
    f = Pairs.Freq_R1(idx);
    ok = isfinite(r)&isfinite(f);
    if nnz(ok)<2, continue; end
    idx = idx(ok); r = r(ok); f = f(ok);

    [~,iL] = min(r);
    [~,iH] = max(r);
    if iL==iH, continue; end
    if r(iH)-r(iL) <= min_abs_diff_spikes, continue; end
    if Pairs.VisitDate(idx(iL)) == Pairs.VisitDate(idx(iH)), continue; end

    Spike_low(k)=r(iL); Spike_high(k)=r(iH);
    Freq_low(k)=f(iL);  Freq_high(k)=f(iH);
    Visit_low(k)=Pairs.VisitDate(idx(iL));
    Visit_high(k)=Pairs.VisitDate(idx(iH));
end

valid = isfinite(Spike_low)&isfinite(Spike_high)&isfinite(Freq_low)&isfinite(Freq_high);
Freq_low = Freq_low(valid);
Freq_high = Freq_high(valid);
nPairs = numel(Freq_low);

if nPairs < 1
    r_effect = NaN; p_signed = NaN;
    return;
end

[p_signed, ~, stats_signed] = signrank(Freq_high, Freq_low, 'method','approx');
r_effect = stats_signed.zval / sqrt(nPairs);

end

%% =====================================================================
%% ======================= SPEARMAN FIG FUNCTION =======================
%% =====================================================================
% (This is your function with only minimal “safety tightening”)

function [SpearmanResults, rs_all, p_all, n_all, rho_lo, rho_hi, subtype_ci] = ...
    spearman_plotting_function(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
                              canonical3, spearman_xLims, spearman_yLims, ...
                              fig_out, ...
                              freqFieldName, labelSuffix, nonZeroOnly)
% spearman_plotting_function
%   Makes 2x2 scatter panels:
%     A) All epilepsy (gray)
%     B) Frontal (yellow)
%     C) Temporal (red)
%     D) General (blue)
%
%   Computes Spearman rho + p (overall and by subtype) + bootstrap CI.
%
%   Assumes helpers exist in your file:
%     - bootstrap_spearman_ci(x,y,nBoot,alpha)
%     - set_log10_ticks(ax, whichAxis, eps_val, axisLims, maxPow)

% ----------------------- params -----------------------
nBoot = 5000;
alpha = 0.05;

fontL = 20;

% ---- Fixed palette (A gray, then yellow, red, blue) ----
COL_all   = [0.45 0.45 0.45];   % gray
COL_front = [0.93 0.69 0.13];   % yellow-ish
COL_temp  = [0.85 0.33 0.10];   % red/orange-ish
COL_gen   = [0.00 0.45 0.74];   % blue

% panel order/title (matches your Fig 2 labeling)
panelOrder = {'Frontal','Temporal','General'};
panelTitle = {'B. Frontal','C. Temporal','D. General'};

% ----------------------- overall stats -----------------------
x_all = double(PatientSpikeSz_All.MeanSpikeRate_perHour);
y_all = double(PatientSpikeSz_All.(freqFieldName));

mask_all = isfinite(x_all) & isfinite(y_all);
if nonZeroOnly
    mask_all = mask_all & (x_all > 0) & (y_all > 0);
end
n_all = sum(mask_all);

if n_all >= 3
    [rs_all, p_all] = corr(x_all(mask_all), y_all(mask_all), ...
        'Type','Spearman','Rows','complete');
    [~, rho_lo, rho_hi] = bootstrap_spearman_ci(x_all(mask_all), y_all(mask_all), nBoot, alpha);
else
    rs_all = NaN; p_all = NaN; rho_lo = NaN; rho_hi = NaN;
end

% ----------------------- subtype stats -----------------------
rowsOut = {};
subtype_ci = table(string(canonical3(:)), nan(numel(canonical3),1), nan(numel(canonical3),1), nan(numel(canonical3),1), ...
    'VariableNames', {'Group','rho','ci_lo','ci_hi'});

for ii = 1:numel(canonical3)
    g = canonical3(ii);

    mBase = (PatientSpikeSz_Typed.EpiType3 == g);
    x = double(PatientSpikeSz_Typed.MeanSpikeRate_perHour(mBase));
    y = double(PatientSpikeSz_Typed.(freqFieldName)(mBase));

    mask = isfinite(x) & isfinite(y);
    if nonZeroOnly
        mask = mask & (x > 0) & (y > 0);
    end
    n = sum(mask);

    if n >= 3
        [rs, p] = corr(x(mask), y(mask), 'Type','Spearman','Rows','complete');
        [rho, lo, hi] = bootstrap_spearman_ci(x(mask), y(mask), nBoot, alpha);
    else
        rs = NaN; p = NaN; rho = NaN; lo = NaN; hi = NaN;
    end

    subtype_ci.rho(ii)   = rho;
    subtype_ci.ci_lo(ii) = lo;
    subtype_ci.ci_hi(ii) = hi;

    rowsOut(end+1,:) = {char(g), n, rs, p}; %#ok<AGROW>
end

SpearmanResults = cell2table(rowsOut, ...
    'VariableNames', {'Group','N','Spearman_r','p_raw'});

k = height(SpearmanResults);
SpearmanResults.p_bonf = min(SpearmanResults.p_raw * k, 1);

% ----------------------- eps guards for log plotting -----------------------
x_used = x_all(mask_all);
y_used = y_all(mask_all);

minpos_rate = min(x_used(x_used > 0));
if isempty(minpos_rate) || ~isfinite(minpos_rate), minpos_rate = 1e-6; end

minpos_sz = min(y_used(y_used > 0));
if isempty(minpos_sz) || ~isfinite(minpos_sz), minpos_sz = 1e-6; end

eps_rate = 0.5 * minpos_rate;
eps_sz   = 0.5 * minpos_sz;

xZero = log10(eps_sz);
yZero = log10(eps_rate);
xLims = spearman_xLims;
yLims = spearman_yLims;

% Build "Tall" for panel A (ONLY rows used in corr)
Tall = table;
Tall.SpikeRate_perHour = x_used;
Tall.SzFreq            = y_used;
Tall.logSpikeRate = log10(Tall.SpikeRate_perHour + (Tall.SpikeRate_perHour <= 0) .* eps_rate);
Tall.logSzFreq    = log10(Tall.SzFreq + (Tall.SzFreq <= 0) .* eps_sz);

isZeroSz   = (Tall.SzFreq == 0);
isZeroRate = (Tall.SpikeRate_perHour == 0);
onlySz     =  isZeroSz & ~isZeroRate;
onlyRate   = ~isZeroSz &  isZeroRate;
bothZero   =  isZeroSz &  isZeroRate;
nonZero    = ~(isZeroSz | isZeroRate);

% ----------------------- figure -----------------------
f2 = figure('Color','w','Position',[60 60 1200 900]);
tiledlayout(f2, 2, 2, 'Padding','compact','TileSpacing','compact');

% ======================= Panel A =======================
axA = nexttile(1); hold(axA,'on'); grid(axA,'on'); box(axA,'off');
xline(axA, xZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
yline(axA, yZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

% points
scatter(axA, Tall.logSzFreq(nonZero), Tall.logSpikeRate(nonZero), 14, COL_all, ...
    'filled', 'MarkerFaceAlpha', 0.25);
plot(axA, Tall.logSzFreq(onlySz),   Tall.logSpikeRate(onlySz),   '*', ...
    'Color', COL_all, 'MarkerSize', 7, 'LineWidth', 1);
plot(axA, Tall.logSzFreq(onlyRate), Tall.logSpikeRate(onlyRate), '*', ...
    'Color', COL_all, 'MarkerSize', 8, 'LineWidth', 1);
plot(axA, Tall.logSzFreq(bothZero), Tall.logSpikeRate(bothZero), '*', ...
    'Color', COL_all, 'MarkerSize', 8, 'LineWidth', 1.2);

% regression line (black)
if n_all >= 3
    X = [ones(n_all,1), Tall.logSzFreq];
    b = X \ Tall.logSpikeRate;
    xgrid = linspace(xLims(1), xLims(2), 300)';
    plot(axA, xgrid, b(1) + b(2)*xgrid, 'k-', 'LineWidth', 2);
end

xlim(axA, xLims); ylim(axA, yLims);
xlabel(axA,'Seizures per month (log scale)','FontSize',fontL);
ylabel(axA,'Spikes per hour (log scale)','FontSize',fontL);
set_log10_ticks(axA, 'x', eps_sz,   xLims);
set_log10_ticks(axA, 'y', eps_rate, yLims);

title(axA, sprintf('A. All epilepsy%s (N=%d)', labelSuffix, n_all), ...
    'FontSize', fontL, 'FontWeight','bold');

txtA = sprintf('\\rho=%.2f [%.2f-%.2f], %s', rs_all, rho_lo, rho_hi, p_label(p_all));

text(axA, 0.98, 0.95, txtA, 'Units','normalized', ...
    'HorizontalAlignment','right', 'VerticalAlignment','top', ...
    'FontSize', fontL-2, 'FontWeight','bold');
set(axA,'FontSize',fontL);

% ======================= Panels B/C/D =======================
for p = 1:3
    ax = nexttile(p+1); hold(ax,'on'); grid(ax,'on'); box(ax,'off');

    gStr = string(panelOrder{p});

    % pick color
    switch gStr
        case "Frontal"
            col = COL_front;
        case "Temporal"
            col = COL_temp;
        case "General"
            col = COL_gen;
        otherwise
            col = COL_all;
    end

    idx = (string(PatientSpikeSz_Typed.EpiType3) == gStr) & ...
          isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perHour) & ...
          isfinite(PatientSpikeSz_Typed.(freqFieldName));

    if nonZeroOnly
        idx = idx & (PatientSpikeSz_Typed.MeanSpikeRate_perHour > 0) & ...
                  (PatientSpikeSz_Typed.(freqFieldName) > 0);
    end

    if nnz(idx) == 0
        axis(ax,'off');
        continue
    end

    x_raw = double(PatientSpikeSz_Typed.(freqFieldName)(idx));
    y_raw = double(PatientSpikeSz_Typed.MeanSpikeRate_perHour(idx));

    logX = log10(x_raw + (x_raw <= 0) .* eps_sz);
    logY = log10(y_raw + (y_raw <= 0) .* eps_rate);

    % classify zeros (for star markers)
    isZx = (x_raw == 0);
    isZy = (y_raw == 0);
    onlySz_g   =  isZx & ~isZy;
    onlyRate_g = ~isZx &  isZy;
    bothZero_g =  isZx &  isZy;
    nonZero_g  = ~(isZx | isZy);

    xline(ax, xZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
    yline(ax, yZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    % points
    scatter(ax, logX(nonZero_g), logY(nonZero_g), 18, col, ...
        'filled', 'MarkerFaceAlpha', 0.35);
    if any(onlySz_g)
        plot(ax, logX(onlySz_g), logY(onlySz_g), '*', ...
            'Color', col, 'MarkerSize', 8, 'LineWidth', 1.1);
    end
    if any(onlyRate_g)
        plot(ax, logX(onlyRate_g), logY(onlyRate_g), '*', ...
            'Color', col, 'MarkerSize', 8, 'LineWidth', 1.1);
    end
    if any(bothZero_g)
        plot(ax, logX(bothZero_g), logY(bothZero_g), '*', ...
            'Color', col, 'MarkerSize', 9, 'LineWidth', 1.2);
    end

    % regression line in subtype color (fit on finite, paired points)
    if nnz(nonZero_g) >= 3
        Xg = [ones(nnz(nonZero_g),1), logX(nonZero_g)];
        bg = Xg \ logY(nonZero_g);
        xg = linspace(xLims(1), xLims(2), 250)';
        plot(ax, xg, bg(1) + bg(2)*xg, '-', 'Color', col, 'LineWidth', 2);
    end

    xlim(ax, xLims); ylim(ax, yLims);
    xlabel(ax,'Seizures per month (log scale)','FontSize',fontL);
    ylabel(ax,'Spikes per hour (log scale)','FontSize',fontL);
    set_log10_ticks(ax, 'x', eps_sz,   xLims);
    set_log10_ticks(ax, 'y', eps_rate, yLims);

    % annotation text from your computed tables
    row = SpearmanResults(strcmp(string(SpearmanResults.Group), gStr), :);
    rowCI = subtype_ci(subtype_ci.Group == gStr, :);

    txt = sprintf('\\rho=%.2f [%.2f-%.2f], p_{bonf}%s', ...
    row.Spearman_r, rowCI.ci_lo, rowCI.ci_hi, ...
    string(regexprep(char(p_label(row.p_bonf)), '^p', '')));

    nNow = nnz(idx);
    title(ax, sprintf('%s%s (N=%d)', panelTitle{p}, labelSuffix, nNow), ...
        'FontSize', fontL, 'FontWeight','bold');
    text(ax, 0.98, 0.95, txt, 'Units','normalized', ...
        'HorizontalAlignment','right', 'VerticalAlignment','top', ...
        'FontSize', fontL-3, 'FontWeight','bold');

    set(ax,'FontSize',fontL);
end

% ----------------------- save -----------------------
if ~exist(fileparts(fig_out),'dir')
    mkdir(fileparts(fig_out));
end
exportgraphics(f2, fig_out, 'Resolution', 300);
fprintf('Saved Spearman figure: %s\n', fig_out);

end

%% =====================================================================
%% ======================= SMALL UTIL HELPERS ==========================
%% =====================================================================

function arr = json_to_string_array(s)
s = strtrim(string(s));
if s=="" || s=="[]" || s=="<missing>"
    arr = strings(0,1); return;
end
dec = jsondecode(char(s));
if iscell(dec)
    arr = strings(numel(dec),1);
    for k=1:numel(dec)
        x = dec{k};
        if ischar(x) || (isstring(x) && isscalar(x)), arr(k)=string(x);
        else, arr(k)="";
        end
    end
elseif ischar(dec) || (isstring(dec)&&isscalar(dec))
    arr = string(dec);
elseif isnumeric(dec)
    arr = string(dec(:));
else
    error('Unsupported JSON string-array type: %s', class(dec));
end
arr = string(arr(:));
end

function arr = json_to_double_array(s)
s = strtrim(string(s));
if s=="" || s=="[]" || s=="<missing>"
    arr = double([]); return;
end
s = regexprep(s,'null','NaN','ignorecase');
dec = jsondecode(char(s));
arr = double(dec(:));
end

function out = max_hasSz(x)
x = x(isfinite(x));
if isempty(x), out = NaN; else, out = max(x); end
end

function out = local_frac_hasSz1(hsVec)
hs = hsVec(isfinite(hsVec));
if isempty(hs), out = NaN; return; end
valid = (hs==0 | hs==1);
if ~any(valid), out = NaN; else, out = nnz(hs(valid)==1)/nnz(valid); end
end

function out = local_min_omitnan(a)
a = a(isfinite(a));
if isempty(a), out = NaN; else, out = min(a); end
end

function out = local_first_nonmissing(s)
s = s(:);
s = s(~ismissing(s) & strlength(s)>0);
if isempty(s), out=""; else, out=s(1); end
end

function ylog = to_log10_per_hour(x_per_hour, eps_rate)
x = double(x_per_hour);
x(~isfinite(x) | x<=0) = eps_rate;
ylog = log10(x);
end

function ylog = to_log10_per_month(freq_per_month, eps_freq)
f = double(freq_per_month);
f(~isfinite(f) | f<=0) = eps_freq;
ylog = log10(f);
end

function Yj = add_y_jitter_eps(Y, Y_ZERO, Y_LIMS, frac)
Yj = Y;
mask = abs(Y - Y_ZERO) < 1e-9;
if any(mask)
    amp = frac * diff(Y_LIMS);
    Yj(mask) = Yj(mask) + (rand(sum(mask),1)-0.5)*amp;
end
end

function add_sigbar(ax, x1, x2, y, ptext)
tick = 0.03 * diff(ax.YLim);
plot(ax, [x1 x1 x2 x2], [y-tick, y, y, y-tick], 'k-', 'LineWidth', 1.3);
text(ax, mean([x1 x2]), y + 0.02*diff(ax.YLim), ptext, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);
end

function pStr = p_label(p)
if isnan(p), pStr="p=NaN"; return; end
if p < 0.001, pStr="p<0.001"; return; end
if p < 0.01,  pStr=sprintf("p=%.2g", p); return; end
pStr = sprintf("p=%.2f", p);
end

function pStr = format_p_html(p)
if isnan(p), pStr='p = NaN'; return; end
if p < 0.001, pStr='p &lt; 0.001'; return; end
if p < 0.01
    s = sprintf('%.2g', p);
    if startsWith(s,'.'), s=['0' s]; end
    pStr = ['p = ' s];
    return;
end
pStr = sprintf('p = %.2f', p);
end

function set_log10_ticks(ax, whichAxis, eps_val, axisLims, maxPow)
if nargin<5 || isempty(maxPow), maxPow=6; end
whichAxis = lower(whichAxis);
decades = 10.^(0:maxPow);
log_dec = log10(decades);
keep = (log_dec>=axisLims(1)) & (log_dec<=axisLims(2));
ticks = log_dec(keep);
labels = string(decades(keep));

log_eps = log10(eps_val);
if log_eps>=axisLims(1) && log_eps<=axisLims(2)
    ticks  = [log_eps; ticks(:)];
    labels = ["0"; labels(:)];
end

if whichAxis=="x"
    ax.XTick = ticks; ax.XTickLabel = labels;
else
    ax.YTick = ticks; ax.YTickLabel = labels;
end
end

function [med, lo, hi] = bootstrap_median_ci(x, nBoot, alpha)
x = double(x(:)); x = x(isfinite(x));
if isempty(x), med=NaN; lo=NaN; hi=NaN; return; end
med = median(x,'omitnan');
n = numel(x);
bootM = nan(nBoot,1);
for b=1:nBoot
    idx = randi(n,n,1);
    bootM(b) = median(x(idx));
end
lo = prctile(bootM, 100*(alpha/2));
hi = prctile(bootM, 100*(1-alpha/2));
end

function add_median_ci_overlay(ax, xpos, med_raw, lo_raw, hi_raw, eps_floor)
yMed = to_log10_per_hour(med_raw, eps_floor);
yLo  = to_log10_per_hour(lo_raw,  eps_floor);
yHi  = to_log10_per_hour(hi_raw,  eps_floor);
plot(ax, [xpos xpos], [yLo yHi], 'k-', 'LineWidth', 3);
plot(ax, xpos, yMed, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6);
end

function add_median_ci_overlay_month(ax, xpos, med_raw, lo_raw, hi_raw, eps_floor)
yMed = to_log10_per_month(med_raw, eps_floor);
yLo  = to_log10_per_month(lo_raw,  eps_floor);
yHi  = to_log10_per_month(hi_raw,  eps_floor);
plot(ax, [xpos xpos], [yLo yHi], 'k-', 'LineWidth', 3);
plot(ax, xpos, yMed, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6);
end

function [rho_hat, lo, hi] = bootstrap_spearman_ci(x, y, nBoot, alpha)
x = double(x(:)); y = double(y(:));
mask = isfinite(x) & isfinite(y);
x = x(mask); y = y(mask);
n = numel(x);
if n<3, rho_hat=NaN; lo=NaN; hi=NaN; return; end
rho_hat = corr(x,y,'Type','Spearman','Rows','complete');
rho_boot = nan(nBoot,1);
for b=1:nBoot
    idx = randi(n,n,1);
    rho_boot(b) = corr(x(idx), y(idx), 'Type','Spearman','Rows','complete');
end
lo = prctile(rho_boot, 100*(alpha/2));
hi = prctile(rho_boot, 100*(1-alpha/2));
end

function d = cliff_delta(x1, x2)
x1 = x1(:); x2 = x2(:);
x1 = x1(isfinite(x1)); x2 = x2(isfinite(x2));
n1=numel(x1); n2=numel(x2);
if n1==0 || n2==0, d=NaN; return; end
[~,~,stats] = ranksum(x1,x2,'method','approx');
R1 = stats.ranksum;
U1 = R1 - n1*(n1+1)/2;
d = (2*U1/(n1*n2)) - 1;
end

function Diag = paired_closest_visit_diag(SessWithDate, Vuniq_R1, EpPatients, maxDaysBefore, maxDaysAfter)
% paired_closest_visit_diag
% For each epilepsy patient:
%   1) For each EEG session, find the closest clinic visit within window:
%        -maxDaysBefore <= (EEG_Date - VisitDate) <= +maxDaysAfter
%   2) Among EEGs with a matched visit, pick the low-spike EEG and high-spike EEG
%   3) Return closest visit date + signed gap + pre/post gaps for each of low/high
%
% Inputs:
%   SessWithDate: table with columns Patient, Session, EEG_Date, SpikesPerHour
%   Vuniq_R1:     table with columns Patient, VisitDate, Freq_R1
%   EpPatients:   vector of patient IDs (double ok)
%   maxDaysBefore, maxDaysAfter: scalars (days)
%
% Output:
%   Diag: one row per patient with both low/high resolved (requires >=2 matched EEGs)

% ---- strict column checks ----
needS = ["Patient","Session","EEG_Date","SpikesPerHour"];
needV = ["Patient","VisitDate","Freq_R1"];
missS = setdiff(needS, string(SessWithDate.Properties.VariableNames));
missV = setdiff(needV, string(Vuniq_R1.Properties.VariableNames));
if ~isempty(missS), error("SessWithDate missing: %s", strjoin(missS,", ")); end
if ~isempty(missV), error("Vuniq_R1 missing: %s", strjoin(missV,", ")); end

EpPatients = unique(double(EpPatients(:)));

% restrict to epilepsy patients
SessWithDate = innerjoin(SessWithDate, table(EpPatients,'VariableNames',{'Patient'}), 'Keys','Patient');
Vuniq_R1     = innerjoin(Vuniq_R1,     table(EpPatients,'VariableNames',{'Patient'}), 'Keys','Patient');

% build per-(patient,session) closest visit match
Pairs = table('Size',[0 9], ...
    'VariableTypes', {'double','double','datetime','double','datetime','double','double','double','double'}, ...
    'VariableNames', {'Patient','Session','EEG_Date','SpikesPerHour','ClosestVisit','GapDays', ...
                      'PreGapDays','PostGapDays','Freq_R1'});

for pid = EpPatients'
    Ssub = SessWithDate(SessWithDate.Patient==pid,:);
    Vsub = Vuniq_R1(Vuniq_R1.Patient==pid,:);
    if isempty(Ssub) || isempty(Vsub), continue; end

    for j = 1:height(Ssub)
        EEGd = Ssub.EEG_Date(j);

        % GapDays = EEG - Visit ( + means visit before EEG; - means visit after EEG )
        gap = days(EEGd - Vsub.VisitDate);

        keep = (gap >= -maxDaysBefore) & (gap <= maxDaysAfter);
        if ~any(keep), continue; end

        kk = find(keep);
        [~,bestRel] = min(abs(gap(kk)));
        bestIdx = kk(bestRel);

        gapBest = gap(bestIdx);
        preGap  = max(gapBest, 0);     % days visit BEFORE EEG
        postGap = max(-gapBest, 0);    % days visit AFTER EEG

        Pairs = [Pairs; { ...
            pid, double(Ssub.Session(j)), EEGd, double(Ssub.SpikesPerHour(j)), ...
            Vsub.VisitDate(bestIdx), gapBest, preGap, postGap, double(Vsub.Freq_R1(bestIdx)) ...
            }]; %#ok<AGROW>
    end
end

% require finite values for spike and freq
Pairs = Pairs(isfinite(Pairs.SpikesPerHour) & isfinite(Pairs.Freq_R1), :);
if isempty(Pairs)
    Diag = table();
    return
end

% collapse duplicates if any (keep closest already; but just in case)
[grpPS, ~] = findgroups(Pairs.Patient, Pairs.Session);
keepRow = false(height(Pairs),1);
for g = 1:max(grpPS)
    idx = find(grpPS==g);
    [~,best] = min(abs(Pairs.GapDays(idx)));
    keepRow(idx(best)) = true;
end
Pairs = Pairs(keepRow,:);

% ---- pick low/high spike EEG per patient ----
[gp, uP] = findgroups(Pairs.Patient);

Diag = table('Size',[0 17], ...
    'VariableTypes', {'double','double', ...
                      'double','datetime','datetime','double','double','double','double', ...
                      'double','datetime','datetime','double','double','double','double','double'}, ...
    'VariableNames', {'Patient','nEEG_matched', ...
                      'Spike_low','EEG_low','Visit_low','Gap_low','PreGap_low','PostGap_low','Freq_low', ...
                      'Spike_high','EEG_high','Visit_high','Gap_high','PreGap_high','PostGap_high','Freq_high', ...
                      'AbsDiffSpike'});

for k = 1:numel(uP)
    pid = uP(k);
    idx = find(gp==k);
    if numel(idx) < 2, continue; end

    r = Pairs.SpikesPerHour(idx);
    [~,iL] = min(r);
    [~,iH] = max(r);
    if iL==iH, continue; end

    rowL = idx(iL);
    rowH = idx(iH);

    Diag = [Diag; { ...
        pid, numel(idx), ...
        Pairs.SpikesPerHour(rowL), Pairs.EEG_Date(rowL), Pairs.ClosestVisit(rowL), ...
        Pairs.GapDays(rowL), Pairs.PreGapDays(rowL), Pairs.PostGapDays(rowL), Pairs.Freq_R1(rowL), ...
        Pairs.SpikesPerHour(rowH), Pairs.EEG_Date(rowH), Pairs.ClosestVisit(rowH), ...
        Pairs.GapDays(rowH), Pairs.PreGapDays(rowH), Pairs.PostGapDays(rowH), Pairs.Freq_R1(rowH), ...
        Pairs.SpikesPerHour(rowH) - Pairs.SpikesPerHour(rowL) ...
        }]; %#ok<AGROW>
end

end
