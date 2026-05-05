function run_spike_sz_pipeline_clean()
%% ===================== INFORMATION ==================
%{
This contains the code for running all analyses associated with the
manuscript "Interictal Spike Rate on Routine Outpatient EEG Is Associated
With Seizure Frequency in a Large Epilepsy Cohort" by Conrad et al., 2026

Requirements:
1. spike_counts.csv and clinical_data_deidentified.csv
     - both are available for download at: https://upenn.box.com/s/yy4o1t6nit7yu35flz59ux6zf54slg9m
2. Matlab (vR2024a) and the Statistics and Machine Learning Toolbox
3. This codebase: https://github.com/erinconrad/seizure_severity
4. Create output and data directories in the paths below, place CSVs in data directory.

To run:
>> run_spike_sz_pipeline_clean
%}

%% ======================= RNG =======================
rng(1);

%% ======================= PATHS =======================
spikeSummaryMultiCsv = '../data/spike_counts.csv';
reportCsv            = '../data/clinical_data_deidentified.csv';

fig1_out    = '../output/Fig1.png';
fig2_out    = '../output/Fig2.png';
figS1_out   = '../output/FigS1.png';
figS2_out   = '../output/FigS2.png';
figMain_out = '../output/FigModel.png';
Table1Csv   = fullfile('..','output','Table1.csv');
tableS1Csv  = fullfile('..','output','TableS1.csv');
resultsHtml = '../output/results_summary.html';

%% ======================= PARAMETERS =======================
MAX_ROUTINE_HOURS = 4;
NESD_LABEL        = "Non-Epileptic Seizure Disorder";
badTypes          = lower(["Uncertain if Epilepsy","Unknown or MRN not found",""]);
canonical3        = ["General","Temporal","Frontal"];

EPS_RATE       = 30e-3;
Y_ZERO         = log10(EPS_RATE);
Y_LIMS         = [-2 4];
TITLE_Y_OFFSET = 0.02;

spearman_xLims = [-3.5, 4];
spearman_yLims = [-1.5, 3];

nBoot    = 5000;
alpha    = 0.05;
countCol = "count_0_46";
durCol   = "Duration_sec";

allowable_visits = [
    "CONSULT VISIT","ESTABLISHED PATIENT VISIT",...
    "FOLLOW-UP PATIENT CLINIC","NEW PATIENT CLINIC","NEW PATIENT VISIT",...
    "NPV MANAGEMENT DURING COVID-19","NPV NEUROLOGY",...
    "RETURN ANNUAL VISIT","RETURN PATIENT EXTENDED","RETURN PATIENT VISIT",...
    "RPV MANAGEMENT DURING COVID-19","TELEHEALTH VIDEO VISIT RETURN"
];

%% ======================= LOAD DATA =======================
SpikeSummaryTable = readtable(spikeSummaryMultiCsv,'TextType','string','VariableNamingRule','preserve');
require_cols(SpikeSummaryTable, ["Patient","Session",countCol,durCol], "SpikeSummaryTable");
SpikeSummaryTable.SpikeRate_perHour = SpikeSummaryTable.(countCol) ./ SpikeSummaryTable.(durCol) * 3600;

ReportTable = readtable(reportCsv,'TextType','string','VariableNamingRule','preserve');
require_cols(ReportTable, ...
    ["patient_id","session_number","acquired_on", ...
     "report_PATIENT_CLASS","jay_in_or_out", ...
     "visit_type","visit_dates_deid","sz_freqs","visit_hasSz", ...
     "epilepsy_type","epilepsy_specific","nlp_gender","deid_birth_date","start_time_deid", ...
     "report_SPORADIC_EPILEPTIFORM_DISCHARGES","jay_focal_epi","jay_multifocal_epi","jay_gen_epi"], ...
    "ReportTable");

%% ======================= FILTER =======================
ReportTable = filter_visit_arrays_by_type(ReportTable, allowable_visits);
[SpikeSummaryTable, ReportTable] = filter_outpatient_routine(SpikeSummaryTable, ReportTable, durCol, MAX_ROUTINE_HOURS);

assert_unique_keys(SpikeSummaryTable, "Patient","Session","SpikeSummaryTable");
assert_unique_keys(ReportTable, "patient_id","session_number","ReportTable");

%% ======================= BUILD COHORT =======================
Vuniq            = build_visit_level_table_R1(ReportTable);
PatientTypingAll = build_patient_typing_from_report(ReportTable, canonical3);
SzFreqPerPatient = build_patient_seizure_metrics(Vuniq);

Views = build_filtered_view(SpikeSummaryTable, ReportTable, PatientTypingAll, ...
    SzFreqPerPatient, NESD_LABEL, badTypes, canonical3);

%% ======================= BUILD PAIR TABLE =======================
% Canonical-subtype patients — for M1-M3
PairTable = build_eeg_visit_pairs(Vuniq, Views.SessionLevelSpikeRates, ...
    Views.ReportForKeptSessions, Views.PatientTypingFiltered);
fprintf('Canonical-subtype pairs: %d, patients: %d\n', ...
    height(PairTable), numel(unique(PairTable.Patient)));

% All epilepsy patients — for M4-M6 (no subtype restriction)
PairTable_Full = build_eeg_visit_pairs_nosubtype(Vuniq, Views.SessionLevelSpikeRates, ...
    Views.ReportForKeptSessions, Views.PatientTyping_AllEpilepsy);
fprintf('All-epilepsy pairs: %d, patients: %d\n', ...
    height(PairTable_Full), numel(unique(PairTable_Full.Patient)));

%% ======================= FIT MODELS =======================
MMR = fit_mixed_effects_models(PairTable, PairTable_Full, nBoot, alpha);

%% ======================= FIGURE 1 =======================
[Fig1, Fig1Stats] = make_fig1_controls(Views, EPS_RATE, Y_ZERO, Y_LIMS, TITLE_Y_OFFSET, nBoot, alpha);
save_fig(Fig1, fig1_out);
fprintf('Saved Fig 1: %s\n', fig1_out);

%% ======================= FIGURE 2 + FIG S1 =======================
[SpearmanResults_main, rs_all_main, p_all_main, n_all_main, rho_lo_main, rho_hi_main, subtype_ci_main] = ...
    spearman_plotting_function(Views.PatientSpikeSz_All, Views.PatientSpikeSz_Typed, ...
        canonical3, spearman_xLims, spearman_yLims, fig2_out, 'MeanSzFreq', '', false);

[SpearmanResults_S1, rs_all_S1, p_all_S1, n_all_S1, rho_lo_S1, rho_hi_S1, subtype_ci_S1] = ...
    spearman_plotting_function(Views.PatientSpikeSz_All, Views.PatientSpikeSz_Typed, ...
        canonical3, spearman_xLims, spearman_yLims, figS1_out, 'MeanSzFreq', ...
        ' (positive spike/seizures only)', true);

fprintf('Saved Fig 2: %s\n', fig2_out);
fprintf('Saved Fig S1: %s\n', figS1_out);

%% Model figs
% Main figure
FigMain = make_model_figure(MMR, figMain_out);

% Supplemental
figSupLag_out = '../output/FigSupLag.png';
FigSupLag = make_figSup_lag(MMR, Vuniq, Views.ReportForKeptSessions, figSupLag_out);

%% ======================= FIG S2 =======================
FigS2 = make_figS2_sz_by_reported_spikes(Views, SzFreqPerPatient, nBoot, alpha, spearman_xLims);
save_fig(FigS2, figS2_out);
fprintf('Saved Fig S2: %s\n', figS2_out);

%% ======================= TABLE 1 =======================
Table1_flat = build_table1_flat(Views, SzFreqPerPatient, Vuniq, EPS_RATE, nBoot, alpha);
if ~exist(fileparts(Table1Csv),'dir'), mkdir(fileparts(Table1Csv)); end
writetable(Table1_flat, Table1Csv);
fprintf('Wrote Table 1: %s\n', Table1Csv);

%% ======================= TABLE S1 =======================
write_tableS1(MMR, tableS1Csv);
fprintf('Wrote Table S1: %s\n', tableS1Csv);

%% ======================= HTML RESULTS =======================
write_results_html(resultsHtml, Views, SzFreqPerPatient, Fig1Stats, ...
    SpearmanResults_main, rs_all_main, p_all_main, n_all_main, ...
    rho_lo_main, rho_hi_main, subtype_ci_main, ...
    SpearmanResults_S1, rs_all_S1, p_all_S1, n_all_S1, ...
    rho_lo_S1, rho_hi_S1, subtype_ci_S1, ...
    Views.ReportForKeptSessions, MMR);
fprintf('Wrote HTML: %s\n', resultsHtml);


end

%% =====================================================================
%% CORE PIPELINE HELPERS
%% =====================================================================

function require_cols(T, cols, name)
missing = setdiff(string(cols), string(T.Properties.VariableNames));
if ~isempty(missing)
    error('%s is missing required columns: %s', name, strjoin(missing, ", "));
end
end

function assert_unique_keys(T, pidCol, sesCol, name)
key = string(T.(pidCol)) + "|" + string(T.(sesCol));
if numel(unique(key)) ~= numel(key)
    error('%s has duplicated (Patient,Session) keys.', name);
end
end

function save_fig(figH, outPath)
if ~exist(fileparts(outPath),'dir'), mkdir(fileparts(outPath)); end
exportgraphics(figH, outPath, 'Resolution', 300);
end

function R = filter_visit_arrays_by_type(R, allowable_visits)
vt_raw_all = strtrim(string(R.visit_type));
mask_null_only = (vt_raw_all == "[null]") | (vt_raw_all == "null");
R.visit_type(mask_null_only)       = "[]";
R.visit_dates_deid(mask_null_only) = "[]";
R.sz_freqs(mask_null_only)         = "[]";
R.visit_hasSz(mask_null_only)      = "[]";

total_before = 0; total_after = 0;
for i = 1:height(R)
    vt_raw    = strtrim(string(R.visit_type(i)));
    dates_raw = strtrim(string(R.visit_dates_deid(i)));
    sz_raw    = strtrim(string(R.sz_freqs(i)));
    hs_raw    = strtrim(string(R.visit_hasSz(i)));

    if vt_raw=="" || vt_raw=="[]" || vt_raw=="<missing>"
        R.visit_type(i)="[]"; R.visit_dates_deid(i)="[]";
        R.sz_freqs(i)="[]";   R.visit_hasSz(i)="[]";
        continue
    end

    vt    = json_to_string_array(vt_raw);
    dates = json_to_string_array(dates_raw);
    sz    = json_to_double_array(sz_raw);
    hs    = json_to_double_array(hs_raw);

    if ~(numel(vt)==numel(dates) && numel(vt)==numel(sz) && numel(vt)==numel(hs))
        error('Row %d: visit arrays have mismatched lengths.', i);
    end

    total_before = total_before + numel(vt);
    keepMask = ismember(string(vt), allowable_visits);

    if ~any(keepMask)
        R.visit_type(i)="[]"; R.visit_dates_deid(i)="[]";
        R.sz_freqs(i)="[]";   R.visit_hasSz(i)="[]";
        continue
    end

    vt_f    = cellstr(string(vt(keepMask)));
    dates_f = cellstr(string(dates(keepMask)));
    sz_f    = sz(keepMask);
    hs_f    = hs(keepMask);
    total_after = total_after + numel(vt_f);

    R.visit_type(i)       = string(jsonencode(vt_f));
    R.visit_dates_deid(i) = string(jsonencode(dates_f));
    R.sz_freqs(i)         = string(jsonencode(sz_f));
    R.visit_hasSz(i)      = string(jsonencode(hs_f));
end
fprintf('[Visit-type filter] %d -> %d visits (kept %.1f%%)\n', ...
    total_before, total_after, 100*total_after/max(1,total_before));
end

function [S, R] = filter_outpatient_routine(S, R, durCol, MAX_ROUTINE_HOURS)
nR0 = height(R); nS0 = height(S);

acqStr   = lower(strtrim(string(R.acquired_on)));
classStr = lower(strtrim(string(R.report_PATIENT_CLASS)));
jayStr   = lower(strtrim(string(R.jay_in_or_out)));

isOutpt = contains(acqStr,"spe") | contains(acqStr,"radnor") | ...
          (classStr == "outpatient") | (jayStr == "out");

OutptKeys = unique(R(isOutpt, {'patient_id','session_number'}));
OutptKeys.Properties.VariableNames = {'Patient','Session'};
if isempty(OutptKeys), error('No outpatient sessions identified.'); end

isRoutine   = isfinite(S.(durCol)) & S.(durCol) <= MAX_ROUTINE_HOURS*3600;
RoutineKeys = unique(S(isRoutine, {'Patient','Session'}));

OutptRoutineKeys = innerjoin(OutptKeys, RoutineKeys, 'Keys', {'Patient','Session'});
S = innerjoin(S, OutptRoutineKeys, 'Keys', {'Patient','Session'});
R = innerjoin(R, OutptRoutineKeys, 'LeftKeys', {'patient_id','session_number'}, ...
    'RightKeys', {'Patient','Session'});

fprintf('[Outpatient+routine] Kept %d/%d spike rows, %d/%d report rows\n', ...
    height(S), nS0, height(R), nR0);
end

function Vuniq = build_visit_level_table_R1(R)
PV = table('Size',[0 4], 'VariableTypes',{'double','datetime','double','double'}, ...
    'VariableNames',{'Patient','VisitDate','Freq','HasSz'});

for j = 1:height(R)
    pid = double(R.patient_id(j));
    ds  = strtrim(string(R.visit_dates_deid(j)));
    if ds=="[]" || ds=="", continue; end
    d = datetime(string(jsondecode(char(ds))),'InputFormat','yyyy-MM-dd');

    s = strtrim(string(R.sz_freqs(j)));
    s = regexprep(s,'null','NaN','ignorecase');
    v = double(jsondecode(char(s)));
    v(~isfinite(v)|v<0) = NaN;

    hs = strtrim(string(R.visit_hasSz(j)));
    h  = double(jsondecode(char(hs)));
    h(h==2) = NaN;

    if ~(numel(d)==numel(v) && numel(d)==numel(h))
        error('Row %d: visit arrays mismatched.', j); end
    if numel(d)==0, continue; end

    PV = [PV; table(repmat(pid,numel(d),1), d(:), v(:), h(:), ...
        'VariableNames', PV.Properties.VariableNames)]; %#ok<AGROW>
end

[gv, pid_keys, date_keys] = findgroups(PV.Patient, PV.VisitDate);
Freq_agg = splitapply(@(x) mean(x(isfinite(x)),'omitnan'), PV.Freq, gv);
Has_agg  = splitapply(@max_hasSz, PV.HasSz, gv);

Vuniq = table(pid_keys, date_keys, Freq_agg, Has_agg, ...
    'VariableNames', {'Patient','VisitDate','Freq','HasSz'});
Vuniq.Freq_R1 = Vuniq.Freq;
mask_rule1 = ~isfinite(Vuniq.Freq_R1) & (Vuniq.HasSz==0);
Vuniq.Freq_R1(mask_rule1) = 0;
end

function SzP = build_patient_seizure_metrics(Vuniq)
[gp, pids] = findgroups(Vuniq.Patient);
MeanSzFreq        = splitapply(@(x) mean(x,'omitnan'), Vuniq.Freq_R1, gp);
FracVisits_HasSz1 = splitapply(@local_frac_hasSz1, Vuniq.HasSz, gp);
SzP = table(pids, MeanSzFreq, FracVisits_HasSz1, ...
    'VariableNames', {'Patient','MeanSzFreq','FracVisits_HasSz1'});
end

function T = build_patient_typing_from_report(R, canonical3)
pid   = double(R.patient_id);
etype = strtrim(string(R.epilepsy_type));
espec = strtrim(string(R.epilepsy_specific));

okT = ~ismissing(etype) & strlength(etype)>0;
okS = ~ismissing(espec) & strlength(espec)>0;

Ttype = table(pid(okT), etype(okT), 'VariableNames', {'Patient','EpilepsyType_raw'});
Ttype = sortrows(Ttype, 'Patient');
[uid, ~, g] = unique(Ttype.Patient, 'stable');
etype_one = strings(numel(uid),1);
for k=1:numel(uid)
    vals = unique(Ttype.EpilepsyType_raw(g==k));
    vals = vals(strlength(vals)>0);
    if isempty(vals), continue; end
    if numel(vals)>1, error('Conflicting epilepsy_type for Patient %g', uid(k)); end
    etype_one(k) = vals(1);
end

Tspec = table(pid(okS), espec(okS), 'VariableNames', {'Patient','EpilepsySpecific_raw'});
Tspec = sortrows(Tspec, 'Patient');
[uidS, ~, gS] = unique(Tspec.Patient, 'stable');
espec_one = strings(numel(uidS),1);
for k=1:numel(uidS)
    vals = unique(Tspec.EpilepsySpecific_raw(gS==k));
    vals = vals(strlength(vals)>0);
    if isempty(vals), continue; end
    if numel(vals)>1, error('Conflicting epilepsy_specific for Patient %g', uidS(k)); end
    espec_one(k) = vals(1);
end

assert(all(uid==uidS))
T = table(uid, etype_one, 'VariableNames', {'Patient','EpilepsyType'});
T = outerjoin(T, table(uidS, espec_one, 'VariableNames', {'Patient','EpilepsySpecific'}), ...
    'Keys','Patient','MergeKeys',true);

spec_norm = lower(strtrim(string(T.EpilepsySpecific)));
type_norm = lower(strtrim(string(T.EpilepsyType)));
E3 = strings(height(T),1);
E3(contains(spec_norm,"temporal"))                   = "Temporal";
E3(contains(spec_norm,"frontal"))                    = "Frontal";
E3((strlength(E3)==0) & strcmp(type_norm,"general")) = "General";
T.EpiType3 = categorical(E3, canonical3);
end

%% =====================================================================
%% VIEW BUNDLE
%% =====================================================================

function Views = build_filtered_view(SessionsFiltered, ReportIn, PatientTypingAll, ...
    SzFreqPerPatient, NESD_LABEL, badTypes, canonical3)

%% Join report to kept sessions
SessKeys = unique(SessionsFiltered(:,{'Patient','Session'}));
ReportForKeptSessions = innerjoin(ReportIn, SessKeys, ...
    'LeftKeys',{'patient_id','session_number'}, 'RightKeys',{'Patient','Session'});
ReportForKeptSessions.Patient = double(ReportForKeptSessions.patient_id);
ReportForKeptSessions.Session = double(ReportForKeptSessions.session_number);

%% Patient-level mean spike rate + typing
PatientsKept   = unique(double(SessionsFiltered.Patient));
TypingFiltered = innerjoin(PatientTypingAll, ...
    table(PatientsKept,'VariableNames',{'Patient'}), 'Keys','Patient');

[gpS, pidsS] = findgroups(double(SessionsFiltered.Patient));
MeanSpikeRate_perHour = splitapply(@(x) mean(x,'omitnan'), SessionsFiltered.SpikeRate_perHour, gpS);
PatientLevelSpikeRates = table(double(pidsS), MeanSpikeRate_perHour, ...
    'VariableNames', {'Patient','MeanSpikeRate_perHour'});
PatientLevelSpikeRates = innerjoin(PatientLevelSpikeRates, ...
    TypingFiltered(:,{'Patient','EpilepsyType','EpilepsySpecific','EpiType3'}), 'Keys','Patient');

%% Epilepsy masks
etype_norm     = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));
IsNESDMask     = (etype_norm == lower(strtrim(NESD_LABEL)));
IsBadType      = ismember(etype_norm, badTypes) | ismissing(etype_norm) | (strlength(etype_norm)==0);
IsEpilepsyMask = ~IsNESDMask & ~IsBadType;

%% Session-level spike rates
SessionLevelSpikeRates = SessionsFiltered(:, {'Patient','Session','SpikeRate_perHour'});
SessionLevelSpikeRates.Properties.VariableNames{'SpikeRate_perHour'} = 'SpikesPerHour';

%% Define cohort: epilepsy + documented seizure frequency
SzFreqFiltered = innerjoin(SzFreqPerPatient, ...
    table(PatientsKept,'VariableNames',{'Patient'}), 'Keys','Patient');
EpPatients     = table(PatientLevelSpikeRates.Patient(IsEpilepsyMask), 'VariableNames',{'Patient'});
SzFreqEpilepsy = innerjoin(SzFreqFiltered, EpPatients, 'Keys','Patient');

CohortPatients = unique(SzFreqEpilepsy.Patient(isfinite(SzFreqEpilepsy.MeanSzFreq)));
CohortTable    = table(CohortPatients, 'VariableNames',{'Patient'});

%% Restrict all tables to cohort
PatientLevelSpikeRates = innerjoin(PatientLevelSpikeRates, CohortTable, 'Keys','Patient');
TypingFiltered         = innerjoin(TypingFiltered,         CohortTable, 'Keys','Patient');
SessionsFiltered       = innerjoin(SessionsFiltered,       CohortTable, 'Keys','Patient');
ReportForKeptSessions  = innerjoin(ReportForKeptSessions,  CohortTable, 'Keys','Patient');
SzFreqEpilepsy         = innerjoin(SzFreqEpilepsy,         CohortTable, 'Keys','Patient');

etype_norm     = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));
IsNESDMask     = (etype_norm == lower(strtrim(NESD_LABEL)));
IsBadType      = ismember(etype_norm, badTypes) | ismissing(etype_norm) | (strlength(etype_norm)==0);
IsEpilepsyMask = ~IsNESDMask & ~IsBadType;
assert(all(IsEpilepsyMask))

fprintf('[Cohort] %d epilepsy patients with documented seizure frequency\n', numel(CohortPatients));

%% All-epilepsy typing table (no subtype restriction) for primary model
PatientTyping_AllEpilepsy = TypingFiltered;

%% Canonical-subtype typing table for Spearman only
keepCanon3    = ~ismissing(PatientLevelSpikeRates.EpiType3) & ...
                ismember(string(PatientLevelSpikeRates.EpiType3), canonical3);
TypedPatients = PatientLevelSpikeRates.Patient(IsEpilepsyMask & keepCanon3);
TypingFiltered_Canonical = TypingFiltered(ismember(TypingFiltered.Patient, TypedPatients), :);

%% Spearman input tables
PatientSpikeSz_All = innerjoin(...
    PatientLevelSpikeRates(IsEpilepsyMask, {'Patient','MeanSpikeRate_perHour'}), ...
    SzFreqEpilepsy, 'Keys','Patient');
keepAll = isfinite(PatientSpikeSz_All.MeanSpikeRate_perHour) & ...
          isfinite(PatientSpikeSz_All.MeanSzFreq);
PatientSpikeSz_All = PatientSpikeSz_All(keepAll,:);

SzFreqCanon = innerjoin(SzFreqEpilepsy, ...
    PatientLevelSpikeRates(IsEpilepsyMask & keepCanon3, {'Patient','EpiType3'}), 'Keys','Patient');
PatientSpikeSz_Typed = innerjoin(...
    PatientLevelSpikeRates(IsEpilepsyMask & ismember(PatientLevelSpikeRates.Patient, TypedPatients), ...
        {'Patient','MeanSpikeRate_perHour'}), ...
    SzFreqCanon, 'Keys','Patient');
keepTyped = isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perHour) & ...
            isfinite(PatientSpikeSz_Typed.MeanSzFreq) & ...
            ~ismissing(PatientSpikeSz_Typed.EpiType3);
PatientSpikeSz_Typed = PatientSpikeSz_Typed(keepTyped,:);

%% Canonical3 group stats for Fig 1B
Canonical3_SubsetTable = PatientLevelSpikeRates(IsEpilepsyMask & keepCanon3, ...
    {'Patient','EpiType3','MeanSpikeRate_perHour'});
Canonical3_SubsetTable.Properties.VariableNames{'EpiType3'} = 'EpiType4';

[g3, cats3] = findgroups(Canonical3_SubsetTable.EpiType4);
medVals = splitapply(@(x) median(x,'omitnan'), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
p25Vals = splitapply(@(x) prctile(x,25),      Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
p75Vals = splitapply(@(x) prctile(x,75),      Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
nVals   = splitapply(@(x) sum(isfinite(x)),   Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
Canonical3_Stats = table(cats3, nVals, medVals, p25Vals, p75Vals, ...
    'VariableNames',{'EpiType4','GroupCount','Median','P25','P75'});

Canonical3_Pairs = ["General","Temporal"; "General","Frontal"; "Temporal","Frontal"];
p_pair = NaN(3,1);
for i = 1:3
    xa = Canonical3_SubsetTable.MeanSpikeRate_perHour(...
        Canonical3_SubsetTable.EpiType4 == categorical(Canonical3_Pairs(i,1), canonical3));
    xb = Canonical3_SubsetTable.MeanSpikeRate_perHour(...
        Canonical3_SubsetTable.EpiType4 == categorical(Canonical3_Pairs(i,2), canonical3));
    if nnz(isfinite(xa))>=3 && nnz(isfinite(xb))>=3
        p_pair(i) = ranksum(xa, xb, 'method','approx');
    end
end

%% Bundle
Views.SessionsForFigures         = SessionsFiltered;
Views.ReportForKeptSessions      = ReportForKeptSessions;
Views.PatientTypingFiltered      = TypingFiltered_Canonical;
Views.PatientTyping_AllEpilepsy  = PatientTyping_AllEpilepsy;
Views.SessionLevelSpikeRates     = SessionLevelSpikeRates;
Views.PatientLevelSpikeRates     = PatientLevelSpikeRates;
Views.PatientSpikeSz_All         = PatientSpikeSz_All;
Views.PatientSpikeSz_Typed       = PatientSpikeSz_Typed;
Views.IsEpilepsyMask             = IsEpilepsyMask;
Views.IsNESDMask                 = IsNESDMask;
Views.Canonical3_SubsetTable     = Canonical3_SubsetTable;
Views.Canonical3_Stats           = Canonical3_Stats;
Views.Canonical3_Pairs           = Canonical3_Pairs;
Views.PvalsPairwise              = p_pair;
Views.PvalsPairwiseBonf          = min(p_pair*3, 1);
end

%% =====================================================================
%% PAIR TABLE BUILDER
%% =====================================================================
function PairTable = build_eeg_visit_pairs(Vuniq, SessionLevelSpikeRates, ...
    ReportForKeptSessions, PatientTyping)
%% Canonical-subtype patients — includes EpiType3

EEG_raw = ReportForKeptSessions.start_time_deid;
if isdatetime(EEG_raw), EEG_dt = EEG_raw;
else, EEG_dt = datetime(strtrim(string(EEG_raw)), 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss");
end

EEG_dates = table(double(ReportForKeptSessions.Patient), ...
    double(ReportForKeptSessions.Session), EEG_dt, ...
    'VariableNames',{'Patient','Session','EEG_Date'});
EEG_dates = EEG_dates(~isnat(EEG_dates.EEG_Date), :);

EEG_tbl = innerjoin(EEG_dates, ...
    SessionLevelSpikeRates(:,{'Patient','Session','SpikesPerHour'}), ...
    'Keys',{'Patient','Session'});

Vtyped = innerjoin(Vuniq(:,{'Patient','VisitDate','Freq_R1','HasSz'}), ...
    PatientTyping(:,{'Patient','EpiType3'}), 'Keys','Patient');

patients  = intersect(unique(EEG_tbl.Patient), unique(Vtyped.Patient));
nEstimate = 0;
for i = 1:numel(patients)
    p = patients(i);
    nEstimate = nEstimate + sum(EEG_tbl.Patient==p) * sum(Vtyped.Patient==p);
end

Patient_out       = nan(nEstimate,1);
Session_out       = nan(nEstimate,1);
VisitDate_out     = NaT(nEstimate,1);
SpikesPerHour_out = nan(nEstimate,1);
SzFreq_out        = nan(nEstimate,1);
HasSz_out         = nan(nEstimate,1);
SignedLag_out     = nan(nEstimate,1);
EpiType3_out      = strings(nEstimate,1);

row = 0;
for i = 1:numel(patients)
    p        = patients(i);
    eeg_rows = EEG_tbl(EEG_tbl.Patient==p, :);
    vis_rows = Vtyped(Vtyped.Patient==p, :);
    for e = 1:height(eeg_rows)
        for v = 1:height(vis_rows)
            row = row + 1;
            Patient_out(row)       = p;
            Session_out(row)       = eeg_rows.Session(e);
            VisitDate_out(row)     = vis_rows.VisitDate(v);
            SpikesPerHour_out(row) = eeg_rows.SpikesPerHour(e);
            SzFreq_out(row)        = vis_rows.Freq_R1(v);
            HasSz_out(row)         = vis_rows.HasSz(v);
            SignedLag_out(row)     = days(vis_rows.VisitDate(v) - eeg_rows.EEG_Date(e));
            EpiType3_out(row)      = string(vis_rows.EpiType3(v));
        end
    end
end

PairTable = table(Patient_out(1:row), Session_out(1:row), VisitDate_out(1:row), ...
    SpikesPerHour_out(1:row), SzFreq_out(1:row), HasSz_out(1:row), ...
    SignedLag_out(1:row), EpiType3_out(1:row), ...
    'VariableNames',{'Patient','Session','VisitDate','SpikesPerHour','SzFreq', ...
                     'HasSz','SignedLag_days','EpiType3'});

PairTable.EEG_ID = categorical(string(PairTable.Patient) + "_" + string(PairTable.Session));

keepMask  = isfinite(PairTable.SpikesPerHour) & isfinite(PairTable.SzFreq) & ...
            isfinite(PairTable.SignedLag_days) & strlength(PairTable.EpiType3)>0;
n_before  = height(PairTable);
PairTable = PairTable(keepMask,:);
fprintf('[build_eeg_visit_pairs] %d patients, %d pairs (%d removed)\n', ...
    numel(patients), height(PairTable), n_before-height(PairTable));

EPS_SPIKE = 1e-3;
PairTable.LogSpikesPerHour = log(PairTable.SpikesPerHour + EPS_SPIKE);
PairTable.SignedLag_years  = PairTable.SignedLag_days / 365.25;
PairTable.PatientID        = categorical(PairTable.Patient);
end

%% =====================================================================
%% MIXED EFFECTS MODELS
%% =====================================================================

function MMR = fit_mixed_effects_models(PairTable, PairTable_Full, nBoot, alpha)
% PairTable      — canonical-subtype patients only (for M1-M3)
% PairTable_Full — all epilepsy patients, no subtype restriction (for M4-M6)
if nargin < 3, nBoot = 0;    end
if nargin < 4, alpha = 0.05; end

%% ---- TABLE A: canonical subtype patients (M1-M3) ----
canonical3 = ["Frontal","General","Temporal"];

keepMask = ismember(string(PairTable.EpiType3), canonical3) & ...
           isfinite(PairTable.LogSpikesPerHour) & ...
           isfinite(PairTable.SignedLag_years)  & ...
           isfinite(PairTable.HasSz);
T = PairTable(keepMask, :);

T.HasSz_bin    = double(T.HasSz == 1);
EPS_SZ         = 1e-3;
T.LogSzFreq    = log(T.SzFreq + EPS_SZ);
T.AbsLag_years = abs(T.SignedLag_years);
T.LagDirection = sign(T.SignedLag_years);
T.LagDirection(T.SignedLag_years == 0) = 1;
T.EpiType3_cat = reordercats(categorical(string(T.EpiType3), canonical3), ...
    ["Temporal","Frontal","General"]);

fprintf('[Model table A - with subtypes] %d pairs, %d patients\n', ...
    height(T), numel(unique(T.Patient)));

%% ---- TABLE B: all epilepsy patients, no subtype restriction (M4-M6) ----
keepMask_full = isfinite(PairTable_Full.LogSpikesPerHour) & ...
                isfinite(PairTable_Full.SignedLag_years)  & ...
                isfinite(PairTable_Full.HasSz);
T_full = PairTable_Full(keepMask_full, :);

T_full.HasSz_bin    = double(T_full.HasSz == 1);
T_full.LogSzFreq    = log(T_full.SzFreq + EPS_SZ);
T_full.AbsLag_years = abs(T_full.SignedLag_years);
T_full.LagDirection = sign(T_full.SignedLag_years);
T_full.LagDirection(T_full.SignedLag_years == 0) = 1;

fprintf('[Model table B - no subtypes] %d pairs, %d patients\n', ...
    height(T_full), numel(unique(T_full.Patient)));

%% ---- MODEL FORMULAS ----
% With subtypes (M1-M3)
formula_M1 = ['HasSz_bin ~ ' ...
    'LogSpikesPerHour * AbsLag_years + ' ...
    'LogSpikesPerHour * LagDirection + ' ...
    'EpiType3_cat + (1|PatientID)'];
formula_M2 = ['HasSz_bin ~ ' ...
    'LogSpikesPerHour + AbsLag_years + LagDirection + ' ...
    'EpiType3_cat + (1|PatientID)'];
formula_M3 = ['LogSzFreq ~ ' ...
    'LogSpikesPerHour * AbsLag_years + ' ...
    'LogSpikesPerHour * LagDirection + ' ...
    'EpiType3_cat + (1|PatientID)'];

% No subtypes (M4-M6)
formula_M4 = ['HasSz_bin ~ ' ...
    'LogSpikesPerHour * AbsLag_years + ' ...
    'LogSpikesPerHour * LagDirection + ' ...
    '(1|PatientID)'];
formula_M5 = ['HasSz_bin ~ ' ...
    'LogSpikesPerHour + AbsLag_years + LagDirection + ' ...
    '(1|PatientID)'];
formula_M6 = ['LogSzFreq ~ ' ...
    'LogSpikesPerHour * AbsLag_years + ' ...
    'LogSpikesPerHour * LagDirection + ' ...
    '(1|PatientID)'];

%% ---- FIT MODELS ----
glme_opts = {'Distribution','Binomial','Link','logit', ...
    'FitMethod','Laplace','CovariancePattern','Diagonal'};

fprintf('\nFitting M1 (logistic + subtypes + interactions)...\n');
try; mdl_M1 = fitglme(T, formula_M1, glme_opts{:}); fprintf('M1 converged.\n'); disp(mdl_M1);
catch ME; fprintf('M1 failed: %s\n', ME.message); mdl_M1 = []; end

fprintf('\nFitting M2 (logistic + subtypes, no interactions)...\n');
try; mdl_M2 = fitglme(T, formula_M2, glme_opts{:}); fprintf('M2 converged.\n'); disp(mdl_M2);
catch ME; fprintf('M2 failed: %s\n', ME.message); mdl_M2 = []; end

fprintf('\nFitting M3 (linear + subtypes + interactions)...\n');
try; mdl_M3 = fitlme(T, formula_M3, 'FitMethod','REML'); fprintf('M3 converged.\n'); disp(mdl_M3);
catch ME; fprintf('M3 failed: %s\n', ME.message); mdl_M3 = []; end

fprintf('\nFitting M4 (logistic, no subtypes + interactions)...\n');
try; mdl_M4 = fitglme(T_full, formula_M4, glme_opts{:}); fprintf('M4 converged.\n'); disp(mdl_M4);
catch ME; fprintf('M4 failed: %s\n', ME.message); mdl_M4 = []; end

fprintf('\nFitting M5 (logistic, no subtypes, no interactions)...\n');
try; mdl_M5 = fitglme(T_full, formula_M5, glme_opts{:}); fprintf('M5 converged.\n'); disp(mdl_M5);
catch ME; fprintf('M5 failed: %s\n', ME.message); mdl_M5 = []; end

fprintf('\nFitting M6 (linear, no subtypes + interactions)...\n');
try; mdl_M6 = fitlme(T_full, formula_M6, 'FitMethod','REML'); fprintf('M6 converged.\n'); disp(mdl_M6);
catch ME; fprintf('M6 failed: %s\n', ME.message); mdl_M6 = []; end

%% ---- LRTs ----
fprintf('\n=== LIKELIHOOD RATIO TESTS ===\n');
lrt_p_A = NaN; lrt_p_B = NaN;

if ~isempty(mdl_M1) && ~isempty(mdl_M2)
    fprintf('\nLRT A: M1 vs M2 (interactions, with-subtype patients)\n');
    lrt_A = compare(mdl_M2, mdl_M1); disp(lrt_A);
    lrt_p_A = lrt_A.pValue(2);
end

if ~isempty(mdl_M4) && ~isempty(mdl_M5)
    fprintf('\nLRT B: M4 vs M5 (interactions, no-subtype restriction)\n');
    lrt_B = compare(mdl_M5, mdl_M4); disp(lrt_B);
    lrt_p_B = lrt_B.pValue(2);
end

%% ---- FIXED EFFECTS TABLES ----
fprintf('\n=== FIXED EFFECTS ===\n');
T_fe1=[]; T_fe2=[]; T_fe3=[]; T_fe4=[]; T_fe5=[]; T_fe6=[];

if ~isempty(mdl_M1)
    [b,bn,s] = fixedEffects(mdl_M1, 'Alpha', alpha);
    T_fe1 = make_fe_table_logistic(bn, b, s);
    fprintf('\nM1 (logistic + subtypes + interactions):\n'); disp(T_fe1);
end
if ~isempty(mdl_M2)
    [b,bn,s] = fixedEffects(mdl_M2, 'Alpha', alpha);
    T_fe2 = make_fe_table_logistic(bn, b, s);
    fprintf('\nM2 (logistic + subtypes, no interactions):\n'); disp(T_fe2);
end
if ~isempty(mdl_M3)
    [b,bn,s] = fixedEffects(mdl_M3, 'Alpha', alpha);
    T_fe3 = table(string(bn.Name), b, s.SE, s.tStat, s.pValue, ...
        'VariableNames',{'Term','Beta','SE','t','p'});
    fprintf('\nM3 (linear + subtypes + interactions):\n'); disp(T_fe3);
end
if ~isempty(mdl_M4)
    [b,bn,s] = fixedEffects(mdl_M4, 'Alpha', alpha);
    T_fe4 = make_fe_table_logistic(bn, b, s);
    fprintf('\nM4 (logistic, no subtypes + interactions):\n'); disp(T_fe4);
end
if ~isempty(mdl_M5)
    [b,bn,s] = fixedEffects(mdl_M5, 'Alpha', alpha);
    T_fe5 = make_fe_table_logistic(bn, b, s);
    fprintf('\nM5 (logistic, no subtypes, no interactions):\n'); disp(T_fe5);
end
if ~isempty(mdl_M6)
    [b,bn,s] = fixedEffects(mdl_M6, 'Alpha', alpha);
    T_fe6 = table(string(bn.Name), b, s.SE, s.tStat, s.pValue, ...
        'VariableNames',{'Term','Beta','SE','t','p'});
    fprintf('\nM6 (linear, no subtypes + interactions):\n'); disp(T_fe6);
end

%% ---- BOOTSTRAP (M1 and M4 only — primary logistic models) ----
[T_boot1, boot_betas1] = run_bootstrap(T,      mdl_M1, formula_M1, nBoot, alpha, 'M1');
[T_boot4, boot_betas4] = run_bootstrap(T_full, mdl_M4, formula_M4, nBoot, alpha, 'M4');
[T_boot3, boot_betas3] = run_bootstrap_lme(T,      mdl_M3, formula_M3, nBoot, alpha, 'M3');
[T_boot6, boot_betas6] = run_bootstrap_lme(T_full, mdl_M6, formula_M6, nBoot, alpha, 'M6');

if ~isempty(boot_betas1) && size(boot_betas1,1) > 10
    plot_bootstrap_diagnostics(boot_betas1, mdl_M1, 'M1');
end
if ~isempty(boot_betas4) && size(boot_betas4,1) > 10
    plot_bootstrap_diagnostics(boot_betas4, mdl_M4, 'M4');
end

%% ---- BUNDLE ----
MMR.ModelTable      = T;
MMR.ModelTable_Full = T_full;

MMR.mdl_M1 = mdl_M1;  MMR.mdl_M2 = mdl_M2;  MMR.mdl_M3 = mdl_M3;
MMR.mdl_M4 = mdl_M4;  MMR.mdl_M5 = mdl_M5;  MMR.mdl_M6 = mdl_M6;

MMR.FE_M1 = T_fe1;  MMR.FE_M2 = T_fe2;  MMR.FE_M3 = T_fe3;
MMR.FE_M4 = T_fe4;  MMR.FE_M5 = T_fe5;  MMR.FE_M6 = T_fe6;

MMR.BootstrapBetas1 = boot_betas1;  MMR.BootstrapTable1 = T_boot1;
MMR.BootstrapBetas4 = boot_betas4;  MMR.BootstrapTable4 = T_boot4;
MMR.BootstrapTable3 = T_boot3;      MMR.BootstrapBetas3 = boot_betas3;
MMR.BootstrapTable6 = T_boot6;      MMR.BootstrapBetas6 = boot_betas6;

MMR.LRT_p_A = lrt_p_A;   % M1 vs M2, chi^2(2)
MMR.LRT_p_B = lrt_p_B;   % M4 vs M5, chi^2(2)

fprintf('\nDone. Primary models: M1 (with subtypes) and M4 (no subtypes).\n');
end


function T_fe = make_fe_table_logistic(bn, b, s)
T_fe = table(string(bn.Name), b, s.SE, s.tStat, s.pValue, ...
    exp(b), exp(b - 1.96*s.SE), exp(b + 1.96*s.SE), ...
    'VariableNames',{'Term','Beta','SE','t','p','OR','OR_lo','OR_hi'});
end

function [T_boot, boot_betas] = run_bootstrap(T, mdl, formula, nBoot, alpha, label)
T_boot     = [];
boot_betas = [];
if isempty(mdl) || nBoot == 0, return; end

fprintf('\nBootstrapping %s (%d iterations)...\n', label, nBoot);
patients   = unique(T.PatientID);
nPat       = numel(patients);
nFixed     = size(fixedEffects(mdl), 1);
boot_betas = nan(nBoot, nFixed);

parfor b = 1:nBoot
    idx      = randi(nPat, nPat, 1);
    bootPats = patients(idx);
    Tboot    = cell(nPat, 1);
    for k = 1:nPat
        Tboot{k} = T(T.PatientID == bootPats(k), :);
        Tboot{k}.PatientID = categorical(repmat(k, height(Tboot{k}), 1));
    end
    Tboot = vertcat(Tboot{:});
    try
        mboot = fitglme(Tboot, formula, 'Distribution','Binomial','Link','logit', ...
            'FitMethod','Laplace','CovariancePattern','Diagonal');
        boot_betas(b,:) = fixedEffects(mboot)';
    catch
    end
end

converged  = all(isfinite(boot_betas), 2);
boot_betas = boot_betas(converged, :);
fprintf('%s bootstrap: %d/%d converged (%.1f%%)\n', ...
    label, sum(converged), nBoot, 100*mean(converged));

ci_lo = prctile(boot_betas, 100*(alpha/2),   1);
ci_hi = prctile(boot_betas, 100*(1-alpha/2), 1);
[beta_obs, betanames_obs] = fixedEffects(mdl);

T_boot = table(string(betanames_obs.Name), beta_obs, ci_lo(:), ci_hi(:), ...
    exp(beta_obs), exp(ci_lo(:)), exp(ci_hi(:)), ...
    'VariableNames',{'Term','Beta','Boot_CI_lo','Boot_CI_hi','OR','OR_CI_lo','OR_CI_hi'});
fprintf('%s bootstrapped ORs:\n', label); disp(T_boot);
end

function plot_bootstrap_diagnostics(boot_betas, mdl, label)
[beta_diag, betanames_diag] = fixedEffects(mdl);
term_names = string(betanames_diag.Name);
nTerms = numel(term_names);
nCols  = ceil(sqrt(nTerms));
nRows  = ceil(nTerms/nCols);
figDiag = figure('Color','w','Position',[100 100 300*nCols 250*nRows], ...
    'Name',sprintf('Bootstrap — %s', label));
tl = tiledlayout(figDiag, nRows, nCols, 'TileSpacing','compact','Padding','loose');
for k = 1:nTerms
    ax = nexttile(tl);
    histogram(ax, boot_betas(:,k), 40, ...
        'FaceColor',[0.2 0.5 0.8],'EdgeColor','none','FaceAlpha',0.7);
    xline(ax, beta_diag(k), 'r-', 'LineWidth', 2);
    xline(ax, prctile(boot_betas(:,k), 2.5),  'k--', 'LineWidth', 1.5);
    xline(ax, prctile(boot_betas(:,k), 97.5), 'k--', 'LineWidth', 1.5);
    title(ax, term_names(k), 'FontSize',9,'Interpreter','none');
    xlabel(ax, '\beta','FontSize',8); box(ax,'off');
end
sgtitle(figDiag, sprintf('Bootstrap — %s (red=observed, dashed=95%% CI)', label), 'FontSize',11);
end

%% =====================================================================
%% FIGURE 1
%% =====================================================================

function [f1, Fig1Stats] = make_fig1_controls(Views, EPS_RATE, Y_ZERO, Y_LIMS, TITLE_Y_OFFSET, nBoot, alpha)
SessionLevel = Views.SessionLevelSpikeRates;
Report       = Views.ReportForKeptSessions;

ReportSlim = resolve_reported_spike_status(Report);
JoinA = innerjoin(SessionLevel(:,{'Patient','Session','SpikesPerHour'}), ReportSlim, ...
    'Keys',{'Patient','Session'});

x_abs = JoinA.SpikesPerHour(JoinA.ReportStatus=="absent");
x_pre = JoinA.SpikesPerHour(JoinA.ReportStatus=="present");
pA      = ranksum(x_abs, x_pre, 'method','approx');
effectA = cliff_delta(x_pre, x_abs);

[med_abs, lo_abs, hi_abs] = bootstrap_median_ci(x_abs, nBoot, alpha);
[med_pre, lo_pre, hi_pre] = bootstrap_median_ci(x_pre, nBoot, alpha);

Y_A = add_y_jitter_eps(to_log10_per_hour([x_abs; x_pre], EPS_RATE), log10(EPS_RATE), Y_LIMS, 0.02);
G_A = [repmat("Absent",numel(x_abs),1); repmat("Present",numel(x_pre),1)];

Sub = Views.PatientSpikeSz_Typed(:,{'EpiType3','MeanSpikeRate_perHour'});
Sub.Properties.VariableNames{'EpiType3'} = 'EpiType4';
Y_C = add_y_jitter_eps(to_log10_per_hour(Sub.MeanSpikeRate_perHour, EPS_RATE), Y_ZERO, Y_LIMS, 0.02);

[p_kw, tbl_kw] = kruskalwallis(Sub.MeanSpikeRate_perHour, Sub.EpiType4, 'off');
eta2_kw = tbl_kw{2,2} / tbl_kw{end,2};

f1 = figure('Color','w','Position',[60 60 950 520]);
tiledlayout(f1,1,2,'TileSpacing','compact','Padding','loose');

axA = nexttile; hold(axA,'on'); box(axA,'off'); grid(axA,'on');
boxchart(axA, categorical(G_A), Y_A, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(axA, categorical(G_A), Y_A, 18, 'filled','MarkerFaceAlpha',0.18);
yline(axA, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axA, Y_LIMS); ylabel(axA,'Spikes/hour (log scale)');
set_log10_ticks(axA,'y',EPS_RATE,Y_LIMS);
add_median_ci_overlay(axA, 1, med_abs, lo_abs, hi_abs, EPS_RATE);
add_median_ci_overlay(axA, 2, med_pre, lo_pre, hi_pre, EPS_RATE);
t = title(axA,'A. Reported presence or absence of spikes');
add_sigbar(axA, 1, 2, Y_LIMS(2)-0.08*range(Y_LIMS), p_label(pA));
set(axA,'FontSize',20);
t.Units='normalized'; t.Position(2) = t.Position(2)+TITLE_Y_OFFSET;
labelsA = string(axA.XTickLabel);
labelsA(labelsA=="Absent")  = sprintf('Absent (N=%d)',  nnz(isfinite(x_abs)));
labelsA(labelsA=="Present") = sprintf('Present (N=%d)', nnz(isfinite(x_pre)));
axA.XTickLabel = labelsA; axA.XTickLabelRotation = 20;

axC = nexttile; hold(axC,'on'); box(axC,'off'); grid(axC,'on');
boxchart(axC, Sub.EpiType4, Y_C, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(axC, Sub.EpiType4, Y_C, 18, 'filled','MarkerFaceAlpha',0.18);
yline(axC, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(axC, Y_LIMS); ylabel(axC,'Spikes/hour (log scale)');
set_log10_ticks(axC,'y',EPS_RATE,Y_LIMS);
cats = categories(Sub.EpiType4);
for k = 1:numel(cats)
    xg = Sub.MeanSpikeRate_perHour(Sub.EpiType4 == cats{k});
    [mg, log, hig] = bootstrap_median_ci(xg, nBoot, alpha);
    add_median_ci_overlay(axC, k, mg, log, hig, EPS_RATE);
end
t = title(axC,'B. Epilepsy subtype');
set(axC,'FontSize',20);
t.Units='normalized'; t.Position(2) = t.Position(2)+TITLE_Y_OFFSET;
SubStats = Views.Canonical3_Stats;
labelsC = string(axC.XTickLabel);
for i = 1:height(SubStats)
    lab = string(SubStats.EpiType4(i));
    labelsC(labelsC==lab) = sprintf('%s (N=%d)', lab, SubStats.GroupCount(i));
end
axC.XTickLabel = labelsC; axC.XTickLabelRotation = 20;

SubtypePairs = Views.Canonical3_Pairs;
p_pair_bonf  = Views.PvalsPairwiseBonf;
cats3 = categorical(string(categories(Sub.EpiType4)));
yTop = Y_LIMS(2); yStep = 0.08*range(Y_LIMS); y0 = yTop - 0.05*range(Y_LIMS);
for i = 1:size(SubtypePairs,1)
    x1  = find(cats3 == categorical(string(SubtypePairs(i,1))));
    x2  = find(cats3 == categorical(string(SubtypePairs(i,2))));
    pval = p_pair_bonf(i);
    if isnan(pval), continue; end
    if pval<1e-3, lab="***"; elseif pval<1e-2, lab="**"; elseif pval<5e-2, lab="*"; else, lab="ns"; end
    add_sigbar(axC, x1, x2, y0-(i-1)*yStep, lab);
end

Fig1Stats.p_rankSum_A   = pA;
Fig1Stats.effectA_cliff = effectA;
Fig1Stats.m_pre = med_pre; Fig1Stats.lo_pre = lo_pre; Fig1Stats.hi_pre = hi_pre;
Fig1Stats.m_abs = med_abs; Fig1Stats.lo_abs = lo_abs; Fig1Stats.hi_abs = hi_abs;
Fig1Stats.p_kw_C    = p_kw;
Fig1Stats.eta2_kw_C = eta2_kw;
Fig1Stats.p_pair_bonf = Views.PvalsPairwiseBonf;

Sub2 = Views.Canonical3_SubsetTable;
subNames = string(categories(Sub2.EpiType4));
subMed = nan(numel(subNames),1); subCI = nan(numel(subNames),3);
for k=1:numel(subNames)
    x = Sub2.MeanSpikeRate_perHour(Sub2.EpiType4==subNames(k));
    [m,lo,hi] = bootstrap_median_ci(x, nBoot, alpha);
    subMed(k)=m; subCI(k,:)=[m lo hi];
end
Fig1Stats.SubtypeStatsTable = table(subNames, subMed, subCI(:,2), subCI(:,3), ...
    'VariableNames',{'Group','Median','CI_lo','CI_hi'});
end

%% =====================================================================
%% MAIN MODEL FIGURE (panels C and D only)
%% =====================================================================

function FigMain = make_model_figure(MMR, outPath)
if nargin < 2, outPath = ''; end

FONT_SIZE = 20;
T   = MMR.ModelTable;
mdl = MMR.mdl_M1;

FigMain = figure('Color','w','Position',[60 60 1300 560]);

% Two panels side by side, more room now that A/B are removed
ax_w     = 0.35;
ax_h     = 0.72;
col1_bot = 0.21;   % C left edge — extra room for long ytick labels
col2     = 0.64;
row1     = 0.12;

axC = axes('Position',[col1_bot, row1, ax_w, ax_h]);
axD = axes('Position',[col2,     row1, ax_w, ax_h]);

%% Panel C: Forest plot
axes(axC); hold(axC,'on'); box(axC,'off'); grid(axC,'on');

[beta_m, betanames_m, stats_m] = fixedEffects(mdl);
raw_names = string(betanames_m.Name);

if ~isempty(MMR.BootstrapTable1)
    BT = MMR.BootstrapTable1; ci_label = '95% Bootstrap CI';
    OR_C    = nan(numel(raw_names),1);
    OR_lo_C = nan(numel(raw_names),1);
    OR_hi_C = nan(numel(raw_names),1);
    for k = 1:numel(raw_names)
        bt_row = BT(string(BT.Term)==raw_names(k),:);
        if ~isempty(bt_row)
            OR_C(k)=bt_row.OR; OR_lo_C(k)=bt_row.OR_CI_lo; OR_hi_C(k)=bt_row.OR_CI_hi;
        else
            OR_C(k)=exp(beta_m(k));
            OR_lo_C(k)=exp(beta_m(k)-1.96*stats_m.SE(k));
            OR_hi_C(k)=exp(beta_m(k)+1.96*stats_m.SE(k));
        end
    end
else
    OR_C    = exp(beta_m);
    OR_lo_C = exp(beta_m-1.96*stats_m.SE);
    OR_hi_C = exp(beta_m+1.96*stats_m.SE);
    ci_label = '95% Laplace CI';
end

disp_names = raw_names;
disp_names(disp_names=="(Intercept)")                   = "Intercept";
disp_names(disp_names=="LogSpikesPerHour")              = "Log spike rate";
disp_names(disp_names=="AbsLag_years")                  = "EEG-visit gap (years)";
disp_names(disp_names=="LagDirection")                  = "Visit after vs before EEG";
disp_names(disp_names=="EpiType3_cat_Frontal")          = "Frontal vs Temporal";
disp_names(disp_names=="EpiType3_cat_General")          = "Generalized vs Temporal";
disp_names(disp_names=="LogSpikesPerHour:AbsLag_years") = "Spike rate effect per year of gap";
disp_names(disp_names=="LogSpikesPerHour:LagDirection") = "Spike rate effect: visit before or after";

isInt       = (disp_names=="Intercept");
OR_C        = OR_C(~isInt);
OR_lo_C     = OR_lo_C(~isInt);
OR_hi_C     = OR_hi_C(~isInt);
pvals_C     = stats_m.pValue(~isInt);
disp_names  = disp_names(~isInt);
nTerms      = numel(OR_C);
plot_order  = nTerms:-1:1;

for k = 1:nTerms
    idx = plot_order(k);
    col = [0.1 0.3 0.7]; if pvals_C(idx)>=0.05, col=[0.6 0.6 0.6]; end
    plot(axC,[OR_lo_C(idx),OR_hi_C(idx)],[k k],'-','Color',col,'LineWidth',2.5);
    scatter(axC,OR_C(idx),k,100,col,'filled');
    p = pvals_C(idx);
    if p<0.001, pstr='p<0.001'; elseif p<0.05, pstr=sprintf('p=%.3f',p); else, pstr=sprintf('p=%.2f',p); end
    text(axC,OR_hi_C(idx)+0.005,k,pstr,'FontSize',FONT_SIZE-5,'VerticalAlignment','middle');
end
xline(axC,1,'k--','LineWidth',1.5);
set(axC,'FontSize',FONT_SIZE);
set(axC,'YTick',1:nTerms,'YTickLabel',disp_names(plot_order),'FontSize',FONT_SIZE-4);
xlabel(axC,sprintf('Odds Ratio (%s)',ci_label),'FontSize',FONT_SIZE);
th = title(axC,{'A. Spike rate, epilepsy type, and EEG-visit gap', 'predict seizure occurrence'}, ...
    'FontSize',FONT_SIZE,'FontWeight','bold');
th.Units = 'normalized'; th.Position(2) = th.Position(2) + 0.02; th.Position(1) = th.Position(1) - 0.13;
all_ors = [OR_lo_C; OR_hi_C];
xlim(axC,[max(0.4,prctile(all_ors,2)), min(2.0,prctile(all_ors,98))]);

%% Panel D: Predicted probability vs spike rate by lag distance
axes(axD); hold(axD,'on'); box(axD,'off'); grid(axD,'on');

spike_grid_raw = linspace(0,50,200);
EPS_SPIKE      = 1e-3;
spike_grid_log = log(spike_grid_raw + EPS_SPIKE);
lag_vals       = [6/12, 2, 4];
lag_labels     = ["6 months","2 years","4 years"];
lag_colors     = [0.05 0.30 0.70;
                  0.15 0.50 0.80;
                  0.40 0.65 0.85];

[beta_pred, betanames_pred] = fixedEffects(mdl);
bnames      = string(betanames_pred.Name);
b_intercept = beta_pred(bnames=="(Intercept)");
b_spike     = beta_pred(bnames=="LogSpikesPerHour");
b_abslag    = beta_pred(bnames=="AbsLag_years");
b_dir       = beta_pred(bnames=="LagDirection");
b_int_prox  = beta_pred(bnames=="LogSpikesPerHour:AbsLag_years");
b_int_dir   = beta_pred(bnames=="LogSpikesPerHour:LagDirection");
dir_val     = 1;   % visit after EEG; temporal epilepsy reference

for k = 1:numel(lag_vals)
    lag = lag_vals(k);
    eta = b_intercept ...
        + b_spike   .* spike_grid_log ...
        + b_abslag  .* lag ...
        + b_dir     .* dir_val ...
        + b_int_prox.* spike_grid_log .* lag ...
        + b_int_dir .* spike_grid_log .* dir_val;
    prob = 1 ./ (1 + exp(-eta));
    plot(axD, spike_grid_raw, prob, '-', ...
        'Color',lag_colors(k,:),'LineWidth',2.5,'DisplayName',lag_labels(k));
end

% OR annotations
x_anno = 22;
for k = 1:numel(lag_vals)
    lag_k  = lag_vals(k);
    or_k   = exp(b_spike + b_int_prox*lag_k);
    eta_k  = b_intercept + b_spike*log(x_anno+EPS_SPIKE) + b_abslag*lag_k + ...
        b_dir*dir_val + b_int_prox*log(x_anno+EPS_SPIKE)*lag_k + ...
        b_int_dir*log(x_anno+EPS_SPIKE)*dir_val;
    prob_k = 1/(1+exp(-eta_k));
    text(axD, x_anno+1.0, prob_k+0.01, sprintf('OR=%.2f',or_k), ...
        'FontSize',FONT_SIZE-7,'Color',lag_colors(k,:), ...
        'VerticalAlignment','middle','HorizontalAlignment','left');
end

xlabel(axD,'Spike rate (spikes/hour)','FontSize',FONT_SIZE);
ylabel(axD,'P(seizure reported at visit)','FontSize',FONT_SIZE);
th = title(axD,{'B. Spike rates are most predictive', ...
    'when EEG is obtained close to the visit'}, ...
    'FontSize',FONT_SIZE,'FontWeight','bold');
th.Units = 'normalized'; th.Position(2) = th.Position(2) + 0.02;
lg = legend(axD,'Location','southeast','FontSize',FONT_SIZE-6);
title(lg,'EEG-visit lag');
xlim(axD,[0 30]); ylim(axD,[0.3 0.61]);
set(axD,'FontSize',FONT_SIZE);

%% Save
if strlength(string(outPath)) > 0
    if ~exist(fileparts(outPath),'dir'), mkdir(fileparts(outPath)); end
    exportgraphics(FigMain, outPath, 'Resolution',300);
    fprintf('Saved main model figure: %s\n', outPath);
end
end


%% =====================================================================
%% SUPPLEMENTAL FIGURE: EEG-visit context panels (old A and B)
%% =====================================================================

function FigSup = make_figSup_lag(MMR, Vuniq, ReportForKeptSessions, outPath)
if nargin < 4, outPath = ''; end

FONT_SIZE = 20;
T       = MMR.ModelTable;
refDate = datetime(2000,1,1);

FigSup = figure('Color','w','Position',[60 60 1300 560]);

ax_w  = 0.38;
ax_h  = 0.72;
col1  = 0.07;
col2  = 0.57;
row1  = 0.12;

axA = axes('Position',[col1, row1, ax_w, ax_h]);
axB = axes('Position',[col2, row1, ax_w, ax_h]);

%% Panel A: Seizure occurrence over time + median seizure frequency (right y-axis)
cohort_patients = unique(T.Patient);
V_cohort = Vuniq(ismember(Vuniq.Patient,cohort_patients) & ...
    (Vuniq.HasSz==0|Vuniq.HasSz==1),:);
V_cohort.YearsSinceFirst = days(V_cohort.VisitDate-refDate)/365.25;
bin_edges_sz   = [0, 1, 2, 3, 4];
bin_centers_sz = (bin_edges_sz(1:end-1) + bin_edges_sz(2:end)) / 2;
nBins_sz       = numel(bin_centers_sz);

bin_prop_sz  = nan(nBins_sz,1); bin_lo_sz   = nan(nBins_sz,1);
bin_hi_sz    = nan(nBins_sz,1); bin_n_sz    = zeros(nBins_sz,1);
bin_med_freq = nan(nBins_sz,1); bin_lo_freq = nan(nBins_sz,1);
bin_hi_freq  = nan(nBins_sz,1);
EPS_FREQ = 1e-3;

for b = 1:nBins_sz
    mask    = V_cohort.YearsSinceFirst >= bin_edges_sz(b) & ...
              V_cohort.YearsSinceFirst <  bin_edges_sz(b+1);
    vals_hz = V_cohort.HasSz(mask);   vals_hz = vals_hz(isfinite(vals_hz));
    vals_fr = V_cohort.Freq_R1(mask); vals_fr = vals_fr(isfinite(vals_fr));
    if numel(vals_hz) >= 10
        bin_prop_sz(b) = mean(vals_hz);
        bin_n_sz(b)    = numel(vals_hz);
        boot_p = nan(500,1);
        for bb = 1:500
            boot_p(bb) = mean(vals_hz(randi(numel(vals_hz),numel(vals_hz),1)));
        end
        bin_lo_sz(b) = prctile(boot_p,2.5);
        bin_hi_sz(b) = prctile(boot_p,97.5);
    end
    if numel(vals_fr) >= 10
        bin_med_freq(b) = median(vals_fr,'omitnan');
        boot_f = nan(500,1);
        for bb = 1:500
            boot_f(bb) = median(vals_fr(randi(numel(vals_fr),numel(vals_fr),1)),'omitnan');
        end
        bin_lo_freq(b) = prctile(boot_f,2.5);
        bin_hi_freq(b) = prctile(boot_f,97.5);
    end
end

validA    = isfinite(bin_prop_sz);
validFreq = isfinite(bin_med_freq);
COL_PROP  = [0.8 0.3 0.1];
COL_FREQ  = [0.1 0.55 0.55];

axes(axA); hold(axA,'on'); box(axA,'off'); grid(axA,'on');
yyaxis(axA,'left');
patch(axA, [bin_centers_sz(validA), fliplr(bin_centers_sz(validA))], ...
    [bin_lo_sz(validA)', fliplr(bin_hi_sz(validA)')], ...
    COL_PROP,'FaceAlpha',0.2,'EdgeColor','none');
plot(axA, bin_centers_sz(validA), bin_prop_sz(validA), 'o-', ...
    'Color',COL_PROP,'LineWidth',2,'MarkerFaceColor',COL_PROP,'MarkerSize',6);
for b = find(validA)'
    scatter(axA, bin_centers_sz(b), bin_prop_sz(b), bin_n_sz(b)/5, ...
        COL_PROP,'filled','MarkerFaceAlpha',0.25);
end
ylim(axA,[0 1]);
ylabel(axA,'Proportion with seizures','FontSize',FONT_SIZE,'Color',COL_PROP);
axA.YAxis(1).Color = COL_PROP;

yyaxis(axA,'right');
freq_log_med = log10(bin_med_freq  + EPS_FREQ);
freq_log_lo  = log10(bin_lo_freq   + EPS_FREQ);
freq_log_hi  = log10(bin_hi_freq   + EPS_FREQ);
patch(axA, [bin_centers_sz(validFreq), fliplr(bin_centers_sz(validFreq))], ...
    [freq_log_lo(validFreq)', fliplr(freq_log_hi(validFreq)')], ...
    COL_FREQ,'FaceAlpha',0.18,'EdgeColor','none');
plot(axA, bin_centers_sz(validFreq), freq_log_med(validFreq), 's--', ...
    'Color',COL_FREQ,'LineWidth',1.8,'MarkerFaceColor',COL_FREQ,'MarkerSize',5);
Y_LIM_FREQ = [-2 2];
ylim(axA, Y_LIM_FREQ);
set_log10_ticks(axA,'y',EPS_FREQ,Y_LIM_FREQ);
ylabel(axA,'Median sz/month (log scale)','FontSize',FONT_SIZE,'Color',COL_FREQ);
axA.YAxis(2).Color = COL_FREQ;

xlabel(axA,'Years since Jan 2000 (calendar anchor)','FontSize',FONT_SIZE);
th = title(axA,'A. Seizure burden tends to decrease over time', ...
    'FontSize',FONT_SIZE,'FontWeight','bold');
th.Units = 'normalized'; th.Position(2) = th.Position(2) + 0.02;
set(axA,'FontSize',FONT_SIZE);

%% Panel B: EEG-visit lag distribution
axes(axB); hold(axB,'on'); box(axB,'off'); grid(axB,'on');

abs_lags = T.AbsLag_years;
histogram(axB, abs_lags, 40,'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.6, ...
    'EdgeColor','none','Normalization','probability');
xline(axB, 1,'k--','LineWidth',1.2);
xlabel(axB,'Absolute EEG-visit gap (years)','FontSize',FONT_SIZE);
ylabel(axB,'Proportion of pairs','FontSize',FONT_SIZE);
th = title(axB,'B. EEG and visit are often separated by years', ...
    'FontSize',FONT_SIZE,'FontWeight','bold');
th.Units = 'normalized'; th.Position(2) = th.Position(2) + 0.02;
set(axB,'FontSize',FONT_SIZE);

%% Save
if strlength(string(outPath)) > 0
    if ~exist(fileparts(outPath),'dir'), mkdir(fileparts(outPath)); end
    exportgraphics(FigSup, outPath, 'Resolution',300);
    fprintf('Saved supplemental lag figure: %s\n', outPath);
end
end

%% =====================================================================
%% FIGURE S2
%% =====================================================================

function fS2 = make_figS2_sz_by_reported_spikes(Views, SzFreqPerPatient, nBoot, alpha, xLims_log10)
RS = resolve_reported_spike_status(Views.ReportForKeptSessions);
[gp, pid] = findgroups(RS.Patient);
hasPresent = splitapply(@(x) any(x=="present"), RS.ReportStatus, gp);
hasAbsent  = splitapply(@(x) any(x=="absent"),  RS.ReportStatus, gp);
RptP = table(pid, hasPresent, hasAbsent, 'VariableNames',{'Patient','HasPresent','HasAbsent'});

S2 = innerjoin(SzFreqPerPatient, RptP, 'Keys','Patient');
EpPatients = Views.PatientLevelSpikeRates.Patient(Views.IsEpilepsyMask);
S2 = innerjoin(S2, table(EpPatients,'VariableNames',{'Patient'}), 'Keys','Patient');

freq_allAbsent  = S2.MeanSzFreq(S2.HasAbsent & ~S2.HasPresent);
freq_anyPresent = S2.MeanSzFreq(S2.HasPresent);
freq_allAbsent  = freq_allAbsent(isfinite(freq_allAbsent));
freq_anyPresent = freq_anyPresent(isfinite(freq_anyPresent));
n_allAbsent  = numel(freq_allAbsent);
n_anyPresent = numel(freq_anyPresent);

p = ranksum(freq_allAbsent, freq_anyPresent, 'method','approx');
assert(p<0.001)

EPS_FREQ = 1e-3; Y_ZERO = log10(EPS_FREQ);
Y = [to_log10_per_month(freq_allAbsent, EPS_FREQ); to_log10_per_month(freq_anyPresent,EPS_FREQ)];
Y = add_y_jitter_eps(Y, Y_ZERO, xLims_log10, 0.02);
G = [repmat("All EEGs: no spikes",numel(freq_allAbsent),1); ...
     repmat("≥1 EEG: spikes present",numel(freq_anyPresent),1)];

[med1,lo1,hi1] = bootstrap_median_ci(freq_allAbsent,  nBoot, alpha);
[med2,lo2,hi2] = bootstrap_median_ci(freq_anyPresent, nBoot, alpha);

fS2 = figure('Color','w','Position',[100 100 800 520]);
ax = axes(fS2); hold(ax,'on'); box(ax,'off'); grid(ax,'on');
boxchart(ax,categorical(G),Y,'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(ax,categorical(G),Y,18,'filled','MarkerFaceAlpha',0.18);
yline(ax,Y_ZERO,':','Color',[0.4 0.4 0.4],'LineWidth',1.2);
ylim(ax,xLims_log10); ylabel(ax,'Seizures/month (log scale)');
set_log10_ticks(ax,'y',EPS_FREQ,xLims_log10);
add_median_ci_overlay_month(ax,1,med1,lo1,hi1,EPS_FREQ);
add_median_ci_overlay_month(ax,2,med2,lo2,hi2,EPS_FREQ);

yl=ylim(ax); yMaxData=max(Y(isfinite(Y)));
yBar=yMaxData+0.06*range(yl);
yNeedTop=yBar+0.10*range(yl);
if yNeedTop>yl(2), ylim(ax,[yl(1) yNeedTop]); end
add_sigbar(ax,1,2,yBar,p_label(p));
t=title(ax,'Mean seizure frequency by reported spikes across EEGs');
t.Units='normalized'; t.Position(2)=1.03;
set(ax,'FontSize',20);

labels = string(ax.XTickLabel);
labels(labels=="All EEGs: no spikes")    = sprintf('All EEGs: no spikes (N=%d)',    n_allAbsent);
labels(labels=="≥1 EEG: spikes present") = sprintf('≥1 EEG: spikes present (N=%d)', n_anyPresent);
ax.XTickLabel = labels; ax.XTickLabelRotation = 20;

fprintf(['\nMedian [95%% CI] seizure frequency: %.2f [%.2f-%.2f] (no spikes) ' ...
    'vs %.2f [%.2f-%.2f] (spikes present) (p<0.001, Cliff''s d=%.2f)\n'], ...
    med1,lo1,hi1,med2,lo2,hi2,cliff_delta(freq_allAbsent,freq_anyPresent));
end

%% =====================================================================
%% TABLE 1
%% =====================================================================

function Table1_flat = build_table1_flat(Views, SzFreqPerPatient, Vuniq, EPS_RATE, nBoot, alpha)
Rk = Views.ReportForKeptSessions;
PL = Views.PatientLevelSpikeRates;
AllPatients = PL.Patient;
N_total = numel(AllPatients);

birth_str = strtrim(string(Rk.deid_birth_date));
isMiss    = birth_str=="" | birth_str=="null" | birth_str=="[null]";
birth_dt  = NaT(size(birth_str));
birth_dt(~isMiss) = datetime(birth_str(~isMiss),'InputFormat','yyyy-MM-dd');
refDate   = datetime(2000,1,1);
age_first = NaN(size(birth_dt));
age_first(~isnat(birth_dt)) = days(refDate - birth_dt(~isnat(birth_dt)))/365.25;
[ga, pidA] = findgroups(Rk.Patient);
AgeFirst = splitapply(@local_min_omitnan, age_first, ga);
AgeTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), ...
    table(pidA,AgeFirst,'VariableNames',{'Patient','AgeFirst'}),'Keys','Patient');
age_vec = AgeTable.AgeFirst;
age_med = median(age_vec,'omitnan'); age_q = prctile(age_vec,[25,75]);

sex_raw = upper(strtrim(string(Rk.nlp_gender)));
[gs, pidS] = findgroups(Rk.Patient);
sex_per = splitapply(@local_first_nonmissing, sex_raw, gs);
SexTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), ...
    table(pidS,sex_per,'VariableNames',{'Patient','SexCode'}),'Keys','Patient');
n_f = nnz(SexTable.SexCode=="F");
n_m = nnz(SexTable.SexCode=="M");
n_u = N_total - n_f - n_m;

E3    = strtrim(string(PL.EpiType3));
espec = strtrim(string(PL.EpilepsySpecific));
isTemp    = (E3=="Temporal") & Views.IsEpilepsyMask;
isFront   = (E3=="Frontal")  & Views.IsEpilepsyMask;
isGen     = (E3=="General")  & Views.IsEpilepsyMask;
isCanon   = isTemp | isFront | isGen;
isUnknown = Views.IsEpilepsyMask & ~isCanon & ...
    (ismissing(espec)|espec==""|espec=="Unclassified or Unspecified"|espec=="Unknown or MRN not found");
isOther   = Views.IsEpilepsyMask & ~isCanon & ~isUnknown;
n_epi     = nnz(Views.IsEpilepsyMask);
n_temp    = nnz(isTemp); n_front = nnz(isFront); n_gen    = nnz(isGen);
n_other   = nnz(isOther); n_subunk = nnz(isUnknown);

[gv, pidV] = findgroups(Vuniq.Patient);
nVisits = splitapply(@(x) numel(unique(x)), Vuniq.VisitDate, gv);
VisitsTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), ...
    table(pidV,nVisits,'VariableNames',{'Patient','NumVisits'}),'Keys','Patient');
vis_med = median(VisitsTable.NumVisits,'omitnan'); vis_q = prctile(VisitsTable.NumVisits,[25,75]);

Sess = Views.SessionsForFigures;
[ge, pidE] = findgroups(Sess.Patient);
nEEG = splitapply(@(x) numel(unique(x)), Sess.Session, ge);
EEGTable = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), ...
    table(pidE,nEEG,'VariableNames',{'Patient','NumEEG'}),'Keys','Patient');
eeg_med = median(EEGTable.NumEEG,'omitnan'); eeg_q = prctile(EEGTable.NumEEG,[25,75]);

SzJ = innerjoin(table(AllPatients,'VariableNames',{'Patient'}), SzFreqPerPatient,'Keys','Patient');
sf_vec = SzJ.MeanSzFreq(isfinite(SzJ.MeanSzFreq));
sf_med = median(sf_vec,'omitnan'); sf_q = prctile(sf_vec,[25,75]);
[~,sf_lo,sf_hi] = bootstrap_median_ci(sf_vec, nBoot, alpha);

sr_vec = PL.MeanSpikeRate_perHour(isfinite(PL.MeanSpikeRate_perHour));
sr_med = median(sr_vec,'omitnan'); sr_q = prctile(sr_vec,[25,75]);
[~,sr_lo,sr_hi] = bootstrap_median_ci(sr_vec, nBoot, alpha);

RS_all = resolve_reported_spike_status(Views.ReportForKeptSessions);
status = string(RS_all.ReportStatus); status(ismissing(RS_all.ReportStatus)) = "unknown";
RS_all.SpikeStatus = categorical(status,["absent","present","unknown"]);
n_rep_pre  = nnz(RS_all.SpikeStatus=="present");
n_rep_abs  = nnz(RS_all.SpikeStatus=="absent");
n_rep_unk  = nnz(RS_all.SpikeStatus=="unknown");
n_eegs_all = height(RS_all);

% Patient-level spike status
[gp_rs, pid_rs] = findgroups(RS_all.Patient);
hasPresent = splitapply(@(x) any(string(x)=="present"), RS_all.SpikeStatus, gp_rs);
hasAbsent  = splitapply(@(x) any(string(x)=="absent"),  RS_all.SpikeStatus, gp_rs);
% Present: any EEG with spikes reported
% Absent: no EEG with spikes reported, but at least one EEG with a report
% Unknown: no EEG with any report
pat_spike_status = repmat("unknown", numel(pid_rs), 1);
pat_spike_status(hasAbsent & ~hasPresent)  = "absent";
pat_spike_status(hasPresent)               = "present";
n_pats_pre = nnz(pat_spike_status=="present");
n_pats_abs = nnz(pat_spike_status=="absent");
n_pats_unk = nnz(pat_spike_status=="unknown");
n_pats_rs  = numel(pid_rs);

OutVar = {}; OutStat = {};
OutVar{end+1,1}="Total N patients"; OutStat{end+1,1}=sprintf('%d',N_total);
OutVar{end+1,1}="Age at first visit (years)"; OutStat{end+1,1}=sprintf('%.1f (%.1f-%.1f)',age_med,age_q(1),age_q(2));
OutVar{end+1,1}="Sex"; OutStat{end+1,1}="";
OutVar{end+1,1}="    Women";         OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_f,100*n_f/N_total);
OutVar{end+1,1}="    Men";           OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_m,100*n_m/N_total);
OutVar{end+1,1}="    Unknown/Other"; OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_u,100*n_u/N_total);
OutVar{end+1,1}="Epilepsy subtype"; OutStat{end+1,1}="";
OutVar{end+1,1}="    Temporal lobe"; OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_temp, 100*n_temp /max(1,n_epi));
OutVar{end+1,1}="    Frontal lobe";  OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_front,100*n_front/max(1,n_epi));
OutVar{end+1,1}="    Generalized";   OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_gen,  100*n_gen  /max(1,n_epi));
OutVar{end+1,1}="    Other";         OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_other,100*n_other/max(1,n_epi));
OutVar{end+1,1}="    Unknown";       OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_subunk,100*n_subunk/max(1,n_epi));
OutVar{end+1,1}="Number of clinic visits"; OutStat{end+1,1}=sprintf('%.1f (%.1f-%.1f)',vis_med,vis_q(1),vis_q(2));
OutVar{end+1,1}="Number of EEGs"; OutStat{end+1,1}=sprintf('%.1f (%.1f-%.1f)',eeg_med,eeg_q(1),eeg_q(2));
OutVar{end+1,1}="Mean seizure frequency (seizures/month)"; OutStat{end+1,1}=sprintf('%.2f (%.2f-%.2f); median CI [%.2f-%.2f]',sf_med,sf_q(1),sf_q(2),sf_lo,sf_hi);
OutVar{end+1,1}="Mean spike rate (spikes/hour)"; OutStat{end+1,1}=sprintf('%.2f (%.2f-%.2f); median CI [%.2f-%.2f]',sr_med,sr_q(1),sr_q(2),sr_lo,sr_hi);
OutVar{end+1,1}="EEGs with reported spikes";  OutStat{end+1,1}="N (% EEGs)";
OutVar{end+1,1}="    Present"; OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_rep_pre,100*n_rep_pre/max(1,n_eegs_all));
OutVar{end+1,1}="    Absent";  OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_rep_abs,100*n_rep_abs/max(1,n_eegs_all));
OutVar{end+1,1}="    Unknown"; OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_rep_unk,100*n_rep_unk/max(1,n_eegs_all));
OutVar{end+1,1}="Patients with reported spikes"; OutStat{end+1,1}="N (% patients)";
OutVar{end+1,1}="    Present"; OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_pats_pre,100*n_pats_pre/max(1,n_pats_rs));
OutVar{end+1,1}="    Absent";  OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_pats_abs,100*n_pats_abs/max(1,n_pats_rs));
OutVar{end+1,1}="    Unknown"; OutStat{end+1,1}=sprintf('%d (%.1f%%)',n_pats_unk,100*n_pats_unk/max(1,n_pats_rs));


Table1_flat = table(string(OutVar), string(OutStat), 'VariableNames',{'Variable','Statistic'});
end

%% =====================================================================
%% TABLE S1
%% =====================================================================

function write_tableS1(MMR, outPath)

function disp = clean_term(raw)
    disp = raw;
    disp = strrep(disp,'(Intercept)',                   'Intercept');
    disp = strrep(disp,'LogSpikesPerHour:AbsLag_years', 'Log spike rate x Absolute lag');
    disp = strrep(disp,'LogSpikesPerHour:LagDirection', 'Log spike rate x Lag direction');
    disp = strrep(disp,'LogSpikesPerHour',              'Log spike rate');
    disp = strrep(disp,'AbsLag_years',                  'Absolute lag (years)');
    disp = strrep(disp,'LagDirection',                  'Lag direction (after vs before)');
    disp = strrep(disp,'EpiType3_cat_Frontal',          'Frontal vs Temporal epilepsy');
    disp = strrep(disp,'EpiType3_cat_General',          'Generalized vs Temporal epilepsy');
end

rows = {};

model_specs = {
    'M1','M1 (logistic, subtypes, interactions)',    'logistic', MMR.FE_M1, MMR.BootstrapTable1;
    'M2','M2 (logistic, subtypes, no interactions)', 'logistic', MMR.FE_M2, [];
    'M3','M3 (linear, subtypes, interactions)',      'linear',   MMR.FE_M3, MMR.BootstrapTable3;
    'M4','M4 (logistic, no subtypes, interactions)', 'logistic', MMR.FE_M4, MMR.BootstrapTable4;
    'M5','M5 (logistic, no subtypes, no interactions)','logistic',MMR.FE_M5, [];
    'M6','M6 (linear, no subtypes, interactions)',   'linear',   MMR.FE_M6, MMR.BootstrapTable6;
};

for mi = 1:size(model_specs,1)
    model_label = model_specs{mi,2};
    is_logistic = strcmp(model_specs{mi,3},'logistic');
    FE  = model_specs{mi,4};
    BT  = model_specs{mi,5};
    if isempty(FE), continue; end

    for i = 1:height(FE)
        term = string(FE.Term(i));
        if term=="(Intercept)", continue; end
        p_val = FE.p(i);

        if is_logistic
            est_pt = FE.OR(i);
            if ~isempty(BT)
                bt_row = BT(string(BT.Term)==term,:);
                if ~isempty(bt_row)
                    ci_lo=bt_row.OR_CI_lo; ci_hi=bt_row.OR_CI_hi; ci_src='Bootstrap';
                else
                    ci_lo=FE.OR_lo(i); ci_hi=FE.OR_hi(i); ci_src='Laplace';
                end
            else
                ci_lo=FE.OR_lo(i); ci_hi=FE.OR_hi(i); ci_src='Laplace';
            end
            est_str=sprintf('%.3f',est_pt); lo_str=sprintf('%.3f',ci_lo); hi_str=sprintf('%.3f',ci_hi);
        else
            est_pt = FE.Beta(i); se = FE.SE(i);
            if ~isempty(BT)
                bt_row = BT(string(BT.Term)==term,:);
                if ~isempty(bt_row)
                    ci_lo=bt_row.Boot_CI_lo; ci_hi=bt_row.Boot_CI_hi; ci_src='Bootstrap';
                else
                    ci_lo=est_pt-1.96*se; ci_hi=est_pt+1.96*se; ci_src='Laplace';
                end
            else
                ci_lo=est_pt-1.96*se; ci_hi=est_pt+1.96*se; ci_src='Laplace';
            end
            est_str=sprintf('%.3f',est_pt); lo_str=sprintf('%.3f',ci_lo); hi_str=sprintf('%.3f',ci_hi);
        end

        if p_val<0.001, p_str='<0.001'; elseif p_val<0.01, p_str=sprintf('%.3f',p_val); else, p_str=sprintf('%.2f',p_val); end
        rows(end+1,:) = {model_label, clean_term(char(term)), est_str, lo_str, hi_str, ci_src, p_str}; %#ok<AGROW>
    end
end

T_out = cell2table(rows, 'VariableNames', ...
    {'Model','Term','Estimate','CI_lower','CI_upper','CI_method','p_value'});
if ~exist(fileparts(outPath),'dir'), mkdir(fileparts(outPath)); end
writetable(T_out, outPath);
end
%% =====================================================================
%% HTML RESULTS
%% =====================================================================

function write_results_html(outPath, Views, SzFreqPerPatient, Fig1Stats, ...
    SpearmanResults_main, rs_all_main, p_all_main, n_all_main, ...
    rho_lo_main, rho_hi_main, subtype_ci_main, ...
    SpearmanResults_S1, rs_all_S1, p_all_S1, n_all_S1, ...
    rho_lo_S1, rho_hi_S1, subtype_ci_S1, ...
    ReportForKeptSessions, MMR)

if ~exist(fileparts(outPath),'dir'), mkdir(fileparts(outPath)); end
fid = fopen(outPath,'w');
if fid==-1, error('Could not open %s', outPath); end

PL         = Views.PatientLevelSpikeRates;
N_total    = numel(PL.Patient);
n_eegs_all = height(ReportForKeptSessions);

Sf     = innerjoin(table(PL.Patient,'VariableNames',{'Patient'}), SzFreqPerPatient,'Keys','Patient');
sf_vec = Sf.MeanSzFreq(isfinite(Sf.MeanSzFreq));
sf_med = median(sf_vec,'omitnan');
[~,sf_ci_lo,sf_ci_hi] = bootstrap_median_ci(sf_vec,5000,0.05);

sr_vec = PL.MeanSpikeRate_perHour(isfinite(PL.MeanSpikeRate_perHour));
sr_med = median(sr_vec,'omitnan');
[~,sr_ci_lo,sr_ci_hi] = bootstrap_median_ci(sr_vec,5000,0.05);

fprintf(fid,'<html><head><meta charset="UTF-8"><title>Results</title></head><body>\n');

%% Cohort
fprintf(fid,'<h2>Cohort summary</h2>\n');
fprintf(fid,'<p>We included %d patients (%d EEGs). Median [95%% CI] monthly seizure frequency was %1.2f [%1.2f-%1.2f], and median spikes/hour was %1.2f [%1.2f-%1.2f] (Table 1).</p>\n', ...
    N_total, n_eegs_all, sf_med, sf_ci_lo, sf_ci_hi, sr_med, sr_ci_lo, sr_ci_hi);

%% Figure 1
fprintf(fid,'<h2>Spike rates by patient groups (Figure 1)</h2>\n');
fprintf(fid,'<p>Spike rates were higher in EEGs with clinically-reported spikes (median %.2f [95%% CI %.2f-%.2f] spikes/hour) than without (%.2f [%.2f-%.2f] spikes/hour) (%s, Cliff''s &delta;=%.2f; Fig. 1A). ', ...
    Fig1Stats.m_pre, Fig1Stats.lo_pre, Fig1Stats.hi_pre, ...
    Fig1Stats.m_abs, Fig1Stats.lo_abs, Fig1Stats.hi_abs, ...
    format_p_html(Fig1Stats.p_rankSum_A), Fig1Stats.effectA_cliff);
fprintf(fid,'Spike rates differed across epilepsy subtypes (Kruskal-Wallis %s, &eta;&sup2;&asymp;%.3f; Fig. 1B).</p>\n', ...
    format_p_html(Fig1Stats.p_kw_C), Fig1Stats.eta2_kw_C);

%% Figure 2
fprintf(fid,'<h2>Spike rate and seizure frequency (Figure 2)</h2>\n');
fprintf(fid,'<p>Spike rate and seizure frequency were positively correlated across all epilepsy patients (N=%d, &rho;=%.2f [95%% CI %.2f-%.2f], %s). ', ...
    n_all_main, rs_all_main, rho_lo_main, rho_hi_main, format_p_html(p_all_main));
fprintf(fid,'Subtype-specific correlations were significant for generalized epilepsy (N=%d, &rho;=%.2f [%.2f-%.2f], Bonferroni-adjusted %s) and temporal lobe epilepsy (N=%d, &rho;=%.2f [%.2f-%.2f], %s), but not frontal lobe epilepsy (N=%d, &rho;=%.2f [%.2f-%.2f], %s; Fig. 2). ', ...
    SpearmanResults_main.N(1), SpearmanResults_main.Spearman_r(1), subtype_ci_main.ci_lo(1), subtype_ci_main.ci_hi(1), format_p_html(SpearmanResults_main.p_bonf(1)), ...
    SpearmanResults_main.N(2), SpearmanResults_main.Spearman_r(2), subtype_ci_main.ci_lo(2), subtype_ci_main.ci_hi(2), format_p_html(SpearmanResults_main.p_bonf(2)), ...
    SpearmanResults_main.N(3), SpearmanResults_main.Spearman_r(3), subtype_ci_main.ci_lo(3), subtype_ci_main.ci_hi(3), format_p_html(SpearmanResults_main.p_bonf(3)));

% S1 Spearman (non-zero only)
canonical_order = ["General","Temporal","Frontal"];
g = canonical_order(2);
row   = SpearmanResults_S1(strcmp(string(SpearmanResults_S1.Group), g), :);
rowCI = subtype_ci_S1(subtype_ci_S1.Group == g, :);
    
if isfinite(rowCI.ci_lo) && isfinite(rowCI.ci_hi)
    ci_str = sprintf('[%.2f-%.2f]', rowCI.ci_lo, rowCI.ci_hi);
else
    ci_str = '';
end
fprintf(fid,['When restricting to patients with non-zero spike rates and seizure frequencies, '...
    'results were similar, although the temporal epilepsy correlation was slightly smaller (&rho;=%.2f %s) '...
    'and no longer significant in this subgroup (Fig. S1). '],row.Spearman_r, ci_str);


fprintf(fid,'Patients with spikes on at least one EEG had higher mean seizure frequencies (Fig. S2).</p>\n');

%% Figure 3
fprintf(fid,'<h2>Mixed effects model (Figure 3)</h2>\n');

if ~isempty(MMR.BootstrapTable1)
    BT = MMR.BootstrapTable1; ci_source = 'bootstrap';
else
    BT = MMR.FE_M1;
    BT.Properties.VariableNames{'OR_lo'} = 'OR_CI_lo';
    BT.Properties.VariableNames{'OR_hi'} = 'OR_CI_hi';
    ci_source = 'Laplace approximation';
end
fprintf(fid,'<p><em>CIs from %s.</em></p>\n', ci_source);

getRow = @(tbl,nm) tbl(string(tbl.Term)==nm,:);
r_spike    = getRow(BT,'LogSpikesPerHour');
r_dir      = getRow(BT,'LagDirection');
r_int_prox = getRow(BT,'LogSpikesPerHour:AbsLag_years');
r_int_dir  = getRow(BT,'LogSpikesPerHour:LagDirection');
r_frontal  = getRow(BT,'EpiType3_cat_Frontal');
r_general  = getRow(BT,'EpiType3_cat_General');
FE1  = MMR.FE_M1;
getP = @(nm) FE1.p(string(FE1.Term)==nm);

n_pairs = height(MMR.ModelTable);
n_pats  = numel(unique(MMR.ModelTable.Patient));

fprintf(fid,['<p>Seizure frequency varies over time within individuals, '...
    'and we hypothesized that spike rates track this variability, ' ...
    'predicting a stronger spike-seizure association for clinic visits close in time to EEGs. '...
    'To test this, we fit logistic mixed effects models on all EEG-visit pairs for patients '...
    'with known epilepsy subtype (N=%d pairs, %d patients), with interaction terms allowing '...
    'the spike-seizure association to vary with the temporal distance between EEG and visit (Supplemental Figure 3). A likelihood ratio '...
    'test confirmed that these interactions '...
    'jointly improved model fit over a model without them (&chi;&sup2;(2), %s). </p>\n'], ...
    n_pairs, n_pats, format_p_html(MMR.LRT_p_A));

fprintf(fid,'<p>Higher spike rates were associated with higher odds of reporting seizures at a clinic visit (OR=%.2f [95%% CI %.2f-%.2f], %s; Fig. 3D). ', ...
    r_spike.OR, r_spike.OR_CI_lo, r_spike.OR_CI_hi, format_p_html(getP('LogSpikesPerHour')));
fprintf(fid,['The spike-seizure association attenuated with greater EEG-visit distance, although the effect was small '...
    '(OR=%.3f [%.3f-%.3f], %s): e.g., the model-predicted OR for spike rate was '], ...
    r_int_prox.OR, r_int_prox.OR_CI_lo, r_int_prox.OR_CI_hi, format_p_html(getP('LogSpikesPerHour:AbsLag_years')));


b_spike_val    = MMR.FE_M1.Beta(string(MMR.FE_M1.Term) == "LogSpikesPerHour");
    b_int_prox_val = MMR.FE_M1.Beta(string(MMR.FE_M1.Term) == "LogSpikesPerHour:AbsLag_years");
    or_at_lag      = exp(b_spike_val + b_int_prox_val .* [0.5, 2, 4]);
fprintf(fid,['%.2f at a 6-month lag, %.2f at 2 years, and %.2f at 4 years ' ...
    '(Fig. 3B). '], ...
    or_at_lag(1), or_at_lag(2), or_at_lag(3));


fprintf(fid,['Spike rates from EEGs obtained before versus after a clinic '...
    'visit were similarly associated with seizure occurrence  (interaction OR=%.3f [%.3f-%.3f], %s). '], ...
    r_int_dir.OR, r_int_dir.OR_CI_lo, r_int_dir.OR_CI_hi, format_p_html(getP('LogSpikesPerHour:LagDirection')));
fprintf(fid,'Visits occurring after the EEG had a lower baseline odds of seizure reporting (OR=%.2f [%.2f-%.2f], %s), consistent with gradual clinical improvement over time (Supplemental Figure 3). ', ...
    r_dir.OR, r_dir.OR_CI_lo, r_dir.OR_CI_hi, format_p_html(getP('LagDirection')));
fprintf(fid,'Compared with temporal lobe epilepsy, generalized epilepsy had a lower baseline odds of seizure reporting (OR=%.2f [%.2f-%.2f], %s), while frontal lobe epilepsy did not differ significantly (OR=%.2f [%.2f-%.2f], %s). ', ...
    r_general.OR, r_general.OR_CI_lo, r_general.OR_CI_hi, format_p_html(getP('EpiType3_cat_General')), ...
    r_frontal.OR, r_frontal.OR_CI_lo, r_frontal.OR_CI_hi, format_p_html(getP('EpiType3_cat_Frontal')));
fprintf(fid,'Results were consistent when using log-transformed continuous seizure frequency as the outcome (Table S1). ');
fprintf(fid,['Together, these results confirm a positive spike rate–seizure association '...
    'and suggest it is strongest when the EEG is obtained close to the clinic visit, '...
    'consistent with spike rates tracking within-individual seizure burden over time.</p>\n']);

fprintf(fid,'</body></html>\n');
fclose(fid);
end

%% =====================================================================
%% REPORT STATUS RESOLUTION
%% =====================================================================

function ReportSlim = resolve_reported_spike_status(ReportForKeptSessions)
any_spikes    = string(ReportForKeptSessions.report_SPORADIC_EPILEPTIFORM_DISCHARGES);
isMainPresent = (any_spikes=="present");
isMainAbsent  = (any_spikes=="absent");

rawF = lower(strtrim(string(ReportForKeptSessions.jay_focal_epi)));
rawM = lower(strtrim(string(ReportForKeptSessions.jay_multifocal_epi)));
rawG = lower(strtrim(string(ReportForKeptSessions.jay_gen_epi)));

isF_P=rawF=="present"; isF_A=rawF=="absent";
isM_P=rawM=="present"; isM_A=rawM=="absent";
isG_P=rawG=="present"; isG_A=rawG=="absent";

presentJay_any = isF_P|isM_P|isG_P;
allJay_absent  = isF_A&isM_A&isG_A;
blankMain      = ~(isMainPresent|isMainAbsent);
blankJayAll    = ~(isF_P|isF_A)&~(isM_P|isM_A)&~(isG_P|isG_A);

if any((allJay_absent&isMainPresent)|(isMainAbsent&isF_P&isM_P&isG_P))
    error('Discordant spike presence between main and jay_* columns.');
end

repCombined = strings(height(ReportForKeptSessions),1);
repCombined(isMainPresent|presentJay_any) = "present";
repCombined(allJay_absent&blankMain)      = "absent";
repCombined(isMainAbsent&blankJayAll)     = "absent";
repCombined(repCombined=="")              = "unknown";

ReportSlim = ReportForKeptSessions(:,{'Patient','Session'});
ReportSlim.ReportStatus = categorical(repCombined,["absent","present","unknown"]);
end

%% =====================================================================
%% SPEARMAN FIGURE
%% =====================================================================

function [SpearmanResults, rs_all, p_all, n_all, rho_lo, rho_hi, subtype_ci] = ...
    spearman_plotting_function(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
    canonical3, spearman_xLims, spearman_yLims, fig_out, freqFieldName, labelSuffix, nonZeroOnly)

nBoot = 5000; alpha = 0.05; fontL = 20;
COL_all   = [0.45 0.45 0.45];
COL_front = [0.93 0.69 0.13];
COL_temp  = [0.85 0.33 0.10];
COL_gen   = [0.00 0.45 0.74];
panelOrder = {'Frontal','Temporal','General'};
panelTitle = {'B. Frontal','C. Temporal','D. General'};

x_all = double(PatientSpikeSz_All.MeanSpikeRate_perHour);
y_all = double(PatientSpikeSz_All.(freqFieldName));
mask_all = isfinite(x_all) & isfinite(y_all);
if nonZeroOnly, mask_all = mask_all & (x_all>0) & (y_all>0); end
n_all = sum(mask_all);

if n_all>=3
    [rs_all,p_all] = corr(x_all(mask_all),y_all(mask_all),'Type','Spearman','Rows','complete');
    [~,rho_lo,rho_hi] = bootstrap_spearman_ci(x_all(mask_all),y_all(mask_all),nBoot,alpha);
else
    rs_all=NaN; p_all=NaN; rho_lo=NaN; rho_hi=NaN;
end

rowsOut = {};
subtype_ci = table(string(canonical3(:)), nan(numel(canonical3),1), nan(numel(canonical3),1), nan(numel(canonical3),1), ...
    'VariableNames',{'Group','rho','ci_lo','ci_hi'});
for ii = 1:numel(canonical3)
    g = canonical3(ii);
    mBase = (PatientSpikeSz_Typed.EpiType3==g);
    x = double(PatientSpikeSz_Typed.MeanSpikeRate_perHour(mBase));
    y = double(PatientSpikeSz_Typed.(freqFieldName)(mBase));
    mask = isfinite(x)&isfinite(y);
    if nonZeroOnly, mask=mask&(x>0)&(y>0); end
    n = sum(mask);
    if n>=3
        [rs,p] = corr(x(mask),y(mask),'Type','Spearman','Rows','complete');
        [rho,lo,hi] = bootstrap_spearman_ci(x(mask),y(mask),nBoot,alpha);
    else
        rs=NaN; p=NaN; rho=NaN; lo=NaN; hi=NaN;
    end
    subtype_ci.rho(ii)=rho; subtype_ci.ci_lo(ii)=lo; subtype_ci.ci_hi(ii)=hi;
    rowsOut(end+1,:) = {char(g),n,rs,p}; %#ok<AGROW>
end
SpearmanResults = cell2table(rowsOut,'VariableNames',{'Group','N','Spearman_r','p_raw'});
SpearmanResults.p_bonf = min(SpearmanResults.p_raw*height(SpearmanResults),1);

x_used = x_all(mask_all); y_used = y_all(mask_all);
minpos_rate = min(x_used(x_used>0)); if isempty(minpos_rate)||~isfinite(minpos_rate), minpos_rate=1e-6; end
minpos_sz   = min(y_used(y_used>0)); if isempty(minpos_sz)||~isfinite(minpos_sz), minpos_sz=1e-6; end
eps_rate = 0.5*minpos_rate; eps_sz = 0.5*minpos_sz;
xZero = log10(eps_sz); yZero = log10(eps_rate);

Tall = table(x_used, y_used, ...
    log10(x_used+(x_used<=0).*eps_rate), log10(y_used+(y_used<=0).*eps_sz), ...
    'VariableNames',{'SpikeRate_perHour','SzFreq','logSpikeRate','logSzFreq'});
isZeroSz  = (Tall.SzFreq==0); isZeroRate = (Tall.SpikeRate_perHour==0);
nonZero   = ~(isZeroSz|isZeroRate);

f2 = figure('Color','w','Position',[60 60 1200 900]);
tiledlayout(f2,2,2,'Padding','compact','TileSpacing','compact');

axA = nexttile(1); hold(axA,'on'); grid(axA,'on'); box(axA,'off');
xline(axA,xZero,':','Color',[0.4 0.4 0.4],'LineWidth',1.2);
yline(axA,yZero,':','Color',[0.4 0.4 0.4],'LineWidth',1.2);
scatter(axA,Tall.logSzFreq(nonZero),Tall.logSpikeRate(nonZero),14,COL_all,'filled','MarkerFaceAlpha',0.25);
plot(axA,Tall.logSzFreq(isZeroSz&~isZeroRate), Tall.logSpikeRate(isZeroSz&~isZeroRate),'*','Color',COL_all,'MarkerSize',7,'LineWidth',1);
plot(axA,Tall.logSzFreq(~isZeroSz&isZeroRate), Tall.logSpikeRate(~isZeroSz&isZeroRate),'*','Color',COL_all,'MarkerSize',8,'LineWidth',1);
plot(axA,Tall.logSzFreq(isZeroSz&isZeroRate),  Tall.logSpikeRate(isZeroSz&isZeroRate), '*','Color',COL_all,'MarkerSize',8,'LineWidth',1.2);
if n_all>=3
    X=[ones(n_all,1),Tall.logSzFreq]; b=X\Tall.logSpikeRate;
    xgrid=linspace(spearman_xLims(1),spearman_xLims(2),300)';
    plot(axA,xgrid,b(1)+b(2)*xgrid,'k-','LineWidth',2);
end
xlim(axA,spearman_xLims); ylim(axA,spearman_yLims);
xlabel(axA,'Seizures per month (log scale)','FontSize',fontL);
ylabel(axA,'Spikes per hour (log scale)','FontSize',fontL);
set_log10_ticks(axA,'x',eps_sz,spearman_xLims);
set_log10_ticks(axA,'y',eps_rate,spearman_yLims);
labs=string(axA.XTickLabel); [~,iMax]=max(axA.XTick); labs(iMax)=""; axA.XTickLabel=labs;
title(axA,sprintf('A. All epilepsy%s (N=%d)',labelSuffix,n_all),'FontSize',fontL,'FontWeight','bold');
text(axA,0.98,0.95,sprintf('\\rho=%.2f [%.2f-%.2f], %s',rs_all,rho_lo,rho_hi,p_label(p_all)),...
    'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',fontL-2,'FontWeight','bold');
set(axA,'FontSize',fontL);

colMap = struct('Frontal',COL_front,'Temporal',COL_temp,'General',COL_gen);
for p = 1:3
    ax = nexttile(p+1); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    gStr = string(panelOrder{p});
    col  = colMap.(gStr);
    idx  = (string(PatientSpikeSz_Typed.EpiType3)==gStr) & ...
           isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perHour) & ...
           isfinite(PatientSpikeSz_Typed.(freqFieldName));
    if nonZeroOnly
        idx = idx & (PatientSpikeSz_Typed.MeanSpikeRate_perHour>0) & ...
                    (PatientSpikeSz_Typed.(freqFieldName)>0);
    end
    if nnz(idx)==0, axis(ax,'off'); continue; end

    x_raw=double(PatientSpikeSz_Typed.(freqFieldName)(idx));
    y_raw=double(PatientSpikeSz_Typed.MeanSpikeRate_perHour(idx));
    logX=log10(x_raw+(x_raw<=0).*eps_sz);
    logY=log10(y_raw+(y_raw<=0).*eps_rate);
    isZx=(x_raw==0); isZy=(y_raw==0);

    xline(ax,xZero,':','Color',[0.4 0.4 0.4],'LineWidth',1.2);
    yline(ax,yZero,':','Color',[0.4 0.4 0.4],'LineWidth',1.2);
    scatter(ax,logX(~isZx&~isZy),logY(~isZx&~isZy),18,col,'filled','MarkerFaceAlpha',0.35);
    if any(isZx&~isZy), plot(ax,logX(isZx&~isZy),logY(isZx&~isZy),'*','Color',col,'MarkerSize',8,'LineWidth',1.1); end
    if any(~isZx&isZy), plot(ax,logX(~isZx&isZy),logY(~isZx&isZy),'*','Color',col,'MarkerSize',8,'LineWidth',1.1); end
    if any(isZx&isZy),  plot(ax,logX(isZx&isZy), logY(isZx&isZy), '*','Color',col,'MarkerSize',9,'LineWidth',1.2); end
    if nnz(~isZx&~isZy)>=3
        Xg=[ones(nnz(~isZx&~isZy),1),logX(~isZx&~isZy)]; bg=Xg\logY(~isZx&~isZy);
        xg=linspace(spearman_xLims(1),spearman_xLims(2),250)';
        plot(ax,xg,bg(1)+bg(2)*xg,'-','Color',col,'LineWidth',2);
    end
    xlim(ax,spearman_xLims); ylim(ax,spearman_yLims);
    xlabel(ax,'Seizures per month (log scale)','FontSize',fontL);
    ylabel(ax,'Spikes per hour (log scale)','FontSize',fontL);
    set_log10_ticks(ax,'x',eps_sz,spearman_xLims);
    set_log10_ticks(ax,'y',eps_rate,spearman_yLims);
    labs=string(ax.XTickLabel); [~,iMax]=max(ax.XTick); labs(iMax)=""; ax.XTickLabel=labs;

    row   = SpearmanResults(strcmp(string(SpearmanResults.Group),gStr),:);
    rowCI = subtype_ci(subtype_ci.Group==gStr,:);
    txt = sprintf('\\rho=%.2f [%.2f-%.2f], p_{bonf}%s', row.Spearman_r, rowCI.ci_lo, rowCI.ci_hi, ...
        string(regexprep(char(p_label(row.p_bonf)),'^p','')));
    title(ax,sprintf('%s%s (N=%d)',panelTitle{p},labelSuffix,nnz(idx)),'FontSize',fontL,'FontWeight','bold');
    text(ax,0.98,0.95,txt,'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',fontL-3,'FontWeight','bold');
    set(ax,'FontSize',fontL);
end

if ~exist(fileparts(fig_out),'dir'), mkdir(fileparts(fig_out)); end
exportgraphics(f2, fig_out,'Resolution',300);
fprintf('Saved Spearman figure: %s\n', fig_out);
end

%% =====================================================================
%% SMALL UTILITIES
%% =====================================================================

function arr = json_to_string_array(s)
s = strtrim(string(s));
if s==""||s=="[]"||s=="<missing>", arr=strings(0,1); return; end
dec = jsondecode(char(s));
if iscell(dec)
    arr = strings(numel(dec),1);
    for k=1:numel(dec)
        x=dec{k};
        if ischar(x)||(isstring(x)&&isscalar(x)), arr(k)=string(x); end
    end
elseif ischar(dec)||(isstring(dec)&&isscalar(dec))
    arr = string(dec);
elseif isnumeric(dec)
    arr = string(dec(:));
else
    error('Unsupported JSON string-array type.');
end
arr = string(arr(:));
end

function arr = json_to_double_array(s)
s = strtrim(string(s));
if s==""||s=="[]"||s=="<missing>", arr=double([]); return; end
s   = regexprep(s,'null','NaN','ignorecase');
dec = jsondecode(char(s));
arr = double(dec(:));
end

function out = max_hasSz(x)
x=x(isfinite(x)); if isempty(x), out=NaN; else, out=max(x); end
end

function out = local_frac_hasSz1(hsVec)
hs=hsVec(isfinite(hsVec));
if isempty(hs), out=NaN; return; end
valid=(hs==0|hs==1);
if ~any(valid), out=NaN; else, out=nnz(hs(valid)==1)/nnz(valid); end
end

function out = local_min_omitnan(a)
a=a(isfinite(a)); if isempty(a), out=NaN; else, out=min(a); end
end

function out = local_first_nonmissing(s)
s=s(:); s=s(~ismissing(s)&strlength(s)>0);
if isempty(s), out=""; else, out=s(1); end
end

function ylog = to_log10_per_hour(x, eps_rate)
x=double(x); x(~isfinite(x)|x<=0)=eps_rate; ylog=log10(x);
end

function ylog = to_log10_per_month(f, eps_freq)
f=double(f); f(~isfinite(f)|f<=0)=eps_freq; ylog=log10(f);
end

function Yj = add_y_jitter_eps(Y, Y_ZERO, Y_LIMS, frac)
Yj=Y; mask=abs(Y-Y_ZERO)<1e-9;
if any(mask), amp=frac*diff(Y_LIMS); Yj(mask)=Yj(mask)+(rand(sum(mask),1)-0.5)*amp; end
end

function add_sigbar(ax, x1, x2, y, ptext)
tick=0.03*diff(ax.YLim);
plot(ax,[x1 x1 x2 x2],[y-tick,y,y,y-tick],'k-','LineWidth',1.3);
if ptext=="**"||ptext=="***", yOff=-0.012*diff(ax.YLim); else, yOff=0.003*diff(ax.YLim); end
text(ax,mean([x1 x2]),y+yOff,ptext,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);
end

function pStr = p_label(p)
if isnan(p), pStr="p=NaN"; return; end
if p<0.001, pStr="p<0.001"; return; end
if p<0.01,  pStr=sprintf("p=%.2g",p); return; end
pStr=sprintf("p=%.2f",p);
end

function pStr = format_p_html(p)
if isnan(p), pStr='p = NaN'; return; end
if p<0.001, pStr='p &lt; 0.001'; return; end
if p<0.01
    s=sprintf('%.2g',p); if startsWith(s,'.'), s=['0' s]; end
    pStr=['p = ' s]; return;
end
pStr=sprintf('p = %.2f',p);
end

function set_log10_ticks(ax, whichAxis, eps_val, axisLims, maxPow)
if nargin<5||isempty(maxPow), maxPow=6; end
whichAxis=lower(whichAxis);
decades=10.^(0:maxPow); log_dec=log10(decades);
keep=(log_dec>=axisLims(1))&(log_dec<=axisLims(2));
ticks=log_dec(keep); labels=string(decades(keep));
log_eps=log10(eps_val);
if log_eps>=axisLims(1)&&log_eps<=axisLims(2)
    ticks=[log_eps;ticks(:)]; labels=["0";labels(:)];
end
if whichAxis=="x", ax.XTick=ticks; ax.XTickLabel=labels;
else,              ax.YTick=ticks; ax.YTickLabel=labels; end
end

function [med,lo,hi] = bootstrap_median_ci(x, nBoot, alpha)
x=double(x(:)); x=x(isfinite(x));
if isempty(x)||nBoot==0, med=median(x,'omitnan'); lo=NaN; hi=NaN; return; end
med=median(x,'omitnan'); n=numel(x);
bootM=nan(nBoot,1);
for b=1:nBoot, idx=randi(n,n,1); bootM(b)=median(x(idx)); end
lo=prctile(bootM,100*(alpha/2)); hi=prctile(bootM,100*(1-alpha/2));
end

function add_median_ci_overlay(ax, xpos, med_raw, lo_raw, hi_raw, eps_floor)
yMed=to_log10_per_hour(med_raw,eps_floor);
yLo =to_log10_per_hour(lo_raw, eps_floor);
yHi =to_log10_per_hour(hi_raw, eps_floor);
plot(ax,[xpos xpos],[yLo yHi],'k-','LineWidth',3);
plot(ax,xpos,yMed,'ko','MarkerFaceColor','k','MarkerSize',6);
end

function add_median_ci_overlay_month(ax, xpos, med_raw, lo_raw, hi_raw, eps_floor)
yMed=to_log10_per_month(med_raw,eps_floor);
yLo =to_log10_per_month(lo_raw, eps_floor);
yHi =to_log10_per_month(hi_raw, eps_floor);
plot(ax,[xpos xpos],[yLo yHi],'k-','LineWidth',3);
plot(ax,xpos,yMed,'ko','MarkerFaceColor','k','MarkerSize',6);
end

function [rho_hat,lo,hi] = bootstrap_spearman_ci(x, y, nBoot, alpha)
x=double(x(:)); y=double(y(:));
mask=isfinite(x)&isfinite(y); x=x(mask); y=y(mask); n=numel(x);
if n<3, rho_hat=NaN; lo=NaN; hi=NaN; return; end
rho_hat=corr(x,y,'Type','Spearman','Rows','complete');
rho_boot=nan(nBoot,1);
for b=1:nBoot, idx=randi(n,n,1); rho_boot(b)=corr(x(idx),y(idx),'Type','Spearman','Rows','complete'); end
lo=prctile(rho_boot,100*(alpha/2)); hi=prctile(rho_boot,100*(1-alpha/2));
end

function d = cliff_delta(x1, x2)
x1=x1(isfinite(x1(:))); x2=x2(isfinite(x2(:)));
n1=numel(x1); n2=numel(x2);
if n1==0||n2==0, d=NaN; return; end
[~,~,stats]=ranksum(x1,x2,'method','approx');
R1=stats.ranksum; U1=R1-n1*(n1+1)/2;
d=(2*U1/(n1*n2))-1;
end


function [T_boot, boot_betas] = run_bootstrap_lme(T, mdl, formula, nBoot, alpha, label)
T_boot     = [];
boot_betas = [];
if isempty(mdl) || nBoot == 0, return; end

fprintf('\nBootstrapping %s (%d iterations)...\n', label, nBoot);
patients   = unique(T.PatientID);
nPat       = numel(patients);
nFixed     = size(fixedEffects(mdl), 1);
boot_betas = nan(nBoot, nFixed);

parfor b = 1:nBoot
    idx      = randi(nPat, nPat, 1);
    bootPats = patients(idx);
    Tboot    = cell(nPat, 1);
    for k = 1:nPat
        Tboot{k} = T(T.PatientID == bootPats(k), :);
        Tboot{k}.PatientID = categorical(repmat(k, height(Tboot{k}), 1));
    end
    Tboot = vertcat(Tboot{:});
    try
        mboot = fitlme(Tboot, formula, 'FitMethod','REML');
        boot_betas(b,:) = fixedEffects(mboot)';
    catch
    end
end

converged  = all(isfinite(boot_betas), 2);
boot_betas = boot_betas(converged, :);
fprintf('%s bootstrap: %d/%d converged (%.1f%%)\n', ...
    label, sum(converged), nBoot, 100*mean(converged));

ci_lo = prctile(boot_betas, 100*(alpha/2),   1);
ci_hi = prctile(boot_betas, 100*(1-alpha/2), 1);
[beta_obs, betanames_obs] = fixedEffects(mdl);

T_boot = table(string(betanames_obs.Name), beta_obs, ci_lo(:), ci_hi(:), ...
    'VariableNames',{'Term','Beta','Boot_CI_lo','Boot_CI_hi'});
fprintf('%s bootstrapped betas:\n', label); disp(T_boot);
end


function PairTable = build_eeg_visit_pairs_nosubtype(Vuniq, SessionLevelSpikeRates, ...
    ReportForKeptSessions, PatientTyping_AllEpilepsy)
%% All-epilepsy patients — no EpiType3 restriction

EEG_raw = ReportForKeptSessions.start_time_deid;
if isdatetime(EEG_raw), EEG_dt = EEG_raw;
else, EEG_dt = datetime(strtrim(string(EEG_raw)), 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss");
end

EEG_dates = table(double(ReportForKeptSessions.Patient), ...
    double(ReportForKeptSessions.Session), EEG_dt, ...
    'VariableNames',{'Patient','Session','EEG_Date'});
EEG_dates = EEG_dates(~isnat(EEG_dates.EEG_Date), :);

EEG_tbl = innerjoin(EEG_dates, ...
    SessionLevelSpikeRates(:,{'Patient','Session','SpikesPerHour'}), ...
    'Keys',{'Patient','Session'});

% Use all epilepsy patients — join only to get Patient list, no EpiType3
AllEpiPatients = unique(PatientTyping_AllEpilepsy.Patient);
Vep = Vuniq(ismember(Vuniq.Patient, AllEpiPatients), {'Patient','VisitDate','Freq_R1','HasSz'});

patients  = intersect(unique(EEG_tbl.Patient), unique(Vep.Patient));
nEstimate = 0;
for i = 1:numel(patients)
    p = patients(i);
    nEstimate = nEstimate + sum(EEG_tbl.Patient==p) * sum(Vep.Patient==p);
end

Patient_out       = nan(nEstimate,1);
Session_out       = nan(nEstimate,1);
VisitDate_out     = NaT(nEstimate,1);
SpikesPerHour_out = nan(nEstimate,1);
SzFreq_out        = nan(nEstimate,1);
HasSz_out         = nan(nEstimate,1);
SignedLag_out     = nan(nEstimate,1);

row = 0;
for i = 1:numel(patients)
    p        = patients(i);
    eeg_rows = EEG_tbl(EEG_tbl.Patient==p, :);
    vis_rows = Vep(Vep.Patient==p, :);
    for e = 1:height(eeg_rows)
        for v = 1:height(vis_rows)
            row = row + 1;
            Patient_out(row)       = p;
            Session_out(row)       = eeg_rows.Session(e);
            VisitDate_out(row)     = vis_rows.VisitDate(v);
            SpikesPerHour_out(row) = eeg_rows.SpikesPerHour(e);
            SzFreq_out(row)        = vis_rows.Freq_R1(v);
            HasSz_out(row)         = vis_rows.HasSz(v);
            SignedLag_out(row)     = days(vis_rows.VisitDate(v) - eeg_rows.EEG_Date(e));
        end
    end
end

PairTable = table(Patient_out(1:row), Session_out(1:row), VisitDate_out(1:row), ...
    SpikesPerHour_out(1:row), SzFreq_out(1:row), HasSz_out(1:row), ...
    SignedLag_out(1:row), ...
    'VariableNames',{'Patient','Session','VisitDate','SpikesPerHour','SzFreq', ...
                     'HasSz','SignedLag_days'});

PairTable.EEG_ID = categorical(string(PairTable.Patient) + "_" + string(PairTable.Session));

keepMask  = isfinite(PairTable.SpikesPerHour) & isfinite(PairTable.SzFreq) & ...
            isfinite(PairTable.SignedLag_days);
n_before  = height(PairTable);
PairTable = PairTable(keepMask,:);
fprintf('[build_eeg_visit_pairs_nosubtype] %d patients, %d pairs (%d removed)\n', ...
    numel(patients), height(PairTable), n_before-height(PairTable));

EPS_SPIKE = 1e-3;
PairTable.LogSpikesPerHour = log(PairTable.SpikesPerHour + EPS_SPIKE);
PairTable.SignedLag_years  = PairTable.SignedLag_days / 365.25;
PairTable.PatientID        = categorical(PairTable.Patient);
end