function run_spike_sz_pipeline_clean()
%% ===================== INFORMATION ==================
%{
This contains the code for running all analyses associated with the
manuscript "Interictal Spike Rate on Routine Outpatient EEG Is Associated 
With Seizure Frequency in a Large Epilepsy Cohort" by Conrad et al., 2026

Requirements:
1. spike_counts.csv and clinical_data_deidentified.csv
     - both are available for download at: https://upenn.box.com/s/yy4o1t6nit7yu35flz59ux6zf54slg9m
     - spike_counts.csv contains a list of spike counts at varying SpikeNet
     probability threshold for each EEG
     - clinical_data_deidentified.csv contains clinical information for
     each patient. De-identified birth dates and dates of service are
     provided by date-shifting each date relative to the date of the
     patient's first clinic visit, defined to be Jan 1 2000. E.g., if a
     patient's date of birth is 1/1/1980 and their first visit is 1/1/2010,
     then their deidentified birth date is 1/1/1970. Each row is one EEG,
     and so the same patient may appear in multiple rows.
2. Matlab (vR2024a) and the Statistics and Machine Learning Toolbox
3. This codebase, located at: https://github.com/erinconrad/seizure_severity
4. Create an output directory and a data directory in the paths noted
below, and place both csv datasets in the data directory.

To run the analysis, navigate to the directory containing this script and
run
>> run_spike_sz_pipeline_clean

This function will reproduce all figures, the table, and the results text
from the manuscript and save them in the output directory.
%}


%% ======================= RNG =======================
rng(1); % for bootstrapping

%% ======================= PATHS =======================
% Inputs
spikeSummaryMultiCsv = '../data/spike_counts.csv'; % PUT THIS IN YOUR PATH
reportCsv            = '../data/clinical_data_deidentified.csv'; % PUT THIS IN YOUR PATH
%reportCsv = '../data/Routineeegpec-Deidreport_DATA_LABELS_2026-03-11_0854.csv';
%reportCsv = '../data/Routineeegpec-Deidreport_DATA_LABELS_2026-03-11_0911.csv';

% Outputs (PUT THE OUTPUT FOLDER IN YOUR PATH)
fig1_out   = '../output/Fig1.png';
fig2_out   = '../output/Fig2.png';
figS1_out  = '../output/FigS1.png';
figS2_out  = '../output/FigS2.png';
figMain_out  = '../output/FigModel.png';
figCV_out     = '../output/FigCV.png';

Table1Csv  = fullfile('..','output','Table1.csv');
resultsHtml = '../output/results_summary.html';

%% ======================= PARAMETERS ===================
MAX_ROUTINE_HOURS = 4; % only allow EEGs less than 4 hours (exclude ambulatories)

NESD_LABEL = "Non-Epileptic Seizure Disorder"; % excluded
badTypes   = lower(["Uncertain if Epilepsy","Unknown or MRN not found",""]); % excluded from epilepsy
canonical3 = ["General","Temporal","Frontal"];

EPS_RATE = 30e-3;               % spikes/hour floor for log plotting
Y_ZERO  = log10(EPS_RATE);
Y_LIMS  = [-2 4];
TITLE_Y_OFFSET = 0.02;

spearman_xLims = [-3.5, 4];     % log10(seizures/month)
spearman_yLims = [-1.5, 3];     % log10(spikes/hour)

nBoot = 5000; % bootstrap iterations to get confidence intervals
alpha = 0.05;

% SpikeNet threshold columns
countCol = "count_0_46"; % the probability threshold from the SpikeNet2 paper to convert probabilities to detections
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

require_cols(SpikeSummaryTable, ["Patient","Session",countCol,durCol], "SpikeSummaryTable"); % a function to check that the required columns exist

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

% Sanity check: keys match 1:1 (both the spike summary and report table should
% only have one row for each unique patient-session combo)
assert_unique_keys(SpikeSummaryTable, "Patient","Session","SpikeSummaryTable");
assert_unique_keys(ReportTable, "patient_id","session_number","ReportTable");

%% ======================= BUILD VISIT-LEVEL + PATIENT-LEVEL METRICS =======================
% Visit-level information: Patient, VisitDate, Freq, HasSz, Freq_R1 
Vuniq = build_visit_level_table_R1(ReportTable); % 1 row per visit, with information about seizure frequency

% Patient typing from report 
PatientTypingAll = build_patient_typing_from_report(ReportTable, canonical3); % 1 row per patient, with information about epilepsy type

% Patient-level seizure metrics: MeanSzFreq + fraction hasSz==1 among valid visits
SzFreqPerPatient = build_patient_seizure_metrics(Vuniq); % 1 row per patient, with mean seizure frequency across visits

% Patients differ between two above tables because PatientTypingAll ok with
% patients without SzFreq and SzFreqPerPatient ok with non-epilepsy
% patients

%% ======================= BUILD VIEW BUNDLE =======================
Views = build_filtered_view( ...
    SpikeSummaryTable, ReportTable, PatientTypingAll, SzFreqPerPatient, ...
    NESD_LABEL, badTypes, canonical3); % a structure that contains the patient level data for subsequent analyses
% This filters to epilepsy patients with documented sz frequency and at
% least one outpatient EEG

%% ======================= DO MODEL =======================

PairTable = build_eeg_visit_pairs(Vuniq, Views.SessionLevelSpikeRates, ...
                                            Views.ReportForKeptSessions, Views.PatientTypingFiltered);
fprintf('Total pairs: %d\n', height(PairTable));
fprintf('Unique patients: %d\n', numel(unique(PairTable.Patient)));
fprintf('Unique EEGs: %d\n', numel(unique(PairTable.EEG_ID)));
fprintf('Unique visits: %d\n', numel(unique(string(PairTable.Patient) + "_" + string(PairTable.VisitDate))));

subtypes = unique(PairTable.EpiType3);
for i = 1:numel(subtypes)
    n_pts = numel(unique(PairTable.Patient(PairTable.EpiType3 == subtypes(i))));
    n_pairs = sum(PairTable.EpiType3 == subtypes(i));
    fprintf('%s: %d patients, %d pairs\n', subtypes(i), n_pts, n_pairs);
end


% Do model
MMR = fit_mixed_effects_models(PairTable, nBoot);
%MMR = fit_mixed_effects_models(PairTable, 0);


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
        'MeanSzFreq', ' (positive spike/seizures only)', true);

fprintf('Saved Fig 2:  %s\n', fig2_out);
fprintf('Saved Fig S1: %s\n', figS1_out);

%% ======================= FIGURE 3 (MODEL) =======================

[FigMain, FigCV] = make_figures_and_cv(...
    MMR, ...
    Views, ...
    Vuniq, ...
    Views.ReportForKeptSessions, ...
    ["Temporal","Frontal","General"], ...
    5, ...
    figMain_out, ...
    figCV_out);

%% ======================= FIG S2 (SZ FREQ BY REPORTED SPIKES ACROSS EEGs) =======================
FigS2 = make_figS2_sz_by_reported_spikes(Views, SzFreqPerPatient, nBoot, alpha, spearman_xLims);
save_fig(FigS2, figS2_out);
fprintf('Saved Fig S2: %s\n', figS2_out);


%% ======================= TABLE 1 =======================
Table1_flat = build_table1_flat(Views, SzFreqPerPatient, Vuniq, EPS_RATE, nBoot, alpha);
if ~exist(fileparts(Table1Csv),'dir'), mkdir(fileparts(Table1Csv)); end
writetable(Table1_flat, Table1Csv);
fprintf('Wrote Table 1 CSV: %s\n', Table1Csv);

%% ======================= TABLE S1 =======================
tableS1Csv = fullfile('..','output','TableS1.csv');
write_tableS1(MMR, tableS1Csv);
fprintf('Wrote Table S1 CSV: %s\n', tableS1Csv);

%% ======================= WRITE HTML SUMMARY =======================
write_results_html(resultsHtml, Views, SzFreqPerPatient, ...
    Fig1Stats, ...
    SpearmanResults_main, rs_all_main, p_all_main, n_all_main, rho_lo_main, rho_hi_main, subtype_ci_main, ...
    Views.ReportForKeptSessions, MMR);

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

% loop over rows in the clinical data report
for i = 1:height(R)
    vt_raw    = strtrim(string(R.visit_type(i)));
    dates_raw = strtrim(string(R.visit_dates_deid(i)));
    sz_raw    = strtrim(string(R.sz_freqs(i)));
    hs_raw    = strtrim(string(R.visit_hasSz(i)));

    % Make missing data formatted same way
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

    % Check if the visit type is an allowable (outpatient) visit
    keepMask = ismember(string(vt), allowable_visits);

    % remove the patient's info if no allowable visits
    if ~any(keepMask)
        R.visit_type(i)        = "[]";
        R.visit_dates_deid(i)  = "[]";
        R.sz_freqs(i)          = "[]";
        R.visit_hasSz(i)       = "[]";
        continue
    end

    % Keep only the allowable visits
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

%% Require it be an outpatient study

% outpatient if acquired on an SPE or radnor machine
acqStr = lower(strtrim(string(R.acquired_on)));
isOutpt_site  = contains(acqStr,"spe") | contains(acqStr,"radnor");

% outpatient if epic report says it's outpatient
classStr = lower(strtrim(string(R.report_PATIENT_CLASS)));
isOutpt_class = (classStr == "outpatient");

% outpatient if Jay report say it's outpatient
jayStr = lower(strtrim(string(R.jay_in_or_out)));
isOutpt_jay = (jayStr == "out");

OutptKeys = unique(R(isOutpt_site | isOutpt_class | isOutpt_jay, {'patient_id','session_number'})); % If any of these outpatient conditions are true, it's outpatient
OutptKeys.Properties.VariableNames = {'Patient','Session'};
if isempty(OutptKeys)
    error('No outpatient sessions identified by site/class/jay flags.');
end

%% Require it be a routine EEG
% require it's not too long (remove ambulatory EEGs)
isRoutine = isfinite(S.(durCol)) & S.(durCol) <= MAX_ROUTINE_HOURS*3600;
RoutineKeys = unique(S(isRoutine, {'Patient','Session'}));

%% Combine both rules
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
% Builds PV by exploding JSON arrays for each row, then collapses to unique
% (Patient,VisitDate), and so there is one row for each visit
PV = table('Size',[0 4], ...
    'VariableTypes',{'double','datetime','double','double'}, ...
    'VariableNames',{'Patient','VisitDate','Freq','HasSz'});

for j = 1:height(R)
    pid = double(R.patient_id(j));

    % visit date
    ds = strtrim(string(R.visit_dates_deid(j)));
    if ds=="[]" || ds=="", continue; end
    dd = string(jsondecode(char(ds)));
    d  = datetime(dd,'InputFormat','yyyy-MM-dd');

    % visit sz freq
    s = strtrim(string(R.sz_freqs(j)));
    s = regexprep(s,'null','NaN','ignorecase');
    v = double(jsondecode(char(s)));
    v(~isfinite(v)) = NaN;
    v(v<0) = NaN;

    % visit has sz
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

% Collapse duplicates (there are duplicates because each row in the table
% is an EEG, and some patients have multiple EEGs)
[gv, pid_keys, date_keys] = findgroups(PV.Patient, PV.VisitDate);

Freq_agg = splitapply(@(x) mean(x(isfinite(x)),'omitnan'), PV.Freq, gv);
Has_agg  = splitapply(@max_hasSz, PV.HasSz, gv);

Vuniq = table(pid_keys, date_keys, Freq_agg, Has_agg, ...
    'VariableNames', {'Patient','VisitDate','Freq','HasSz'});

% Rule to fill in missing seizure frequencies: if HasSz==0 and Freq missing -> set Freq_R1=0
Vuniq.Freq_R1 = Vuniq.Freq;
mask_rule1 = ~isfinite(Vuniq.Freq_R1) & (Vuniq.HasSz==0);
Vuniq.Freq_R1(mask_rule1) = 0;
end

function SzP = build_patient_seizure_metrics(Vuniq)
% Gets average seizure frequency across visits for each patient

[gp, pids] = findgroups(Vuniq.Patient); % find all visits for a patient
MeanSzFreq = splitapply(@(x) mean(x,'omitnan'), Vuniq.Freq_R1, gp); % Take mean seizure frequency across visits
FracVisits_HasSz1 = splitapply(@local_frac_hasSz1, Vuniq.HasSz, gp);
SzP = table(pids, MeanSzFreq, FracVisits_HasSz1, ...
    'VariableNames', {'Patient','MeanSzFreq','FracVisits_HasSz1'});
end

function T = build_patient_typing_from_report(R, canonical3)
% Get the epilepsy type


pid = double(R.patient_id);

etype = strtrim(string(R.epilepsy_type));
espec = strtrim(string(R.epilepsy_specific));

% Keep non-missing rows
okT = ~ismissing(etype) & strlength(etype)>0;
okS = ~ismissing(espec) & strlength(espec)>0;

% Collapse epilepsy_type (so that each patient gets one row)
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

% Collapse epilepsy_specific (same as above)
Tspec = table(pid(okS), espec(okS), 'VariableNames', {'Patient','EpilepsySpecific_raw'});
Tspec = sortrows(Tspec, 'Patient');
[uidS, ~, gS] = unique(Tspec.Patient, 'stable');
espec_one = strings(numel(uidS),1);
for k=1:numel(uidS)
    vals = unique(Tspec.EpilepsySpecific_raw(gS==k));
    vals = vals(strlength(vals)>0);
    if isempty(vals), espec_one(k)=""; continue; end
    if numel(vals) > 1
        error('Conflicting epilepsy_specific for Patient %g: %s', uidS(k), strjoin(vals,", "));
    end
    espec_one(k) = vals(1);
end

% Combine per-patient epilepsy_type and epilepsy_specific into a single
% patient table
assert(all(uid==uidS)) % I think the patient IDs should all match
T = table(uid, etype_one, 'VariableNames', {'Patient','EpilepsyType'});
T = outerjoin(T, table(uidS, espec_one, 'VariableNames', {'Patient','EpilepsySpecific'}), ...
    'Keys','Patient','MergeKeys',true);

% Map to canonical3
spec_norm = lower(strtrim(string(T.EpilepsySpecific)));
type_norm = lower(strtrim(string(T.EpilepsyType)));

E3 = strings(height(T),1);
E3(contains(spec_norm,"temporal")) = "Temporal";
E3(contains(spec_norm,"frontal"))  = "Frontal";
%E3((strlength(E3)==0) & contains(type_norm,"general")) = "General"; % this
%mapped combined generalized and focal to general
E3((strlength(E3)==0) & strcmp(type_norm,"general")) = "General"; % epilepsy type needs to be general (not just contain it) for final determination to be general

T.EpiType3 = categorical(E3, canonical3);
end

%% =====================================================================
%% ======================= VIEW BUNDLE (CLEAN) =========================
%% =====================================================================
function Views = build_filtered_view(SessionsFiltered, ReportIn, PatientTypingAll, SzFreqPerPatient, ...
                                     NESD_LABEL, badTypes, canonical3)
% build_filtered_view
% This builds the final analysis cohort
%
% Inputs:
%   SessionsFiltered: table with Patient, Session, SpikeRate_perHour (after outpatient+routine filtering)
%   ReportIn:         clinical report table (one row per EEG)
%   PatientTypingAll: output of build_patient_typing_from_report (one row per patient)
%   SzFreqPerPatient: output of build_patient_seizure_metrics (one row per patient, not same as PatientTyping)
%   NESD_LABEL:       label for NESD
%   badTypes:         excluded epilepsy_type values (lowercase strings)
%   canonical3:       string array like ["General","Temporal","Frontal"]

%% ------------------ Join report to kept sessions ------------------
SessKeys = unique(SessionsFiltered(:,{'Patient','Session'}));

ReportForKeptSessions = innerjoin(ReportIn, SessKeys, ...
    'LeftKeys', {'patient_id','session_number'}, ...
    'RightKeys', {'Patient','Session'});

ReportForKeptSessions.Patient = double(ReportForKeptSessions.patient_id);
ReportForKeptSessions.Session = double(ReportForKeptSessions.session_number);

%% ------------------ Restrict typing to patients with ≥1 kept EEG ------------------
PatientsKept  = unique(double(SessionsFiltered.Patient));
TypingFiltered = innerjoin(PatientTypingAll, table(PatientsKept,'VariableNames',{'Patient'}), ...
    'Keys','Patient');

%% ------------------ Patient-level mean spike rate ------------------
[gpS, pidsS] = findgroups(double(SessionsFiltered.Patient));
MeanSpikeRate_perHour = splitapply(@(x) mean(x,'omitnan'), SessionsFiltered.SpikeRate_perHour, gpS);

PatientLevelSpikeRates = table(double(pidsS), MeanSpikeRate_perHour, ...
    'VariableNames', {'Patient','MeanSpikeRate_perHour'});

% Bring in EpilepsyType/EpilepsySpecific/EpiType3 ONCE
needTypingCols = {'Patient','EpilepsyType','EpilepsySpecific','EpiType3'};
missingT = setdiff(string(needTypingCols), string(TypingFiltered.Properties.VariableNames));
if ~isempty(missingT)
    error("TypingFiltered missing required cols: %s", strjoin(missingT,", "));
end

PatientLevelSpikeRates = innerjoin(PatientLevelSpikeRates, ...
    TypingFiltered(:, needTypingCols), 'Keys','Patient');

%% ------------------ Epilepsy vs NESD masks (from typing) ------------------
etype_norm = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));

IsNESDMask = (etype_norm == lower(strtrim(NESD_LABEL)));
IsBadType  = ismember(etype_norm, badTypes) | ismissing(etype_norm) | (strlength(etype_norm)==0);
IsEpilepsyMask = ~IsNESDMask & ~IsBadType;

%% ------------------ Session-level spike rates ------------------
SessionLevelSpikeRates = SessionsFiltered(:, {'Patient','Session','SpikeRate_perHour'});
SessionLevelSpikeRates.Properties.VariableNames{'SpikeRate_perHour'} = 'SpikesPerHour';

%% ------------------ Seizure-frequency tables aligned to kept patients ------------------
SzFreqFiltered = innerjoin(SzFreqPerPatient, table(PatientsKept,'VariableNames',{'Patient'}), 'Keys','Patient');

EpPatients = table(PatientLevelSpikeRates.Patient(IsEpilepsyMask), 'VariableNames', {'Patient'});
SzFreqEpilepsy = innerjoin(SzFreqFiltered, EpPatients, 'Keys','Patient'); % only keep those with epilepsy

%% ===================== DEFINE MAIN ANALYSIS COHORT =====================
% Epilepsy + documented seizure frequency
CohortPatients = SzFreqEpilepsy.Patient(isfinite(SzFreqEpilepsy.MeanSzFreq));
CohortPatients = unique(CohortPatients);
CohortTable = table(CohortPatients, 'VariableNames', {'Patient'});

%% ===================== RESTRICT ALL DOWNSTREAM TABLES ==================
% So all analyses are restricted to epilepsy patients with documented
% seizure frequency and at least one eeg
PatientLevelSpikeRates = innerjoin(PatientLevelSpikeRates, CohortTable, 'Keys','Patient'); % one row per patient
TypingFiltered         = innerjoin(TypingFiltered,         CohortTable, 'Keys','Patient');

SessionsFiltered      = innerjoin(SessionsFiltered,      CohortTable, 'Keys','Patient'); % one row per eeg (spike rates)
ReportForKeptSessions = innerjoin(ReportForKeptSessions, CohortTable, 'Keys','Patient'); % one row per eeg (report table)

SzFreqEpilepsy = innerjoin(SzFreqEpilepsy, CohortTable, 'Keys','Patient');

% Recompute masks after restriction (aligned to PatientLevelSpikeRates rows)
etype_norm = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));
IsNESDMask = (etype_norm == lower(strtrim(NESD_LABEL)));
IsBadType  = ismember(etype_norm, badTypes) | ismissing(etype_norm) | (strlength(etype_norm)==0);
IsEpilepsyMask = ~IsNESDMask & ~IsBadType;
assert(all(IsEpilepsyMask)) % Now everyone should have epilepsy

fprintf('[Cohort restriction] Using %d epilepsy patients with documented seizure frequency\n', ...
    numel(CohortPatients));

%% ------------------ Build Spearman input tables ------------------
PatientSpikeSz_All = innerjoin( ...
    PatientLevelSpikeRates(IsEpilepsyMask, {'Patient','MeanSpikeRate_perHour'}), ...
    SzFreqEpilepsy, 'Keys','Patient'); % Has spike rate and seizure frequency

keepAll = isfinite(PatientSpikeSz_All.MeanSpikeRate_perHour) & isfinite(PatientSpikeSz_All.MeanSzFreq);
PatientSpikeSz_All = PatientSpikeSz_All(keepAll,:);

% Typed subset: require EpiType3 in canonical3 and not missing
keepCanon3 = ~ismissing(PatientLevelSpikeRates.EpiType3) & ismember(string(PatientLevelSpikeRates.EpiType3), canonical3);
TypedPatients = PatientLevelSpikeRates.Patient(IsEpilepsyMask & keepCanon3);

SzFreqCanon = innerjoin(SzFreqEpilepsy, ...
    PatientLevelSpikeRates(IsEpilepsyMask & keepCanon3, {'Patient','EpiType3'}), ...
    'Keys','Patient');

PatientSpikeSz_Typed = innerjoin( ...
    PatientLevelSpikeRates(IsEpilepsyMask & ismember(PatientLevelSpikeRates.Patient, TypedPatients), ...
        {'Patient','MeanSpikeRate_perHour'}), ...
    SzFreqCanon, 'Keys','Patient');

keepTyped = isfinite(PatientSpikeSz_Typed.MeanSpikeRate_perHour) & ...
            isfinite(PatientSpikeSz_Typed.MeanSzFreq) & ...
            ~ismissing(PatientSpikeSz_Typed.EpiType3);
PatientSpikeSz_Typed = PatientSpikeSz_Typed(keepTyped,:);

%% ------------------ Canonical3 group stats (for Fig 1 subtype panel) ------------------
Canonical3_SubsetTable = PatientLevelSpikeRates(IsEpilepsyMask & keepCanon3, ...
    {'Patient','EpiType3','MeanSpikeRate_perHour'});

% Rename for downstream compatibility if you want to keep existing plotting code:
Canonical3_SubsetTable.Properties.VariableNames{'EpiType3'} = 'EpiType4';

[g3, cats3] = findgroups(Canonical3_SubsetTable.EpiType4);
medVals = splitapply(@(x) median(x,'omitnan'), Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
p25Vals = splitapply(@(x) prctile(x,25),      Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
p75Vals = splitapply(@(x) prctile(x,75),      Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);
nVals   = splitapply(@(x) sum(isfinite(x)),   Canonical3_SubsetTable.MeanSpikeRate_perHour, g3);

Canonical3_Stats = table(cats3, nVals, medVals, p25Vals, p75Vals, ...
    'VariableNames', {'EpiType4','GroupCount','Median','P25','P75'});

Canonical3_Pairs = ["General","Temporal";
                    "General","Frontal";
                    "Temporal","Frontal"];

p_pair = NaN(3,1);
for i = 1:3
    A = Canonical3_Pairs(i,1);
    B = Canonical3_Pairs(i,2);
    xa = Canonical3_SubsetTable.MeanSpikeRate_perHour( ...
        Canonical3_SubsetTable.EpiType4 == categorical(A, canonical3));
    xb = Canonical3_SubsetTable.MeanSpikeRate_perHour( ...
        Canonical3_SubsetTable.EpiType4 == categorical(B, canonical3));
    if nnz(isfinite(xa))>=3 && nnz(isfinite(xb))>=3
        p_pair(i) = ranksum(xa, xb, 'method','approx');
    end
end
p_pair_bonf = min(p_pair*3, 1);

%% ------------------ Bundle outputs ------------------
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

% Panel A: reported spikes present/absent 
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

% Panel B: subtype
Sub = Views.PatientSpikeSz_Typed(:,{'EpiType3','MeanSpikeRate_perHour'});
Sub.Properties.VariableNames{'EpiType3'} = 'EpiType4'; % to minimize downstream edits

Y_C = add_y_jitter_eps(to_log10_per_hour(Sub.MeanSpikeRate_perHour, EPS_RATE), Y_ZERO, Y_LIMS, 0.02);

[p_kw, tbl_kw] = kruskalwallis(Sub.MeanSpikeRate_perHour, Sub.EpiType4, 'off');
SS_total = tbl_kw{end,2}; SS_group = tbl_kw{2,2};
eta2_kw = SS_group/SS_total;

% Draw
f1 = figure('Color','w','Position',[60 60 950 520]);
tiledlayout(f1,1,2,'TileSpacing','compact','Padding','loose');

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


t = title(axC, sprintf('B. Epilepsy subtype'));
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


% --- Pairwise significance bars for Fig 1B (Bonferroni-corrected) ---

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

% Panel B omnibus
Fig1Stats.p_kw_C = p_kw;
Fig1Stats.eta2_kw_C = eta2_kw;

% Panel B pairwise (already in Views)
Fig1Stats.p_pair_bonf = Views.PvalsPairwiseBonf;

% Panel B subtype medians + CI (bootstrap medians per group)
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
assert(p<0.001)

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

yMaxData = max(Y(isfinite(Y)));
pad = 0.10 * range(xLims_log10);
ySig = yMaxData + pad;

if ySig > xLims_log10(2) - 0.02*range(xLims_log10)
    ylim(ax, [xLims_log10(1), ySig + 0.05*range(xLims_log10)]);
end

yl = ylim(ax);
yMaxData = max(Y(isfinite(Y)));

% put bar a bit above the max point
yBar = yMaxData + 0.06*range(yl);

% ensure there is room for the bar + its text
yNeedTop = yBar + 0.10*range(yl);
if yNeedTop > yl(2)
    ylim(ax, [yl(1) yNeedTop]);
    yl = ylim(ax);
end

add_sigbar(ax, 1, 2, yBar, p_label(p));
t = title(ax,'Mean seizure frequency by reported spikes across EEGs');
t.Units = 'normalized';
t.Position(2) = 1.03;   % small bump


set(ax,'FontSize',20);

% ---- Add Ns to x tick labels ----
labels = string(ax.XTickLabel);

labels(labels=="All EEGs: no spikes") = ...
    sprintf('All EEGs: no spikes (N=%d)', n_allAbsent);

labels(labels=="≥1 EEG: spikes present") = ...
    sprintf('≥1 EEG: spikes present (N=%d)', n_anyPresent);

ax.XTickLabel = labels;
ax.XTickLabelRotation = 20;

%% Print stats in the command window so I can put it in the supplemental figure legend
fprintf(['\nThe median [95%% CI] seizure frequency was %1.2f [%1.2f-%1.2f] '...
    'for patients whose EEGs had no reported spikes '...
    ' and %1.2f [%1.2f-%1.2f] for patients whose EEGs had reported spikes '...
    '(p < 0.001, Cliff''s d = %1.2f).\n'], ...
    med1,lo1,hi1,med2,lo2,hi2,cliff_delta(freq_allAbsent, freq_anyPresent))

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

% --- Subtype (among epilepsy) ---
E3 = strtrim(string(PL.EpiType3));
espec = strtrim(string(PL.EpilepsySpecific));

isTemp  = (E3=="Temporal") & Views.IsEpilepsyMask;
isFront = (E3=="Frontal")  & Views.IsEpilepsyMask;
isGen   = (E3=="General")  & Views.IsEpilepsyMask;

isCanon = isTemp | isFront | isGen;

% Unknown: no canonical subtype AND explicitly unclassified/unspecified (or missing)
isUnknown = Views.IsEpilepsyMask & ~isCanon & ( ...
    ismissing(espec) | espec=="" | ...
    espec=="Unclassified or Unspecified" | ...
    espec=="Unknown or MRN not found");

% Other: epilepsy with some non-canonical localization/type info
isOther = Views.IsEpilepsyMask & ~isCanon & ~isUnknown;

% Counts
n_temp    = nnz(isTemp);
n_front   = nnz(isFront);
n_gen     = nnz(isGen);
n_other   = nnz(isOther);
n_subunk  = nnz(isUnknown);



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

%{
OutVar{end+1,1} = "Epilepsy diagnosis"; OutStat{end+1,1} = "";
OutVar{end+1,1} = "    Epilepsy"; OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_epi, 100*n_epi/N_total);
OutVar{end+1,1} = "    PNES";     OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_pnes, 100*n_pnes/N_total);
OutVar{end+1,1} = "    Unknown";  OutStat{end+1,1} = sprintf('%d (%.1f%%)', n_dx_unk, 100*n_dx_unk/N_total);
%}

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


function write_results_html(outPath, Views, SzFreqPerPatient, ...
    Fig1Stats, ...
    SpearmanResults_main, rs_all_main, p_all_main, n_all_main, rho_lo_main, rho_hi_main, subtype_ci_main, ...
    ReportForKeptSessions, MMR)

if ~exist(fileparts(outPath),'dir'), mkdir(fileparts(outPath)); end
fid = fopen(outPath,'w');
if fid==-1, error('Could not open %s', outPath); end

PL = Views.PatientLevelSpikeRates;
AllPatients = PL.Patient;
N_total = numel(AllPatients);

% EEG count (after outpatient+routine filtering)
n_eegs_all = height(ReportForKeptSessions);

% cohort medians + CI
cohortPatients = Views.PatientLevelSpikeRates.Patient;
Sf = innerjoin(table(cohortPatients,'VariableNames',{'Patient'}), SzFreqPerPatient, 'Keys','Patient');
sf_vec = Sf.MeanSzFreq;
sf_vec = sf_vec(isfinite(sf_vec));
sf_med = median(sf_vec,'omitnan');
[~, sf_ci_lo, sf_ci_hi] = bootstrap_median_ci(sf_vec, 5000, 0.05);

sr_vec = PL.MeanSpikeRate_perHour;
sr_vec = sr_vec(isfinite(sr_vec));
sr_med = median(sr_vec,'omitnan');
[~, sr_ci_lo, sr_ci_hi] = bootstrap_median_ci(sr_vec, 5000, 0.05);

fprintf(fid,'<html><head><meta charset="UTF-8"><title>Results</title></head><body>\n');

%% ---- Cohort summary ----
fprintf(fid, '<h2>Cohort summary</h2>\n');
fprintf(fid,['<p>We included %d patients (%d EEGs). Median [95%% CI] monthly seizure ' ...
    'frequency was %1.2f [%1.2f-%1.2f], and median [95%% CI] spikes/hour across EEGs ' ...
    'was %1.2f [%1.2f-%1.2f] (Table 1).</p>'], ...
    N_total, n_eegs_all, ...
    sf_med, sf_ci_lo, sf_ci_hi, ...
    sr_med, sr_ci_lo, sr_ci_hi);

%% ---- Figure 1 ----
fprintf(fid, '<h2>Spike rates by patient groups</h2>\n');

pA_str = format_p_html(Fig1Stats.p_rankSum_A);
fprintf(fid,['<p>Automatically detected spike rates were higher in EEGs with ' ...
    'clinically-reported spikes (median spike rate %.2f [95%% CI %.2f-%.2f] spikes/hour) ' ...
    'than in EEGs without reported spikes (%.2f [%.2f-%.2f] spikes/hour) ' ...
    '(%s, Cliff''s &delta; = %.2f; Fig. 1A), suggesting a low false-positive detection rate. '], ...
    Fig1Stats.m_pre, Fig1Stats.lo_pre, Fig1Stats.hi_pre, ...
    Fig1Stats.m_abs, Fig1Stats.lo_abs, Fig1Stats.hi_abs, ...
    pA_str, Fig1Stats.effectA_cliff);

fprintf(fid,['Spike rates differed across epilepsy subtypes (Kruskal-Wallis %s, ' ...
    '&eta;&sup2; &asymp; %.3f), with higher rates in generalized epilepsy than ' ...
    'temporal or frontal lobe epilepsy (Fig. 1B).</p>'], ...
    format_p_html(Fig1Stats.p_kw_C), Fig1Stats.eta2_kw_C);

%% ---- Figure 2 ----
fprintf(fid, '<h2>Relationship between spike rate and seizure frequency (Figure 2)</h2>\n');

fprintf(fid,['<p>Spike rate and seizure frequency were positively correlated ' ...
    '(N = %d, &rho; = %.2f [95%% CI %.2f-%.2f], %s). '], ...
    n_all_main, rs_all_main, rho_lo_main, rho_hi_main, format_p_html(p_all_main));

fprintf(fid,['Subtype-specific correlations were significant for generalized epilepsy ' ...
    '(N = %d, &rho; = %.2f [%.2f-%.2f], Bonferroni-adjusted %s) and temporal lobe epilepsy ' ...
    '(N = %d, &rho; = %.2f [%.2f-%.2f], %s), but not frontal lobe epilepsy ' ...
    '(N = %d, &rho; = %.2f [%.2f-%.2f], %s; Fig. 2A-D). '], ...
    SpearmanResults_main.N(1), SpearmanResults_main.Spearman_r(1), ...
    subtype_ci_main.ci_lo(1), subtype_ci_main.ci_hi(1), ...
    format_p_html(SpearmanResults_main.p_bonf(1)), ...
    SpearmanResults_main.N(2), SpearmanResults_main.Spearman_r(2), ...
    subtype_ci_main.ci_lo(2), subtype_ci_main.ci_hi(2), ...
    format_p_html(SpearmanResults_main.p_bonf(2)), ...
    SpearmanResults_main.N(3), SpearmanResults_main.Spearman_r(3), ...
    subtype_ci_main.ci_lo(3), subtype_ci_main.ci_hi(3), ...
    format_p_html(SpearmanResults_main.p_bonf(3)));

fprintf(fid,['Results were similar when restricting analyses to patients with detectable ' ...
    'spikes and non-zero seizure frequencies (Fig. S1). ']);
fprintf(fid,['Patients with spikes clinically-reported on at least one EEG had higher ' ...
    'mean seizure frequencies than those without spikes (Fig. S2).</p>']);

%% ---- Figure 3: Mixed effects model ----
fprintf(fid, '<h2>Mixed effects model results (Figure 3)</h2>\n');

% Pull results from MMR — use bootstrap CIs if available, else Laplace
if ~isempty(MMR.BootstrapTable4)
    BT = MMR.BootstrapTable4;
    ci_source = 'bootstrap';
else
    % Fall back to Laplace CIs from fixed effects table
    BT = MMR.FE_logistic_proxonly;
    BT.Properties.VariableNames{'OR_lo'} = 'OR_CI_lo';
    BT.Properties.VariableNames{'OR_hi'} = 'OR_CI_hi';
    ci_source = 'Laplace approximation';
end

fprintf(fid, '<p><em>CIs from %s.</em></p>\n', ci_source);

% Helper to extract a row by term name
getRow = @(tbl, nm) tbl(string(tbl.Term) == nm, :);

% Key terms
r_spike   = getRow(BT, 'LogSpikesPerHour');
r_abslag  = getRow(BT, 'AbsLag_years');
r_dir     = getRow(BT, 'LagDirection');
r_int     = getRow(BT, 'LogSpikesPerHour:AbsLag_years');
r_frontal = getRow(BT, 'EpiType3_cat_Frontal');
r_general = getRow(BT, 'EpiType3_cat_General');

% Pull p-values from fixed effects table (not bootstrap table)
FE4 = MMR.FE_logistic_proxonly;
getP = @(nm) FE4.p(string(FE4.Term) == nm);

% Pair table info
T_mod = MMR.ModelTable;
n_pairs = height(T_mod);
n_pats  = numel(unique(T_mod.Patient));

fprintf(fid,['<p>To formally test the association between spike rate and seizure ' ...
    'occurrence while accounting for repeated EEG-visit observations within patients, ' ...
    'we constructed a table of all EEG-visit pairs (N = %d pairs from %d patients) ' ...
    'and fit a logistic mixed effects model with spike rate, absolute EEG-visit lag, ' ...
    'lag direction, and epilepsy subtype as fixed effects and a random intercept for ' ...
    'patient. '], n_pairs, n_pats);

fprintf(fid,['Higher log spike rates were associated with significantly higher odds of ' ...
    'reporting seizures at a clinic visit (OR = %.2f [95%% CI %.2f-%.2f], %s; Fig. 3A). '], ...
    r_spike.OR, r_spike.OR_CI_lo, r_spike.OR_CI_hi, ...
    format_p_html(getP('LogSpikesPerHour')));

fprintf(fid,['The spike-seizure association was significantly attenuated at greater ' ...
    'temporal distances between the EEG and clinic visit (log spike rate &times; ' ...
    'absolute lag interaction: OR = %.3f [%.3f-%.3f], %s; Fig. 3B), indicating that ' ...
    'spike rate is most informative when measured close in time to the clinical assessment. '], ...
    r_int.OR, r_int.OR_CI_lo, r_int.OR_CI_hi, ...
    format_p_html(getP('LogSpikesPerHour:AbsLag_years')));

% LRT results
fprintf(fid,['A likelihood ratio test confirmed that the proximity interaction ' ...
    'significantly improved model fit over a model without interactions ' ...
    '(&chi;&sup2;(1), %s). '], ...
    format_p_html(8.44e-7)); % LRT 3 p-value — hardcoded for now, see note below

fprintf(fid,['Clinic visits occurring after the EEG were associated with lower ' ...
    'baseline odds of seizure reporting than visits before the EEG ' ...
    '(OR = %.2f [%.2f-%.2f], %s), consistent with patients improving over time ' ...
    'following their initial EEG evaluation. '], ...
    r_dir.OR, r_dir.OR_CI_lo, r_dir.OR_CI_hi, ...
    format_p_html(getP('LagDirection')));

fprintf(fid,['Compared with temporal lobe epilepsy, generalized epilepsy was ' ...
    'associated with lower baseline odds of seizure reporting at the same spike rate ' ...
    '(OR = %.2f [%.2f-%.2f], %s), while frontal lobe epilepsy did not differ ' ...
    'significantly (OR = %.2f [%.2f-%.2f], %s). '], ...
    r_general.OR, r_general.OR_CI_lo, r_general.OR_CI_hi, ...
    format_p_html(getP('EpiType3_cat_General')), ...
    r_frontal.OR, r_frontal.OR_CI_lo, r_frontal.OR_CI_hi, ...
    format_p_html(getP('EpiType3_cat_Frontal')));

fprintf(fid,['Results were consistent when using log-transformed continuous seizure ' ...
    'frequency as the outcome in a linear mixed effects model (Table S1).</p>\n']);

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
%% ======================= SPEARMAN FIG FUNCTION =======================
%% =====================================================================

function [SpearmanResults, rs_all, p_all, n_all, rho_lo, rho_hi, subtype_ci] = ...
    spearman_plotting_function(PatientSpikeSz_All, PatientSpikeSz_Typed, ...
                              canonical3, spearman_xLims, spearman_yLims, ...
                              fig_out, ...
                              freqFieldName, labelSuffix, nonZeroOnly)

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
    mask_all = mask_all & (x_all > 0) & (y_all > 0); % exclude zeros for both spike rate AND seizure frequency
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
        mask = mask & (x > 0) & (y > 0); % require non zero
    %spikes and non zero szs
        %mask = mask & (x > 0); % only require non-zero spikes, allow zero szs
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

% Remove the label for the maximum x tick (keep the tick itself)
labs = string(axA.XTickLabel);
[~, iMax] = max(axA.XTick);
labs(iMax) = "";
axA.XTickLabel = labs;


title(axA, sprintf('A. All epilepsy%s (N=%d)', labelSuffix, n_all), ...
    'FontSize', fontL, 'FontWeight','bold');
%{
title(axA, { ...
    'A. All epilepsy', ...
    labelSuffix, ...
    sprintf('(N=%d)', n_all) ...
    }, 'FontSize', fontL, 'FontWeight','bold');

%}
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


        idx = idx & (PatientSpikeSz_Typed.MeanSpikeRate_perHour > 0); % only require non zero spike rate
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

    % Remove the label for the maximum x tick (keep the tick itself)
    labs = string(ax.XTickLabel);
    [~, iMax] = max(ax.XTick);
    labs(iMax) = "";
    ax.XTickLabel = labs;

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

% --- label-specific vertical offset ---
if ptext == "**" || ptext == "***"
    yOff = -0.012 * diff(ax.YLim);   
else
    yOff = 0.003 * diff(ax.YLim);   
end

text(ax, mean([x1 x2]), y + yOff, ptext, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize',20);
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


function PairTable = build_eeg_visit_pairs(Vuniq, SessionLevelSpikeRates, ...
                                            ReportForKeptSessions, PatientTyping)
%% build_eeg_visit_pairs
% Builds a table with one row per (Patient, EEG, Visit) combination.
%
% Inputs:
%   Vuniq                  : visit-level table with Patient, VisitDate, Freq_R1, HasSz
%   SessionLevelSpikeRates : EEG-level table with Patient, Session, SpikesPerHour
%   ReportForKeptSessions  : EEG report table with Patient, Session, start_time_deid
%   PatientTyping          : patient-level table with Patient, EpiType3
%
% Output:
%   PairTable : one row per (Patient, EEG Session, Visit) with:
%       Patient, Session, VisitDate, SpikesPerHour,
%       SzFreq, HasSz, SignedLag_days, EpiType3

%% ---- Parse EEG datetimes ----
EEG_raw = ReportForKeptSessions.start_time_deid;
if isdatetime(EEG_raw)
    EEG_dt = EEG_raw;
else
    EEG_dt = datetime(strtrim(string(EEG_raw)), ...
        'InputFormat', "yyyy-MM-dd'T'HH:mm:ss");
end

% Build slim EEG date table (one row per EEG)
EEG_dates = table(...
    double(ReportForKeptSessions.Patient), ...
    double(ReportForKeptSessions.Session), ...
    EEG_dt, ...
    'VariableNames', {'Patient','Session','EEG_Date'});

% Remove EEGs with unparseable dates
EEG_dates = EEG_dates(~isnat(EEG_dates.EEG_Date), :);

%% ---- Join spike rates to EEG dates ----
% SessionLevelSpikeRates has Patient, Session, SpikesPerHour
EEG_tbl = innerjoin(...
    EEG_dates, ...
    SessionLevelSpikeRates(:, {'Patient','Session','SpikesPerHour'}), ...
    'Keys', {'Patient','Session'});

%% ---- Join epilepsy type to visits ----
% PatientTyping has Patient, EpiType3
Vtyped = innerjoin(...
    Vuniq(:, {'Patient','VisitDate','Freq_R1','HasSz'}), ...
    PatientTyping(:, {'Patient','EpiType3'}), ...
    'Keys', 'Patient');

%% ---- Cross-join EEGs x Visits within each patient ----
patients = unique(EEG_tbl.Patient);
patients = intersect(patients, unique(Vtyped.Patient)); % only patients in both

nEstimate = 0; % preallocate estimate
for i = 1:numel(patients)
    p = patients(i);
    nE = sum(EEG_tbl.Patient == p);
    nV = sum(Vtyped.Patient == p);
    nEstimate = nEstimate + nE * nV;
end

% Preallocate output arrays
Patient_out      = nan(nEstimate, 1);
Session_out      = nan(nEstimate, 1);
VisitDate_out    = NaT(nEstimate, 1);
SpikesPerHour_out = nan(nEstimate, 1);
SzFreq_out       = nan(nEstimate, 1);
HasSz_out        = nan(nEstimate, 1);
SignedLag_out    = nan(nEstimate, 1);
EpiType3_out     = strings(nEstimate, 1);

row = 0;

for i = 1:numel(patients)
    p = patients(i);

    % All EEGs for this patient
    eeg_rows = EEG_tbl(EEG_tbl.Patient == p, :);
    % All visits for this patient
    vis_rows = Vtyped(Vtyped.Patient == p, :);

    nE = height(eeg_rows);
    nV = height(vis_rows);

    % Cross join: all EEG x visit combinations
    for e = 1:nE
        for v = 1:nV
            row = row + 1;

            % Signed lag: positive = visit AFTER EEG, negative = visit BEFORE EEG
            lag_days = days(vis_rows.VisitDate(v) - eeg_rows.EEG_Date(e));

            Patient_out(row)       = p;
            Session_out(row)       = eeg_rows.Session(e);
            VisitDate_out(row)     = vis_rows.VisitDate(v);
            SpikesPerHour_out(row) = eeg_rows.SpikesPerHour(e);
            SzFreq_out(row)        = vis_rows.Freq_R1(v);
            HasSz_out(row)         = vis_rows.HasSz(v);
            SignedLag_out(row)     = lag_days;
            EpiType3_out(row)      = string(vis_rows.EpiType3(v));
        end
    end
end

% Trim to actual rows used
Patient_out       = Patient_out(1:row);
Session_out       = Session_out(1:row);
VisitDate_out     = VisitDate_out(1:row);
SpikesPerHour_out = SpikesPerHour_out(1:row);
SzFreq_out        = SzFreq_out(1:row);
HasSz_out         = HasSz_out(1:row);
SignedLag_out     = SignedLag_out(1:row);
EpiType3_out      = EpiType3_out(1:row);

%% ---- Assemble output table ----
PairTable = table(...
    Patient_out, ...
    Session_out, ...
    VisitDate_out, ...
    SpikesPerHour_out, ...
    SzFreq_out, ...
    HasSz_out, ...
    SignedLag_out, ...
    EpiType3_out, ...
    'VariableNames', {...
        'Patient', 'Session', 'VisitDate', ...
        'SpikesPerHour', 'SzFreq', 'HasSz', ...
        'SignedLag_days', 'EpiType3'});

%% ---- Create globally unique EEG identifier ----
PairTable.EEG_ID = categorical(...
    string(PairTable.Patient) + "_" + string(PairTable.Session));

%% ---- Remove rows with missing key variables ----
keepMask = isfinite(PairTable.SpikesPerHour) & ...
           isfinite(PairTable.SzFreq) & ...
           isfinite(PairTable.SignedLag_days) & ...
           strlength(PairTable.EpiType3) > 0;

n_before = height(PairTable);
PairTable = PairTable(keepMask, :);
n_after  = height(PairTable);

fprintf('[build_eeg_visit_pairs] %d patients, %d EEG-visit pairs (%d removed for missing data)\n', ...
    numel(patients), n_after, n_before - n_after);

%% ---- Normalize spike rate for model stability ----
% Log-transform spike rate (adding small offset for zeros)
EPS_SPIKE = 1e-3;
PairTable.LogSpikesPerHour = log(PairTable.SpikesPerHour + EPS_SPIKE);

% Scale signed lag to years (easier to interpret coefficients)
PairTable.SignedLag_years = PairTable.SignedLag_days / 365.25;

% Convert Patient to categorical for use as random effect grouping variable
PairTable.PatientID = categorical(PairTable.Patient);

end

function MixedModelResults = fit_mixed_effects_models(PairTable, nBoot, alpha)
%% fit_mixed_effects_models
% Fits linear and logistic mixed effects models for the association between
% spike rate and seizure frequency/occurrence, with absolute time lag,
% lag direction, and epilepsy subtype as fixed effects, and patient as
% random effect.
%
% Inputs:
%   PairTable : output of build_eeg_visit_pairs
%   nBoot     : number of bootstrap iterations (0 to skip, 5000 for final)
%   alpha     : significance level (e.g. 0.05)
%
% Output:
%   MixedModelResults : struct with model objects, summary tables, and
%                       bootstrap CIs

if nargin < 2, nBoot = 0;    end
if nargin < 3, alpha = 0.05; end

%% ================================================================
%% PREPARE MODEL TABLE
%% ================================================================
canonical3 = ["Frontal","General","Temporal"];

keepMask = ismember(string(PairTable.EpiType3), canonical3) & ...
           isfinite(PairTable.LogSpikesPerHour) & ...
           isfinite(PairTable.SignedLag_years)  & ...
           isfinite(PairTable.HasSz);

T = PairTable(keepMask, :);

% Epilepsy type: Temporal as reference (largest group)
T.EpiType3_cat = categorical(string(T.EpiType3), canonical3);
T.EpiType3_cat = reordercats(T.EpiType3_cat, ["Temporal","Frontal","General"]);

% Binary outcome
T.HasSz_bin = double(T.HasSz == 1);

% Continuous outcome (log-transformed with epsilon offset)
EPS_SZ      = 1e-3;
T.LogSzFreq = log(T.SzFreq + EPS_SZ);

% Absolute lag: proximity regardless of direction
T.AbsLag_years = abs(T.SignedLag_years);

% Direction: +1 if visit after EEG, -1 if before
% Same-day visits treated as after (direction = +1)
T.LagDirection                          = sign(T.SignedLag_years);
T.LagDirection(T.SignedLag_years == 0)  = 1;

% Report lag statistics
n_after  = sum(T.LagDirection ==  1);
n_before = sum(T.LagDirection == -1);
fprintf('[Lag direction] After EEG: %d (%.1f%%), Before EEG: %d (%.1f%%)\n', ...
    n_after,  100*n_after /height(T), ...
    n_before, 100*n_before/height(T));
fprintf('[Abs lag] Median %.2f years, IQR [%.2f, %.2f], range [%.2f, %.2f]\n', ...
    median(T.AbsLag_years), ...
    prctile(T.AbsLag_years, 25), ...
    prctile(T.AbsLag_years, 75), ...
    min(T.AbsLag_years), ...
    max(T.AbsLag_years));
fprintf('[fit_mixed_effects_models] Fitting on %d pairs, %d patients\n', ...
    height(T), numel(unique(T.Patient)));

%% ================================================================
%% MODEL FORMULAS
%% ================================================================

% Model 1: Full model — proximity + direction interactions
formula_logistic = ['HasSz_bin ~ ' ...
    'LogSpikesPerHour * AbsLag_years + ' ...
    'LogSpikesPerHour * LagDirection + ' ...
    'EpiType3_cat + (1|PatientID)'];

% Model 2: Linear sensitivity — same fixed effects as Model 1
formula_linear = ['LogSzFreq ~ ' ...
    'LogSpikesPerHour * AbsLag_years + ' ...
    'LogSpikesPerHour * LagDirection + ' ...
    'EpiType3_cat + (1|PatientID)'];

% Model 3: No interactions (for LRT baseline)
formula_logistic_noint = ['HasSz_bin ~ ' ...
    'LogSpikesPerHour + AbsLag_years + LagDirection + ' ...
    'EpiType3_cat + (1|PatientID)'];

% Model 4: Proximity interaction only — PRIMARY MODEL
formula_logistic_proxonly = ['HasSz_bin ~ ' ...
    'LogSpikesPerHour * AbsLag_years + ' ...
    'LagDirection + ' ...
    'EpiType3_cat + (1|PatientID)'];

%% ================================================================
%% FIT MODEL 1: Logistic — proximity + direction interactions
%% ================================================================
fprintf('\nFitting Model 1: Logistic — proximity + direction interactions...\n');
try
    mdl_logistic = fitglme(T, formula_logistic, ...
        'Distribution',      'Binomial', ...
        'Link',              'logit', ...
        'FitMethod',         'Laplace', ...
        'CovariancePattern', 'Diagonal');
    fprintf('Model 1 converged successfully.\n');
    disp(mdl_logistic);
catch ME
    disp(['Model 1 failed: ' ME.message]);
    mdl_logistic = [];
end

%% ================================================================
%% FIT MODEL 2: Linear — proximity + direction interactions
%% ================================================================
fprintf('\nFitting Model 2: Linear — proximity + direction interactions...\n');
try
    mdl_linear = fitlme(T, formula_linear, ...
        'FitMethod', 'REML');
    fprintf('Model 2 converged successfully.\n');
    disp(mdl_linear);
catch ME
    disp(['Model 2 failed: ' ME.message]);
    mdl_linear = [];
end

%% ================================================================
%% FIT MODEL 3: Logistic — main effects only (LRT baseline)
%% ================================================================
fprintf('\nFitting Model 3: Logistic — main effects only...\n');
try
    mdl_logistic_noint = fitglme(T, formula_logistic_noint, ...
        'Distribution',      'Binomial', ...
        'Link',              'logit', ...
        'FitMethod',         'Laplace', ...
        'CovariancePattern', 'Diagonal');
    fprintf('Model 3 converged successfully.\n');
    disp(mdl_logistic_noint);
catch ME
    disp(['Model 3 failed: ' ME.message]);
    mdl_logistic_noint = [];
end

%% ================================================================
%% FIT MODEL 4: Logistic — proximity interaction only (PRIMARY)
%% ================================================================
fprintf('\nFitting Model 4 (PRIMARY): Logistic — proximity interaction only...\n');
try
    mdl_logistic_proxonly = fitglme(T, formula_logistic_proxonly, ...
        'Distribution',      'Binomial', ...
        'Link',              'logit', ...
        'FitMethod',         'Laplace', ...
        'CovariancePattern', 'Diagonal');
    fprintf('Model 4 converged successfully.\n');
    disp(mdl_logistic_proxonly);
catch ME
    disp(['Model 4 failed: ' ME.message]);
    mdl_logistic_proxonly = [];
end

%% ================================================================
%% LIKELIHOOD RATIO TESTS
%% ================================================================
fprintf('\n========================================\n');
fprintf('LIKELIHOOD RATIO TESTS\n');
fprintf('========================================\n');

% LRT 1: Both interactions vs no interactions
if ~isempty(mdl_logistic) && ~isempty(mdl_logistic_noint)
    fprintf('\nLRT 1: Both interactions vs no interactions\n');
    lrt1 = compare(mdl_logistic_noint, mdl_logistic);
    disp(lrt1);
end

% LRT 2: Both interactions vs proximity only (tests direction interaction)
if ~isempty(mdl_logistic) && ~isempty(mdl_logistic_proxonly)
    fprintf('\nLRT 2: Both interactions vs proximity interaction only\n');
    lrt2 = compare(mdl_logistic_proxonly, mdl_logistic);
    disp(lrt2);
end

% LRT 3: Proximity interaction vs no interactions (tests proximity interaction)
if ~isempty(mdl_logistic_proxonly) && ~isempty(mdl_logistic_noint)
    fprintf('\nLRT 3: Proximity interaction only vs no interactions\n');
    lrt3 = compare(mdl_logistic_noint, mdl_logistic_proxonly);
    disp(lrt3);
end

%% ================================================================
%% FIXED EFFECTS SUMMARY — ALL MODELS
%% ================================================================
fprintf('\n========================================\n');
fprintf('FIXED EFFECTS SUMMARY\n');
fprintf('========================================\n');

T_fe1 = []; T_fe2 = []; T_fe3 = []; T_fe4 = [];

if ~isempty(mdl_logistic)
    fprintf('\nModel 1 (Logistic, proximity + direction interactions):\n');
    [beta, betanames, stats] = fixedEffects(mdl_logistic, 'Alpha', alpha);
    T_fe1 = table(...
        string(betanames.Name), beta, stats.SE, stats.tStat, stats.pValue, ...
        exp(beta), exp(beta - 1.96*stats.SE), exp(beta + 1.96*stats.SE), ...
        'VariableNames', {'Term','Beta','SE','t','p','OR','OR_lo','OR_hi'});
    disp(T_fe1);
end

if ~isempty(mdl_linear)
    fprintf('\nModel 2 (Linear, proximity + direction interactions):\n');
    [beta2, betanames2, stats2] = fixedEffects(mdl_linear, 'Alpha', alpha);
    T_fe2 = table(...
        string(betanames2.Name), beta2, stats2.SE, stats2.tStat, stats2.pValue, ...
        'VariableNames', {'Term','Beta','SE','t','p'});
    disp(T_fe2);
end

if ~isempty(mdl_logistic_noint)
    fprintf('\nModel 3 (Logistic, no interactions):\n');
    [beta3, betanames3, stats3] = fixedEffects(mdl_logistic_noint, 'Alpha', alpha);
    T_fe3 = table(...
        string(betanames3.Name), beta3, stats3.SE, stats3.tStat, stats3.pValue, ...
        exp(beta3), exp(beta3 - 1.96*stats3.SE), exp(beta3 + 1.96*stats3.SE), ...
        'VariableNames', {'Term','Beta','SE','t','p','OR','OR_lo','OR_hi'});
    disp(T_fe3);
end

if ~isempty(mdl_logistic_proxonly)
    fprintf('\nModel 4 (PRIMARY — Logistic, proximity interaction only):\n');
    [beta4, betanames4, stats4] = fixedEffects(mdl_logistic_proxonly, 'Alpha', alpha);
    T_fe4 = table(...
        string(betanames4.Name), beta4, stats4.SE, stats4.tStat, stats4.pValue, ...
        exp(beta4), exp(beta4 - 1.96*stats4.SE), exp(beta4 + 1.96*stats4.SE), ...
        'VariableNames', {'Term','Beta','SE','t','p','OR','OR_lo','OR_hi'});
    disp(T_fe4);
end

%% ================================================================
%% BOOTSTRAP CIs — MODEL 4 (PRIMARY)
%% ================================================================
T_boot4     = [];
boot_betas4 = [];

if ~isempty(mdl_logistic_proxonly) && nBoot > 0
    fprintf('\nBootstrapping Model 4 (primary) CIs (%d iterations)...\n', nBoot);

    patients    = unique(T.PatientID);
    nPat        = numel(patients);
    nFixed4     = size(fixedEffects(mdl_logistic_proxonly), 1);
    boot_betas4 = nan(nBoot, nFixed4);

    parfor b = 1:nBoot
        idx      = randi(nPat, nPat, 1);
        bootPats = patients(idx);

        Tboot = cell(nPat, 1);
        for k = 1:nPat
            Tboot{k} = T(T.PatientID == bootPats(k), :);
            Tboot{k}.PatientID = categorical(repmat(k, height(Tboot{k}), 1));
        end
        Tboot = vertcat(Tboot{:});

        try
            mboot = fitglme(Tboot, formula_logistic_proxonly, ...
                'Distribution',      'Binomial', ...
                'Link',              'logit', ...
                'FitMethod',         'Laplace', ...
                'CovariancePattern', 'Diagonal');
            boot_betas4(b,:) = fixedEffects(mboot)';
        catch
            % leave as NaN — excluded below
        end
    end

    % Remove failed iterations
    converged4  = all(isfinite(boot_betas4), 2);
    boot_betas4 = boot_betas4(converged4, :);
    fprintf('Model 4 bootstrap: %d/%d iterations converged (%.1f%%)\n', ...
        sum(converged4), nBoot, 100*mean(converged4));

    ci_lo4 = prctile(boot_betas4, 100*(alpha/2),   1);
    ci_hi4 = prctile(boot_betas4, 100*(1-alpha/2), 1);

    [beta_obs4, betanames_obs4] = fixedEffects(mdl_logistic_proxonly);

    T_boot4 = table(...
        string(betanames_obs4.Name), ...
        beta_obs4, ...
        ci_lo4(:), ci_hi4(:), ...
        exp(beta_obs4), exp(ci_lo4(:)), exp(ci_hi4(:)), ...
        'VariableNames', {'Term','Beta','Boot_CI_lo','Boot_CI_hi', ...
                          'OR','OR_CI_lo','OR_CI_hi'});

    fprintf('\nModel 4 bootstrapped ORs and CIs:\n');
    disp(T_boot4);
end

%% ================================================================
%% BOOTSTRAP CIs — MODEL 1 (SUPPLEMENTAL)
%% ================================================================
T_boot1     = [];
boot_betas1 = [];

if ~isempty(mdl_logistic) && nBoot > 0
    fprintf('\nBootstrapping Model 1 (supplemental) CIs (%d iterations)...\n', nBoot);

    patients    = unique(T.PatientID);
    nPat        = numel(patients);
    nFixed1     = size(fixedEffects(mdl_logistic), 1);
    boot_betas1 = nan(nBoot, nFixed1);

    parfor b = 1:nBoot
        idx      = randi(nPat, nPat, 1);
        bootPats = patients(idx);

        Tboot = cell(nPat, 1);
        for k = 1:nPat
            Tboot{k} = T(T.PatientID == bootPats(k), :);
            Tboot{k}.PatientID = categorical(repmat(k, height(Tboot{k}), 1));
        end
        Tboot = vertcat(Tboot{:});

        try
            mboot = fitglme(Tboot, formula_logistic, ...
                'Distribution',      'Binomial', ...
                'Link',              'logit', ...
                'FitMethod',         'Laplace', ...
                'CovariancePattern', 'Diagonal');
            boot_betas1(b,:) = fixedEffects(mboot)';
        catch
            % leave as NaN
        end
    end

    converged1  = all(isfinite(boot_betas1), 2);
    boot_betas1 = boot_betas1(converged1, :);
    fprintf('Model 1 bootstrap: %d/%d iterations converged (%.1f%%)\n', ...
        sum(converged1), nBoot, 100*mean(converged1));

    ci_lo1 = prctile(boot_betas1, 100*(alpha/2),   1);
    ci_hi1 = prctile(boot_betas1, 100*(1-alpha/2), 1);

    [beta_obs1, betanames_obs1] = fixedEffects(mdl_logistic);

    T_boot1 = table(...
        string(betanames_obs1.Name), ...
        beta_obs1, ...
        ci_lo1(:), ci_hi1(:), ...
        exp(beta_obs1), exp(ci_lo1(:)), exp(ci_hi1(:)), ...
        'VariableNames', {'Term','Beta','Boot_CI_lo','Boot_CI_hi', ...
                          'OR','OR_CI_lo','OR_CI_hi'});

    fprintf('\nModel 1 bootstrapped ORs and CIs:\n');
    disp(T_boot1);
end

%% ================================================================
%% DIAGNOSTIC PLOTS FOR BOOTSTRAP DISTRIBUTIONS
%% ================================================================
if ~isempty(boot_betas4) && size(boot_betas4,1) > 10
    [~, betanames_diag] = fixedEffects(mdl_logistic_proxonly);
    term_names = string(betanames_diag.Name);
    nTerms     = numel(term_names);

    nCols = ceil(sqrt(nTerms));
    nRows = ceil(nTerms / nCols);

    figDiag = figure('Color','w', ...
        'Position',[100 100 300*nCols 250*nRows], ...
        'Name','Bootstrap distributions — Model 4');
    tl_diag = tiledlayout(figDiag, nRows, nCols, ...
        'TileSpacing','compact','Padding','loose');

    [beta_diag, ~, stats_diag] = fixedEffects(mdl_logistic_proxonly);

    for k = 1:nTerms
        ax = nexttile(tl_diag);
        histogram(ax, boot_betas4(:,k), 40, ...
            'FaceColor',[0.2 0.5 0.8],'EdgeColor','none','FaceAlpha',0.7);
        xline(ax, beta_diag(k), 'r-',  'LineWidth', 2);
        xline(ax, prctile(boot_betas4(:,k), 2.5),  'k--', 'LineWidth', 1.5);
        xline(ax, prctile(boot_betas4(:,k), 97.5), 'k--', 'LineWidth', 1.5);
        title(ax, term_names(k), 'FontSize', 9, 'Interpreter','none');
        xlabel(ax, '\beta', 'FontSize', 8);
        box(ax,'off');
    end

    sgtitle(figDiag, 'Bootstrap distributions — Model 4 (red=observed, dashed=95% CI)', ...
        'FontSize', 11);
end

%% ================================================================
%% BUNDLE OUTPUT
%% ================================================================
MixedModelResults.ModelTable             = T;

% Model objects
MixedModelResults.mdl_logistic           = mdl_logistic;
MixedModelResults.mdl_linear             = mdl_linear;
MixedModelResults.mdl_logistic_noint     = mdl_logistic_noint;
MixedModelResults.mdl_logistic_proxonly  = mdl_logistic_proxonly;

% Fixed effects tables (Laplace CIs)
MixedModelResults.FE_logistic            = T_fe1;
MixedModelResults.FE_linear              = T_fe2;
MixedModelResults.FE_logistic_noint      = T_fe3;
MixedModelResults.FE_logistic_proxonly   = T_fe4;  % primary

% Bootstrap results — Model 4 (primary)
MixedModelResults.BootstrapBetas4        = boot_betas4;
MixedModelResults.BootstrapTable4        = T_boot4;

% Bootstrap results — Model 1 (supplemental)
MixedModelResults.BootstrapBetas1        = boot_betas1;
MixedModelResults.BootstrapTable1        = T_boot1;

fprintf('\nDone. Primary model is Model 4 (mdl_logistic_proxonly).\n');
if ~isempty(T_boot4)
    fprintf('Bootstrap CIs available in MMR.BootstrapTable4\n');
else
    fprintf('No bootstrap CIs computed (nBoot=0). Rerun with nBoot>0 for final results.\n');
end

end


function [FigMain, FigCV] = make_figures_and_cv(MMR, Views, Vuniq, ...
                                                  ReportForKeptSessions, ...
                                                  canonical3, nFolds, ...
                                                  outPath_main, outPath_cv)
%% make_figures_and_cv
% Creates:
%   1. FigMain : Combined context + model figure (4 panels A-D)
%   2. FigCV   : Cross-validated classifier figure (ROC + PR curves)
%
% Inputs:
%   MMR                   : output of fit_mixed_effects_models
%   Views                 : output of build_filtered_view
%   Vuniq                 : visit-level table (Patient, VisitDate, HasSz, Freq_R1)
%   ReportForKeptSessions : EEG report table (Patient, start_time_deid)
%   canonical3            : string array e.g. ["Temporal","Frontal","General"]
%   nFolds                : number of CV folds (e.g. 5)
%   outPath_main          : output path for main figure
%   outPath_cv            : output path for CV figure

if nargin < 6, nFolds      = 5;  end
if nargin < 7, outPath_main = ''; end
if nargin < 8, outPath_cv   = ''; end

FONT_SIZE = 20;
T         = MMR.ModelTable;
mdl       = MMR.mdl_logistic_proxonly;

%% ================================================================
%% FIGURE 1: COMBINED CONTEXT + MODEL (4 panels)
%% A: EEG and visit timing
%% B: Seizure occurrence over time
%% C: EEG-visit lag distribution
%% D: Forest plot
%% ================================================================

FigMain = figure('Color','w','Position',[60 60 1100 900]);
tl = tiledlayout(FigMain, 2, 2, 'TileSpacing','compact', 'Padding','loose');

%% ----- Parse EEG dates -----
EEG_raw = ReportForKeptSessions.start_time_deid;
if isdatetime(EEG_raw)
    EEG_dt = EEG_raw;
else
    EEG_dt = datetime(strtrim(string(EEG_raw)), ...
        'InputFormat', "yyyy-MM-dd'T'HH:mm:ss");
end

refDate   = datetime(2000, 1, 1);
EEG_yrs   = days(EEG_dt - refDate) / 365.25;
EEG_yrs   = EEG_yrs(isfinite(EEG_yrs));
Visit_yrs = days(Vuniq.VisitDate - refDate) / 365.25;
Visit_yrs = Visit_yrs(isfinite(Visit_yrs));

%% ----- Panel A: EEG timing vs visit timing -----
axA = nexttile(tl, 1);
hold(axA,'on'); box(axA,'off'); grid(axA,'on');

bin_edges_t = 0:1:15;

histogram(axA, EEG_yrs, bin_edges_t, ...
    'FaceColor',[0.2 0.5 0.8],'FaceAlpha',0.6,'EdgeColor','none', ...
    'Normalization','probability','DisplayName','EEG recordings');
histogram(axA, Visit_yrs, bin_edges_t, ...
    'FaceColor',[0.8 0.3 0.1],'FaceAlpha',0.6,'EdgeColor','none', ...
    'Normalization','probability','DisplayName','Clinic visits');

med_eeg   = median(EEG_yrs,   'omitnan');
med_visit = median(Visit_yrs, 'omitnan');
yl = ylim(axA);

xline(axA, med_eeg,   '--','Color',[0.2 0.5 0.8],'LineWidth',1.5, ...
    'HandleVisibility','off');
xline(axA, med_visit, '--','Color',[0.8 0.3 0.1],'LineWidth',1.5, ...
    'HandleVisibility','off');

text(axA, med_eeg+0.15,   yl(2)*0.90, sprintf('Median\n%.1fy', med_eeg), ...
    'Color',[0.2 0.5 0.8],'FontSize',FONT_SIZE-6);
text(axA, med_visit+0.15, yl(2)*0.72, sprintf('Median\n%.1fy', med_visit), ...
    'Color',[0.8 0.3 0.1],'FontSize',FONT_SIZE-6);

xlabel(axA, 'Years since first visit', 'FontSize', FONT_SIZE);
ylabel(axA, 'Proportion',              'FontSize', FONT_SIZE);
title(axA,  'A. EEG and visit timing', 'FontSize', FONT_SIZE, 'FontWeight','bold');
legend(axA, 'Location','northeast',    'FontSize', FONT_SIZE-6);
set(axA, 'FontSize', FONT_SIZE);

%% ----- Panel B: Proportion of visits with HasSz == 1 over time -----
axB = nexttile(tl, 2);
hold(axB,'on'); box(axB,'off'); grid(axB,'on');

cohort_patients = unique(T.Patient);
V_cohort = Vuniq(ismember(Vuniq.Patient, cohort_patients) & ...
                 (Vuniq.HasSz == 0 | Vuniq.HasSz == 1), :);
V_cohort.YearsSinceFirst = days(V_cohort.VisitDate - refDate) / 365.25;

bin_edges_sz   = 0:1:14;
bin_centers_sz = bin_edges_sz(1:end-1) + 0.5;
nBins_sz       = numel(bin_centers_sz);
bin_prop_sz    = nan(nBins_sz, 1);
bin_lo_sz      = nan(nBins_sz, 1);
bin_hi_sz      = nan(nBins_sz, 1);
bin_n_sz       = zeros(nBins_sz, 1);

for b = 1:nBins_sz
    mask = V_cohort.YearsSinceFirst >= bin_edges_sz(b) & ...
           V_cohort.YearsSinceFirst <  bin_edges_sz(b+1);
    vals = V_cohort.HasSz(mask);
    vals = vals(isfinite(vals));
    if numel(vals) >= 10
        bin_prop_sz(b) = mean(vals);
        bin_n_sz(b)    = numel(vals);
        boot_p = nan(500,1);
        for bb = 1:500
            boot_p(bb) = mean(vals(randi(numel(vals), numel(vals), 1)));
        end
        bin_lo_sz(b) = prctile(boot_p, 2.5);
        bin_hi_sz(b) = prctile(boot_p, 97.5);
    end
end

validB = isfinite(bin_prop_sz);

patch(axB, ...
    [bin_centers_sz(validB), fliplr(bin_centers_sz(validB))], ...
    [bin_lo_sz(validB)', fliplr(bin_hi_sz(validB)')], ...
    [0.8 0.3 0.1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

plot(axB, bin_centers_sz(validB), bin_prop_sz(validB), 'o-', ...
    'Color',[0.8 0.3 0.1],'LineWidth',2, ...
    'MarkerFaceColor',[0.8 0.3 0.1],'MarkerSize',6);

for b = find(validB)'
    scatter(axB, bin_centers_sz(b), bin_prop_sz(b), ...
        bin_n_sz(b)/5, [0.8 0.3 0.1],'filled','MarkerFaceAlpha',0.25);
end

yline(axB, 0.5, 'k--', 'LineWidth', 1.2, ...
    'Label','50%','LabelHorizontalAlignment','right', ...
    'FontSize', FONT_SIZE-6);

overall_prop = mean(V_cohort.HasSz);
text(axB, 0.98, 0.95, sprintf('Overall: %.0f%%', 100*overall_prop), ...
    'Units','normalized','HorizontalAlignment','right', ...
    'FontSize',FONT_SIZE-6,'Color',[0.5 0.5 0.5]);

ylim(axB, [0 1]);
xlabel(axB, 'Years since first visit',             'FontSize', FONT_SIZE);
ylabel(axB, 'Proportion of visits with seizures',  'FontSize', FONT_SIZE);
title(axB,  'B. Seizure occurrence over time',     'FontSize', FONT_SIZE, 'FontWeight','bold');
set(axB, 'FontSize', FONT_SIZE);

%% ----- Panel C: Absolute lag distribution -----
axC = nexttile(tl, 3);
hold(axC,'on'); box(axC,'off'); grid(axC,'on');

abs_lags = T.AbsLag_years;

histogram(axC, abs_lags, 40, ...
    'FaceColor',[0.3 0.3 0.3],'FaceAlpha',0.6,'EdgeColor','none', ...
    'Normalization','probability');

med_lag = median(abs_lags);
q1_lag  = prctile(abs_lags, 25);
q3_lag  = prctile(abs_lags, 75);
yl      = ylim(axC);

xline(axC, med_lag, 'k-',  'LineWidth', 2.0);
xline(axC, q1_lag,  'k--', 'LineWidth', 1.5);
xline(axC, q3_lag,  'k--', 'LineWidth', 1.5);

text(axC, med_lag+0.1, yl(2)*0.93, sprintf('Median: %.1fy', med_lag), ...
    'FontSize', FONT_SIZE-6);
text(axC, q1_lag+0.1,  yl(2)*0.80, sprintf('Q1: %.1fy', q1_lag), ...
    'FontSize', FONT_SIZE-6);
text(axC, q3_lag+0.1,  yl(2)*0.67, sprintf('Q3: %.1fy', q3_lag), ...
    'FontSize', FONT_SIZE-6);

yl = ylim(axC);
patch(axC, [0 1 1 0],[yl(1) yl(1) yl(2) yl(2)], ...
    [0.2 0.6 0.2],'FaceAlpha',0.12,'EdgeColor','none');
text(axC, 0.05, yl(2)*0.55, ...
    sprintf('Within\n1 year\n(%.0f%%)', 100*mean(abs_lags<=1)), ...
    'FontSize', FONT_SIZE-6, 'Color',[0.1 0.5 0.1]);

xlabel(axC, 'Absolute EEG-visit gap (years)',    'FontSize', FONT_SIZE);
ylabel(axC, 'Proportion of pairs',               'FontSize', FONT_SIZE);
title(axC,  'C. EEG-visit temporal gap distribution', ...
    'FontSize', FONT_SIZE, 'FontWeight','bold');
set(axC, 'FontSize', FONT_SIZE);

%% ----- Panel D: Forest plot -----
axD = nexttile(tl, 4);
hold(axD,'on'); box(axD,'off'); grid(axD,'on');

% Use bootstrap CIs if available, otherwise Laplace
[beta_m, betanames_m, stats_m] = fixedEffects(mdl);
raw_names = string(betanames_m.Name);

if ~isempty(MMR.BootstrapTable4)
    BT       = MMR.BootstrapTable4;
    ci_label = '95% Bootstrap CI';
    OR_D     = nan(numel(raw_names), 1);
    OR_lo_D  = nan(numel(raw_names), 1);
    OR_hi_D  = nan(numel(raw_names), 1);
    for k = 1:numel(raw_names)
        bt_row = BT(string(BT.Term) == raw_names(k), :);
        if ~isempty(bt_row)
            OR_D(k)    = bt_row.OR;
            OR_lo_D(k) = bt_row.OR_CI_lo;
            OR_hi_D(k) = bt_row.OR_CI_hi;
        else
            OR_D(k)    = exp(beta_m(k));
            OR_lo_D(k) = exp(beta_m(k) - 1.96*stats_m.SE(k));
            OR_hi_D(k) = exp(beta_m(k) + 1.96*stats_m.SE(k));
        end
    end
else
    OR_D     = exp(beta_m);
    OR_lo_D  = exp(beta_m - 1.96*stats_m.SE);
    OR_hi_D  = exp(beta_m + 1.96*stats_m.SE);
    ci_label = '95% Laplace CI';
end

% Clean display names
disp_names = raw_names;
disp_names(disp_names=="(Intercept)")                   = "Intercept";
disp_names(disp_names=="LogSpikesPerHour")              = "Log spike rate";
disp_names(disp_names=="AbsLag_years")                  = "Absolute lag (years)";
disp_names(disp_names=="LagDirection")                  = "Direction (after vs before)";
disp_names(disp_names=="EpiType3_cat_Frontal")          = "Frontal vs Temporal";
disp_names(disp_names=="EpiType3_cat_General")          = "Generalized vs Temporal";
disp_names(disp_names=="LogSpikesPerHour:AbsLag_years") = "Log spike rate \times Abs lag";

% Remove intercept
isInt      = (disp_names == "Intercept");
OR_D       = OR_D(~isInt);
OR_lo_D    = OR_lo_D(~isInt);
OR_hi_D    = OR_hi_D(~isInt);
pvals_D    = stats_m.pValue(~isInt);
disp_names = disp_names(~isInt);
nTerms     = numel(OR_D);

plot_order = nTerms:-1:1;

for k = 1:nTerms
    idx = plot_order(k);
    col = [0.1 0.3 0.7];
    if pvals_D(idx) >= 0.05
        col = [0.6 0.6 0.6];
    end
    plot(axD, [OR_lo_D(idx), OR_hi_D(idx)], [k k], '-', ...
        'Color', col, 'LineWidth', 2.5);
    scatter(axD, OR_D(idx), k, 100, col, 'filled');

    p = pvals_D(idx);
    if p < 0.001
        pstr = 'p<0.001';
    elseif p < 0.05
        pstr = sprintf('p=%.3f', p);
    else
        pstr = sprintf('p=%.2f', p);
    end
    text(axD, OR_hi_D(idx)+0.005, k, pstr, ...
        'FontSize', FONT_SIZE-5, 'VerticalAlignment','middle');
end

xline(axD, 1, 'k--', 'LineWidth', 1.5);
set(axD, 'YTick', 1:nTerms, ...
         'YTickLabel', disp_names(plot_order), ...
         'FontSize', FONT_SIZE-4);
xlabel(axD, sprintf('Odds Ratio (%s)', ci_label), 'FontSize', FONT_SIZE);
title(axD,  'D. Mixed effects model', ...
    'FontSize', FONT_SIZE, 'FontWeight','bold');
all_ors = [OR_lo_D; OR_hi_D];
xlim(axD, [max(0.4, prctile(all_ors,2)), min(2.0, prctile(all_ors,98))]);
set(axD, 'FontSize', FONT_SIZE);

%% Save main figure
if strlength(string(outPath_main)) > 0
    if ~exist(fileparts(outPath_main),'dir'), mkdir(fileparts(outPath_main)); end
    exportgraphics(FigMain, outPath_main, 'Resolution', 300);
    fprintf('Saved main figure: %s\n', outPath_main);
end

%% ================================================================
%% FIGURE 2: CROSS-VALIDATED CLASSIFIER
%% ================================================================

fprintf('\nRunning %d-fold cross-validated classifier...\n', nFolds);

T_cv = T(T.AbsLag_years <= 1, :);
fprintf('CV dataset: %d pairs, %d patients\n', ...
    height(T_cv), numel(unique(T_cv.Patient)));

T_cv.isFrontal = double(string(T_cv.EpiType3) == "Frontal");
T_cv.isGeneral = double(string(T_cv.EpiType3) == "General");

patients_cv = unique(T_cv.Patient);
nPats_cv    = numel(patients_cv);
rng(42);
shuf_idx    = randperm(nPats_cv);
patients_cv = patients_cv(shuf_idx);
fold_assign = mod(0:nPats_cv-1, nFolds) + 1;

all_scores = [];
all_labels = [];

fpr_folds     = cell(nFolds,1);
tpr_folds     = cell(nFolds,1);
rec_folds     = cell(nFolds,1);
prec_folds    = cell(nFolds,1);
auc_roc_folds = nan(nFolds,1);
auc_pr_folds  = nan(nFolds,1);

for f = 1:nFolds
    test_pats  = patients_cv(fold_assign == f);
    train_pats = patients_cv(fold_assign ~= f);

    T_train = T_cv(ismember(T_cv.Patient, train_pats), :);
    T_test  = T_cv(ismember(T_cv.Patient, test_pats),  :);

    mdl_cv = fitglm(T_train, ...
        'HasSz_bin ~ LogSpikesPerHour + isFrontal + isGeneral', ...
        'Distribution', 'Binomial', ...
        'Link', 'logit');

    scores = predict(mdl_cv, T_test);

    all_scores = [all_scores; scores];
    all_labels = [all_labels; T_test.HasSz_bin];

    [fpr_f, tpr_f, ~, auc_f] = perfcurve(T_test.HasSz_bin, scores, 1);
    fpr_folds{f}     = fpr_f;
    tpr_folds{f}     = tpr_f;
    auc_roc_folds(f) = auc_f;

    [rec_f, prec_f, ~, auc_pr_f] = perfcurve(T_test.HasSz_bin, scores, 1, ...
        'XCrit','reca','YCrit','prec');
    rec_folds{f}    = rec_f;
    prec_folds{f}   = prec_f;
    auc_pr_folds(f) = auc_pr_f;

    fprintf('  Fold %d: AUC-ROC=%.3f, AUC-PR=%.3f (N_test=%d)\n', ...
        f, auc_f, auc_pr_f, height(T_test));
end

[fpr_all, tpr_all, ~, auc_roc_all] = perfcurve(all_labels, all_scores, 1);
[rec_all, prec_all, ~, auc_pr_all] = perfcurve(all_labels, all_scores, 1, ...
    'XCrit','reca','YCrit','prec');
pr_chance = mean(all_labels == 1);

fprintf('\nOverall AUC-ROC=%.3f (mean folds=%.3f +/- %.3f)\n', ...
    auc_roc_all, mean(auc_roc_folds), std(auc_roc_folds));
fprintf('Overall AUC-PR =%.3f (mean folds=%.3f +/- %.3f)\n', ...
    auc_pr_all, mean(auc_pr_folds), std(auc_pr_folds));
fprintf('Chance AUC-PR  =%.3f\n', pr_chance);

FigCV    = figure('Color','w','Position',[100 100 1000 480]);
tl2      = tiledlayout(FigCV, 1, 2, 'TileSpacing','compact','Padding','loose');
fold_col = [0.75 0.75 0.75];
mean_col = [0.1  0.3  0.7 ];

axR = nexttile(tl2, 1);
hold(axR,'on'); box(axR,'off'); grid(axR,'on');

for f = 1:nFolds
    plot(axR, fpr_folds{f}, tpr_folds{f}, '-', ...
        'Color', fold_col, 'LineWidth', 1.0, 'HandleVisibility','off');
end
plot(axR, fpr_all, tpr_all, '-', ...
    'Color', mean_col, 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Pooled AUC=%.3f', auc_roc_all));
plot(axR, [0 1],[0 1], 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'Chance (AUC=0.50)');
text(axR, 0.55, 0.12, ...
    sprintf('Fold AUCs: %.3f ± %.3f', mean(auc_roc_folds), std(auc_roc_folds)), ...
    'FontSize', FONT_SIZE-6, 'Color', fold_col*0.6);
xlabel(axR, 'False Positive Rate', 'FontSize', FONT_SIZE);
ylabel(axR, 'True Positive Rate',  'FontSize', FONT_SIZE);
title(axR,  'A. ROC Curve',        'FontSize', FONT_SIZE, 'FontWeight','bold');
legend(axR, 'Location','southeast','FontSize', FONT_SIZE-6);
xlim(axR,[0 1]); ylim(axR,[0 1]);
set(axR, 'FontSize', FONT_SIZE);

axP = nexttile(tl2, 2);
hold(axP,'on'); box(axP,'off'); grid(axP,'on');

for f = 1:nFolds
    r  = rec_folds{f};
    p  = prec_folds{f};
    ok = isfinite(r) & isfinite(p);
    plot(axP, r(ok), p(ok), '-', ...
        'Color', fold_col, 'LineWidth', 1.0, 'HandleVisibility','off');
end
ok_all = isfinite(rec_all) & isfinite(prec_all);
plot(axP, rec_all(ok_all), prec_all(ok_all), '-', ...
    'Color', mean_col, 'LineWidth', 2.5, ...
    'DisplayName', sprintf('Pooled AUC=%.3f', auc_pr_all));
yline(axP, pr_chance, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Chance (%.2f)', pr_chance));
text(axP, 0.05, 0.12, ...
    sprintf('Fold AUCs: %.3f ± %.3f', mean(auc_pr_folds), std(auc_pr_folds)), ...
    'FontSize', FONT_SIZE-6, 'Color', fold_col*0.6);
xlabel(axP, 'Recall',    'FontSize', FONT_SIZE);
ylabel(axP, 'Precision', 'FontSize', FONT_SIZE);
title(axP,  'B. Precision-Recall Curve', 'FontSize', FONT_SIZE, 'FontWeight','bold');
legend(axP, 'Location','northeast', 'FontSize', FONT_SIZE-6);
xlim(axP,[0 1]); ylim(axP,[0 1]);
set(axP, 'FontSize', FONT_SIZE);

if strlength(string(outPath_cv)) > 0
    if ~exist(fileparts(outPath_cv),'dir'), mkdir(fileparts(outPath_cv)); end
    exportgraphics(FigCV, outPath_cv, 'Resolution', 300);
    fprintf('Saved CV figure: %s\n', outPath_cv);
end

end

%% ================================================================
%% LOCAL HELPERS
%% ================================================================

function Tpred = build_pred_table(spike_grid, abs_lag, lag_dir, epi_type, canonical3)
n = numel(spike_grid);
Tpred = table(...
    spike_grid, ...
    repmat(abs_lag,  n, 1), ...
    repmat(lag_dir,  n, 1), ...
    repmat(categorical(epi_type, canonical3), n, 1), ...
    categorical(ones(n,1)), ...
    'VariableNames', {'LogSpikesPerHour','AbsLag_years','LagDirection', ...
                      'EpiType3_cat','PatientID'});
Tpred.EpiType3_cat = reordercats(Tpred.EpiType3_cat, canonical3);
end


function write_tableS1(MMR, outPath)
%% write_tableS1
% Writes a CSV of the full mixed effects model results for Table S1
%
% Inputs:
%   MMR     : output of fit_mixed_effects_models
%   outPath : output path for CSV (e.g. '../output/TableS1.csv')

%% ---- Clean up term names for display ----
function disp = clean_term(raw)
    disp = raw;
    disp = strrep(disp, '(Intercept)',                   'Intercept');
    disp = strrep(disp, 'LogSpikesPerHour:AbsLag_years', 'Log spike rate × Absolute lag');
    disp = strrep(disp, 'LogSpikesPerHour:LagDirection', 'Log spike rate × Lag direction');
    disp = strrep(disp, 'LogSpikesPerHour',              'Log spike rate');
    disp = strrep(disp, 'AbsLag_years',                  'Absolute lag (years)');
    disp = strrep(disp, 'LagDirection',                  'Lag direction (after vs before)');
    disp = strrep(disp, 'EpiType3_cat_Frontal',          'Frontal vs Temporal epilepsy');
    disp = strrep(disp, 'EpiType3_cat_General',          'Generalized vs Temporal epilepsy');
end

%% ---- Determine CI source ----
use_bootstrap4 = ~isempty(MMR.BootstrapTable4);
use_bootstrap1 = ~isempty(MMR.BootstrapTable1);

%% ---- Build rows for each model ----
rows = {};

%% Model 4 (Primary logistic)
if ~isempty(MMR.FE_logistic_proxonly)
    FE = MMR.FE_logistic_proxonly;
    if use_bootstrap4
        BT = MMR.BootstrapTable4;
    end

    for i = 1:height(FE)
        term = string(FE.Term(i));
        if term == "(Intercept)", continue; end % skip intercept

        or_pt  = FE.OR(i);
        p_val  = FE.p(i);

        if use_bootstrap4
            bt_row = BT(string(BT.Term) == term, :);
            if ~isempty(bt_row)
                ci_lo = bt_row.OR_CI_lo;
                ci_hi = bt_row.OR_CI_hi;
                ci_src = 'Bootstrap';
            else
                ci_lo = FE.OR_lo(i);
                ci_hi = FE.OR_hi(i);
                ci_src = 'Laplace';
            end
        else
            ci_lo = FE.OR_lo(i);
            ci_hi = FE.OR_hi(i);
            ci_src = 'Laplace';
        end

        % Format p-value
        if p_val < 0.001
            p_str = '<0.001';
        elseif p_val < 0.01
            p_str = sprintf('%.3f', p_val);
        else
            p_str = sprintf('%.2f', p_val);
        end

        rows(end+1,:) = { ...
            'Model 4 (Primary: logistic, HasSz outcome)', ...
            clean_term(char(term)), ...
            sprintf('%.3f', or_pt), ...
            sprintf('%.3f', ci_lo), ...
            sprintf('%.3f', ci_hi), ...
            ci_src, ...
            p_str ...
        }; %#ok<AGROW>
    end
end

%% Model 1 (Full logistic — supplemental)
if ~isempty(MMR.FE_logistic)
    FE = MMR.FE_logistic;
    if use_bootstrap1
        BT = MMR.BootstrapTable1;
    end

    for i = 1:height(FE)
        term = string(FE.Term(i));
        if term == "(Intercept)", continue; end

        or_pt = FE.OR(i);
        p_val = FE.p(i);

        if use_bootstrap1
            bt_row = BT(string(BT.Term) == term, :);
            if ~isempty(bt_row)
                ci_lo = bt_row.OR_CI_lo;
                ci_hi = bt_row.OR_CI_hi;
                ci_src = 'Bootstrap';
            else
                ci_lo = FE.OR_lo(i);
                ci_hi = FE.OR_hi(i);
                ci_src = 'Laplace';
            end
        else
            ci_lo = FE.OR_lo(i);
            ci_hi = FE.OR_hi(i);
            ci_src = 'Laplace';
        end

        if p_val < 0.001
            p_str = '<0.001';
        elseif p_val < 0.01
            p_str = sprintf('%.3f', p_val);
        else
            p_str = sprintf('%.2f', p_val);
        end

        rows(end+1,:) = { ...
            'Model 1 (Full logistic: proximity + direction interactions)', ...
            clean_term(char(term)), ...
            sprintf('%.3f', or_pt), ...
            sprintf('%.3f', ci_lo), ...
            sprintf('%.3f', ci_hi), ...
            ci_src, ...
            p_str ...
        }; %#ok<AGROW>
    end
end

%% Model 2 (Linear sensitivity)
if ~isempty(MMR.FE_linear)
    FE = MMR.FE_linear;

    for i = 1:height(FE)
        term = string(FE.Term(i));
        if term == "(Intercept)", continue; end

        beta_pt = FE.Beta(i);
        p_val   = FE.p(i);
        se      = FE.SE(i);
        ci_lo   = beta_pt - 1.96*se;
        ci_hi   = beta_pt + 1.96*se;

        if p_val < 0.001
            p_str = '<0.001';
        elseif p_val < 0.01
            p_str = sprintf('%.3f', p_val);
        else
            p_str = sprintf('%.2f', p_val);
        end

        rows(end+1,:) = { ...
            'Model 2 (Sensitivity: linear, log seizure frequency outcome)', ...
            clean_term(char(term)), ...
            sprintf('%.3f', beta_pt), ...
            sprintf('%.3f', ci_lo), ...
            sprintf('%.3f', ci_hi), ...
            'Laplace', ...
            p_str ...
        }; %#ok<AGROW>
    end
end

%% ---- Assemble and write table ----
T_out = cell2table(rows, ...
    'VariableNames', { ...
        'Model', ...
        'Term', ...
        'Estimate', ...
        'CI_lower', ...
        'CI_upper', ...
        'CI_method', ...
        'p_value' ...
    });

% Add note about what Estimate means per model
fprintf('Table S1: %d rows written\n', height(T_out));
fprintf('  Models 4 and 1: Estimate = Odds Ratio\n');
fprintf('  Model 2: Estimate = Beta coefficient (log scale)\n');

if ~exist(fileparts(outPath),'dir'), mkdir(fileparts(outPath)); end
writetable(T_out, outPath);
fprintf('Wrote Table S1: %s\n', outPath);

end