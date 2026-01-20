function run_spike_sz_pipeline_clean()
%% ===================== INFORMATION ==================
%{
This contains the code for running all analyses associated with the
manuscript "Automated estimation of seizure and spike burden and their 
association in a large epilepsy cohort" by Conrad et al.

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
2. Matlab (R2024a) and the Statistics and Machine Learning Toolbox
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
spikeSummaryMultiCsv = '../data/spike_counts.csv';
reportCsv            = '../data/clinical_data_deidentified.csv';

% Outputs
fig1_out   = '../output/Fig1.png';
fig2_out   = '../output/Fig2.png';
figS1_out  = '../output/FigS1.png';
figS2_out  = '../output/FigS2.png';
figS3_out  = '../output/FigS3.png';
figS4_out  = '../output/FigS4.png'; % Sex comparison: spike rate in women vs men
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

nBoot = 5000; % bootstrap iterations to get confidence intervals
alpha = 0.05;

% SpikeNet threshold columns
countCol = "count_0_46"; % the probability threshold from the SpikeNet2 paper
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


%% Stats for colin gender grant (not used in the manuscript)
%{
% Sex (first nonmissing)
sex_raw = upper(strtrim(string(ReportTable.nlp_gender)));
[gs, pidS] = findgroups(ReportTable.patient_id);
sex_per = splitapply(@local_first_nonmissing, sex_raw, gs);
SexTableColin = innerjoin(table(unique(ReportTable.patient_id),'VariableNames',{'Patient'}), ...
                     table(pidS,sex_per,'VariableNames',{'Patient','SexCode'}), 'Keys','Patient');
n_f_colin = nnz(SexTableColin.SexCode=="F");
%}

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
Vuniq = build_visit_level_table_R1(ReportTable); % 1 row per visit, with information about seizure frequency

% Patient typing from report 
PatientTypingAll = build_patient_typing_from_report(ReportTable, canonical3); % 1 row per patient, with information about epilepsy type

% restrict to visits within +/- 12 months of an EEG (per patient)
Vuniq_12mo = restrict_visits_to_eeg_window(Vuniq, ReportTable, 365);

% Patient seizure metrics: MeanSzFreq + fraction hasSz==1 among valid visits
SzFreqPerPatient = build_patient_seizure_metrics(Vuniq); % 1 row per patient, with mean seizure frequency across visits

% recompute patient seizure metrics using ONLY +/- 12 month visits
SzFreqPerPatient_12mo = build_patient_seizure_metrics(Vuniq_12mo);


%% ======================= BUILD VIEW BUNDLE =======================
Views = build_filtered_view( ...
    SpikeSummaryTable, ReportTable, PatientTypingAll, SzFreqPerPatient, ...
    NESD_LABEL, badTypes, canonical3); % a structure that contains the patient level data for subsequent analyses


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
        'MeanSzFreq', ' (positive spike rates only)', true);

fprintf('Saved Fig 2:  %s\n', fig2_out);
fprintf('Saved Fig S1: %s\n', figS1_out);

%% ======================= FIG S2 (SZ FREQ BY REPORTED SPIKES ACROSS EEGs) =======================
FigS2 = make_figS2_sz_by_reported_spikes(Views, SzFreqPerPatient, nBoot, alpha, spearman_xLims);
save_fig(FigS2, figS2_out);
fprintf('Saved Fig S2: %s\n', figS2_out);


%% ======================= FIG S3 (NEAR VS FAR CORRELATION (PAIRED BOOTSTRAP)) ====================
plot_delta_rho_histogram( ...
    Views, Vuniq, Views.ReportForKeptSessions, ...
    0.333, 0.667, nBoot, alpha, ...
    figS3_out);

%% ======================= FIG S4 (SPIKE RATE BY SEX, not used in manuscript) =======================
%{
[FigS4, FigS4Stats] = make_figS4_spike_by_sex(Views, EPS_RATE, Y_ZERO, Y_LIMS, TITLE_Y_OFFSET, nBoot, alpha);
save_fig(FigS4, figS4_out);
fprintf('Saved Fig S4: %s\n', figS4_out);

% Optional: print a concise stat line
fprintf('[Fig S4] ranksum women vs men: %s\n', p_label(FigS4Stats.p_rankSum));
%}

%% ======================= TABLE 1 =======================
Table1_flat = build_table1_flat(Views, SzFreqPerPatient, Vuniq, EPS_RATE, nBoot, alpha);
if ~exist(fileparts(Table1Csv),'dir'), mkdir(fileparts(Table1Csv)); end
writetable(Table1_flat, Table1Csv);
fprintf('Wrote Table 1 CSV: %s\n', Table1Csv);

%% ======================= WRITE HTML SUMMARY =======================
write_results_html(resultsHtml, Views, SzFreqPerPatient, ...
    Fig1Stats, ...
    SpearmanResults_main, rs_all_main, p_all_main, n_all_main, rho_lo_main, rho_hi_main, subtype_ci_main, ...
    Views.ReportForKeptSessions);

fprintf('Wrote HTML summary: %s\n', resultsHtml);



end

%% =====================================================================
%% ======================= CORE PIPELINE HELPERS =======================
%% =====================================================================

function Vuniq_f = restrict_visits_to_eeg_window(Vuniq, ReportForKeptSessions, windowDays)
% Keep only clinic visits that are within +/- windowDays of ANY EEG for that patient.
%
% Inputs:
%   Vuniq: table with columns Patient (double), VisitDate (datetime), Freq_R1 (double), ...
%   ReportForKeptSessions: table with columns patient_id and start_time_deid (or Patient + start_time_deid)
%   windowDays: scalar (e.g., 183 for ~6 months)
%
% Output:
%   Vuniq_f: subset of Vuniq

% ---- strict checks ----
needV = ["Patient","VisitDate","Freq_R1"];
missV = setdiff(needV, string(Vuniq.Properties.VariableNames));
if ~isempty(missV)
    error("Vuniq missing required columns: %s", strjoin(missV,", "));
end

% Accept either (patient_id) or (Patient) in report
if ismember("Patient", string(ReportForKeptSessions.Properties.VariableNames))
    pid = double(ReportForKeptSessions.Patient);
elseif ismember("patient_id", string(ReportForKeptSessions.Properties.VariableNames))
    pid = double(ReportForKeptSessions.patient_id);
else
    error("ReportForKeptSessions must contain Patient or patient_id.");
end

if ~ismember("start_time_deid", string(ReportForKeptSessions.Properties.VariableNames))
    error("ReportForKeptSessions missing start_time_deid.");
end

% ---- parse EEG datetimes ----
EEG_raw = ReportForKeptSessions.start_time_deid;
if isdatetime(EEG_raw)
    EEG_dt = EEG_raw;
else
    % handle strings like 2000-01-01T12:34:56
    EEG_dt = datetime(strtrim(string(EEG_raw)), 'InputFormat',"yyyy-MM-dd'T'HH:mm:ss");
end

okEEG = ~isnat(EEG_dt) & isfinite(pid);
EEG_tbl = table(pid(okEEG), EEG_dt(okEEG), 'VariableNames', {'Patient','EEG_Date'});
if isempty(EEG_tbl)
    error("No valid EEG dates found (start_time_deid could not be parsed).");
end

% ---- compute per-visit min distance to ANY EEG for that patient ----
Vuniq_f = Vuniq;
keep = false(height(Vuniq_f),1);

[gv, pidV] = findgroups(Vuniq_f.Patient);

for k = 1:numel(pidV)
    p = pidV(k);

    idxV = find(gv == k);
    vDates = Vuniq_f.VisitDate(idxV);

    eegDates = EEG_tbl.EEG_Date(EEG_tbl.Patient == p);
    if isempty(eegDates)
        % no EEG dates -> drop all visits for this patient
        continue
    end

    % min abs gap to any EEG for each visit
    minAbsGap = inf(numel(vDates),1);
    for j = 1:numel(vDates)
        minAbsGap(j) = min(abs(days(eegDates - vDates(j))));
    end

    keep(idxV) = (minAbsGap <= windowDays);
end

n_before = height(Vuniq);
n_after  = nnz(keep);
fprintf('[Visit-EEG window] Kept %d/%d visits (%.1f%%) within +/- %d days of an EEG\n', ...
    n_after, n_before, 100*n_after/max(1,n_before), windowDays);

Vuniq_f = Vuniq_f(keep,:);
end


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
PatientLevelSpikeRates.EpiType4 = TypingFiltered.EpiType3;


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

% ===================== DEFINE MAIN ANALYSIS COHORT =====================
% Epilepsy + documented seizure frequency
CohortPatients = SzFreqEpilepsy.Patient(isfinite(SzFreqEpilepsy.MeanSzFreq));
CohortPatients = unique(CohortPatients);
CohortTable = table(CohortPatients, 'VariableNames', {'Patient'});

% ===================== RESTRICT ALL DOWNSTREAM TABLES ==================
% Patient-level tables
PatientLevelSpikeRates = innerjoin(PatientLevelSpikeRates, CohortTable, 'Keys','Patient');
TypingFiltered         = innerjoin(TypingFiltered,         CohortTable, 'Keys','Patient');

% Session-level + report tables (Fig 1A, Fig S2)
SessionsFiltered      = innerjoin(SessionsFiltered,      CohortTable, 'Keys','Patient');
ReportForKeptSessions = innerjoin(ReportForKeptSessions, CohortTable, 'Keys','Patient');

% Seizure-frequency table (keep consistent)
SzFreqEpilepsy = innerjoin(SzFreqEpilepsy, CohortTable, 'Keys','Patient');

% Recompute epilepsy masks AFTER restriction
etype_norm = lower(strtrim(string(PatientLevelSpikeRates.EpilepsyType)));
IsNESDMask = (etype_norm == lower(strtrim(NESD_LABEL)));
IsBadType  = ismember(etype_norm, badTypes) | ismissing(etype_norm) | (strlength(etype_norm)==0);
IsEpilepsyMask = ~IsNESDMask & ~IsBadType;

fprintf('[Cohort restriction] Using %d epilepsy patients with documented seizure frequency\n', ...
    numel(CohortPatients));
% end of new cohort restriction


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

% Use the SAME subtype definition as Fig 2 (EpiType3)
InCanonical3Mask = ismember(TypingFiltered.EpiType3, canonical3);

Canonical3_SubsetTable = innerjoin( ...
    TypingFiltered(InCanonical3Mask, {'Patient','EpiType3'}), ...
    PatientLevelSpikeRates(:, {'Patient','MeanSpikeRate_perHour'}), ...
    'Keys','Patient');

Canonical3_SubsetTable.Properties.VariableNames{'EpiType3'} = 'EpiType4';

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

% Panel B: subtype
%Sub = Views.Canonical3_SubsetTable;
% Use the same typed cohort as Fig 2
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
    ReportForKeptSessions)

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
cohortPatients = Views.PatientLevelSpikeRates.Patient;  % already cohort-restricted
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

%% Table 1 general overview
fprintf(fid, '<h2>Cohort summary</h2>\n');
%{
fprintf(fid,['<p>We included %d patients with at least one outpatient routine EEG '...
    '(%d EEGs). %d patients (%1.1f%%) carried a '...
    'diagnosis of epilepsy. Median [95%% CI] monthly seizure frequency was '...
    '%1.2f [%1.2f-%1.2f], and median [95%% CI] spikes/hour across EEGs '...
    'was %1.2f [%1.2f-%1.2f] (Table 1).</p>'],...
    N_total, n_eegs_all, n_epi, 100*n_epi/max(1,N_total), ...
    sf_med, sf_ci_lo, sf_ci_hi, ...
    sr_med, sr_ci_lo, sr_ci_hi);
%}
fprintf(fid,['<p>We included %d patients '...
    '(%d EEGs). Median [95%% CI] monthly seizure frequency was '...
    '%1.2f [%1.2f-%1.2f], and median [95%% CI] spikes/hour across EEGs '...
    'was %1.2f [%1.2f-%1.2f] (Table 1).</p>'],...
    N_total, n_eegs_all, ...
    sf_med, sf_ci_lo, sf_ci_hi, ...
    sr_med, sr_ci_lo, sr_ci_hi);

%% ---- Figure 1: Control Panels ----
fprintf(fid, '<h2>Spike rates by patient groups</h2>\n');

% Panel A
pA_str = format_p_html(Fig1Stats.p_rankSum_A);
fprintf(fid,['<p>Automatically detected spike rates were higher in EEGs with clinically-reported spikes '...
            '(median spike rate %.2f [95%% CI %.2f-%.2f] spikes/hour) than '...
            'in EEGs without reported spikes (%.2f [%.2f-%.2f] '...
            'spikes/hour) (%s, Cliff''s &delta; = %.2f; Fig. 1A), '...
            'suggesting a low false-positive detection rate. '],...
            Fig1Stats.m_pre, Fig1Stats.lo_pre, Fig1Stats.hi_pre, ...
            Fig1Stats.m_abs, Fig1Stats.lo_abs, Fig1Stats.hi_abs, ...
            pA_str, Fig1Stats.effectA_cliff);

% Panel B omnibus
fprintf(fid,['Spike rates differed across epilepsy subtypes (Kruskal–Wallis '...
            '%s, η² ≈ %.3f), '], format_p_html(Fig1Stats.p_kw_C), Fig1Stats.eta2_kw_C);

% Panel B pairwise narrative (assumes General, Temporal, Frontal exist)
Tsub = Fig1Stats.SubtypeStatsTable;
getRow = @(nm) Tsub(strcmp(string(Tsub.Group), string(nm)), :);

rG = getRow("General");
rT = getRow("Temporal");
rF = getRow("Frontal");

p_pair_bonf = Fig1Stats.p_pair_bonf; % [G vs T; G vs F; T vs F]
fprintf(fid,['with higher rates in generalized epilepsy than '...
    'temporal or frontal lobe epilepsy (Fig. 1B).'])

fprintf(['Generalized epilepsy demonstrated higher spike rates '...
            '(%.2f [%.2f-%.2f]) than temporal '...
            '(%.2f [%.2f-%.2f]; Bonferroni-adjusted '...
            '%s) and frontal '...
            '(%.2f [%.2f-%.2f]; '...
            '%s). The temporal versus frontal comparison '...
            'was not significant (%s; Fig. 1B).</p>'],...
            rG.Median, rG.CI_lo, rG.CI_hi, ...
            rT.Median, rT.CI_lo, rT.CI_hi, format_p_html(p_pair_bonf(1)), ...
            rF.Median, rF.CI_lo, rF.CI_hi, format_p_html(p_pair_bonf(2)), ...
            format_p_html(p_pair_bonf(3)));
%}

%% ---- Figure 2: Spearman correlations ----
fprintf(fid, '<h2>Relationship between spike rate and seizure frequency</h2>\n');

fprintf(fid,['<p>Spike rate '...
    'and seizure frequency were positively correlated '...
    '(N = %d, &rho; = %.2f [95%% CI %.2f-%.2f], %s).'], ...
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

fprintf(fid,['Results were similar when restricting analyses to patients with detectable spikes, '...
    'although the correlation for frontal lobe epilepsy was larger and significant in this secondary analysis (Fig. S1). ']);
fprintf(fid,['Patients with spikes clinically-reported on at least one EEG also had '...
    'higher mean seizure frequencies than those without spikes (Fig. S2). ']);
fprintf(fid,['Spike-seizure correlations were stronger when clinic '...
    'visits occurred closer in time to EEG acquisition (Fig. S3), '...
    'suggesting that spike rates may track seizure burden over time.</p>'])
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
%% OLD FUNCTION, UNDERPOWERED ANALYSIS, NOT USING ANYMORE
%{
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
%}

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
    %mask_all = mask_all & (x_all > 0) & (y_all > 0); % require non zero
    %spikes and non zero szs
    mask_all = mask_all & (x_all > 0); % only require non-zero spikes, allow zero szs
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
        %mask = mask & (x > 0) & (y > 0); % require non zero
    %spikes and non zero szs
        mask = mask & (x > 0); % only require non-zero spikes, allow zero szs
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
        %{
        idx = idx & (PatientSpikeSz_Typed.MeanSpikeRate_perHour > 0) & ...
                  (PatientSpikeSz_Typed.(freqFieldName) > 0); % require non zero
    %spikes and non zero szs
        %}

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

    % annotation text from your computed tables
    row = SpearmanResults(strcmp(string(SpearmanResults.Group), gStr), :);
    rowCI = subtype_ci(subtype_ci.Group == gStr, :);

    txt = sprintf('\\rho=%.2f [%.2f-%.2f], p_{bonf}%s', ...
    row.Spearman_r, rowCI.ci_lo, rowCI.ci_hi, ...
    string(regexprep(char(p_label(row.p_bonf)), '^p', '')));

    nNow = nnz(idx);
    
    title(ax, sprintf('%s%s (N=%d)', panelTitle{p}, labelSuffix, nNow), ...
        'FontSize', fontL, 'FontWeight','bold');
    %{
    if isempty(labelSuffix)
        titleText = { ...
            panelTitle{p}, ...
            sprintf('(N=%d)', nNow) ...
        };
    else
        titleText = { ...
            panelTitle{p}, ...
            labelSuffix, ...
            sprintf('(N=%d)', nNow) ...
        };
    end

    title(ax, titleText, 'FontSize', fontL, 'FontWeight','bold');
    %}
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
%{
function add_sigbar(ax, x1, x2, y, ptext)
tick = 0.03 * diff(ax.YLim);
plot(ax, [x1 x1 x2 x2], [y-tick, y, y, y-tick], 'k-', 'LineWidth', 1.3);
text(ax, mean([x1 x2]), y +0.001*diff(ax.YLim), ptext, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);
end
%}

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



function GapTable = compute_visit_eeg_gaps(Vuniq, ReportForKeptSessions)
% Returns one row per visit with min |visit - EEG| in days

% ---- checks ----
require_cols(Vuniq, ["Patient","VisitDate"], "Vuniq");

if ismember("Patient", string(ReportForKeptSessions.Properties.VariableNames))
    pid = double(ReportForKeptSessions.Patient);
elseif ismember("patient_id", string(ReportForKeptSessions.Properties.VariableNames))
    pid = double(ReportForKeptSessions.patient_id);
else
    error("ReportForKeptSessions must have Patient or patient_id");
end

require_cols(ReportForKeptSessions, ["start_time_deid"], "ReportForKeptSessions");

% ---- EEG datetimes ----
EEG_raw = ReportForKeptSessions.start_time_deid;
if isdatetime(EEG_raw)
    EEG_dt = EEG_raw;
else
    EEG_dt = datetime(strtrim(string(EEG_raw)), ...
        'InputFormat',"yyyy-MM-dd'T'HH:mm:ss");
end

EEG_tbl = table(pid, EEG_dt, 'VariableNames', {'Patient','EEG_Date'});
EEG_tbl = EEG_tbl(~isnat(EEG_tbl.EEG_Date), :);

% ---- compute per-visit min gap ----
GapTable = Vuniq(:, {'Patient','VisitDate'});
GapTable.MinAbsGap_days = NaN(height(GapTable),1);

[gv, pidV] = findgroups(GapTable.Patient);

for k = 1:numel(pidV)
    p = pidV(k);
    idxV = find(gv == k);

    eegDates = EEG_tbl.EEG_Date(EEG_tbl.Patient == p);
    if isempty(eegDates)
        continue
    end

    for j = idxV'
        GapTable.MinAbsGap_days(j) = ...
            min(abs(days(eegDates - GapTable.VisitDate(j))));
    end
end

GapTable.MinAbsGap_years = GapTable.MinAbsGap_days / 365.25;
end




function V_f = restrict_visits_by_min_abs_gap(Vuniq, ReportForKeptSessions, minDays, maxDays)
% Keep visits for which min |VisitDate - EEG_Date| across that patient's EEGs is in [minDays, maxDays].
%
% minDays can be 0, maxDays can be Inf

needV = ["Patient","VisitDate","Freq_R1"];
missV = setdiff(needV, string(Vuniq.Properties.VariableNames));
if ~isempty(missV)
    error("Vuniq missing required columns: %s", strjoin(missV,", "));
end

% Accept either Patient or patient_id
if ismember("Patient", string(ReportForKeptSessions.Properties.VariableNames))
    pid = double(ReportForKeptSessions.Patient);
elseif ismember("patient_id", string(ReportForKeptSessions.Properties.VariableNames))
    pid = double(ReportForKeptSessions.patient_id);
else
    error("ReportForKeptSessions must contain Patient or patient_id.");
end

if ~ismember("start_time_deid", string(ReportForKeptSessions.Properties.VariableNames))
    error("ReportForKeptSessions missing start_time_deid.");
end

% Parse EEG datetimes
EEG_raw = ReportForKeptSessions.start_time_deid;
if isdatetime(EEG_raw)
    EEG_dt = EEG_raw;
else
    EEG_dt = datetime(strtrim(string(EEG_raw)), 'InputFormat',"yyyy-MM-dd'T'HH:mm:ss");
end

okEEG = ~isnat(EEG_dt) & isfinite(pid);
EEG_tbl = table(pid(okEEG), EEG_dt(okEEG), 'VariableNames', {'Patient','EEG_Date'});
if isempty(EEG_tbl)
    error("No valid EEG dates found (start_time_deid could not be parsed).");
end

V_f = Vuniq;
keep = false(height(V_f),1);

[gv, pidV] = findgroups(V_f.Patient);

for k = 1:numel(pidV)
    p = pidV(k);
    idxV = find(gv == k);
    vDates = V_f.VisitDate(idxV);

    eegDates = EEG_tbl.EEG_Date(EEG_tbl.Patient == p);
    if isempty(eegDates)
        continue
    end

    minAbsGap = inf(numel(vDates),1);
    for j = 1:numel(vDates)
        minAbsGap(j) = min(abs(days(eegDates - vDates(j))));
    end

    keep(idxV) = (minAbsGap >= minDays) & (minAbsGap <= maxDays);
end

fprintf('[Visit-EEG distance] Kept %d/%d visits with min|gap| in [%g, %g] days\n', ...
    nnz(keep), height(Vuniq), minDays, maxDays);

V_f = V_f(keep,:);
end

function NearFarStats = plot_delta_rho_histogram(Views, Vuniq, ReportForKeptSessions, ...
                                                 nearQ, farQ, nBoot, alpha, outPng)
% plot_delta_rho_histogram
% Plot (top) distribution of min Visit–EEG gaps with tercile shading + quantile cutoffs,
% and (bottom) bootstrap distribution of Δρ = ρ_near − ρ_far (Spearman),
% with Δρ x-axis centered around 0.
%
% Inputs:
%   nearQ, farQ : quantiles in [0,1], e.g., 0.333 and 0.667
%   nBoot       : bootstrap iterations (e.g., 5000)
%   alpha       : CI alpha (e.g., 0.05)
%   outPng      : output file path ("" or '' to skip saving)
%
% NOTE: This function DEFINES near/far by quantiles of MinAbsGap_days computed across
%       visits in Vuniq.

if nargin < 8, outPng = ""; end
if nargin < 4 || isempty(nearQ), nearQ = 0.333; end
if nargin < 5 || isempty(farQ),  farQ  = 0.667; end
if nargin < 6 || isempty(nBoot), nBoot = 5000; end
if nargin < 7 || isempty(alpha), alpha = 0.05; end

% ------------------ Base cohort patients (match main Spearman cohort) ------------------
if isfield(Views, 'PatientSpikeSz_All') && ~isempty(Views.PatientSpikeSz_All)
    basePatients = unique(double(Views.PatientSpikeSz_All.Patient));
else
    basePatients = unique(double(Views.PatientLevelSpikeRates.Patient));
end

SpikeTbl = Views.PatientLevelSpikeRates(:, {'Patient','MeanSpikeRate_perHour'});
SpikeTbl.Patient = double(SpikeTbl.Patient);
SpikeTbl = innerjoin(SpikeTbl, table(basePatients,'VariableNames',{'Patient'}), 'Keys','Patient');

% ------------------ Compute min gaps for ALL visits ------------------
GapTable = compute_visit_eeg_gaps(Vuniq, ReportForKeptSessions);
gaps = double(GapTable.MinAbsGap_days);
gaps = gaps(isfinite(gaps));

if isempty(gaps)
    error("No finite MinAbsGap_days found.");
end

nearDays = quantile(gaps, nearQ);
farDays  = quantile(gaps, farQ);

% ------------------ Build NEAR/FAR visit subsets ------------------
V_near = restrict_visits_by_min_abs_gap(Vuniq, ReportForKeptSessions, 0, nearDays);
V_far  = restrict_visits_by_min_abs_gap(Vuniq, ReportForKeptSessions, farDays, Inf);

Sz_near = build_patient_seizure_metrics(V_near);  % Patient, MeanSzFreq, ...
Sz_far  = build_patient_seizure_metrics(V_far);

% Rename BEFORE join so names are guaranteed
Sz_near2 = Sz_near(:, {'Patient','MeanSzFreq'});
Sz_far2  = Sz_far(:,  {'Patient','MeanSzFreq'});

Sz_near2 = renamevars(Sz_near2, "MeanSzFreq", "MeanSzFreq_near");
Sz_far2  = renamevars(Sz_far2,  "MeanSzFreq", "MeanSzFreq_far");

% Join to spikes, require both near and far seizure metrics
J_near = innerjoin(SpikeTbl, Sz_near2, 'Keys','Patient');
J_far  = innerjoin(SpikeTbl, Sz_far2,  'Keys','Patient');
J = innerjoin(J_near, J_far(:,{'Patient','MeanSzFreq_far'}), 'Keys','Patient');

% Finite mask
mask = isfinite(J.MeanSpikeRate_perHour) & isfinite(J.MeanSzFreq_near) & isfinite(J.MeanSzFreq_far);
J = J(mask,:);

n = height(J);
if n < 3
    error("Not enough patients with BOTH near and far seizure metrics (n=%d).", n);
end

x  = double(J.MeanSpikeRate_perHour);
yN = double(J.MeanSzFreq_near);
yF = double(J.MeanSzFreq_far);

% ------------------ Observed correlations ------------------
rho_near = corr(x, yN, 'Type','Spearman', 'Rows','complete');
rho_far  = corr(x, yF, 'Type','Spearman', 'Rows','complete');
delta_obs = rho_near - rho_far;

% ------------------ Patient-level bootstrap ------------------
delta = nan(nBoot,1);
for b = 1:nBoot
    idx = randi(n, n, 1);
    rn = corr(x(idx), yN(idx), 'Type','Spearman', 'Rows','complete');
    rf = corr(x(idx), yF(idx), 'Type','Spearman', 'Rows','complete');
    delta(b) = rn - rf;
end

ci_lo = prctile(delta, 100*(alpha/2));
ci_hi = prctile(delta, 100*(1-alpha/2));
delta_med = median(delta,'omitnan');

p_one = mean(delta <= 0);                           % H1: near > far
p_two = 2*min(mean(delta <= 0), mean(delta >= 0));  % two-sided bootstrap p

% ------------------ Plot: gaps (top) + delta rho (bottom) ------------------
fig = figure('Color','w','Position',[120 80 950 780]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

% ===================== TOP: Gap distribution + tercile shading =====================
ax1 = nexttile(tl,1); hold(ax1,'on'); box(ax1,'off'); grid(ax1,'on');

% Histogram first so y-lims are defined
hG = histogram(ax1, gaps, 60, 'EdgeColor','none');

xlabel(ax1, '|Visit - EEG| gap (days)','fontsize',20);
ylabel(ax1, 'Visit count','fontsize',20);

% Determine x-range for shading
xL = min(hG.BinEdges);
xU = max(hG.BinEdges);

% Clamp cutoffs to plotting range
nearX = max(min(nearDays, xU), xL);
farX  = max(min(farDays,  xU), xL);

% Use current y-limits
yl = get(ax1,'YLim');
yl_new = [yl(1) yl(1) + 1.3*(yl(2)-yl(1))];
set(ax1,'YLim', yl_new);

% Helper for shaded region [xa, xb]
makeShade = @(xa, xb, a) patch(ax1, ...
    [xa xb xb xa], [yl(1) yl(1) yl(2) yl(2)], ...
    [0 0 0], 'FaceAlpha', a, 'EdgeColor','none');

% Shade: Near / Middle / Far (slightly stronger on near/far)
pNear   = makeShade(xL,   nearX, 0.12);
pMiddle = makeShade(nearX, farX, 0.06);
pFar    = makeShade(farX,  xU,   0.12);

% Put shading behind bars
uistack([pNear pMiddle pFar], 'bottom');
uistack(hG, 'top');

% Cutoff lines
xline(ax1, nearDays, 'k--', 'LineWidth',2);
xline(ax1, farDays,  'k--', 'LineWidth',2);

% Title 
title(ax1, 'A. Visit–EEG gap distribution with lower and upper third cutoffs', ...
    'FontSize', 20, 'Interpreter','none');

% Region labels 
yl = get(ax1,'YLim');
yText = yl(2) * 0.99;

text(ax1, mean([xL nearX]), yText, sprintf("Short gap\n(lower third)"), ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', ...
    'FontSize',20, 'Interpreter','tex');

text(ax1, mean([farX xU]), yText, sprintf("Long gap\n(upper third)"), ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', ...
    'FontSize',20, 'Interpreter','tex');


% ===================== BOTTOM: Δρ distribution (centered at 0) =====================
ax2 = nexttile(tl,2); hold(ax2,'on'); box(ax2,'off'); grid(ax2,'on');

histogram(ax2, delta, 40, 'EdgeColor','none');
xline(ax2, 0,         'k--', 'LineWidth',2);
xline(ax2, delta_med, 'k-',  'LineWidth',2);
%xline(ax2, ci_lo,     'k:',  'LineWidth',2);
%xline(ax2, ci_hi,     'k:',  'LineWidth',2);

% symmetric x-lims around 0 (robust even if delta all positive)
maxAbs = max(abs(delta));
if ~isfinite(maxAbs) || maxAbs==0, maxAbs = 1e-3; end
pad = 0.08 * maxAbs;
xlim(ax2, [-maxAbs-pad, maxAbs+pad]);

xlabel(ax2, '\Delta\rho = \rho_{short gap} - \rho_{long gap}');
ylabel(ax2, 'Bootstrap count');

% Title text WITH interpreter explicitly set to avoid warnings
t2 = sprintf(['B. Distribution of differences in spike-seizure correlation\nbetween short and long visit-EEG gaps\n' ...
              '95%% CI [%.3f, %.3f], p = %.3g'], ...
              ci_lo, ci_hi, p_one);

title(ax2, t2, 'FontSize', 20, 'Interpreter','tex');

set([ax1 ax2],'FontSize',20);

if strlength(string(outPng)) > 0
    if ~exist(fileparts(outPng),'dir'), mkdir(fileparts(outPng)); end
    exportgraphics(fig, outPng, 'Resolution', 300);
end

fprintf(['\nFig S3 analysis:\n'...
    'N patients: %d\n'...
    'Median rho short gap: %1.2f\n'...
    'Median rho long gap: %1.2f\n'...
    'Median [95%% CI] difference in rho: %1.3f [%1.2f-%1.2f]\n'...
    'p = %1.4f.\n'],...
    n,rho_near,rho_far,delta_obs,ci_lo,ci_hi,p_one);

% ------------------ Bundle output ------------------
NearFarStats = struct();
NearFarStats.nPatients = n;

NearFarStats.nearQ = nearQ;
NearFarStats.farQ  = farQ;
NearFarStats.nearDays = nearDays;
NearFarStats.farDays  = farDays;

NearFarStats.rho_near = rho_near;
NearFarStats.rho_far  = rho_far;
NearFarStats.delta_obs = delta_obs;

NearFarStats.delta_boot = delta;
NearFarStats.delta_median = delta_med;
NearFarStats.delta_ci_lo = ci_lo;
NearFarStats.delta_ci_hi = ci_hi;

NearFarStats.p_one_sided = p_one;
NearFarStats.p_two_sided = p_two;

NearFarStats.tableUsed = J;
NearFarStats.gapsUsed = gaps;

end

function [fS4, SexStats] = make_figS4_spike_by_sex(Views, EPS_RATE, Y_ZERO, Y_LIMS, TITLE_Y_OFFSET, nBoot, alpha)
% make_figS4_spike_by_sex
% Plot patient-level mean spike rate (spikes/hour) for women vs men WITH EPILEPSY
% and compare with Wilcoxon rank-sum test.

PL = Views.PatientLevelSpikeRates;

% Restrict to epilepsy patients using the mask aligned to PL rows
if ~isfield(Views, 'IsEpilepsyMask')
    error('Views is missing IsEpilepsyMask.');
end
PL_epi = PL(Views.IsEpilepsyMask, :);

% Pull sex per patient from report (first nonmissing nlp_gender)
Rk = Views.ReportForKeptSessions;
require_cols(Rk, ["Patient","nlp_gender"], "ReportForKeptSessions");

sex_raw = upper(strtrim(string(Rk.nlp_gender)));
[gs, pidS] = findgroups(double(Rk.Patient));
sex_per = splitapply(@local_first_nonmissing, sex_raw, gs);

SexPerPatient = table(double(pidS), sex_per, 'VariableNames', {'Patient','SexCode'});

% Join to epilepsy cohort
T = innerjoin(PL_epi(:, {'Patient','MeanSpikeRate_perHour'}), SexPerPatient, 'Keys','Patient');

% Keep only M/F
isF = (T.SexCode == "F");
isM = (T.SexCode == "M");
T = T(isF | isM, :);

xF = T.MeanSpikeRate_perHour(T.SexCode=="F");
xM = T.MeanSpikeRate_perHour(T.SexCode=="M");

xF = xF(isfinite(xF));
xM = xM(isfinite(xM));

nF = numel(xF);
nM = numel(xM);

% Rank-sum test (if enough data)
if nF >= 3 && nM >= 3
    p = ranksum(xF, xM, 'method','approx');
else
    p = NaN;
end

% Effect size (Cliff's delta; optional but often useful)
d = cliff_delta(xF, xM);

% Bootstrap medians + CI
[medF, loF, hiF] = bootstrap_median_ci(xF, nBoot, alpha);
[medM, loM, hiM] = bootstrap_median_ci(xM, nBoot, alpha);

% Log plotting values (same convention as Fig 1)
Y = to_log10_per_hour([xF; xM], EPS_RATE);
Y = add_y_jitter_eps(Y, Y_ZERO, Y_LIMS, 0.02);
G = categorical( ...
    [repmat("Women", nF, 1); repmat("Men", nM, 1)], ...
    ["Men","Women"], ...      % <-- explicit axis order
    'Ordinal', true);

% ---- draw ----
fS4 = figure('Color','w','Position',[110 110 800 520]);
ax = axes(fS4); hold(ax,'on'); box(ax,'off'); grid(ax,'on');

boxchart(ax, categorical(G), Y, 'BoxFaceAlpha',0.25,'MarkerStyle','none');
swarmchart(ax, categorical(G), Y, 18, 'filled','MarkerFaceAlpha',0.18);

yline(ax, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(ax, Y_LIMS);
ylabel(ax,'Spikes/hour (log scale)');
set_log10_ticks(ax,'y',EPS_RATE,Y_LIMS);

% Median + CI overlays (men at x=1, wommen at x=2)
add_median_ci_overlay(ax, 1, medM, loM, hiM, EPS_RATE);
add_median_ci_overlay(ax, 2, medF, loF, hiF, EPS_RATE);

% Significance bar
yl = ylim(ax);
yMaxData = max(Y(isfinite(Y)));
yBar = yMaxData + 0.06*range(yl);
yNeedTop = yBar + 0.10*range(yl);
if yNeedTop > yl(2)
    ylim(ax, [yl(1) yNeedTop]);
    yl = ylim(ax);
    yBar = yMaxData + 0.06*range(yl);
end
add_sigbar(ax, 1, 2, yBar, p_label(p));

t = title(ax, 'Spike rate by sex (epilepsy only)');
t.Units='normalized';
t.Position(2) = t.Position(2) + TITLE_Y_OFFSET;

set(ax,'FontSize',20);

% Add Ns in labels
labels = string(ax.XTickLabel);
labels(labels=="Women") = sprintf('Women (N=%d)', nF);
labels(labels=="Men")   = sprintf('Men (N=%d)', nM);
ax.XTickLabel = labels;
ax.XTickLabelRotation = 20;

% ---- stats bundle ----
SexStats = struct();
SexStats.nWomen = nF;
SexStats.nMen   = nM;
SexStats.p_rankSum = p;
SexStats.cliffs_delta = d;
SexStats.medWomen = medF; SexStats.loWomen = loF; SexStats.hiWomen = hiF;
SexStats.medMen   = medM; SexStats.loMen   = loM; SexStats.hiMen   = hiM;

% Optional command-window printout (nice for legends)
fprintf(['\nSex analysis (epilepsy only):\n' ...
    'Women: median [95%% CI] = %.2f [%.2f-%.2f] spikes/hour (N=%d)\n' ...
    'Men:   median [95%% CI] = %.2f [%.2f-%.2f] spikes/hour (N=%d)\n' ...
    '%s, Cliff''s δ = %.2f\n'], ...
    medF, loF, hiF, nF, ...
    medM, loM, hiM, nM, ...
    p_label(p), d);

end
