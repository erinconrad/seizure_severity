function paired_spikerate_by_hasSz(which_runs)
% paired_spikerate_by_hasSz(which_runs)
%   0 (default) = whole file
%   1           = first run (~1h)
%   2           = first 24 runs (~24h)
%
% Purpose:
%   For each patient, compare the mean spike rate across all matched visits
%   with HasSz==0 vs HasSz==1 (paired per patient), for ALL valid epilepsy
%   patients (no subtype restriction in the analysis).
%
% Matching:
%   - Each EEG is mapped to the nearest clinic visit within asymmetric windows:
%       -preGapDays <= (EEG_Date - Visit_Date) <= +postGapDays
%   - Among eligible visits, the closest (by |gap|) is chosen.
%   - Per patient, paired test compares the mean of EEGs matched to HasSz==0
%     vs the mean of EEGs matched to HasSz==1.
%
% Outpatient definition (when only_outpt==1):
%   outpatient if and only if:
%     (1) acquired_on contains "spe" or "radnor" (case-insensitive)
%         OR
%     (2) report_PATIENT_CLASS == "Outpatient"
%         OR
%     (3) jay_in_or_out == "out"
%
% Outputs:
%   - Single spaghetti plot of per-patient paired means (overall only),
%     plotted as log10(spike rate in spikes/hour).
%   - Wilcoxon sign-rank on RAW spike rates (spikes/hour).
%   - Console summary.
%
% Erin 2025-11-21

if nargin < 1, which_runs = 0; end

%% ======================= CONFIG =======================
% Paths
outCsv     = '../data/SN_counts/spike_counts_summary.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-19_1356.csv';
pairsOut   = '../output/eeg_visit_pairs_pairedmeans.csv';
saveAudit  = true;   % write the matched-pairs audit CSV

% Filters
only_amb   = 2;   % 1 = ONLY ambulatory (Duration_sec >= 12h); 2 = ONLY routine (<=12h); 0 = all
only_outpt = 1;   % 1 = restrict to outpatient EEGs using outpatient rule (see header)

% Asymmetric windows (days) for matching EEG to visit
preGapDays  = 365;   % EEG may be up to this many days BEFORE the visit
postGapDays = 30;    % EEG may be up to this many days AFTER  the visit

% Patient-level exclusion list (applies to EpilepsyType string)
badTypes = lower(["Non-Epileptic Seizure Disorder","Uncertain if Epilepsy","Unknown or MRN not found",""]);

% For plotting log10
EPS_PLOT = 100e-3;   % spikes/hour guard for zeros

%% ===== Select spike-count/duration columns per which_runs =====
switch which_runs
    case 1
        colCount = "FirstRun_Spikes";
        colDur   = "FirstRun_Duration_sec";
        segLabel = 'First run (~1h)';
    case 2
        colCount = "First24Runs_Spikes";
        colDur   = "First24Runs_Duration_sec";
        segLabel = 'First 24 runs (~24h)';
    otherwise
        colCount = "Total_Spikes";
        colDur   = "Duration_sec";
        segLabel = 'Whole file';
end
segLabel = string(segLabel);

%% ================== LOAD SPIKE SUMMARY ==================
S = readtable(outCsv, 'TextType','string', 'VariableNamingRule','preserve');
reqS = unique({'EEG_Name','Patient','Session','Total_Spikes','Duration_sec', char(colCount), char(colDur)});
if ~all(ismember(reqS, S.Properties.VariableNames))
    error('Spike summary missing required columns: %s', strjoin(setdiff(reqS,S.Properties.VariableNames), ', '));
end
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end

% Segment spike rate (spikes/hour) for this analysis
S.SpikeRate_perHour = double(S.(colCount)) ./ (double(S.(colDur))/3600);

% Ambulatory/routine filter (by FULL duration)
if only_amb == 1
    S(S.Duration_sec < 3600*12,:) = [];
elseif only_amb == 2
    S(S.Duration_sec > 3600*12,:) = [];
end

%% ================== LOAD REPORT SHEET ===================
R = readtable(reportFile, 'TextType','string', 'VariableNamingRule','preserve');
reqR = {'patient_id','session_number','start_time_deid','visit_dates_deid','visit_hasSz', ...
        'epilepsy_type','epilepsy_specific'};
assert(all(ismember(reqR, R.Properties.VariableNames)), 'Report file missing required fields.');

% normalize IDs
if ~isnumeric(R.patient_id),     R.patient_id     = double(str2double(string(R.patient_id))); end
if ~isnumeric(R.session_number), R.session_number = double(str2double(string(R.session_number))); end

% EEG date
rawStart = string(R.start_time_deid);
okStart  = ~ismissing(rawStart) & strlength(strtrim(rawStart))>0;
R.EEG_Date = NaT(height(R),1);
try
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart)), 'InputFormat','M/d/yy HH:mm');
catch
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart)));
end

% Parse visit arrays: VisitDates + HasSzVec
R.VisitDates = cell(height(R),1);
R.HasSzVec   = cell(height(R),1);
for i = 1:height(R)
    % Dates
    vstr = strtrim(string(R.visit_dates_deid(i)));
    if strlength(vstr)>0 && vstr~="[]"
        try
            vcell = jsondecode(vstr);
            vdt   = datetime(vcell, 'InputFormat','yyyy-MM-dd');
        catch
            toks = regexp(vstr, '\d{4}-\d{2}-\d{2}', 'match');
            vdt  = datetime(toks, 'InputFormat','yyyy-MM-dd');
        end
        R.VisitDates{i} = vdt(:);
    else
        R.VisitDates{i} = NaT(0,1);
    end

    % HasSz
    hstr = strtrim(string(R.visit_hasSz(i)));
    if strlength(hstr)>0 && hstr~="[]"
        try
            hv = jsondecode(char(hstr));
            R.HasSzVec{i} = double(hv(:));
        catch
            nums = regexp(hstr,'[-]?\d+\.?\d*','match');
            R.HasSzVec{i} = double(str2double(string(nums)));
        end
    else
        R.HasSzVec{i} = nan(0,1);
    end
end

%% ============== Patient-level inclusion by epilepsy_type (NO subtype requirement) ==============
PerPat = derive_patient_types_(R);   % gives Patient, EpilepsyType, EpiType3

% Exclude NESD/uncertain/unknown/empty
et = string(PerPat.EpilepsyType);
isEmpty = ismissing(et) | strlength(strtrim(et))==0;
isBad   = ismember(lower(strtrim(et)), badTypes) | isEmpty;

validP = double(PerPat.Patient(~isBad));

% Restrict S/R to valid epilepsy patients
S = S(ismember(S.Patient, validP), :);
R = R(ismember(R.patient_id, validP), :);

%% ============== OUTPATIENT FILTER (Spe/Radnor OR PATIENT_CLASS OR jay_in_or_out) ==============
if only_outpt == 1
    fprintf('\n=== BUILDING OUTPATIENT KEYS (site rule OR report_PATIENT_CLASS OR jay_in_or_out) ===\n');

    % required outpatient columns
    reqCols = {'report_PATIENT_CLASS','jay_in_or_out'};
    missingReq = setdiff(reqCols, R.Properties.VariableNames);
    if ~isempty(missingReq)
        error('Report table lacks required outpatient columns: %s', strjoin(missingReq, ', '));
    end

    hasAcqCol = ismember('acquired_on', R.Properties.VariableNames);
    if ~hasAcqCol
        warning('Report lacks "acquired_on"; Spe/Radnor outpatient-by-site rule will be unavailable.');
        acqStr_report = strings(height(R),1);
    else
        acqStr_report = lower(strtrim(string(R.acquired_on)));
    end

    % ---- 1) Site-based (Spe/Radnor) ----
    OutptBySite = table('Size',[0 2], ...
                         'VariableTypes',{'double','double'}, ...
                         'VariableNames',{'Patient','Session'});
    if hasAcqCol
        isOutpt_site = ~ismissing(acqStr_report) & ...
                       (contains(acqStr_report,"spe") | contains(acqStr_report,"radnor"));
        if any(isOutpt_site)
            OutptBySite = unique(R(isOutpt_site, {'patient_id','session_number'}));
            OutptBySite.Properties.VariableNames = {'Patient','Session'};
            if ~isnumeric(OutptBySite.Patient)
                OutptBySite.Patient = double(str2double(string(OutptBySite.Patient)));
            end
            if ~isnumeric(OutptBySite.Session)
                OutptBySite.Session = double(str2double(string(OutptBySite.Session)));
            end
        end
        fprintf('[Outpatient rule] Site-based (Spe/Radnor) outpatients: %d sessions.\n', ...
            height(OutptBySite));
    else
        warning('No acquired_on column: cannot use Spe/Radnor site rule.');
    end

    % ---- 2) report_PATIENT_CLASS == "Outpatient" ----
    OutptByClass = table('Size',[0 2], ...
                         'VariableTypes',{'double','double'}, ...
                         'VariableNames',{'Patient','Session'});
    classStr = lower(strtrim(string(R.report_PATIENT_CLASS)));
    isOutpt_class = (classStr == "outpatient");
    if any(isOutpt_class)
        OutptByClass = unique(R(isOutpt_class, {'patient_id','session_number'}));
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

    % ---- 3) jay_in_or_out == "out" ----
    OutptByJay = table('Size',[0 2], ...
                       'VariableTypes',{'double','double'}, ...
                       'VariableNames',{'Patient','Session'});
    jayStr = lower(strtrim(string(R.jay_in_or_out)));
    isOutpt_jay = (jayStr == "out");
    if any(isOutpt_jay)
        OutptByJay = unique(R(isOutpt_jay, {'patient_id','session_number'}));
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

    % ---- 4) union of all three ----
    OutptKeys = [OutptBySite; OutptByClass; OutptByJay];
    if ~isempty(OutptKeys)
        OutptKeys = unique(OutptKeys);
    else
        error(['No outpatient sessions identified from any rule:\n' ...
               '  Spe/Radnor site, report_PATIENT_CLASS=="Outpatient", or jay_in_or_out=="out".']);
    end

    % ---- 5) apply filter to S and R ----
    S = innerjoin(S, OutptKeys, 'Keys', {'Patient','Session'});
    R = innerjoin(R, OutptKeys, ...
        'LeftKeys', {'patient_id','session_number'}, ...
        'RightKeys', {'Patient','Session'});

    fprintf(['\n[Outpatient filter] Outpatients defined as ANY of:\n' ...
         '  - acquired_on containing "Spe" or "Radnor" (case-insensitive) OR\n' ...
         '  - report_PATIENT_CLASS == "Outpatient" OR\n' ...
         '  - jay_in_or_out == "out"\n']);
    fprintf('[Outpatient filter] Kept %d spike rows and %d report rows.\n', ...
        height(S), height(R));
end

if isempty(S)
    warning('No EEGs left after filters; exiting.');
    return;
end

%% ================= JOIN + MATCH (asymmetric windows) =================
rightVars = {'EEG_Date','VisitDates','HasSzVec','epilepsy_type'};
if ~ismember('acquired_on',          R.Properties.VariableNames), R.acquired_on = strings(height(R),1); end
if ~ismember('report_PATIENT_CLASS', R.Properties.VariableNames), R.report_PATIENT_CLASS = strings(height(R),1); end
rightVars{end+1} = 'acquired_on';
rightVars{end+1} = 'report_PATIENT_CLASS';

JR = innerjoin(S, R, ...
    'LeftKeys',{'Patient','Session'}, ...
    'RightKeys',{'patient_id','session_number'}, ...
    'RightVariables', rightVars);

% Keep rows with valid EEG date and at least one visit date
JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) & any(~isnat(v)), JR.VisitDates), :);

pairs = table('Size',[0 9], ...
    'VariableTypes', {'double','double','string','datetime','double','datetime','double','double','string'}, ...
    'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour', ...
                      'Visit_Date','HasSz','GapDays_abs','EpilepsyType'});

for i = 1:height(JR)
    eegDate = JR.EEG_Date(i);
    vdates  = JR.VisitDates{i};
    hasV    = JR.HasSzVec{i};
    pid     = JR.Patient(i);

    if isempty(vdates) || isempty(hasV), continue; end
    nAlign = min(numel(vdates), numel(hasV));
    if nAlign==0, continue; end
    vdates = vdates(1:nAlign);
    hasV   = hasV(1:nAlign);

    dSigned = days(eegDate - vdates);
    elig = (dSigned >= -preGapDays) & (dSigned <= postGapDays);
    if ~any(elig), continue; end

    [~, idxRel] = min(abs(dSigned(elig)));
    idxVec = find(elig);
    j = idxVec(idxRel);

    thisHas = hasV(j);
    if ~(isfinite(thisHas) && (thisHas==0 || thisHas==1))
        continue
    end

    pairs = [pairs; {pid, JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perHour(i), ...
                     vdates(j), double(thisHas), abs(dSigned(j)), string(JR.epilepsy_type(i))}]; %#ok<AGROW>
end

fprintf('Matched %d EEG–visit pairs under [%d days pre, %d days post].\n', height(pairs), preGapDays, postGapDays);
if isempty(pairs)
    warning('No pairs constructed.'); 
    if saveAudit, writetable(pairs, pairsOut); end
    return;
end

% Require ≥2 pairs per patient
gc = groupcounts(pairs,'Patient');
keepP2 = gc.Patient(gc.GroupCount>=2);
pairs = pairs(ismember(pairs.Patient, keepP2), :);
if isempty(pairs)
    warning('Nothing left after ≥2-pair filter.');
    if saveAudit, writetable(pairs, pairsOut); end
    return;
end

%% ================= Per-patient MEANS: HasSz==0 vs HasSz==1 =================
% RAW spike rates (spikes/hour) for stats,
% log10(spikes/hour) for plotting.

[G, pid] = findgroups(double(pairs.Patient));

mu0_raw = splitapply(@(y,hs) mean(y(hs==0), 'omitnan'), pairs.SpikeRate_perHour, pairs.HasSz, G);
mu1_raw = splitapply(@(y,hs) mean(y(hs==1), 'omitnan'), pairs.SpikeRate_perHour, pairs.HasSz, G);
n0      = splitapply(@(hs) sum(hs==0),       pairs.HasSz, G);
n1      = splitapply(@(hs) sum(hs==1),       pairs.HasSz, G);

keep = isfinite(mu0_raw) & isfinite(mu1_raw) & (n0>0) & (n1>0);
yNoSz_raw     = mu0_raw(keep);       % HasSz==0 mean spike rate (spikes/hour)
ySz_raw       = mu1_raw(keep);       % HasSz==1 mean spike rate (spikes/hour)
patients_kept = pid(keep);

Npairs = numel(ySz_raw);

if Npairs < 3
    warning('Fewer than 3 patients with both HasSz==0 and HasSz==1; no stats/plot.');
    if saveAudit, writetable(pairs, pairsOut); end
    return;
end

% Wilcoxon on RAW spike rates
diffs_raw   = ySz_raw - yNoSz_raw;
p_signrank  = signrank(ySz_raw, yNoSz_raw, 'method','approx');
medDeltaRaw = median(diffs_raw, 'omitnan');
med0        = median(yNoSz_raw, 'omitnan');
med1        = median(ySz_raw,   'omitnan');

% For plotting: log10(spikes/hour) with EPS_PLOT guard
logNoSz = log10(max(yNoSz_raw, EPS_PLOT));
logSz   = log10(max(ySz_raw,   EPS_PLOT));

%% ================= Plot (log10 spike rate in spikes/min) =================

% Convert RAW means (spikes/hour) → spikes/min
yNoSz_perMin = yNoSz_raw / 60;
ySz_perMin   = ySz_raw   / 60;

% For plotting: log10(spikes/min) with EPS guard
logNoSz = log10(max(yNoSz_perMin, EPS_PLOT/60));  % convert EPS_PLOT to per-minute scale
logSz   = log10(max(ySz_perMin,   EPS_PLOT/60));

jitterWidth = 0.08;   % adjust as needed (0.05–0.12 looks good)

% Jittered x-positions for points (but NOT for lines)
x0_jit = -jitterWidth/2 + jitterWidth*rand(Npairs,1);      % around 0
x1_jit =  1 - jitterWidth/2 + jitterWidth*rand(Npairs,1);  % around 1

figure('Color','w','Position',[200 200 750 580]);
hold on; grid on; box off; set(gca,'FontSize',16);

% --- Paired lines (NO jitter) ---
for i = 1:Npairs
    plot([0 1], [logNoSz(i) logSz(i)], '-', 'Color',[0.7 0.7 0.7]);
end

% --- Jittered scatter points ---
scatter(x0_jit, logNoSz, 50, 'filled', 'MarkerFaceAlpha',0.8);
scatter(x1_jit, logSz,   50, 'filled', 'MarkerFaceAlpha',0.8);

% Aesthetics
xlim([-0.25 1.25]);
set(gca,'XTick',[0 1], 'XTickLabel',{'HasSz=0 mean','HasSz=1 mean'});
ylabel('log_{10} spike rate (spikes/min)');
title(sprintf('All epilepsy patients — %s', segLabel), 'FontSize',18,'FontWeight','bold');

subtitle(sprintf('n=%d  median HasSz=0=%.3f/hr, HasSz=1=%.3f/hr, median Δ=%.3f/hr,  p=%.3g', ...
    Npairs, med0, med1, medDeltaRaw, p_signrank));


%% ================= Console summary =================
fprintf('\nPaired sign-rank on per-patient mean spike rates (HasSz=1 minus HasSz=0) — %s:\n', segLabel);
fprintf('  n patients with both states = %d\n', Npairs);
fprintf('  Median HasSz=0 mean spike rate      = %.3f spikes/hour\n', med0);
fprintf('  Median HasSz=1 mean spike rate      = %.3f spikes/hour\n', med1);
fprintf('  Median difference (HasSz=1 - HasSz=0) = %.3f spikes/hour\n', medDeltaRaw);
fprintf('  Wilcoxon sign-rank p = %.3g\n', p_signrank);

%% =================== Save audit (optional) ===================
if saveAudit
    writetable(pairs, pairsOut);
    fprintf('Saved EEG–visit pairs audit to: %s\n', pairsOut);
end
end % main


%% =================== HELPER ===================
function PerPat = derive_patient_types_(R)
% Returns table with Patient, EpilepsyType, EpiType3 (General/Temporal/Frontal/Other)
pid   = double(R.patient_id);
etype = string(R.epilepsy_type);
espec = string(R.epilepsy_specific);

T = table(pid, etype, espec, 'VariableNames',{'Patient','EpilepsyType','EpilepsySpecific'});
T = T(~ismissing(T.Patient),:);
T = sortrows(T, 'Patient');
[up, ia] = unique(T.Patient,'stable');
PerPat = T(ia, {'Patient','EpilepsyType','EpilepsySpecific'});
PerPat.Patient = double(PerPat.Patient);

% bucket to display types (not used in this simplified analysis, but kept for compatibility)
etLower = lower(strtrim(PerPat.EpilepsyType));
esLower = lower(strtrim(PerPat.EpilepsySpecific));
isGen   = contains(etLower, 'general');
isTemp  = contains(esLower, 'temporal');
isFront = contains(esLower, 'frontal');

E3 = strings(height(PerPat),1);
E3(isGen)                      = "General";
E3(~isGen & isTemp)            = "Temporal";
E3(~isGen & ~isTemp & isFront) = "Frontal";
E3(E3=="")                     = "Other";
PerPat.EpiType3 = E3;
end
