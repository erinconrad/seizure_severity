function paired_spikerate_by_hasSz(which_runs)
% paired_spikerate_by_hasSz(which_runs)
%   0 (default) = whole file
%   1           = first run (~1h)
%   2           = first 24 runs (~24h)
%
% Purpose:
%   For each patient, compare the mean spike rate across all matched visits
%   with HasSz==0 vs HasSz==1 (paired per patient). Overall includes ALL
%   valid epilepsy patients (not limited to General/Temporal/Frontal).
%   Panels for General/Temporal/Frontal are shown if present.
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
%     (2) report_patient_class == "Outpatient"
%
% Outputs:
%   - 2x2 spaghetti plot of per-patient paired means
%   - Console summary with raw and Bonferroni-corrected p-values
%   - CSV audit of matched pairs (optional toggle)

if nargin < 1, which_runs = 0; end

%% ======================= CONFIG =======================
% Paths
outCsv     = '../data/SN_counts/spike_counts_summary.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-19_1231.csv';
pairsOut   = '../output/eeg_visit_pairs_pairedmeans.csv';
saveAudit  = true;   % write the matched-pairs audit CSV

% Filters
only_amb     = 2;          % 1 = ONLY ambulatory (Duration_sec >= 12h); 2 = ONLY routine (<=12h); 0 = all
only_outpt   = 1;          % 1 = restrict to outpatient EEGs using outpatient rule (see header)
outpt_tokens = ["spe","radnor"]; % lowercase tokens to match in 'acquired_on'

% Asymmetric windows (days) for matching EEG to visit
preGapDays  = 365;   % EEG may be up to this many days BEFORE the visit
postGapDays = 30;    % EEG may be up to this many days AFTER  the visit

% Transform (visualization only)
useLogTransform = false;  % log1p(spikes/hour) before averaging/plotting

% Panels to display (not used for inclusion)
displayTypes3 = ["General","Temporal","Frontal"];

% Patient-level exclusion list (applies to EpilepsyType string)
badTypes = lower(["Non-Epileptic Seizure Disorder","Uncertain if Epilepsy","Unknown or MRN not found",""]);

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

% Segment spike rate (spikes/hour)
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

hasAcq   = ismember('acquired_on',          R.Properties.VariableNames); % optional
hasClass = ismember('report_patient_class', R.Properties.VariableNames); % optional

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

% Parse visit arrays
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

%% ============== Patient-level inclusion by epilepsy_type (NO type requirement) ==============
% Build per-patient type mapping (EpilepsyType/EpiType3) without restricting to G/T/F.
PerPat = derive_patient_types_(R);

% Exclude patients with missing/empty epilepsy_type or in badTypes list
et = string(PerPat.EpilepsyType);
isEmpty = ismissing(et) | strlength(strtrim(et))==0;
isBad   = ismember(lower(strtrim(et)), badTypes) | isEmpty;

validP = double(PerPat.Patient(~isBad));

% Restrict S/R to valid epilepsy patients (no requirement for General/Temporal/Frontal)
S = S(ismember(S.Patient, validP), :);
R = R(ismember(R.patient_id, validP), :);

%% ============== Optional outpatient restriction (new rule) ==============
if only_outpt == 1
    if ~(hasAcq || hasClass)
        warning('Report lacks both "acquired_on" and "report_patient_class"; outpatient filter skipped.');
    else
        % Build a small join table with whichever outpatient fields exist
        keepCols = {'patient_id','session_number'};
        rightVars = {};
        if hasAcq
            keepCols{end+1} = 'acquired_on';
            rightVars{end+1} = 'acquired_on';
        end
        if hasClass
            keepCols{end+1} = 'report_patient_class';
            rightVars{end+1} = 'report_patient_class';
        end

        RJ = R(:, keepCols);
        SJ = innerjoin(S, RJ, ...
            'LeftKeys',  {'Patient','Session'}, ...
            'RightKeys', {'patient_id','session_number'}, ...
            'RightVariables', rightVars);

        % Outpatient rule:
        %   (acquired_on contains "spe"/"radnor") OR (report_patient_class == "Outpatient")
        nSJ = height(SJ);
        isSpeRad   = false(nSJ,1);
        isOutClass = false(nSJ,1);

        if hasAcq
            acq = strtrim(lower(string(SJ.acquired_on)));
            isSpeRad = ~ismissing(acq) & contains_any_token_(acq, outpt_tokens);
        end
        if hasClass
            cls = strtrim(lower(string(SJ.report_patient_class)));
            isOutClass = ~ismissing(cls) & strcmp(cls, 'outpatient');
        end

        keep_out = isSpeRad | isOutClass;
        SJ = SJ(keep_out, :);

        if isempty(SJ)
            warning('Outpatient filter removed all EEGs; continuing with empty S.');
            S = S([],:);  % will bail later when no pairs
        else
            % Keep original S variables, plus SpikeRate_perHour and outpatient info if present
            keepVars = intersect(S.Properties.VariableNames, SJ.Properties.VariableNames, 'stable');
            keepVars = unique([keepVars, {'SpikeRate_perHour'} , rightVars], 'stable');
            S = SJ(:, keepVars);
        end
    end
end

%% ================= JOIN + MATCH (asymmetric windows) =================
rightVars = {'EEG_Date','VisitDates','HasSzVec','epilepsy_type','epilepsy_specific'};
if hasAcq,   rightVars{end+1} = 'acquired_on'; end
if hasClass, rightVars{end+1} = 'report_patient_class'; end

JR = innerjoin(S, R, ...
    'LeftKeys',{'Patient','Session'}, ...
    'RightKeys',{'patient_id','session_number'}, ...
    'RightVariables', rightVars);

% If 'acquired_on' absent, create placeholder for audit
if ~ismember('acquired_on', JR.Properties.VariableNames)
    JR.acquired_on = strings(height(JR),1);
end
if ~ismember('report_patient_class', JR.Properties.VariableNames)
    JR.report_patient_class = strings(height(JR),1);
end

% Keep rows with valid EEG date and at least one visit date
JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) & any(~isnat(v)), JR.VisitDates), :);

% Fast map for display type bucket (optional)
type3Map = containers.Map(double(PerPat.Patient), cellstr(PerPat.EpiType3));

pairs = table('Size',[0 13], ...
    'VariableTypes', {'double','double','string','datetime','double','datetime','double','double','string','string','double','double','string'}, ...
    'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour', ...
                      'Visit_Date','HasSz','GapDays_abs','EpilepsyType','EpiType3', ...
                      'Meanable0','Meanable1','Segment'});

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

    % annotate display bucket (may be "Other")
    e3 = ""; if isKey(type3Map, pid), e3 = string(type3Map(pid)); end

    pairs = [pairs; {pid, JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perHour(i), ...
                     vdates(j), double(thisHas), abs(dSigned(j)), string(JR.epilepsy_type(i)), e3, NaN, NaN, segLabel}]; %#ok<AGROW>
end

fprintf('Matched %d EEG–visit pairs under [%d days pre, %d days post].\n', height(pairs), preGapDays, postGapDays);
if isempty(pairs)
    warning('No pairs constructed.'); if saveAudit, writetable(pairs, pairsOut); end; return
end

% Require ≥2 pairs per patient
gc = groupcounts(pairs,'Patient');
keepP2 = gc.Patient(gc.GroupCount>=2);
pairs = pairs(ismember(pairs.Patient, keepP2), :);
if isempty(pairs)
    warning('Nothing left after ≥2-pair filter.'); if saveAudit, writetable(pairs, pairsOut); end; return
end

%% ================= Per-patient MEANS: HasSz==0 vs HasSz==1 =================
if useLogTransform
    pairs.SpikeResp = log1p(pairs.SpikeRate_perHour);
    yLabel = 'log(1 + spike rate, spikes/hour)';
else
    pairs.SpikeResp = pairs.SpikeRate_perHour;
    yLabel = 'spike rate (spikes/hour)';
end

% Overall includes ALL valid patients (including "Other")
sets = struct();
sets.Overall  = pairs;
sets.General  = pairs(pairs.EpiType3=="General",:);
sets.Temporal = pairs(pairs.EpiType3=="Temporal",:);
sets.Frontal  = pairs(pairs.EpiType3=="Frontal",:);

order = {'Overall','General','Temporal','Frontal'};
pvals = nan(4,1); Npairs = nan(4,1); medDelta = nan(4,1);

figure('Color','w'); tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for g = 1:4
    name = order{g};
    T = sets.(name);
    if isempty(T)
        nexttile; axis off; title(sprintf('%s — %s', name, segLabel));
        continue
    end

    [G, pid] = findgroups(double(T.Patient));
    mu0 = splitapply(@(y,hs) mean(y(hs==0), 'omitnan'), T.SpikeResp, T.HasSz, G);
    mu1 = splitapply(@(y,hs) mean(y(hs==1), 'omitnan'), T.SpikeResp, T.HasSz, G);
    n0  = splitapply(@(hs) sum(hs==0), T.HasSz, G);
    n1  = splitapply(@(hs) sum(hs==1), T.HasSz, G);

    keep = isfinite(mu0) & isfinite(mu1) & (n0>0) & (n1>0);
    yNoSz = mu0(keep);
    ySz   = mu1(keep);

    Npairs(g) = numel(ySz);

    nexttile; hold on; grid on; box off; set(gca,'FontSize',14);
    if Npairs(g) >= 3
        diffs = ySz - yNoSz;
        p = signrank(ySz, yNoSz, 'method','approx');
        pvals(g) = p;
        medDelta(g) = median(diffs,'omitnan');

        % spaghetti plot
        for i = 1:Npairs(g)
            plot([0 1],[yNoSz(i) ySz(i)],'-','Color',[0.7 0.7 0.7])
        end
        scatter(zeros(Npairs(g),1), yNoSz, 28, 'filled')
        scatter(ones(Npairs(g),1),  ySz,   28, 'filled')
        xlim([-0.25 1.25]); set(gca,'XTick',[0 1],'XTickLabel',{'HasSz=0 mean','HasSz=1 mean'})
        ylabel(yLabel);
        title(sprintf('%s — %s', name, segLabel));

        if useLogTransform
            pct = 100*(exp(medDelta(g)) - 1);
            subtitle(sprintf('n=%d  median Δ=%.3f (%+0.0f%%)  p=%.3g', Npairs(g), medDelta(g), pct, pvals(g)));
        else
            subtitle(sprintf('n=%d  median Δ=%.3f  p=%.3g', Npairs(g), medDelta(g), pvals(g)));
        end
    else
        title(sprintf('%s — %s', name, segLabel));
        text(0.5,0.5,'Insufficient patients with both states','HorizontalAlignment','center'); axis off
    end
end

% Bonferroni across the 4 panels we display
pvals_bonf = min(pvals * 4, 1);
fprintf('\nPaired sign-rank on per-patient means (HasSz=1 minus HasSz=0) — %s:\n', segLabel);
for g = 1:4
    nm = order{g};
    if ~isnan(pvals(g))
        if useLogTransform
            pct = 100*(exp(medDelta(g)) - 1);
            fprintf('  %-9s  n=%-3d  median Δ=%.3f (%+0.0f%%)  p=%.3g  p_bonf=%.3g\n', ...
                nm, Npairs(g), medDelta(g), pct, pvals(g), pvals_bonf(g));
        else
            fprintf('  %-9s  n=%-3d  median Δ=%.3f          p=%.3g  p_bonf=%.3g\n', ...
                nm, Npairs(g), medDelta(g), pvals(g), pvals_bonf(g));
        end
    else
        fprintf('  %-9s  insufficient pairs\n', nm);
    end
end

%% =================== Save audit (optional) ===================
if saveAudit
    writetable(pairs, pairsOut);
    fprintf('Saved EEG–visit pairs audit to: %s\n', pairsOut);
end
end % main

%% =================== HELPERS ===================
function tf = contains_any_token_(strs, tokens)
% strs: string array; tokens: string array, already lowercase
tf = false(size(strs));
for t = 1:numel(tokens)
    tf = tf | contains(strs, tokens(t));
end
end

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

% bucket to display types (used only for subpanels; NOT for inclusion)
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
