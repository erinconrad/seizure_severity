function model_szfreq_vs_spikerate(which_runs)
% model_szfreq_vs_spikerate(which_runs)
%   0 (default) = whole file
%   1           = first run (~1h)
%   2           = first 24 runs (~24h)
%
% Models: log10(seizures/month) ~ log10(spikes/min) + (1|Patient)
% Matching:
%   - Pairs each EEG with the nearest visit under asymmetric bounds:
%       -preGapDays <= (EEG_Date - Visit_Date) <= +postGapDays
%   - Uses sz_freqs at the matched visit when present.
%   - If sz_freqs is missing/NaN/null at that visit AND visit_hasSz==0 → uses 0.
%   - Otherwise, the pair is dropped (cannot infer seizure frequency).
%
% Output:
%   - Figures: Overall scatter with fit; 3 panel facets (General/Temporal/Frontal)
%   - CSV audit of all pairs and selected seizure frequency
%
% Assumptions:
%   - sz_freqs are reported as seizures/month in R.visit-level arrays.
%   - Spike rate computed from selected segment; converted to spikes/min.
%
% Erin-friendly knobs are grouped under CONFIG.

if nargin < 1, which_runs = 0; end

%% ======================= CONFIG =======================
outCsv      = '../data/SN_counts/spike_counts_summary.csv';
reportFile  = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-10_1443.csv';
pairsOut    = '../output/eeg_visit_pairs_szfreq_vs_spikerate.csv';

only_amb    = 0;    % 1 = ONLY ambulatory (Duration_sec >= 12h by FULL duration rule)
                    % 2 = ONLY routine  (Duration_sec <= 12h by FULL duration rule)
                    % 0 = all
only_outpt  = 1;    % 1 = filter to outpatient by 'acquired_on' tokens below
outpt_tokens = ["spe","radnor"];   % lowercase tokens to match for outpatient

% Asymmetric windows for matching EEG to visits
preGapDays  = 180;  % EEG may be up to this many days BEFORE the visit
postGapDays = 30;  % EEG may be up to this many days AFTER  the visit

% Log transforms (base-10) with small epsilons to accommodate zeros
EPS_SPIKES_PER_MIN = 1e-3;  % spikes/min floor for log10
EPS_SZ_PER_MON     = 1e-3;  % seizures/month floor for log10

% Optional epilepsy type handling
wantTypes3 = ["General","Temporal","Frontal"];  % facet order
keep_only_these_types = true;   % restrict to these 3 canonical buckets

% Figure axis limits (applied to all panels for strict comparability)
FIX_AXES = true;
X_LIMS = [-3.5, 4.0];  % log10(spikes/min)
Y_LIMS = [-3.0, 2.0];  % log10(seizures/month)

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
    missing = setdiff(reqS, S.Properties.VariableNames);
    error('Spike summary missing required columns: %s', strjoin(missing, ', '));
end
% numericize keys
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end

% Segment spike rate → spikes/min
S.SpikeRate_perMin = double(S.(colCount)) ./ max(double(S.(colDur)),eps) * 60;

% Ambulatory/routine filter using FULL duration column (unchanged rule)
if only_amb == 1
    S(S.Duration_sec < 3600*12, :) = [];
elseif only_amb == 2
    S(S.Duration_sec > 3600*12, :) = [];
end

%% ================== LOAD REPORT SHEET (R) ===================
R = readtable(reportFile, 'TextType','string', 'VariableNamingRule','preserve');
reqR = {'patient_id','session_number','start_time_deid','visit_dates_deid', ...
        'visit_hasSz','sz_freqs','epilepsy_type','epilepsy_specific','acquired_on'};
assert(all(ismember(reqR, R.Properties.VariableNames)), 'Report file missing required fields.');

% Keys / types
if ~isnumeric(R.patient_id),     R.patient_id     = double(str2double(string(R.patient_id))); end
if ~isnumeric(R.session_number), R.session_number = double(str2double(string(R.session_number))); end

% EEG date
rawStart = string(R.start_time_deid);
okStart  = ~ismissing(rawStart) & strlength(strtrim(rawStart))>0;
R.EEG_Date = NaT(height(R),1);
try
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart)), 'InputFormat','M/d/yy HH:mm');
catch
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart))); % fallback autodetect
end

% Parse visit arrays: dates, hasSz, sz_freqs
R.VisitDates = cell(height(R),1);
R.HasSzVec   = cell(height(R),1);
R.SzFreqVec  = cell(height(R),1);

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

    % hasSz
    hstr = strtrim(string(R.visit_hasSz(i)));
    if strlength(hstr)>0 && hstr~="[]"
        try
            hv = jsondecode(char(hstr));
            R.HasSzVec{i} = double(hv(:));
        catch
            nums = regexp(hstr, '[-]?\d+\.?\d*', 'match'); % allow negatives just in case
            R.HasSzVec{i} = double(str2double(string(nums)));
        end
    else
        R.HasSzVec{i} = nan(0,1);
    end

    % sz_freqs (seizures/month), may contain nulls
    fstr = strtrim(string(R.sz_freqs(i)));
    if strlength(fstr)>0 && fstr~="[]"
        try
            fv = jsondecode(char(fstr));
            % JSON null becomes []/NaN after processing — coerce to double with NaN for empties
            if iscell(fv)
                tmp = nan(numel(fv),1);
                for k = 1:numel(fv)
                    if isempty(fv{k})
                        tmp(k) = NaN;
                    else
                        tmp(k) = double(fv{k});
                    end
                end
                R.SzFreqVec{i} = tmp;
            else
                R.SzFreqVec{i} = double(fv(:));
            end
        catch
            % fallback: scrape numbers and 'null'
            parts = regexp(fstr,'null|[-+]?\d*\.?\d+','match','ignorecase');
            tmp = nan(numel(parts),1);
            for k = 1:numel(parts)
                if strcmpi(parts{k},'null')
                    tmp(k) = NaN;
                else
                    tmp(k) = str2double(parts{k});
                end
            end
            R.SzFreqVec{i} = tmp;
        end
    else
        R.SzFreqVec{i} = nan(0,1);
    end
end

% Per-patient epilepsy buckets → General/Temporal/Frontal
PerPat = derive_patient_types_(R);
if keep_only_these_types
    PerPat = PerPat(ismember(PerPat.EpiType3, wantTypes3), :);
end

% Subset S/R to valid patients and outpatient if requested
validP = PerPat.Patient;
S = S(ismember(S.Patient, validP), :);
R = R(ismember(R.patient_id, validP), :);

if only_outpt == 1
    if ismember('acquired_on', R.Properties.VariableNames)
        RJ = R(:, {'patient_id','session_number','acquired_on'});
        SJ = innerjoin(S, RJ, 'LeftKeys',{'Patient','Session'}, 'RightKeys',{'patient_id','session_number'}, 'RightVariables','acquired_on');
        acq = strtrim(lower(string(SJ.acquired_on)));
        keep_out = ~ismissing(acq) & contains_any_token_(acq, outpt_tokens);
        SJ = SJ(keep_out, :);
        keepVars = intersect(S.Properties.VariableNames, SJ.Properties.VariableNames, 'stable');
        keepVars = unique([keepVars, {'SpikeRate_perMin'}], 'stable');
        S = SJ(:, keepVars);
    else
        warning('Report lacks "acquired_on"; only_outpt flag ignored.');
    end
end

%% ============ JOIN S+R, BUILD PAIRS WITH sz_freq RESOLUTION ============
JR = innerjoin(S, R, ...
    'LeftKeys',{'Patient','Session'}, ...
    'RightKeys',{'patient_id','session_number'}, ...
    'RightVariables',{'EEG_Date','VisitDates','HasSzVec','SzFreqVec','epilepsy_type','epilepsy_specific'});

% Keep rows with valid EEG date and at least one visit date
JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) & any(~isnat(v)), JR.VisitDates), :);

pairs = table('Size',[0 17], ...
    'VariableTypes', {'double','double','string','datetime','double','datetime','double','double','double','logical','string','string','double','double','double','double','string'}, ...
    'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perMin', ...
                      'Visit_Date','HasSz','SzFreq_perMon','SzFreq_SourceIsZeroFill', ...
                      'UsedZeroFill','EpiType','EpiType3','GapDays_abs','GapDays_signed', ...
                      'log10_SzPerMon','log10_SpikesPerMin','Segment'});

% fast maps
typeMap  = containers.Map(double(PerPat.Patient), cellstr(PerPat.EpilepsyType));
type3Map = containers.Map(double(PerPat.Patient), cellstr(PerPat.EpiType3));

for i = 1:height(JR)
    eegDate = JR.EEG_Date(i);
    vdates  = JR.VisitDates{i};
    hasV    = JR.HasSzVec{i};
    freqV   = JR.SzFreqVec{i};
    pid     = JR.Patient(i);

    if isempty(vdates) || isempty(hasV), continue; end
    nAlign = min([numel(vdates), numel(hasV), numel(freqV)]);
    if nAlign==0, continue; end
    vdates = vdates(1:nAlign);
    hasV   = hasV(1:nAlign);
    freqV  = freqV(1:nAlign);

    dSigned = days(eegDate - vdates);
    elig = (dSigned >= -preGapDays) & (dSigned <= postGapDays);
    if ~any(elig), continue; end

    [~, idxRel] = min(abs(dSigned(elig)));
    idxVec = find(elig);
    j = idxVec(idxRel);

    % Resolve seizure frequency for that visit:
    % Rule: if freq is missing/NaN and HasSz==0 → use 0; else if missing and HasSz~=0 → drop
    thisHas = hasV(j);
    thisFreq = freqV(j);
    usedZeroFill = false;  % true if we set freq to 0 per rule
    if ~(isfinite(thisFreq))
        if isfinite(thisHas) && thisHas==0
            thisFreq = 0;
            usedZeroFill = true;
        else
            % cannot determine frequency for this matched visit
            continue
        end
    end

    etype  = ""; if isKey(typeMap, pid),  etype  = string(typeMap(pid));  end
    etype3 = ""; if isKey(type3Map, pid), etype3 = string(type3Map(pid)); end
    if keep_only_these_types && ~ismember(etype3, wantTypes3), continue; end

    gap_signed = dSigned(j);
    gap_abs    = abs(gap_signed);

    % logs
    logSpike = log10(max(JR.SpikeRate_perMin(i), EPS_SPIKES_PER_MIN));
    logSz    = log10(max(thisFreq, EPS_SZ_PER_MON));

    pairs = [pairs; {pid, JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perMin(i), ...
                     vdates(j), double(thisHas), double(thisFreq), double(usedZeroFill), usedZeroFill, ...
                     etype, etype3, gap_abs, gap_signed, logSz, logSpike, segLabel}]; %#ok<AGROW>
end

fprintf('Built %d EEG–visit pairs under [%d days pre, %d days post].\n', height(pairs), preGapDays, postGapDays);
if isempty(pairs)
    warning('No pairs were constructed.'); writetable(pairs, pairsOut); return
end

% ≥2 pairs per patient (stability)
gc = groupcounts(pairs,'Patient');
keepP = gc.Patient(gc.GroupCount>=2);
pairs = pairs(ismember(pairs.Patient, keepP), :);

% Remove patients whose EEGs all matched the same visit date
[G, pidg] = findgroups(pairs.Patient);
nUniqueVisits = splitapply(@(d) numel(unique(d)), pairs.Visit_Date, G);
keepPatients = pidg(nUniqueVisits>=2);
pairs = pairs(ismember(pairs.Patient, keepPatients), :);

if isempty(pairs)
    warning('Nothing left after patient/visit filters.'); writetable(pairs, pairsOut); return
end

%% =================== MODEL: log10(sz/month) ~ log10(spikes/min) + (1|Patient) ===================
tbl = pairs;
tbl.Patient = categorical(tbl.Patient);

% Basic LME
lme = fitlme(tbl, 'log10_SzPerMon ~ log10_SpikesPerMin + (1|Patient)', 'DummyVarCoding','reference');
disp(lme);

% Optional: per-type slopes (quietly skip if level missing)
try
    lme_type = fitlme(tbl, 'log10_SzPerMon ~ log10_SpikesPerMin * EpiType3 + (1|Patient)', 'DummyVarCoding','reference');
    disp(lme_type);
catch
    lme_type = [];
end

%% =================== PLOTS ===================
% Panel A: Overall scatter + fitted line (conditional on RE=0)
figure('Color','w'); hold on; box off; grid on; set(gca,'FontSize',16);
scatter(tbl.log10_SpikesPerMin, tbl.log10_SzPerMon, 18, 'filled', 'MarkerFaceAlpha',0.6);
% Fit line from lme (population-level)
xline_grid = linspace(min(tbl.log10_SpikesPerMin), max(tbl.log10_SpikesPerMin), 200)';
newT = table(xline_grid, repmat(tbl.Patient(1),numel(xline_grid),1), 'VariableNames',{'log10_SpikesPerMin','Patient'});
[yhat, yCI] = predict(lme, newT, 'Conditional', false);
plot(xline_grid, yhat, 'LineWidth',2);
for i = 1:numel(xline_grid)
    line([xline_grid(i) xline_grid(i)], [yCI(i,1) yCI(i,2)], 'LineWidth',1);
end
xlabel('log_{10} (spikes/min)'); ylabel('log_{10} (seizures/month)');
title(sprintf('Overall: %s', segLabel));
if FIX_AXES, xlim(X_LIMS); ylim(Y_LIMS); end

% Panel B-D: Facets by EpiType3
types = wantTypes3;
figure('Color','w'); tiledlayout(1, numel(types), 'TileSpacing','compact','Padding','compact');
for k = 1:numel(types)
    nexttile; hold on; box off; grid on; set(gca,'FontSize',14);
    Tk = tbl(tbl.EpiType3==types(k), :);
    if isempty(Tk)
        text(0.5,0.5,'No data','HorizontalAlignment','center'); axis off; continue
    end
    scatter(Tk.log10_SpikesPerMin, Tk.log10_SzPerMon, 18, 'filled', 'MarkerFaceAlpha',0.6);

    % simple OLS fit in this stratum (visual guide; not mixed model)
    X = [ones(height(Tk),1) Tk.log10_SpikesPerMin];
    if rank(X)>=2
        b = X\Tk.log10_SzPerMon;
        xx = linspace(min(Tk.log10_SpikesPerMin), max(Tk.log10_SpikesPerMin), 100);
        yy = b(1) + b(2)*xx;
        plot(xx, yy, '-', 'LineWidth',2);
    end
    xlabel('log_{10} (spikes/min)');
    if k==1, ylabel('log_{10} (seizures/month)'); end
    title(sprintf('%s — %s', types(k), segLabel));
    if FIX_AXES, xlim(X_LIMS); ylim(Y_LIMS); end
end

%% =================== CONSOLE SUMMARY ===================
coef = dataset2table(lme.Coefficients); % for nicer printing if available
b_row = strcmp(lme.Coefficients.Name,'log10_SpikesPerMin');
if any(b_row)
    est = lme.Coefficients.Estimate(b_row);
    se  = lme.Coefficients.SE(b_row);
    p   = lme.Coefficients.pValue(b_row);
    fprintf('\nOverall slope (log10 Sz/mon vs log10 Spikes/min): β = %.3f ± %.3f, p = %.3g\n', est, se, p);
end

%% =================== SAVE AUDIT ===================
writetable(pairs, pairsOut);
fprintf('Saved EEG–visit audit (with sz freq resolution) to: %s\n', pairsOut);

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
% Returns table with Patient, EpilepsyType, EpiType3 ∈ {General, Temporal, Frontal, (others)}
pid = double(R.patient_id);
etype = string(R.epilepsy_type);
espec = string(R.epilepsy_specific);

% take first non-missing per patient
T = table(pid, etype, espec, 'VariableNames',{'Patient','EpilepsyType','EpilepsySpecific'});
T = T(~ismissing(T.Patient),:);
T = sortrows(T, 'Patient');
[up, ia] = unique(T.Patient,'stable');
PerPat = T(ia, {'Patient','EpilepsyType','EpilepsySpecific'});
PerPat.Patient = double(PerPat.Patient);

% bucket
etLower = lower(strtrim(PerPat.EpilepsyType));
esLower = lower(strtrim(PerPat.EpilepsySpecific));
isGen   = contains(etLower, 'general');
isTemp  = contains(esLower, 'temporal');
isFront = contains(esLower, 'frontal');

E3 = strings(height(PerPat),1);
E3(isGen)                        = "General";
E3(~isGen & isTemp)              = "Temporal";
E3(~isGen & ~isTemp & isFront)   = "Frontal";
E3(E3=="")                       = "Other";

PerPat.EpiType3 = E3;
end
