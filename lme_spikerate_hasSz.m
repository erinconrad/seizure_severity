function lme_spikerate_hasSz(which_runs)
% lme_spikerate_hasSz(which_runs)
%   0 (default) = whole file
%   1           = first run (~1h)
%   2           = first 24 runs (~24h)
%
% Model: log10(spikes/min) ~ HasSz + (1 | Patient)
% Matching: nearest clinic visit within ±365 days (HasSz from visit)

if nargin < 1, which_runs = 0; end

%% ======================= CONFIG =======================
outCsv     = '../data/SN_counts/spike_counts_summary.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-10_1443.csv';
pairsOut   = '../output/eeg_visit_pairs_for_LME.csv';
saveAudit  = true;

only_amb     = 0;                 % 1=only ambulatory, 2=only routine, 0=all
only_outpt   = 1;                 % honored only if 'acquired_on' exists
outpt_tokens = ["spe","radnor"];  % lower-case tokens
gapDays      = 300;               % ± days window
wantTypes3   = ["General","Temporal","Frontal"];

%% ===== Select segment =====
switch which_runs
    case 1
        colCount = "FirstRun_Spikes";     colDur = "FirstRun_Duration_sec";     segLabel = 'First run (~1h)';
    case 2
        colCount = "First24Runs_Spikes";  colDur = "First24Runs_Duration_sec";  segLabel = 'First 24 runs (~24h)';
    otherwise
        colCount = "Total_Spikes";        colDur = "Duration_sec";              segLabel = 'Whole file';
end
segLabel = string(segLabel);

%% ================== LOAD SPIKE SUMMARY ==================
S = readtable(outCsv,'TextType','string','VariableNamingRule','preserve');
reqS = unique({'EEG_Name','Patient','Session','Total_Spikes','Duration_sec', char(colCount), char(colDur)});
if ~all(ismember(reqS, S.Properties.VariableNames))
    error('Spike summary missing required columns: %s', strjoin(setdiff(reqS,S.Properties.VariableNames), ', '));
end
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end

dur = double(S.(colDur));
cnt = double(S.(colCount));
validDur = isfinite(dur) & dur > 0;
S.SpikeRate_perMin = nan(height(S),1);
S.SpikeRate_perMin(validDur) = (cnt(validDur) ./ dur(validDur)) * 60;

if only_amb == 1
    S(S.Duration_sec < 3600*12,:) = [];
elseif only_amb == 2
    S(S.Duration_sec > 3600*12,:) = [];
end

%% ================== LOAD REPORT SHEET ===================
R = readtable(reportFile,'TextType','string','VariableNamingRule','preserve');
reqR = {'patient_id','session_number','start_time_deid','visit_dates_deid','visit_hasSz','epilepsy_type','epilepsy_specific'};
if ~all(ismember(reqR, R.Properties.VariableNames))
    error('Report file missing required fields: %s', strjoin(setdiff(reqR,R.Properties.VariableNames), ', '));
end
% 'acquired_on' is optional; handle below.

if ~isnumeric(R.patient_id),     R.patient_id     = double(str2double(string(R.patient_id))); end
if ~isnumeric(R.session_number), R.session_number = double(str2double(string(R.session_number))); end

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
        try, vcell = jsondecode(vstr);  vdt = datetime(vcell,'InputFormat','yyyy-MM-dd');
        catch, toks = regexp(vstr,'\d{4}-\d{2}-\d{2}','match'); vdt = datetime(toks,'InputFormat','yyyy-MM-dd');
        end
        R.VisitDates{i} = vdt(:);
    else
        R.VisitDates{i} = NaT(0,1);
    end

    % HasSz
    hstr = strtrim(string(R.visit_hasSz(i)));
    if strlength(hstr)>0 && hstr~="[]"
        try, hv = jsondecode(char(hstr)); R.HasSzVec{i} = double(hv(:));
        catch, nums = regexp(hstr,'[-]?\d+\.?\d*','match'); R.HasSzVec{i} = double(str2double(string(nums)));
        end
    else
        R.HasSzVec{i} = nan(0,1);
    end
end

%% ============== Patient types (optional; do NOT filter) ==============
% Keep the mapping (for audit/optional covariates) but do not restrict cohort
PerPat = derive_patient_types_(R);           % may contain General/Temporal/Frontal/Other/blank
validP = unique(double(R.patient_id(~isnan(R.patient_id))));  % include everyone present in R
S = S(ismember(S.Patient, validP), :);
R = R(ismember(R.patient_id, validP), :);


%% ============== Optional outpatient restriction ==============
hasAcq = ismember('acquired_on', R.Properties.VariableNames);
if only_outpt == 1 && hasAcq
    RJ = R(:, {'patient_id','session_number','acquired_on'});
    SJ = innerjoin(S, RJ, 'LeftKeys',{'Patient','Session'}, 'RightKeys',{'patient_id','session_number'}, 'RightVariables','acquired_on');
    acq = strtrim(lower(string(SJ.acquired_on)));
    keep_out = ~ismissing(acq) & contains_any_token_(acq, outpt_tokens);
    SJ = SJ(keep_out, :);
    keepVars = intersect(S.Properties.VariableNames, SJ.Properties.VariableNames, 'stable');
    keepVars = unique([keepVars, {'SpikeRate_perMin','acquired_on'}], 'stable');
    S = SJ(:, keepVars);
elseif only_outpt == 1 && ~hasAcq
    warning('Report lacks "acquired_on"; outpatient filter skipped.');
end

%% ================= JOIN + MATCH (±365 d) =================
rightVars = {'EEG_Date','VisitDates','HasSzVec','epilepsy_type','epilepsy_specific'};
if hasAcq, rightVars{end+1} = 'acquired_on'; end

JR = innerjoin(S, R, ...
    'LeftKeys', {'Patient','Session'}, ...
    'RightKeys',{'patient_id','session_number'}, ...
    'RightVariables', rightVars);

% If 'acquired_on' absent, make a placeholder so the audit always works
if ~ismember('acquired_on', JR.Properties.VariableNames)
    JR.acquired_on = strings(height(JR),1);
end

JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) & any(~isnat(v)), JR.VisitDates), :);

typeMap3 = containers.Map(double(PerPat.Patient), cellstr(PerPat.EpiType3));

matches = table('Size',[0 10], ...
    'VariableTypes', {'double','double','string','datetime','double','datetime','double','double','string','string'}, ...
    'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perMin', ...
                      'Visit_Date','HasSz','GapDays_abs','EpiType3','acquired_on'});

for i = 1:height(JR)
    eegDate = JR.EEG_Date(i);
    vdates  = JR.VisitDates{i};
    hasV    = JR.HasSzVec{i};
    pid     = JR.Patient(i);

    if isempty(vdates) || isempty(hasV), continue; end
    nAlign = min(numel(vdates), numel(hasV));
    if nAlign==0, continue; end
    vdates = vdates(1:nAlign); hasV = hasV(1:nAlign);

    dSigned = days(eegDate - vdates);
    elig = abs(dSigned) <= gapDays;
    if ~any(elig), continue; end

    [~, idxRel] = min(abs(dSigned(elig)));
    j = find(elig); j = j(idxRel);

    thisHas = hasV(j);
    if ~(isfinite(thisHas) && (thisHas==0 || thisHas==1)), continue; end

    e3 = ""; if isKey(typeMap3, pid), e3 = string(typeMap3(pid)); end
    matches = [matches; { ...
        pid, JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perMin(i), ...
        vdates(j), double(thisHas), abs(dSigned(j)), e3, string(JR.acquired_on(i))}]; %#ok<AGROW>
end

fprintf('Matched %d EEG–visit rows within ±%d days. Segment: %s\n', height(matches), gapDays, segLabel);
if isempty(matches), if saveAudit, writetable(matches, pairsOut); end, return; end

okRate = isfinite(matches.SpikeRate_perMin) & matches.SpikeRate_perMin >= 0;
matches = matches(okRate, :);
if isempty(matches), if saveAudit, writetable(matches, pairsOut); end, return; end

%% ================= Transform to log10(spikes/min) =================
pos = matches.SpikeRate_perMin(matches.SpikeRate_perMin > 0);
eps_rate = iff(isempty(pos), 1e-6, 0.5*min(pos));
matches.Log10Rate = log10(matches.SpikeRate_perMin + (matches.SpikeRate_perMin<=0).*eps_rate);

matches.HasSz  = categorical(matches.HasSz, [0 1], {'0','1'}); % 0 as reference
matches.Patient = categorical(matches.Patient);

%% ================= Fit LME =================
lme = fitlme(matches, 'Log10Rate ~ HasSz + (1|Patient)', 'FitMethod','REML');

fprintf('\n===== LME: Log10(spikes/min) ~ HasSz + (1|Patient) =====\n');
disp(lme)

% --- make Coef a table regardless of MATLAB version ---
Coef = lme.Coefficients;
if isa(Coef,'dataset')
    % Older MATLAB: convert dataset -> table
    try
        Coef = dataset2table(Coef);   % preferred
    catch
        % Fallback if dataset2table is unavailable in your release
        Coef = struct2table(struct(Coef));
    end
end

% now Coef is a table, proceed safely
vn = string(Coef.Properties.VariableNames);
rowHas = strcmp(Coef.Name, 'HasSz_1');

if any(rowHas)
    est  = Coef.Estimate(rowHas);
    se   = Coef.SE(rowHas);
    pval = Coef.pValue(rowHas);

    hasCI = all(ismember({'Lower','Upper'}, vn));
    if hasCI
        lo = Coef.Lower(rowHas); 
        hi = Coef.Upper(rowHas);
    else
        z = 1.96; 
        lo = est - z*se; 
        hi = est + z*se;
    end

    fprintf('HasSz (1 vs 0): beta = %.3f  [%.3f, %.3f],  SE = %.3f,  p = %.3g\n', ...
        est, lo, hi, se, pval);
end

%% ================= Save audit =================
if saveAudit
    writetable(matches, pairsOut);
    fprintf('Saved matched audit to: %s\n', pairsOut);
end
end

%% ================= Helpers =================
function tf = contains_any_token_(strs, tokens)
tf = false(size(strs));
for t = 1:numel(tokens)
    tf = tf | contains(strs, tokens(t));
end
end

function out = iff(cond, a, b), if cond, out=a; else, out=b; end, end

function PerPat = derive_patient_types_(R)
pid   = double(R.patient_id);
etype = string(R.epilepsy_type);
espec = string(R.epilepsy_specific);

T = table(pid, etype, espec, 'VariableNames',{'Patient','EpilepsyType','EpilepsySpecific'});
T = T(~ismissing(T.Patient),:);
T = sortrows(T,'Patient');
[~, ia] = unique(T.Patient,'stable');
PerPat = T(ia, {'Patient','EpilepsyType','EpilepsySpecific'});
PerPat.Patient = double(PerPat.Patient);

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
