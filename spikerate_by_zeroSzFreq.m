function spikerate_by_zeroSzFreq(which_runs)
% spikerate_by_zeroSzFreq(which_runs)
%   0 = whole file (default)
%   1 = first run (~1h)
%   2 = first 24 runs (~24h)
%
% Compares log10(spikes/min) for patients with MeanSzFreq==0 vs >0
% - Outpatient mirroring (only_outpt) so spike-rate and sz-freq cohorts match
% - Excludes non-epilepsy types (NESD / uncertain / unknown)
% - Fixed y-limits; ranksum p-value + significance bar

if nargin < 1, which_runs = 0; end

%% ===== Config =====
only_amb   = 0; % 1 = ONLY ambulatory (>=12h); 2 = ONLY routine (<=12h); 0 = all
only_outpt = 1; % 1 = keep EEGs where acquired_on contains "spe" or "radnor" (case-insens.)
Y_LIMS     = [-3, 2];
EPS_PER_MIN = 1e-6; % epsilon for log10(spikes/min) when spikes/min <= 0

% Paths
outCsv     = '../data/SN_counts/spike_counts_summary.csv';
reportFile = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-10_1443.csv';
figOut     = '../figures/spikerate_by_zeroSzFreq.png';

%% ===== Load spike counts; choose segment =====
S = readtable(outCsv,'TextType','string','VariableNamingRule','preserve');

switch which_runs
    case 1, colCount="FirstRun_Spikes";    colDur="FirstRun_Duration_sec";    segLabel='First run (~1h)';
    case 2, colCount="First24Runs_Spikes"; colDur="First24Runs_Duration_sec"; segLabel='First 24 runs (~24h)';
    otherwise, colCount="Total_Spikes";    colDur="Duration_sec";             segLabel='Whole file';
end

% Ambulatory/routine filter via FULL duration
if only_amb == 1
    S(S.Duration_sec < 12*3600,:) = [];
elseif only_amb == 2
    S(S.Duration_sec > 12*3600,:) = [];
end

% Keys & rates
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end
if ~isnumeric(S.(colCount)), S.(colCount) = double(S.(colCount)); end
if ~isnumeric(S.(colDur)),   S.(colDur)   = double(S.(colDur));   end

S.SpikeRate_perMin = nan(height(S),1);
okDur = S.(colDur) > 0;
S.SpikeRate_perMin(okDur) = (S.(colCount)(okDur) ./ S.(colDur)(okDur)) * (60);

% Guard against tiny negatives
tinyNeg = S.SpikeRate_perMin < 0 & S.SpikeRate_perMin > -1e-12;
S.SpikeRate_perMin(tinyNeg) = 0;
S.SpikeRate_perMin(S.SpikeRate_perMin < -1e-12) = NaN;

%% ===== Load report =====
R = readtable(reportFile,'TextType','string','VariableNamingRule','preserve');
if ~isnumeric(R.patient_id),     R.patient_id     = double(str2double(R.patient_id)); end
if ~isnumeric(R.session_number), R.session_number = double(str2double(R.session_number)); end
if ~isstring(R.epilepsy_type),   R.epilepsy_type  = string(R.epilepsy_type); end
assert(ismember('sz_freqs',R.Properties.VariableNames),'Column "sz_freqs" not found.');

%% ===== OUTPATIENT mirroring (align S & R on same (Patient,Session)) =====
if ismember('acquired_on', R.Properties.VariableNames)
    RJ = R(:, {'patient_id','session_number','acquired_on'});
    SJ = innerjoin(S, RJ, ...
        'LeftKeys', {'Patient','Session'}, ...
        'RightKeys', {'patient_id','session_number'}, ...
        'RightVariables','acquired_on');
    if only_outpt == 1
        acq = strtrim(lower(string(SJ.acquired_on)));
        keep_out = ~ismissing(acq) & (contains(acq,"spe") | contains(acq,"radnor"));
        SJ = SJ(keep_out, :);
    end
    keepVars = intersect(S.Properties.VariableNames, SJ.Properties.VariableNames, 'stable');
    keepVars = unique([keepVars, {'SpikeRate_perMin','acquired_on'}], 'stable');
    S = SJ(:, keepVars);
else
    if only_outpt==1
        warning('Report lacks "acquired_on"; only_outpt flag ignored.');
    end
end

RSess_keep = unique(S(:, {'Patient','Session'}));
RSess_keep.Properties.VariableNames = {'Patient','Session'};
RSess_keep.Patient = double(RSess_keep.Patient);
RSess_keep.Session = double(RSess_keep.Session);

R = innerjoin(R, RSess_keep, ...
    'LeftKeys', {'patient_id','session_number'}, ...
    'RightKeys', {'Patient','Session'});

%% ===== Exclude non-epilepsy types =====
badTypes = lower(["Non-Epileptic Seizure Disorder","Uncertain if Epilepsy","Unknown or MRN not found",""]);
isEmptyType = ismissing(R.epilepsy_type) | strlength(strtrim(R.epilepsy_type))==0;
Rt = sortrows(R(~isEmptyType, {'patient_id','epilepsy_type'}), 'patient_id');
[up, ia] = unique(Rt.patient_id, 'stable');
PerPatType = table(up, Rt.epilepsy_type(ia), 'VariableNames', {'Patient','EpilepsyType'});

etype = lower(strtrim(PerPatType.EpilepsyType));
isBad = ismember(etype, badTypes);
validPatients = PerPatType.Patient(~isBad);
R = R(ismember(R.patient_id, validPatients), :);
% Drop rows still missing per-row type
if ~isstring(R.epilepsy_type), R.epilepsy_type = string(R.epilepsy_type); end
R = R(~(ismissing(R.epilepsy_type) | strlength(strtrim(R.epilepsy_type))==0), :);

%% ===== Build per-visit table & per-patient MeanSzFreq =====
% Collect all visits across rows for a patient, dedup by date
if ~isstring(R.sz_freqs),           R.sz_freqs = string(R.sz_freqs); end
if ~isstring(R.visit_dates_deid),   R.visit_dates_deid = string(R.visit_dates_deid); end
if ~ismember('visit_hasSz', R.Properties.VariableNames)
    R.visit_hasSz = strings(height(R),1);
elseif ~isstring(R.visit_hasSz)
    R.visit_hasSz = string(R.visit_hasSz);
end

PV = table('Size',[0 4], 'VariableTypes',{'double','datetime','double','double'}, ...
           'VariableNames',{'Patient','VisitDate','Freq','HasSz'});

for j = 1:height(R)
    pid = double(R.patient_id(j));

    % -- Dates
    dates_raw = strtrim(R.visit_dates_deid(j));
    d = datetime.empty(0,1);
    if strlength(dates_raw)>0 && dates_raw~="[]" && dates_raw~=""
        try, dd = jsondecode(char(dates_raw)); dd = string(dd(:));
        catch, dd = string(regexp(dates_raw,'\d{4}-\d{2}-\d{2}','match')); end
        try, d = datetime(dd,'InputFormat','yyyy-MM-dd'); catch, d=datetime(dd); end
    end

    % -- Freqs
    freq_raw = strtrim(R.sz_freqs(j));
    v = [];
    if strlength(freq_raw)>0 && freq_raw~="[]" && freq_raw~=""
        s = regexprep(freq_raw, 'null', 'NaN', 'ignorecase');
        try, v = jsondecode(char(s)); v = double(v(:));
        catch
            nums = regexp(s,'[-+]?\d+(\.\d+)?([eE][-+]?\d+)?','match');
            if ~isempty(nums), v = double(str2double(string(nums(:)))); end
        end
    end
    if isempty(v), v = nan(0,1); end
    % sanitize: nonfinite->NaN, negatives->NaN
    v(~isfinite(v)) = NaN; v(v<0) = NaN;

    % -- HasSz
    has_raw = strtrim(R.visit_hasSz(j));
    h = nan(0,1);
    if strlength(has_raw)>0 && has_raw~="[]" && has_raw~=""
        try, h = jsondecode(char(has_raw)); h = double(h(:));
        catch
            numsH = regexp(has_raw,'\d+','match');
            if ~isempty(numsH), h = double(str2double(string(numsH(:)))); end
        end
    end
    if isempty(h), h = nan(0,1); end

    % Align
    nA = min([numel(d), numel(v), numel(h)]);
    if nA==0, continue; end
    PV = [PV; table(repmat(pid,nA,1), d(1:nA), v(1:nA), h(1:nA), ...
            'VariableNames', PV.Properties.VariableNames)]; %#ok<AGROW>
end

% Dedup per (Patient, VisitDate)
[gv, pidk, datek] = findgroups(PV.Patient, PV.VisitDate);
Freq_agg = splitapply(@(x) mean(x(isfinite(x)),'omitnan'), PV.Freq, gv);
Has_agg  = splitapply(@(x) max(x(isfinite(x))), PV.HasSz, gv);
Vuniq = table(pidk, datek, Freq_agg, Has_agg, ...
    'VariableNames', {'Patient','VisitDate','Freq','HasSz'});

% Substitute zeros where freq is missing but HasSz==0 (your chosen rule)
fillZero = ~isfinite(Vuniq.Freq) & (Vuniq.HasSz==0);
Vuniq.Freq(fillZero) = 0;

% Per-patient mean seizure frequency
[gpV, pidsV] = findgroups(Vuniq.Patient);
MeanSzFreq = splitapply(@(x) mean(x,'omitnan'), Vuniq.Freq, gpV);
Rg = table(pidsV, MeanSzFreq, 'VariableNames', {'Patient','MeanSzFreq'});

% Guard small negatives
tinyNegSz = Rg.MeanSzFreq < 0 & Rg.MeanSzFreq > -1e-12;
Rg.MeanSzFreq(tinyNegSz) = 0;
Rg.MeanSzFreq(Rg.MeanSzFreq < -1e-12) = NaN;

%% ===== Per-patient mean spike rate (spikes/min) =====
[gpS, pidsS] = findgroups(S.Patient);
MeanSpikeRate_perMin = splitapply(@(x) mean(x,'omitnan'), S.SpikeRate_perMin, gpS);
Sg = table(pidsS, MeanSpikeRate_perMin, 'VariableNames', {'Patient','MeanSpikeRate_perMin'});

%% ===== Merge and split groups =====
P = innerjoin(Sg, Rg, 'Keys','Patient');
P = P(isfinite(P.MeanSpikeRate_perMin) & isfinite(P.MeanSzFreq), :);

grpZero  = (P.MeanSzFreq==0);
grpNon   = (P.MeanSzFreq>0);

Y_zero   = log10(max(P.MeanSpikeRate_perMin(grpZero), EPS_PER_MIN));
Y_non    = log10(max(P.MeanSpikeRate_perMin(grpNon),  EPS_PER_MIN));

n0 = numel(Y_zero);
n1 = numel(Y_non);

p_rs = NaN;
if n0>=3 && n1>=3
    p_rs = ranksum(P.MeanSpikeRate_perMin(grpZero), P.MeanSpikeRate_perMin(grpNon), 'method','approx');
end

m0 = median(P.MeanSpikeRate_perMin(grpZero), 'omitnan');
m1 = median(P.MeanSpikeRate_perMin(grpNon),  'omitnan');
iqr0 = iqr(P.MeanSpikeRate_perMin(grpZero));
iqr1 = iqr(P.MeanSpikeRate_perMin(grpNon));

%% ===== Plot =====
f = figure('Color','w','Position',[80 80 720 540]);
ax = axes(f); hold(ax,'on'); box(ax,'off'); grid(ax,'on');

G = [repmat("MeanSzFreq = 0", n0, 1); repmat("MeanSzFreq > 0", n1, 1)];
Y = [Y_zero; Y_non];

if ~isempty(Y)
    boxchart(ax, categorical(G), Y, 'BoxFaceAlpha',0.25);
    try
        swarmchart(ax, categorical(G), Y, 18, 'filled','MarkerFaceAlpha',0.45);
    catch
        scatter(ax, double(categorical(G)) + 0.05*randn(size(Y)), Y, 18, 'filled','MarkerFaceAlpha',0.45);
    end
end

yline(ax, log10(EPS_PER_MIN), ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
ylim(ax, Y_LIMS);
ylabel(ax, 'log_{10}(spikes/min)');
title(ax, sprintf('Mean spike rate by seizure-frequency group — %s', segLabel), 'FontSize', 18);

% significance bar
if isfinite(p_rs)
    add_sigbar(ax, 1, 2, Y_LIMS(2) - 0.08*range(Y_LIMS), sprintf('p = %.3g', p_rs));
else
    text(ax, 0.5,0.97,'Insufficient N','Units','normalized','HorizontalAlignment','center','FontSize',16);
end
set(ax,'FontSize',16);

% Save
if ~exist(fileparts(figOut),'dir'), mkdir(fileparts(figOut)); end
exportgraphics(f, figOut, 'Resolution', 300);
fprintf('Saved figure: %s\n', figOut);

%% ===== Console summary =====
fprintf('\n=== Mean spike rate (spikes/min) by MeanSzFreq group — %s ===\n', segLabel);
fprintf('MeanSzFreq = 0   : n=%d  median=%.3g  IQR=%.3g\n', n0, m0, iqr0);
fprintf('MeanSzFreq > 0   : n=%d  median=%.3g  IQR=%.3g\n', n1, m1, iqr1);
if isfinite(p_rs), fprintf('Wilcoxon rank-sum p=%.3g\n', p_rs); else, fprintf('Insufficient N for rank-sum\n'); end

end % function

%% ===== Helper =====
function add_sigbar(ax, x1, x2, y, ptext)
tick = 0.03 * diff(ax.YLim);
plot(ax, [x1 x1 x2 x2], [y- tick, y, y, y- tick], 'k-', 'LineWidth', 1.3);
text(ax, mean([x1 x2]), y + 0.02*diff(ax.YLim), ptext, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16);
end
