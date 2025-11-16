%% ============================================================
%  Spearman correlations: spike rate vs seizure frequency
%  Display: log10(spikes/min) vs log10(seizures/month)
%  Panels:
%    A) All epilepsy (overall) — r, p (unadjusted)
%    B) Frontal — r, p_bonf    C) Temporal — r, p_bonf    D) General — r, p_bonf
%  Style: log10–log10 scatter, dotted zero guards, zero-origin stars
%  Axes fixed across panels: X [-3.5, 4], Y [-3, 2]
%  ------------------------------------------------------------

% clear; clc;

%% ===== Config =====
only_amb   = 0; % 1 = ONLY ambulatories; 2 = NO ambulatories (routine only); 0 = everything
which_runs = 0; % 0 = whole file, 1 = first run (~1h), 2 = first 24 runs (~24h)
only_outpt = 0; % 1 = keep only EEGs with acquired_on containing "spe" or "radnor" (case-insensitive)

% Paths
outCsv     = '../../data/SN_counts/spike_counts_summary.csv';
reportFile = '../../data/Routineeegpec-Deidreport_DATA_LABELS_2025-11-10_1443.csv';
outFigure  = '../../figures/spearman_spikerate_szfreq_log10_min_month_fixedaxes.png';
resultsCsv = '../../output/spearman_spikerate_szfreq_log10_min_month_fixedaxes.csv';   % [] to skip saving

%% ===== Load spike counts; choose segment by which_runs =====
S = readtable(outCsv,'TextType','string','VariableNamingRule','preserve');

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

% Ambulatory/routine filter via FULL duration
if only_amb == 1
    S(S.Duration_sec < 3600*12,:) = [];
elseif only_amb == 2
    S(S.Duration_sec > 3600*12,:) = [];
end

% Keys
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end

% Rates
if ~isnumeric(S.(colCount)), S.(colCount) = double(S.(colCount)); end
if ~isnumeric(S.(colDur)),   S.(colDur)   = double(S.(colDur));   end
S.SpikeRate_perHour = nan(height(S),1);
validDur = S.(colDur) > 0;
S.SpikeRate_perHour(validDur) = (S.(colCount)(validDur) ./ S.(colDur)(validDur)) * 3600;
S.SpikeRate_perMin = nan(height(S),1);
S.SpikeRate_perMin(validDur) = S.SpikeRate_perHour(validDur) / 60;

% Guard: spikes should never be negative
tinyNeg = S.SpikeRate_perMin < 0 & S.SpikeRate_perMin > -1e-12;
S.SpikeRate_perMin(tinyNeg) = 0;
S.SpikeRate_perMin(S.SpikeRate_perMin < -1e-12) = NaN;

%% ===== Load report =====
R = readtable(reportFile,'TextType','string','VariableNamingRule','preserve');
if ~isnumeric(R.patient_id),     R.patient_id     = double(str2double(R.patient_id)); end
if ~isnumeric(R.session_number), R.session_number = double(str2double(R.session_number)); end
if ~isstring(R.epilepsy_type),   R.epilepsy_type  = string(R.epilepsy_type); end
assert(ismember('sz_freqs',R.Properties.VariableNames),'Column "sz_freqs" not found.');

%% ===== OUTPATIENT/session alignment (mirror filter onto S and R) =====
if ~ismember('acquired_on', R.Properties.VariableNames)
    if only_outpt==1
        warning('Report lacks "acquired_on"; only_outpt flag ignored.');
    end
else
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

    % Keep only S columns plus matched acquired_on (for audit)
    keepVars = intersect(S.Properties.VariableNames, SJ.Properties.VariableNames, 'stable');
    keepVars = unique([keepVars, {'SpikeRate_perHour','SpikeRate_perMin','acquired_on'}], 'stable');
    S = SJ(:, keepVars);
end

% Mirror the (Patient,Session) set onto R (ensures seizure freq uses same cohort)
RSess_keep = unique(S(:, {'Patient','Session'}));
RSess_keep.Properties.VariableNames = {'Patient','Session'};
RSess_keep.Patient = double(RSess_keep.Patient);
RSess_keep.Session = double(RSess_keep.Session);

R = innerjoin(R, RSess_keep, ...
    'LeftKeys', {'patient_id','session_number'}, ...
    'RightKeys', {'Patient','Session'});

%% ===== Exclude non-epilepsy types (patient-level), and drop missing per-row types =====
badTypes = lower(["Non-Epileptic Seizure Disorder","Uncertain if Epilepsy","Unknown or MRN not found",""]);
isEmptyType  = ismissing(R.epilepsy_type) | strlength(strtrim(R.epilepsy_type))==0;
Rt = sortrows(R(~isEmptyType, {'patient_id','epilepsy_type'}), 'patient_id');
[uniq_pid, ia] = unique(Rt.patient_id, 'stable');
PerPatType = table(uniq_pid, Rt.epilepsy_type(ia), 'VariableNames', {'Patient','EpilepsyType'});

epi_norm = lower(strtrim(PerPatType.EpilepsyType));
isBad    = ismember(epi_norm, badTypes);
validPatients = PerPatType.Patient(~isBad);

% Keep only valid epilepsy patients
R = R(ismember(R.patient_id, validPatients), :);

% Drop rows with missing/empty epilepsy_type
if ~isstring(R.epilepsy_type), R.epilepsy_type = string(R.epilepsy_type); end
maskMissingType = ismissing(R.epilepsy_type) | strlength(strtrim(R.epilepsy_type))==0;
R = R(~maskMissingType, :);

fprintf('Total patients in report (post filters): %d\n', numel(unique(R.patient_id(~isnan(R.patient_id)))));
fprintf('Patients with ANY epilepsy type kept:    %d\n', numel(unique(R.patient_id)));
fprintf('Excluded by epilepsy_type:               %d\n\n', sum(isBad));

%% ===== Build patient-level visit table (dedup by date) and MeanSzFreq =====
% This avoids double-counting visits when a patient has multiple EEG rows.
if ~isstring(R.sz_freqs),           R.sz_freqs           = string(R.sz_freqs);           end
if ~isstring(R.visit_dates_deid),   R.visit_dates_deid   = string(R.visit_dates_deid);   end
if ~ismember('visit_hasSz', R.Properties.VariableNames)
    R.visit_hasSz = strings(height(R),1);
elseif ~isstring(R.visit_hasSz)
    R.visit_hasSz = string(R.visit_hasSz);
end

PV = table('Size',[0 4], 'VariableTypes',{'double','datetime','double','double'}, ...
           'VariableNames',{'Patient','VisitDate','Freq','HasSz'});

negDrop_cnt = 0; % diagnostics

for j = 1:height(R)
    pid = double(R.patient_id(j));

    % ---- Dates ----
    dates_raw = strtrim(R.visit_dates_deid(j));
    d = datetime.empty(0,1);
    if strlength(dates_raw)>0 && dates_raw~="[]" && dates_raw~=""
        try
            dd = jsondecode(char(dates_raw)); dd = string(dd(:));
        catch
            dd = string(regexp(dates_raw, '\d{4}-\d{2}-\d{2}', 'match'));
        end
        try
            d = datetime(dd, 'InputFormat','yyyy-MM-dd');
        catch
            d = datetime(dd);
        end
    end

    % ---- Frequencies ----
    freq_raw = strtrim(R.sz_freqs(j));
    v = [];
    if strlength(freq_raw)>0 && freq_raw~="[]" && freq_raw~=""
        s = regexprep(freq_raw, 'null', 'NaN', 'ignorecase');
        try
            v = jsondecode(char(s)); v = double(v(:));
        catch
            nums = regexp(s,'[-+]?\d+(\.\d+)?([eE][-+]?\d+)?','match');
            if ~isempty(nums), v = double(str2double(string(nums(:)))); end
        end
    end
    if isempty(v), v = nan(0,1); end

    % *** SANITIZE visit freqs: nonfinite -> NaN; negatives -> NaN
    if ~isempty(v)
        negDrop_cnt = negDrop_cnt + nnz(isfinite(v) & v < 0);
        v(~isfinite(v)) = NaN;
        v(v < 0)        = NaN;
    end

    % ---- HasSz ----
    has_raw = strtrim(R.visit_hasSz(j));
    h = nan(0,1);
    if strlength(has_raw)>0 && has_raw~="[]" && has_raw~=""
        try
            h = jsondecode(char(has_raw)); h = double(h(:));
        catch
            numsH = regexp(has_raw,'\d+','match');
            if ~isempty(numsH), h = double(str2double(string(numsH(:)))); end
        end
    end
    if isempty(h), h = nan(0,1); end

    % Align lengths
    nAlign = min([numel(d), numel(v), numel(h)]);
    if nAlign==0, continue; end
    d = d(1:nAlign); v = v(1:nAlign); h = h(1:nAlign);

    PV = [PV; table(repmat(pid,nAlign,1), d, v, h, ...
           'VariableNames', PV.Properties.VariableNames)]; %#ok<AGROW>
end

% Deduplicate per (Patient, VisitDate)
[gv, pid_keys, date_keys] = findgroups(PV.Patient, PV.VisitDate);
Freq_agg = splitapply(@(x) mean(x(isfinite(x)),'omitnan'), PV.Freq, gv);
Has_agg  = splitapply(@(x) max(x(isfinite(x))), PV.HasSz, gv); % max of finite flags
Vuniq = table(pid_keys, date_keys, Freq_agg, Has_agg, ...
              'VariableNames', {'Patient','VisitDate','Freq','HasSz'});

% *** SANITIZE after dedup: any negative remnants -> NaN
negMaskV = Vuniq.Freq < 0;
if any(negMaskV), Vuniq.Freq(negMaskV) = NaN; end

% Your original rule: substitute zeros where freq is missing but HasSz==0
fillZero = ~isfinite(Vuniq.Freq) & (Vuniq.HasSz==0);
Vuniq.Freq(fillZero) = 0;

% Per-patient mean seizure frequency
[gpV, pidsV] = findgroups(Vuniq.Patient);
MeanSzFreq = splitapply(@(x) mean(x,'omitnan'), Vuniq.Freq, gpV);
%MeanSzFreq = splitapply(@(x) median(x,'omitnan'), Vuniq.Freq, gpV);
Rg = table(pidsV, MeanSzFreq, 'VariableNames', {'Patient','MeanSzFreq'});

% *** Guard at patient level: tiny negative noise -> 0; true negatives -> NaN
tinyNegSz = Rg.MeanSzFreq < 0 & Rg.MeanSzFreq > -1e-12;
Rg.MeanSzFreq(tinyNegSz) = 0;
Rg.MeanSzFreq(Rg.MeanSzFreq < -1e-12) = NaN;

fprintf('Visits with negative sz_freq dropped during parsing: %d\n', negDrop_cnt);

%% ===== Attach EpilepsySpecific and map to 3 groups =====
if ismember('epilepsy_specific', R.Properties.VariableNames)
    if ~isstring(R.epilepsy_specific), R.epilepsy_specific = string(R.epilepsy_specific); end
    spec_ok = ~ismissing(R.epilepsy_specific) & strlength(strtrim(R.epilepsy_specific))>0;
    RtSpec  = sortrows(R(spec_ok, {'patient_id','epilepsy_specific'}), 'patient_id');
    [uniq_pid_s, ia_s] = unique(RtSpec.patient_id, 'stable');
    PerPatSpecific = table(uniq_pid_s, RtSpec.epilepsy_specific(ia_s), ...
        'VariableNames', {'Patient','EpilepsySpecific'});
    Rg = outerjoin(Rg, PerPatSpecific, 'Keys','Patient', 'MergeKeys', true);
else
    Rg.EpilepsySpecific = strings(height(Rg),1);
end

% Normalize names
if ~isstring(Rg.EpilepsySpecific), Rg.EpilepsySpecific = string(Rg.EpilepsySpecific); end
Rg.EpilepsySpecific = strtrim(Rg.EpilepsySpecific);

% Bring in per-patient epilepsy_type for 'General'
Rg = outerjoin(Rg, PerPatType(~isBad,:), 'Keys','Patient','MergeKeys',true);
if ~isstring(Rg.EpilepsyType), Rg.EpilepsyType = string(Rg.EpilepsyType); end
Rg.EpilepsyType = strtrim(Rg.EpilepsyType);

spec_norm = lower(Rg.EpilepsySpecific);
type_norm = lower(Rg.EpilepsyType);
isTemporal    = contains(spec_norm,"temporal");
isFrontal     = contains(spec_norm,"frontal");
isGeneralType = contains(type_norm,"general");

SpecCanon = strings(size(spec_norm));
SpecCanon(isTemporal) = "Temporal Lobe";
SpecCanon(isFrontal)  = "Frontal Lobe";
SpecCanon((strlength(SpecCanon)==0) & isGeneralType) = "General";
Rg.EpilepsySpecific = SpecCanon;

%% ===== Merge with per-patient spike rate; keep complete cases =====
[grpP, pidsP] = findgroups(S.Patient);
MeanSpikeRate_perMin = splitapply(@(x) mean(x,'omitnan'), S.SpikeRate_perMin, grpP);
%MeanSpikeRate_perMin = splitapply(@(x) median(x,'omitnan'), S.SpikeRate_perMin, grpP);
Sg = table(pidsP, MeanSpikeRate_perMin, 'VariableNames', {'Patient','MeanSpikeRate_perMin'});

P_all = innerjoin(Sg, Rg, 'Keys','Patient');
P_all = P_all(isfinite(P_all.MeanSzFreq) & isfinite(P_all.MeanSpikeRate_perMin), :);

keep3 = ismember(P_all.EpilepsySpecific, ["General","Temporal Lobe","Frontal Lobe"]);
P     = P_all(keep3, :);

fprintf('Segment: %s\n', segLabel);
fprintf('Patients with both measures (all epilepsy): %d\n', height(P_all));
fprintf('Patients in modeled 3 groups: %d\n\n', height(P));

%% ===== Spearman correlations =====
x_all = P_all.MeanSpikeRate_perMin;  % spikes/min
y_all = P_all.MeanSzFreq;            % seizures/month
mask_all = isfinite(x_all) & isfinite(y_all);
[rs_all, p_all] = corr(x_all(mask_all), y_all(mask_all), 'Type','Spearman','Rows','complete');
n_all = sum(mask_all);

groupsWanted = ["Frontal Lobe","Temporal Lobe","General"];
rows = {};
for g = groupsWanted
    m = (P.EpilepsySpecific == g);
    x = P.MeanSpikeRate_perMin(m);
    y = P.MeanSzFreq(m);
    mask = isfinite(x) & isfinite(y);
    n = sum(mask);
    if n >= 3
        [rs, p] = corr(x(mask), y(mask), 'Type','Spearman','Rows','complete');
    else
        rs = NaN; p = NaN;
    end
    rows(end+1,:) = {char(g), n, rs, p}; %#ok<SAGROW>
end
Results = cell2table(rows, 'VariableNames', {'Group','N','Spearman_r','p_raw'});
k = height(Results);
Results.p_bonf = min(Results.p_raw * k, 1);

disp('=== Spearman correlations: SpikeRate_perMin vs MeanSzFreq (per month) ===');
disp([table("Overall (all epilepsy)", n_all, rs_all, p_all, NaN, ...
      'VariableNames', {'Group','N','Spearman_r','p_raw','p_bonf'}); Results])

%% ===== Plotting prep (log10) =====
% Epsilon set from positive mins, avoiding log of <=0
minpos_rate = min(P_all.MeanSpikeRate_perMin(P_all.MeanSpikeRate_perMin>0)); if isempty(minpos_rate), minpos_rate=1e-6; end
minpos_sz   = min(P_all.MeanSzFreq(P_all.MeanSzFreq>0));                     if isempty(minpos_sz),   minpos_sz=1e-6; end
eps_rate = 0.5*minpos_rate;
eps_sz   = 0.5*minpos_sz;

Tall = table;
Tall.logSpikeRate = log10(P_all.MeanSpikeRate_perMin + (P_all.MeanSpikeRate_perMin<=0).*eps_rate);
Tall.logSzFreq    = log10(P_all.MeanSzFreq           + (P_all.MeanSzFreq<=0).*eps_sz);

T = table;
T.EpilepsySpecific = categorical(P.EpilepsySpecific);
T.EpilepsySpecific = reordercats(removecats(T.EpilepsySpecific), {'Frontal Lobe','Temporal Lobe','General'});
T.logSpikeRate = log10(P.MeanSpikeRate_perMin + (P.MeanSpikeRate_perMin<=0).*eps_rate);
T.logSzFreq    = log10(P.MeanSzFreq           + (P.MeanSzFreq<=0).*eps_sz);
ok = isfinite(T.logSpikeRate) & isfinite(T.logSzFreq) & ~ismissing(T.EpilepsySpecific);
T = T(ok,:);
present = categories(T.EpilepsySpecific);

% Quick sanity: confirm logs are real
if any(~isreal(Tall.logSzFreq)) || any(~isreal(Tall.logSpikeRate))
    warning('Complex values detected pre-plot (Tall).');
end
if any(~isreal(T.logSzFreq)) || any(~isreal(T.logSpikeRate))
    warning('Complex values detected pre-plot (T).');
end

% Markers for zero locations (pre-log)
isZeroSz_T   = (P.MeanSzFreq==0);              isZeroSz_T   = isZeroSz_T(ok);
isZeroRate_T = (P.MeanSpikeRate_perMin==0);    isZeroRate_T = isZeroRate_T(ok);

isZeroSz_all   = (P_all.MeanSzFreq==0);
isZeroRate_all = (P_all.MeanSpikeRate_perMin==0);
onlySz_all     =  isZeroSz_all & ~isZeroRate_all;
onlyRate_all   = ~isZeroSz_all &  isZeroRate_all;
bothZero_all   =  isZeroSz_all &  isZeroRate_all;
nonZero_all    = ~(isZeroSz_all | isZeroRate_all);

xZero = log10(eps_sz);
yZero = log10(eps_rate);

%% ===== GLOBAL AXIS LIMITS (consistent across all panels) =====
xLims = [-3.5, 4];     % log10(seizures/month)
yLims = [-3, 2];       % log10(spikes/min)
xrange = diff(xLims);
yrange = diff(yLims);
dx = 0.008 * xrange;   % small nudge for marking true zero cells
dy = 0.008 * yrange;

%% ===== Figure (2x2) =====
f = figure('Color','w','Position',[60 60 1200 900]);
tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');
pal   = lines(3);
fontL = 20;

% A. Overall
axA = nexttile(1); hold(axA,'on'); grid(axA,'on'); box(axA,'off');
xline(axA, xZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
yline(axA, yZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

baseColor = [0.4 0.4 0.4];
scatter(axA, Tall.logSzFreq(nonZero_all), Tall.logSpikeRate(nonZero_all), ...
        14, baseColor, 'filled', 'MarkerFaceAlpha', 0.25);
plot(axA, Tall.logSzFreq(onlySz_all)+dx,   Tall.logSpikeRate(onlySz_all),   '*','Color',baseColor,'MarkerSize',7,'LineWidth',1);
plot(axA, Tall.logSzFreq(onlyRate_all),    Tall.logSpikeRate(onlyRate_all)+dy,'*','Color',baseColor,'MarkerSize',7,'LineWidth',1);
plot(axA, Tall.logSzFreq(bothZero_all)+dx, Tall.logSpikeRate(bothZero_all)+dy,'*','Color',baseColor,'MarkerSize',8,'LineWidth',1.2);

X = [ones(sum(nonZero_all),1), Tall.logSzFreq(nonZero_all)];
b = X \ Tall.logSpikeRate(nonZero_all);
xgrid = linspace(xLims(1), xLims(2), 300)';
plot(axA, xgrid, b(1) + b(2)*xgrid, 'k-', 'LineWidth', 2);

xlim(axA, xLims); ylim(axA, yLims);
xl = xlim(axA); yl = ylim(axA);
text(axA, xZero + 0.03*(xl(2)-xl(1)), mean(yl), 'Seizures/month = 0', ...
     'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
     'FontSize',fontL-2, 'Color',[0.25 0.25 0.25]);
text(axA, mean(xl), yZero + 0.06*(yl(2)-yl(1)), 'Spikes/min = 0', ...
     'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
     'FontSize',fontL-2, 'Color',[0.25 0.25 0.25]);

xlabel(axA,'log_{10} Seizures per month','FontSize',fontL);
ylabel(axA,'log_{10} Spikes per minute','FontSize',fontL);
title(axA, sprintf('A. All epilepsy  (N=%d) — %s', height(Tall), segLabel), 'FontSize',fontL, 'FontWeight','bold');
txtA = sprintf('Spearman r=%.3f, p=%.3g', rs_all, p_all);
text(axA, 0.98, 0.95, txtA, 'Units','normalized', 'HorizontalAlignment','right', 'VerticalAlignment','top', 'FontSize',fontL-2, 'FontWeight','bold');
set(axA,'FontSize',fontL);

% B/C/D
order      = {'Frontal Lobe','Temporal Lobe','General'};
panelTitle = {'B. Frontal','C. Temporal','D. General'};

for p = 1:3
    ax = nexttile(p+1); hold(ax,'on'); grid(ax,'on'); box(ax,'off');

    if ~ismember(order{p}, present)
        axis(ax,'off'); continue;
    end

    g   = order{p};
    gi  = find(strcmp(present, g), 1);
    col = pal(gi,:);
    m   = (T.EpilepsySpecific == g);

    onlySz_g   =  m &  isZeroSz_T & ~isZeroRate_T;
    onlyRate_g =  m & ~isZeroSz_T &  isZeroRate_T;
    bothZero_g =  m &  isZeroSz_T &  isZeroRate_T;
    nonZero_g  =  m & ~(isZeroSz_T | isZeroRate_T);

    xline(ax, xZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
    yline(ax, yZero, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    scatter(ax, T.logSzFreq(nonZero_g), T.logSpikeRate(nonZero_g), 18, col, 'filled', 'MarkerFaceAlpha', 0.35);
    if any(onlySz_g),   plot(ax, T.logSzFreq(onlySz_g)+dx,   T.logSpikeRate(onlySz_g),   '*','Color',col,'MarkerSize',8,'LineWidth',1.1); end
    if any(onlyRate_g), plot(ax, T.logSzFreq(onlyRate_g),    T.logSpikeRate(onlyRate_g)+dy,'*','Color',col,'MarkerSize',8,'LineWidth',1.1); end
    if any(bothZero_g), plot(ax, T.logSzFreq(bothZero_g)+dx, T.logSpikeRate(bothZero_g)+dy,'*','Color',col,'MarkerSize',9,'LineWidth',1.2); end

    if nnz(nonZero_g) >= 3
        Xg = [ones(nnz(nonZero_g),1), T.logSzFreq(nonZero_g)];
        bg = Xg \ T.logSpikeRate(nonZero_g);
        xg = linspace(xLims(1), xLims(2), 250)';
        plot(ax, xg, bg(1)+bg(2)*xg, '-', 'Color', col, 'LineWidth', 2);
    end

    xlim(ax, xLims); ylim(ax, yLims);
    xl = xlim(ax); yl = ylim(ax);
    text(ax, xZero + 0.03*(xl(2)-xl(1)), mean(yl), 'Seizures/month = 0', ...
         'Rotation',90, 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
         'FontSize',fontL-2, 'Color',[0.25 0.25 0.25]);
    text(ax, mean(xl), yZero + 0.06*(yl(2)-yl(1)), 'Spikes/min = 0', ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
         'FontSize',fontL-2, 'Color',[0.25 0.25 0.25]);

    % Panel titles & stats
    row = Results(strcmp(Results.Group, string(g)), :);
    if ~isempty(row) && row.N >= 3 && isfinite(row.Spearman_r)
        txt = sprintf('Spearman r=%.3f, p_{bonf}=%.3g', row.Spearman_r, row.p_bonf);
    else
        txt = 'Insufficient data';
    end

    nNow = sum(T.EpilepsySpecific==g);
    title(ax, sprintf('%s (N=%d)', panelTitle{p}, nNow), 'FontSize', fontL, 'FontWeight','bold');
    text(ax, 0.98, 0.95, txt, 'Units','normalized', ...
         'HorizontalAlignment','right', 'VerticalAlignment','top', ...
         'FontSize', fontL-3, 'FontWeight','bold');

    set(ax,'FontSize',fontL);
end

%% ===== Save figure =====
if ~isempty(outFigure)
    if ~exist(fileparts(outFigure),'dir'), mkdir(fileparts(outFigure)); end
    exportgraphics(f, outFigure, 'Resolution', 300);
    fprintf('Saved figure: %s\n', outFigure);
end

%% ===== Save results (optional) =====
if ~isempty(resultsCsv)
    if ~exist(fileparts(resultsCsv),'dir'), mkdir(fileparts(resultsCsv)); end
    ResultsOut = [table("Overall (all epilepsy)", n_all, rs_all, p_all, NaN, ...
        'VariableNames', {'Group','N','Spearman_r','p_raw','p_bonf'}); Results];
    writetable(ResultsOut, resultsCsv);
    fprintf('Saved results to: %s\n', resultsCsv);
end
