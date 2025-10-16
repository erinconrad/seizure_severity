function within_sz_paired_analysis
%% ======================= CONFIG =======================
outCsv      = '../data/SN_counts/spike_counts_summary.csv';                     % spike counts with Duration_sec
reportFile  = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-07_1611.csv';
pairsOut    = '../output/eeg_visit_pairs.csv';                                  % audit output (built pairs)

maxGapDays        = 365;     % visit must be within this many days of EEG (≈ 6 months; use 365 for a year)
NEG_POLICY        = 'nan';   % how to handle negative sz_freq: 'nan' (drop) or 'zero' (set to 0)
removeConstSzFreq = true;    % drop patients whose sz frequency never changes across their pairs
plotLogRate       = false;   % if true, plot log1p(spike rate); test is run on raw rates either way

%% ================== LOAD & PREP SPIKE SUMMARY ==================
S = readtable(outCsv, 'TextType','string', 'VariableNamingRule','preserve');  
reqS = {'EEG_Name','Patient','Session','Total_Spikes','Left_Spikes','Right_Spikes','Duration_sec'};
assert(all(ismember(reqS, S.Properties.VariableNames)), 'Spike summary missing required columns.');
if ~isnumeric(S.Patient), S.Patient = double(str2double(string(S.Patient))); end
if ~isnumeric(S.Session), S.Session = double(str2double(string(S.Session))); end
S.SpikeRate_perHour = S.Total_Spikes ./ (S.Duration_sec/3600);

%% ================== LOAD & PREP REPORT SHEET ===================
R = readtable(reportFile, 'TextType','string', 'VariableNamingRule','preserve');
reqR = {'patient_id','session_number','start_time_deid','visit_dates_deid','sz_freqs'};
assert(all(ismember(reqR, R.Properties.VariableNames)), ...
    'Report file missing required columns. Found: %s', strjoin(R.Properties.VariableNames, ', '));
if ~isnumeric(R.('patient_id')),     R.('patient_id')     = double(str2double(string(R.('patient_id')))); end
if ~isnumeric(R.('session_number')), R.('session_number') = double(str2double(string(R.('session_number')))); end

% EEG datetime (handles "M/d/yy HH:mm" and falls back to autodetect)
rawStart = string(R.('start_time_deid'));
okStart  = ~ismissing(rawStart) & strlength(strtrim(rawStart)) > 0;
R.EEG_Date = NaT(height(R),1);
try
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart)), 'InputFormat','M/d/yy HH:mm');
catch
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart))); % autodetect
end

% Parse visit dates (ISO) and sz freqs (JSON with possible "null")
rawVisits = string(R.('visit_dates_deid'));
rawSz     = string(R.('sz_freqs'));
R.VisitDates = cell(height(R),1);
R.SzFreqs    = cell(height(R),1);

for i = 1:height(R)
    % Visit dates
    vstr = strtrim(rawVisits(i));
    if strlength(vstr) > 0 && vstr ~= "[]"
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
    % Sz freqs
    sstr = strtrim(rawSz(i));
    if strlength(sstr) > 0 && sstr ~= "[]"
        sfix = regexprep(sstr, 'null', 'NaN', 'ignorecase');
        try
            vv = jsondecode(char(sfix)); 
            R.SzFreqs{i} = vv(:);
        catch
            nums = regexp(sfix, '[-+]?\d+(\.\d+)?([eE][-+]?\d+)?', 'match');
            R.SzFreqs{i} = str2double(string(nums));
        end
    else
        R.SzFreqs{i} = nan(0,1);
    end
end

%% ====== JOIN spike summary with report to get EEG date & visits ======
JR = innerjoin(S, R, ...
    'LeftKeys',  {'Patient','Session'}, ...
    'RightKeys', {'patient_id','session_number'}, ...
    'RightVariables', {'EEG_Date','VisitDates','SzFreqs'});

% Keep rows with EEG date and at least one valid visit date
JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) && all(~isnat(v)), JR.VisitDates), :);

%% ====== BUILD EEG–VISIT PAIRS (closest visit within maxGapDays) ======
pairs = table('Size',[0 9], ...
    'VariableTypes', {'double','double','string','datetime','double','datetime','double','double','double'}, ...
    'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour','Visit_Date','Sz_Freq','DaysGap','EEG_Index'});

for i = 1:height(JR)
    eegDate = JR.EEG_Date(i);
    vdates  = JR.VisitDates{i};
    szf     = JR.SzFreqs{i};

    if isempty(vdates) || isempty(szf), continue; end
    nAlign = min(numel(vdates), numel(szf));
    vdates = vdates(1:nAlign);
    szf    = szf(1:nAlign);

    [dmin, idx] = min(abs(days(vdates - eegDate)));
    if ~isempty(dmin) && dmin <= maxGapDays && ~isnan(szf(idx))
        pairs = [pairs; {JR.Patient(i), JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perHour(i), ...
                         vdates(idx), szf(idx), dmin, i}]; %#ok<AGROW>
    end
end

fprintf('Built %d EEG–visit pairs (within %d days).\n', height(pairs), maxGapDays);
if isempty(pairs), warning('No pairs constructed.'); return; end

% ≥2 pairs per patient
counts   = groupcounts(pairs, 'Patient');
eligible = counts.Patient(counts.GroupCount >= 2);
pairs    = pairs(ismember(pairs.Patient, eligible), :);
fprintf('Eligible patients with ≥2 pairs: %d\n', numel(eligible));
if isempty(pairs), warning('Nothing remains after ≥2-pairs filter.'); return; end

% Remove patients whose EEGs all matched the same visit date
[Gv, pidv] = findgroups(pairs.Patient);
nUniqueVisits = splitapply(@(d) numel(unique(d)), pairs.Visit_Date, Gv);
keepPatients = pidv(nUniqueVisits >= 2);
pairs = pairs(ismember(pairs.Patient, keepPatients), :);
if isempty(pairs), warning('Nothing left after visit-date variation filter.'); return; end

%% ===== Clean sz freq (negatives) and drop non-finite =====
pairs.SpikeRate_perHour = double(pairs.SpikeRate_perHour);
pairs.Sz_Freq           = double(pairs.Sz_Freq);
negMask = isfinite(pairs.Sz_Freq) & pairs.Sz_Freq < 0;
switch lower(NEG_POLICY)
    case 'nan',  pairs.Sz_Freq(negMask) = NaN;
    case 'zero', pairs.Sz_Freq(negMask) = 0;
    otherwise,   error('Unknown NEG_POLICY: %s', NEG_POLICY);
end
pairs = pairs(isfinite(pairs.SpikeRate_perHour) & isfinite(pairs.Sz_Freq), :);
if isempty(pairs), warning('All pairs dropped after sz-freq cleaning.'); return; end

% Optionally remove patients with constant sz frequency (no within-patient variation)
if removeConstSzFreq
    [Gc, pids] = findgroups(pairs.Patient);
    varSz = splitapply(@(v) nanvar(v), pairs.Sz_Freq, Gc);
    keepP2 = pids(varSz > 0);
    pairs = pairs(ismember(pairs.Patient, keepP2), :);
    if isempty(pairs), warning('Nothing left after constant sz-freq filter.'); return; end
end

%% ====== BUILD PER-PATIENT LOW vs HIGH PAIRS (tie-break by closest gap) ======
[Gp, pid] = findgroups(pairs.Patient);

lowRates   = nan(numel(pid),1);
highRates  = nan(numel(pid),1);
lowSz      = nan(numel(pid),1);
highSz     = nan(numel(pid),1);
lowInfo    = strings(numel(pid),1);
highInfo   = strings(numel(pid),1);

for i = 1:numel(pid)
    rows = find(Gp == i);
    szf  = pairs.Sz_Freq(rows);
    rate = pairs.SpikeRate_perHour(rows);
    gap  = pairs.DaysGap(rows);

    if numel(rows) < 2 || all(szf == szf(1))
        continue; % skip constant sz-freq or only one pair
    end

    % LOWEST sz freq (if ties, pick smallest DaysGap; if still ties, first)
    minSz = min(szf);
    idxL  = find(szf == minSz);
    if numel(idxL) > 1
        [~, j] = min(gap(idxL));
        idxL = idxL(j);
    end

    % HIGHEST sz freq (if ties, pick smallest DaysGap)
    maxSz = max(szf);
    idxH  = find(szf == maxSz);
    if numel(idxH) > 1
        [~, j] = min(gap(idxH));
        idxH = idxH(j);
    end

    lowRates(i)  = rate(idxL);
    highRates(i) = rate(idxH);
    lowSz(i)     = szf(idxL);
    highSz(i)    = szf(idxH);

    lowInfo(i)   = pairs.EEG_Name(rows(idxL));
    highInfo(i)  = pairs.EEG_Name(rows(idxH));
end

% Keep informative patients
ok = isfinite(lowRates) & isfinite(highRates) & (lowSz ~= highSz);
lowRates  = lowRates(ok);
highRates = highRates(ok);
lowSz     = lowSz(ok);
highSz    = highSz(ok);

fprintf('Paired analysis: %d patients (low vs high sz-freq).\n', numel(lowRates));
if isempty(lowRates)
    warning('No informative patients for paired analysis.'); 
    return
end

%% ====== WILCOXON SIGNED-RANK (high vs low spike rate) ======
[pp, h, stats] = signrank(highRates, lowRates);   % nonparametric within-patient
zval = stats.zval;
rEff = zval / sqrt(numel(lowRates));              % effect size r

fprintf('Wilcoxon signed-rank (high vs low spike rate): z = %.3f, p = %.3g, r = %.3f\n', zval, pp, rEff);
fprintf('Median low rate = %.3f, median high rate = %.3f (spikes/hour)\n', median(lowRates), median(highRates));

%% ====== PAIRED PLOT ======
figure('Color','w'); hold on;
if plotLogRate
    yLow  = log1p(lowRates);
    yHigh = log1p(highRates);
    ylab  = 'log(1 + spike rate), spikes/hour';
else
    yLow  = lowRates;
    yHigh = highRates;
    ylab  = 'Spike rate (spikes/hour)';
end

for i = 1:numel(yLow)
    plot([1 2], [yLow(i) yHigh(i)], '-o', 'Color',[0.6 0.6 0.6], 'MarkerFaceColor',[0.6 0.6 0.6]);
end
xlim([0.5 2.5]);
set(gca, 'XTick',[1 2], 'XTickLabel',{'Lowest sz freq','Highest sz freq'});
ylabel(ylab);
title(sprintf('Paired spike rates by seizure frequency (n=%d), p=%.3g', numel(yLow), pp));
grid on; box off;

%% ====== SAVE PAIRED DATA (optional) & Pairs audit ======
writetable(pairs, pairsOut);
fprintf('Saved EEG–visit pairs to: %s\n', pairsOut);

% Optional: save the low/high pair summary
% pairedOut = '../output/paired_low_vs_high.csv';
% Tpaired = table(lowRates, highRates, lowSz, highSz, 'VariableNames', ...
%    {'LowSpikeRate_perHour','HighSpikeRate_perHour','LowSzFreq','HighSzFreq'});
% writetable(Tpaired, pairedOut);
% fprintf('Saved paired summary to: %s\n', pairedOut);

end
