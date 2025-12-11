%% ===== Config =====
mainDir   = '/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/output/spikenet/';
outCsv    = '../data/SN_counts/spike_counts_summary.csv';
OVERWRITE = false;   % TRUE = rebuild from scratch; FALSE = append/skip existing

THRESHOLD = 0.43;
FS        = 16;          % Hz (0.0625 s per sample)
WIN_SEC   = 1.0;
WIN_SAMP  = round(WIN_SEC * FS);
HALF_WIN  = floor(WIN_SAMP/2);

%% ===== Find all probability CSVs (two patterns) =====
f1 = dir(fullfile(mainDir, 'sub-Penn*', 'ses-*', 'sub-Penn*_ses-*_task-EEG_eeg_*.csv'));
f2 = dir(fullfile(mainDir, 'sub-Penn*', 'ses-*', 'sub-Penn*_ses-*_combined_*.csv'));
filesAll = [f1; f2];

if isempty(filesAll)
    fprintf('No CSVs found under: %s\n', mainDir);
end
filePaths = string(fullfile({filesAll.folder}, {filesAll.name}));
[~, uniqIdx] = unique(filePaths, 'stable');
files = filesAll(uniqIdx);

%% ===== Schema =====
varNames = { ...
    'EEG_Name','Patient','Session', ...
    'Total_Spikes','Left_Spikes','Right_Spikes','Duration_sec','Longest_Spike_Duration_sec', ...
    'FirstRun_Spikes','FirstRun_Left_Spikes','FirstRun_Right_Spikes','FirstRun_Duration_sec', ...
    'First24Runs_Spikes','First24Runs_Left_Spikes','First24Runs_Right_Spikes','First24Runs_Duration_sec' ...
    };
varTypes = {'string','double','double', ...
    'double','double','double','double','double', ...
    'double','double','double','double', ...
    'double','double','double','double'};

Summary  = table('Size',[0 numel(varNames)], 'VariableTypes', varTypes, 'VariableNames', varNames);

%% ===== Load existing (and build donePairs by Patient-Session) =====
% Column-empty so vertical concatenation is safe
donePairs = string.empty(0,1);

if ~OVERWRITE && isfile(outCsv)
    try
        Existing = readtable(outCsv, 'TextType','string');

        % Add missing columns for backward compatibility
        for vn = varNames
            if ~ismember(vn{1}, Existing.Properties.VariableNames)
                switch vn{1}
                    case 'EEG_Name', Existing.(vn{1}) = strings(height(Existing),1);
                    otherwise,       Existing.(vn{1}) = nan(height(Existing),1);
                end
            end
        end

        % Coerce types and order
        Existing = Existing(:, varNames);
        Existing.EEG_Name = string(Existing.EEG_Name);
        for vn = varNames(2:end)
            Existing.(vn{1}) = double(Existing.(vn{1}));
        end

        % De-duplicate by (Patient, Session) if any accidental duplicates exist
        key = strcat(string(Existing.Patient), "_", string(Existing.Session));
        [~, keepIdx] = unique(key, 'stable');
        if numel(keepIdx) < height(Existing)
            fprintf('Note: found %d duplicate (Patient,Session) rows; keeping first instances.\n', height(Existing) - numel(keepIdx));
            Existing = Existing(keepIdx, :);
        end

        Summary   = Existing;
        donePairs = strcat(string(Summary.Patient), "_", string(Summary.Session));
        donePairs = donePairs(:);  % ensure column vector
        fprintf('Loaded existing summary (%d rows). Using (Patient,Session) to skip.\n', height(Summary));
    catch ME
        warning('Could not read existing summary (%s). Proceeding as empty.\n%s', outCsv, ME.message);
    end
else
    if OVERWRITE
        fprintf('OVERWRITE=true: rebuilding summary from scratch.\n');
    else
        fprintf('No existing summary found. Starting a new one.\n');
    end
end

%% ===== Process each CSV =====
for k = 1:numel(files)
    fpath  = fullfile(files(k).folder, files(k).name);
    fname  = files(k).name;

    % --- Parse Patient/Session robustly from the PATH ---
    tok = regexp(fpath, 'sub-Penn(\d+).+?ses-(\d+)', 'tokens', 'once');
    patientNum = NaN; sessionNum = NaN;
    if ~isempty(tok)
        patientNum = str2double(tok{1});
        sessionNum = str2double(tok{2});
    else
        tok2 = regexp(fname, 'sub-Penn(\d+)_ses-(\d+)', 'tokens', 'once');
        if ~isempty(tok2)
            patientNum = str2double(tok2{1});
            sessionNum = str2double(tok2{2});
        end
    end

    if ~isfinite(patientNum) || ~isfinite(sessionNum)
        warning('Skipping %s (could not parse Patient/Session).', fname);
        continue
    end

    pairKey = strcat(string(patientNum), "_", string(sessionNum));

    % --- Skip if this (Patient,Session) already summarized and not overwriting ---
    if ~OVERWRITE && any(donePairs == pairKey)
        fprintf('Skipping (Patient=%g, Session=%g) already summarized. (%s)\n', patientNum, sessionNum, fname);
        continue
    end

    % --- Read CSV ---
    try
        T = readtable(fpath, 'TextType','string');
    catch ME
        warning('Failed to read %s. Skipping. (%s)', fname, ME.message);
        continue
    end

    required = {'SN2','SN2_zero_right','SN2_zero_left'};
    if ~all(ismember(required, T.Properties.VariableNames))
        warning('Skipping %s: required columns not found.', fname);
        continue
    end

    SN2 = double(T.SN2(:));
    ZR  = double(T.SN2_zero_right(:));
    ZL  = double(T.SN2_zero_left(:));

    % Handle NaNs defensively: treat NaN as below threshold
    SN2(~isfinite(SN2)) = -inf;
    ZR(~isfinite(ZR))   = -inf;
    ZL(~isfinite(ZL))   = -inf;

    N = numel(SN2);
    duration_sec = N / FS;

    % ===== Longest continuous above-threshold duration (whole file) =====
    mask = SN2 > THRESHOLD;
    if any(mask)
        d = diff([false; mask; false]);   % rising/falling edges
        run_starts  = find(d ==  1);
        run_ends    = find(d == -1) - 1;
        run_lengths = run_ends - run_starts + 1; % samples
        longest_run_sec = max(run_lengths) / FS;
    else
        longest_run_sec = 0;
    end

    % ===== Whole-file clustered detections (existing behavior) =====
    [totalCount, leftCount, rightCount] = count_spikes_segment(SN2, ZR, ZL, 1:N, THRESHOLD, WIN_SAMP);

    % ===== Identify a run/record column if present =====
    runVar = "";
    candidates = {'run','Run','RUN','run_idx','run_index','RunIndex','record','Record','REC','run_number','RunNumber'};
    for c = 1:numel(candidates)
        if ismember(candidates{c}, T.Properties.VariableNames)
            runVar = candidates{c};
            break
        end
    end

    % ===== First-run + First-24-runs metrics =====
    firstRun_total = totalCount; firstRun_L = leftCount; firstRun_R = rightCount; firstRun_dur = duration_sec;
    first24_total  = totalCount; first24_L  = leftCount; first24_R  = rightCount; first24_dur  = duration_sec;

    if runVar ~= ""
        runNums = double(T.(runVar)(:));
        runNums(~isfinite(runNums)) = nan;

        % Unique runs in file order
        [uniqRuns, ~, ~] = unique(runNums, 'stable');
        uniqRuns = uniqRuns(isfinite(uniqRuns));

        if ~isempty(uniqRuns)
            % ----- First run only -----
            r1 = uniqRuns(1);
            idxFirst = find(runNums == r1);
            [firstRun_total, firstRun_L, firstRun_R] = count_spikes_segment(SN2, ZR, ZL, idxFirst, THRESHOLD, WIN_SAMP);
            firstRun_dur = numel(idxFirst) / FS;

            % ----- First 24 runs (or entire file if <24 runs, per request) -----
            if numel(uniqRuns) >= 24
                takeRuns = uniqRuns(1:24);
                first24_total = 0; first24_L = 0; first24_R = 0; first24_dur = 0;
                for r = 1:numel(takeRuns)
                    idxr = find(runNums == takeRuns(r));
                    [t_, L_, R_] = count_spikes_segment(SN2, ZR, ZL, idxr, THRESHOLD, WIN_SAMP);
                    first24_total = first24_total + t_;
                    first24_L     = first24_L + L_;
                    first24_R     = first24_R + R_;
                    first24_dur   = first24_dur + numel(idxr)/FS;
                end
            else
                % Fewer than 24 runs => use *full-file* stats
                first24_total = totalCount;
                first24_L     = leftCount;
                first24_R     = rightCount;
                first24_dur   = duration_sec;
            end
        end
    end

    % ===== Append summary row =====
    newRow = { ...
        string(fname), patientNum, sessionNum, ...
        totalCount, leftCount, rightCount, duration_sec, longest_run_sec, ...
        firstRun_total, firstRun_L, firstRun_R, firstRun_dur, ...
        first24_total, first24_L, first24_R, first24_dur ...
        };
    Summary = [Summary; cell2table(newRow, 'VariableNames', varNames)]; %#ok<AGROW>

    % Track the pair as done so duplicates in the same run are skipped
    donePairs(end+1,1) = pairKey; %#ok<AGROW>

    fprintf(['Processed (Patient=%g, Session=%g) %-40s  total=%4d  L=%4d  R=%4d  dur=%.1fs  longest=%.3fs  | ' ...
             '1stRun: total=%4d L=%4d R=%4d dur=%.1fs  | ' ...
             '1st24: total=%4d L=%4d R=%4d dur=%.1fs\n'], ...
        patientNum, sessionNum, fname, ...
        totalCount, leftCount, rightCount, duration_sec, longest_run_sec, ...
        firstRun_total, firstRun_L, firstRun_R, firstRun_dur, ...
        first24_total, first24_L, first24_R, first24_dur);
end

%% ===== Optional: sort for readability =====
if height(Summary) > 0
    try
        Summary = sortrows(Summary, {'Patient','Session','EEG_Name'});
    catch
        % ignore sort errors
    end
end

%% ===== Save summary =====
writetable(Summary, outCsv);
fprintf('Saved summary to: %s  (rows=%d)\n', outCsv, height(Summary));

%% ===== Local function =====
function [totalCount, leftCount, rightCount] = count_spikes_segment(SN2, ZR, ZL, idxs, THRESHOLD, WIN_SAMP)
% Count clustered detections within an index subset (idxs).
% This applies the same 1-second skip rule, restricted to the provided indices.
    totalCount = 0; leftCount = 0; rightCount = 0;
    if isempty(idxs); return; end

    % Work in the local index space
    sn = SN2(idxs);
    zr = ZR(idxs);
    zl = ZL(idxs);

    Nloc = numel(sn);
    i = 1;
    half_win = floor(WIN_SAMP/2);

    while i <= Nloc
        if sn(i) > THRESHOLD
            wStart = max(1, i - half_win);
            wEnd   = min(Nloc, i + half_win);

            maxZR = max(zr(wStart:wEnd));
            maxZL = max(zl(wStart:wEnd));

            if maxZR > maxZL
                leftCount  = leftCount + 1;
            else
                rightCount = rightCount + 1;
            end
            totalCount = totalCount + 1;
            i = i + WIN_SAMP;  % skip ahead within this segment
        else
            i = i + 1;
        end
    end
end
