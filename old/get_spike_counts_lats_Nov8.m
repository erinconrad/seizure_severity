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

% If you prefer to catch *any* candidate CSV under ses-* (use with care):
% f3 = dir(fullfile(mainDir, 'sub-Penn*', 'ses-*', '*.csv'));
% filesAll = [filesAll; f3];

if isempty(filesAll)
    fprintf('No CSVs found under: %s\n', mainDir);
end
filePaths = string(fullfile({filesAll.folder}, {filesAll.name}));
[~, uniqIdx] = unique(filePaths, 'stable');
files = filesAll(uniqIdx);

%% ===== Schema =====
varNames = {'EEG_Name','Patient','Session','Total_Spikes','Left_Spikes','Right_Spikes','Duration_sec','Longest_Spike_Duration_sec'};
varTypes = {'string','double','double','double','double','double','double','double'};
Summary  = table('Size',[0 numel(varNames)], 'VariableTypes', varTypes, 'VariableNames', varNames);

%% ===== Load existing (and build donePairs by Patient-Session) =====
donePairs = string.empty(1,0);

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
    % Matches: .../sub-Penn1147/.../ses-2/...
    tok = regexp(fpath, 'sub-Penn(\d+).+?ses-(\d+)', 'tokens', 'once');
    patientNum = NaN; sessionNum = NaN;
    if ~isempty(tok)
        patientNum = str2double(tok{1});
        sessionNum = str2double(tok{2});
    else
        % Fallback to filename pattern (older style)
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

    % ===== Longest continuous above-threshold duration =====
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

    % ===== Cluster detections =====
    i = max(1, HALF_WIN + 1);
    totalCount = 0; leftCount = 0; rightCount = 0;

    while i <= N
        if SN2(i) > THRESHOLD
            wStart = max(1, i - HALF_WIN);
            wEnd   = min(N, i + HALF_WIN);

            maxZR = max(ZR(wStart:wEnd));
            maxZL = max(ZL(wStart:wEnd));

            if maxZR > maxZL
                leftCount  = leftCount + 1;
            else
                rightCount = rightCount + 1;
            end
            totalCount = totalCount + 1;
            i = i + WIN_SAMP;   % skip ahead
        else
            i = i + 1;
        end
    end

    % ===== Append summary row =====
    newRow = {string(fname), patientNum, sessionNum, totalCount, leftCount, rightCount, duration_sec, longest_run_sec};
    Summary = [Summary; cell2table(newRow, 'VariableNames', varNames)]; %#ok<AGROW>

    % Track the pair as done so duplicates in the same run are skipped
    donePairs = [donePairs; pairKey]; %#ok<AGROW>

    fprintf('Processed (Patient=%g, Session=%g) %-40s  total=%4d  L=%4d  R=%4d  dur=%.1fs  longest=%.3fs\n', ...
        patientNum, sessionNum, fname, totalCount, leftCount, rightCount, duration_sec, longest_run_sec);
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
