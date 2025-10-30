%% ===== Config =====
mainDir   = '/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/output/spikenet/';   % <- set this to your top-level folder
outCsv    = '../data/SN_counts/spike_counts_summary.csv';
OVERWRITE = false;   % <--- set TRUE to rebuild from scratch; FALSE to append/skip existing

THRESHOLD = 0.43;
FS        = 16;          % Hz (0.0625 s per sample)
WIN_SEC   = 1.0;         % cluster window ~1 second
WIN_SAMP  = round(WIN_SEC * FS);
HALF_WIN  = floor(WIN_SAMP/2);

%% ===== Find all probability CSVs in nested structure =====
files = dir(fullfile(mainDir, 'sub-Penn*', 'ses-*', 'sub-Penn*_ses-*_task-EEG_eeg_*.csv'));
if isempty(files)
    fprintf('No CSVs found under: %s\n', mainDir);
end

%% ===== Prepare (load existing summary if not overwriting) =====
varNames = {'EEG_Name','Patient','Session','Total_Spikes','Left_Spikes','Right_Spikes','Duration_sec','Longest_Spike_Duration_sec'};
varTypes = {'string','double','double','double','double','double','double','double'};

Summary = table('Size',[0 numel(varNames)], 'VariableTypes', varTypes, 'VariableNames', varNames);
doneSet = string.empty(1,0);

if ~OVERWRITE && isfile(outCsv)
    try
        Existing = readtable(outCsv, 'TextType','string');
        % Ensure required columns exist (add if missing for backward compatibility)
        for vn = varNames
            if ~ismember(vn{1}, Existing.Properties.VariableNames)
                % Add missing var with sensible default
                switch vn{1}
                    case 'EEG_Name', Existing.(vn{1}) = strings(height(Existing),1);
                    otherwise,       Existing.(vn{1}) = nan(height(Existing),1);
                end
            end
        end
        % Reorder to standard schema and coerce types
        Existing = Existing(:, varNames);
        % Coerce types (important if older CSV stored chars/cells)
        Existing.EEG_Name = string(Existing.EEG_Name);
        for vn = varNames(2:end)
            Existing.(vn{1}) = double(Existing.(vn{1}));
        end

        Summary = Existing;
        doneSet = unique(Summary.EEG_Name);
        fprintf('Loaded existing summary (%d rows). Will skip already processed EEGs.\n', height(Summary));
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
    fpath = fullfile(files(k).folder, files(k).name);
    fname = files(k).name;
    eegName = string(fname);

    % Skip if already processed and not overwriting
    if ~OVERWRITE && any(doneSet == eegName)
        fprintf('Skipping already summarized: %s\n', fname);
        continue
    end

    % Parse patient/session
    patientNum = NaN; sessionNum = NaN;
    tok = regexp(fname, 'sub-Penn(\d+)_ses-(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        patientNum = str2double(tok{1});
        sessionNum = str2double(tok{2});
    else
        tok2 = regexp(fpath, 'sub-Penn(\d+)[/\\]ses-(\d+)', 'tokens', 'once');
        if ~isempty(tok2)
            patientNum = str2double(tok2{1});
            sessionNum = str2double(tok2{2});
        end
    end

    % Read CSV
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

    % ===== Longest continuous above-threshold duration (independent of clustering) =====
    mask = SN2 > THRESHOLD;
    if any(mask)
        d = diff([false; mask; false]);         % rising/falling edges
        run_starts  = find(d ==  1);
        run_ends    = find(d == -1) - 1;
        run_lengths = run_ends - run_starts + 1; % in samples
        longest_run_sec = max(run_lengths) / FS;
    else
        longest_run_sec = 0;
    end

    % ===== Cluster detections (your original logic) =====
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

    % Append summary row
    newRow = {eegName, patientNum, sessionNum, totalCount, leftCount, rightCount, duration_sec, longest_run_sec};
    Summary = [Summary; cell2table(newRow, 'VariableNames', varNames)]; %#ok<AGROW>

    fprintf('Processed %-60s  total=%4d  L=%4d  R=%4d  dur=%.1fs  longest=%.3fs\n', ...
        fname, totalCount, leftCount, rightCount, duration_sec, longest_run_sec);
end

%% ===== Optional: sort by Patient/Session for readability =====
if height(Summary) > 0
    try
        Summary = sortrows(Summary, {'Patient','Session','EEG_Name'});
    catch
        % ignore sort errors if types differ
    end
end

%% ===== Save summary =====
writetable(Summary, outCsv);
fprintf('Saved summary to: %s  (rows=%d)\n', outCsv, height(Summary));
