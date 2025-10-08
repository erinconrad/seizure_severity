%% ===== Config =====
mainDir  = '/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/output/spikenet/';   % <- set this to your top-level folder
outCsv   = '../data/SN_counts/spike_counts_summary.csv';

THRESHOLD = 0.43;
FS        = 16;          % Hz (0.0625 s)
WIN_SEC   = 1.0;         % cluster window ~1 second
WIN_SAMP  = round(WIN_SEC * FS);
HALF_WIN  = floor(WIN_SAMP/2);

%% ===== Find all probability CSVs in nested structure =====
% Matches: [main]/sub-PennXXX/ses-Y/sub-PennXXX_ses-Y_task-EEG_eeg_*.csv
files = dir(fullfile(mainDir, 'sub-Penn*', 'ses-*', 'sub-Penn*_ses-*_task-EEG_eeg_*.csv'));
if isempty(files)
    fprintf('No CSVs found under: %s\n', mainDir);
end

Summary = table('Size',[0 6], 'VariableTypes', ...
    {'string','double','double','double','double','double'}, ...
    'VariableNames', {'EEG_Name','Patient','Session','Total_Spikes','Left_Spikes','Right_Spikes'});

%% ===== Process each CSV =====
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    fname = files(k).name;

    % Parse patient/session from either filename or folder path
    % Works for: sub-Penn102/ses-2/sub-Penn102_ses-2_task-EEG_eeg_0_1.csv
    patientNum = NaN; sessionNum = NaN;

    tok = regexp(fname, 'sub-Penn(\d+)_ses-(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        patientNum = str2double(tok{1});
        sessionNum = str2double(tok{2});
    else
        % Fallback: parse folder names if filename format differs
        tok2 = regexp(fpath, 'sub-Penn(\d+)[/\\]ses-(\d+)', 'tokens', 'once');
        if ~isempty(tok2)
            patientNum = str2double(tok2{1});
            sessionNum = str2double(tok2{2});
        end
    end

    % Read CSV (expects columns: SN2, SN2_zero_right, SN2_zero_left)
    T = readtable(fpath, 'TextType','string');
    required = {'SN2','SN2_zero_right','SN2_zero_left'};
    if ~all(ismember(required, T.Properties.VariableNames))
        warning('Skipping %s: required columns not found.', fname);
        continue
    end

    SN2 = double(T.SN2(:));
    ZR  = double(T.SN2_zero_right(:));  % zero right -> LEFT contribution remains
    ZL  = double(T.SN2_zero_left(:));   % zero left  -> RIGHT contribution remains

    N = numel(SN2);
    i = max(1, HALF_WIN + 1);
    totalCount = 0; leftCount = 0; rightCount = 0;

    % Cluster detections: 1-s non-overlapping windows
    while i <= N
        if SN2(i) > THRESHOLD
            wStart = max(1, i - HALF_WIN);
            wEnd   = min(N, i + HALF_WIN);

            % Decide side using windowed max (matches prior logic)
            maxZR = max(ZR(wStart:wEnd));
            maxZL = max(ZL(wStart:wEnd));

            if maxZR > maxZL
                leftCount  = leftCount + 1;   % higher ZR => LEFT
            else
                rightCount = rightCount + 1;
            end
            totalCount = totalCount + 1;

            i = i + WIN_SAMP;   % skip ahead ~1 s
        else
            i = i + 1;
        end
    end

    % Append summary row
    Summary = [Summary; {string(fname), patientNum, sessionNum, totalCount, leftCount, rightCount}]; %#ok<AGROW>
    fprintf('Processed %-60s  total=%4d  L=%4d  R=%4d\n', fname, totalCount, leftCount, rightCount);
end

%% ===== Save summary =====
writetable(Summary, outCsv);
fprintf('Saved summary to: %s\n', outCsv);
