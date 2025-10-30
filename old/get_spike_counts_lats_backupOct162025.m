%% ===== Config =====
mainDir  = '/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/output/spikenet/';   % <- set this to your top-level folder
outCsv   = '../data/SN_counts/spike_counts_summary.csv';
%mainDir = '../data/SN_examples/';

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

Summary = table('Size',[0 7], 'VariableTypes', ...
    {'string','double','double','double','double','double','double'}, ...
    'VariableNames', {'EEG_Name','Patient','Session','Total_Spikes','Left_Spikes','Right_Spikes','Duration_sec'});

%% ===== Process each CSV =====
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    fname = files(k).name;

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
    T = readtable(fpath, 'TextType','string');
    required = {'SN2','SN2_zero_right','SN2_zero_left'};
    if ~all(ismember(required, T.Properties.VariableNames))
        warning('Skipping %s: required columns not found.', fname);
        continue
    end

    SN2 = double(T.SN2(:));
    ZR  = double(T.SN2_zero_right(:));
    ZL  = double(T.SN2_zero_left(:));

    N = numel(SN2);
    duration_sec = N / FS;

    i = max(1, HALF_WIN + 1);
    totalCount = 0; leftCount = 0; rightCount = 0;

    % Cluster detections
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
    Summary = [Summary; {string(fname), patientNum, sessionNum, totalCount, leftCount, rightCount, duration_sec}]; %#ok<AGROW>
    fprintf('Processed %-60s  total=%4d  L=%4d  R=%4d  dur=%.1fs\n', ...
        fname, totalCount, leftCount, rightCount, duration_sec);
end

%% ===== Save summary =====
writetable(Summary, outCsv);
fprintf('Saved summary to: %s\n', outCsv);
