%% ===== Config =====
mainDir   = '/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/output/spikenet/';
outCsv    = '../data/SN_counts/spike_counts_summary_multiThresh.csv';
OVERWRITE = false;   % TRUE = rebuild from scratch; FALSE = append/skip existing

FS        = 16;          % Hz (0.0625 s per sample)
WIN_SEC   = 1.0;
WIN_SAMP  = round(WIN_SEC * FS);
HALF_WIN  = floor(WIN_SAMP/2);

%% ===== Threshold list =====
% Base thresholds: 0.05, 0.10, ..., 0.95
baseThresh = 0.05:0.05:0.95;

% Add the original threshold explicitly
extraThresh = 0.43;

% Combine and ensure unique, sorted list
threshList = unique([baseThresh, extraThresh]);

% Number of thresholds
nThresh = numel(threshList);

% Create valid MATLAB variable names
threshVarNames = strings(1, nThresh);
for iT = 1:nThresh
    thrStr = sprintf('%.2f', threshList(iT));   % e.g. '0.43'
    thrStr(thrStr == '.') = '_';               % becomes '0_43'
    threshVarNames(iT) = "count_" + thrStr;    % 'count_0_43'
end


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
varNames = [{ ...
    'EEG_Name','Patient','Session','Duration_sec' ...
    }, cellstr(threshVarNames)];

varTypes = [{'string','double','double','double'}, repmat({'double'}, 1, nThresh)];

Summary  = table('Size',[0 numel(varNames)], ...
                 'VariableTypes', varTypes, ...
                 'VariableNames', varNames);

%% ===== Load existing (and build donePairs by Patient-Session) =====
donePairs = string.empty(0,1);  % Column-empty so vertical concatenation is safe

if ~OVERWRITE && isfile(outCsv)
    try
        Existing = readtable(outCsv, 'TextType','string');

        % Add missing columns for backward compatibility
        for vn = varNames
            if ~ismember(vn{1}, Existing.Properties.VariableNames)
                switch vn{1}
                    case 'EEG_Name'
                        Existing.(vn{1}) = strings(height(Existing),1);
                    otherwise
                        Existing.(vn{1}) = nan(height(Existing),1);
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
            fprintf('Note: found %d duplicate (Patient,Session) rows; keeping first instances.\n', ...
                    height(Existing) - numel(keepIdx));
            Existing = Existing(keepIdx, :);
        end

        Summary   = Existing;
        donePairs = strcat(string(Summary.Patient), "_", string(Summary.Session));
        donePairs = donePairs(:);  % ensure column vector
        fprintf('Loaded existing summary (%d rows). Using (Patient,Session) to skip.\n', ...
                height(Summary));
    catch ME
        warning('Could not read existing summary (%s). Proceeding as empty.\n%s', ...
                outCsv, ME.message);
    end
else
    if OVERWRITE
        fprintf('OVERWRITE=true: rebuilding summary from scratch.\n');
    else
        fprintf('No existing summary found. Starting a new one.\n');
    end
end

%% ===== Process each CSV =====
for k = 1:20%numel(files)
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
        fprintf('Skipping (Patient=%g, Session=%g) already summarized. (%s)\n', ...
                patientNum, sessionNum, fname);
        continue
    end

    % --- Read CSV ---
    try
        T = readtable(fpath, 'TextType','string');
    catch ME
        warning('Failed to read %s. Skipping. (%s)', fname, ME.message);
        continue
    end

    % Only need SN2 now
    required = {'SN2'};
    if ~all(ismember(required, T.Properties.VariableNames))
        warning('Skipping %s: required column SN2 not found.', fname);
        continue
    end

    SN2 = double(T.SN2(:));

    % Handle NaNs defensively: treat NaN as below threshold
    SN2(~isfinite(SN2)) = -inf;

    N = numel(SN2);
    duration_sec = N / FS;

    % ===== Whole-file clustered detections for all thresholds =====
    counts = zeros(1, nThresh);
    for iT = 1:nThresh
        counts(iT) = count_spikes_wholefile(SN2, threshList(iT), WIN_SAMP);
    end

    % ===== Append summary row =====
    newRow = [ ...
        {string(fname), patientNum, sessionNum, duration_sec}, ...
        num2cell(counts) ...
        ];
    Summary = [Summary; cell2table(newRow, 'VariableNames', varNames)]; %#ok<AGROW>

    % Track the pair as done so duplicates in the same run are skipped
    donePairs(end+1,1) = pairKey; %#ok<AGROW>

    % Nicely formatted printout (show a few thresholds)
    fprintf('Processed (Patient=%g, Session=%g) %-40s  dur=%.1fs  count@0.10=%.0f  count@0.50=%.0f  count@0.90=%.0f\n', ...
        patientNum, sessionNum, fname, duration_sec, ...
        counts(threshList==0.10), counts(threshList==0.50), counts(threshList==0.90));
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
function totalCount = count_spikes_wholefile(SN2, THRESHOLD, WIN_SAMP)
% Count clustered detections over the WHOLE file for a given threshold.
% Uses the same 1-second skip rule as before, but without left/right labels.

    totalCount = 0;
    Nloc = numel(SN2);
    i = 1;
    half_win = floor(WIN_SAMP/2);

    while i <= Nloc
        if SN2(i) > THRESHOLD
            % Define a local window (for future extensions; currently we just
            % treat any above-threshold sample as a detection)
            wStart = max(1, i - half_win); %#ok<NASGU>
            wEnd   = min(Nloc, i + half_win); %#ok<NASGU>

            totalCount = totalCount + 1;
            i = i + WIN_SAMP;  % skip ahead ~1 second
        else
            i = i + 1;
        end
    end
end
