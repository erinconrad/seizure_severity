%% check_outpt_labels.m
% Scans a BIDS EEG tree:
%   sub-Penn[Patient_ID]/ses-[SESSION_ID]/eeg/sub-Penn[Patient_ID]_ses-[SESSION_ID]_task-EEG_eeg.edf
% For each EDF, reads header via edfinfo and tests whether both labels
% "LOC" and "ROC" are present in info.SignalLabels (exact label match after
% trimming + uppercasing).
%
% Output CSV: outpt_label.csv with columns:
%   filename, patient_id, session_id, outpt_elecs (0/1)

%% ----------------------- CONFIG -----------------------
% Point this to the directory that contains the sub-Penn*/ses-* folders
rootDir = "/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/data/routine_eegs/";   % <-- EDIT ME

outCsv  = fullfile("/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/output/", "outpt_label.csv");

% File pattern (strict BIDS path/filename as you described)
filePattern = fullfile(rootDir, "sub-Penn*", "ses-*", "eeg", ...
                       "sub-Penn*_ses-*_task-EEG_eeg.edf");

%% --------------------- FIND FILES ---------------------
edfFiles = dir(filePattern);

if isempty(edfFiles)
    warning("No EDF files found using pattern:\n  %s", filePattern);
end

% Preallocate containers
n = numel(edfFiles);
filenames   = strings(n,1);
patient_ids = strings(n,1);
session_ids = strings(n,1);
outpt_flags = zeros(n,1);  % 0/1

%% -------------------- MAIN LOOP -----------------------
needLabels = ["LOC","ROC"];  % required labels (case-insensitive, trimmed)

for k = 1:10%:n
    f = edfFiles(k);
    fpath = string(fullfile(f.folder, f.name));
    filenames(k) = fpath;

    % --- Extract patient_id and session_id from path (robust to minor naming variance) ---
    % Prefer folder names "sub-PennXXXX" and "ses-YY" if present
    parts = split(string(f.folder), filesep);
    % Look for parts that start with "sub-Penn" and "ses-"
    subIdx = find(startsWith(parts, "sub-Penn"), 1, 'last');
    sesIdx = find(startsWith(parts, "ses-"), 1, 'last');

    if ~isempty(subIdx), patient_ids(k) = extractAfter(parts(subIdx), "sub-"); end
    if ~isempty(sesIdx), session_ids(k) = extractAfter(parts(sesIdx), "ses-"); end

    % If parsing from folder failed, try filename regex as a fallback
    if patient_ids(k) == "" || session_ids(k) == ""
        % Example name: sub-Penn1147_ses-2_task-EEG_eeg.edf
        tok = regexp(f.name, "sub-(?<pid>Penn[0-9]+)_ses-(?<sid>[^_]+)_task-EEG_eeg\.edf$", "names");
        if ~isempty(tok)
            if patient_ids(k) == "", patient_ids(k) = string(tok.pid); end
            if session_ids(k) == "", session_ids(k) = string(tok.sid); end
        end
    end

    % --- Read EDF header and test for labels ---
    try
        info = edfinfo(fpath);  %#ok<NASGU>  % keep for MATLAB versions without output warnings
        labels = string(info.SignalLabels);    % [N x 1 string]
        labels = upper(strtrim(labels));       % normalize
        hasBoth = all(ismember(needLabels, labels));
        outpt_flags(k) = double(hasBoth);
    catch ME
        % If header read fails, mark as 0 and report
        outpt_flags(k) = 0;
        warning("edfinfo failed for:\n  %s\nReason: %s", fpath, ME.message);
    end
end

%% -------------------- WRITE OUTPUT --------------------
T = table( ...
    filenames, ...
    patient_ids, ...
    session_ids, ...
    outpt_flags, ...
    'VariableNames', {'filename','patient_id','session_id','outpt_elecs'} ...
);

writetable(T, outCsv);
fprintf('Wrote %d rows to: %s\n', height(T), outCsv);

% Optional: quick summary in command window
nYes = sum(T.outpt_elecs == 1);
nNo  = sum(T.outpt_elecs == 0);
fprintf('Contains both LOC & ROC: %d  |  Missing either: %d\n', nYes, nNo);
