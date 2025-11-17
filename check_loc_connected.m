%% check_loc_connected.m
% Scans a BIDS EEG tree:
%   sub-Penn[Patient_ID]/ses-[SESSION_ID]/eeg/sub-Penn[Patient_ID]_ses-[SESSION_ID]_task-EEG_eeg.edf
%
% For each EDF:
%   1) Reads header via edfinfo
%   2) Checks whether BOTH labels "LOC" and "ROC" are present in info.SignalLabels
%      (case-insensitive, trimmed).
%   3) If either label is missing:
%        - outpt_elecs = 0
%        - loc_ratio, roc_ratio = NaN
%      and moves on (no edfread).
%   4) If BOTH labels present:
%        - outpt_elecs = 1
%        - Loads data via edfread
%        - Extracts LOC and ROC channels
%        - Computes:
%             fs       = info.NumSamples(1);
%             p60      = bandpower(chan, fs, [58 62]);
%             pall     = bandpower(chan, fs, [0.5 70]);
%             ratio    = p60 / pall;
%          for both LOC and ROC → loc_ratio, roc_ratio
%
% Output CSV: outpt_info.csv with columns:
%   filename, patient_id, session_id, outpt_elecs, loc_ratio, roc_ratio

%% ----------------------- CONFIG -----------------------
% Root directory containing sub-Penn*/ses-* folders
rootDir = "/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/data/routine_eegs/";
%rootDir = '/Users/erinconrad/Desktop/research/edf_pipeline/test_edf';

outCsv  = fullfile("/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/output/", ...
                   "outpt_info.csv");

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
outpt_flags = zeros(n,1);     % 0/1: both LOC & ROC labels present?
loc_ratio   = nan(n,1);       % p60/pall for LOC
roc_ratio   = nan(n,1);       % p60/pall for ROC

%% -------------------- MAIN LOOP -----------------------
needLabels = ["LOC","ROC"];  % required labels (case-insensitive, trimmed)

for k = 1:n
    f = edfFiles(k);
    fpath = string(fullfile(f.folder, f.name));
    filenames(k) = fpath;

    % --- Extract patient_id and session_id from path (robust to minor naming variance) ---
    parts = split(string(f.folder), filesep);
    % Look for parts that start with "sub-Penn" and "ses-"
    subIdx = find(startsWith(parts, "sub-Penn"), 1, 'last');
    sesIdx = find(startsWith(parts, "ses-"), 1, 'last');

    if ~isempty(subIdx), patient_ids(k) = extractAfter(parts(subIdx), "sub-"); end
    if ~isempty(sesIdx), session_ids(k) = extractAfter(parts(sesIdx), "ses-"); end

    % If parsing from folder failed, try filename regex as a fallback
    if patient_ids(k) == "" || session_ids(k) == ""
        % Example name: sub-Penn1147_ses-2_task-EEG_eeg.edf
        tok = regexp(f.name, ...
            "sub-(?<pid>Penn[0-9]+)_ses-(?<sid>[^_]+)_task-EEG_eeg\.edf$", ...
            "names");
        if ~isempty(tok)
            if patient_ids(k) == "", patient_ids(k) = string(tok.pid); end
            if session_ids(k) == "", session_ids(k) = string(tok.sid); end
        end
    end

    % --- Read EDF header and test for labels ---
    try
        info = edfinfo(fpath);

        labels = string(info.SignalLabels); % [N x 1 string]
        labels = upper(strtrim(labels));    % normalize
        hasBoth = all(ismember(needLabels, labels));

        if ~hasBoth
            % Missing LOC or ROC → outpt_elecs = 0, ratios stay NaN, skip EDF data
            outpt_flags(k) = 0;
            continue;
        end

        % If we get here, BOTH LOC & ROC labels are present in the header
        outpt_flags(k) = 1;

        % ------------------ LOAD DATA & COMPUTE RATIOS ------------------
        try
            TT = edfread(fpath);  % timetable

            % Variable names in timetable may differ in case; match case-insensitively
            varNames = TT.Properties.VariableNames;
            locIdx   = find(strcmpi(varNames, "LOC"), 1);
            rocIdx   = find(strcmpi(varNames, "ROC"), 1);

            if isempty(locIdx) || isempty(rocIdx)
                warning("edfread missing LOC/ROC variable(s) for:\n  %s\nLeaving ratios as NaN.", fpath);
                continue;
            end

            locVarName = varNames{locIdx};
            rocVarName = varNames{rocIdx};

            loc = TT.(locVarName);
            roc = TT.(rocVarName);

            % If variables are cell arrays, convert to numeric
            if iscell(loc), loc = cell2mat(loc); end
            if iscell(roc), roc = cell2mat(roc); end

            % Sampling rate (per your spec)
            fs = info.NumSamples(1);

            % LOC ratios
            p60_loc  = bandpower(loc, fs, [58 62]);
            pall_loc = bandpower(loc, fs, [0.5 70]);
            loc_ratio(k) = p60_loc / pall_loc;

            % ROC ratios
            p60_roc  = bandpower(roc, fs, [58 62]);
            pall_roc = bandpower(roc, fs, [0.5 70]);
            roc_ratio(k) = p60_roc / pall_roc;

        catch ME2
            warning("edfread/bandpower failed for:\n  %s\nReason: %s", fpath, ME2.message);
            % outpt_elecs remains 1 (labels present), ratios stay NaN
        end

    catch ME
        % If header read fails, mark as 0 and leave ratios NaN
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
    loc_ratio, ...
    roc_ratio, ...
    'VariableNames', {'filename','patient_id','session_id','outpt_elecs','loc_ratio','roc_ratio'} ...
);

writetable(T, outCsv);
fprintf('Wrote %d rows to: %s\n', height(T), outCsv);

% Optional: quick summary in command window
nBoth = sum(T.outpt_elecs == 1);
nMiss = sum(T.outpt_elecs == 0);
fprintf('Header has both LOC & ROC: %d  |  Missing either: %d\n', nBoth, nMiss);
