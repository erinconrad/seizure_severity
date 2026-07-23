%% count_spikes_deadtime.m
% For each SpikeNet2 probability CSV, split detections at a single threshold
% into those falling in the leading "dead time" of the EDF vs. those falling
% in the true recording period.
%
% Dead time is defined as:  deadtime_sec = (EDF duration) - (clinical duration_hms)
% and is assumed to sit at the BEGINNING of the file.

%% ===== Config =====
mainDir    = '/mnt/sauce/littlab/users/erinconr/projects/routine_eeg/output/spikenet/';
clinicalCsv = '../data/clinical_data_deidentified.csv';   % <-- EDIT: path to the clinical spreadsheet
outCsv     = '../data/SN_counts/spike_counts_deadtime_0p46.csv';

THRESH   = 0.46;
FS       = 16;                      % Hz
WIN_SEC  = 1.0;
WIN_SAMP = round(WIN_SEC * FS);

% If sub-PennNNN numbering is offset relative to patient_id, set it here
% (Patient_from_path = patient_id + PATIENT_ID_OFFSET).
PATIENT_ID_OFFSET = 0;

% Tolerance (sec) for "negative dead time" before we call it an error rather
% than rounding noise.
NEG_DEAD_TOL = 5;

% Set to Inf for the full run, or e.g. 100 for a quick test pass.
MAX_FILES     = 100;
PRINT_EVERY   = 100;   % progress update interval (files)

%% ===== Load clinical durations =====
assert(isfile(clinicalCsv), 'Clinical CSV not found: %s', clinicalCsv);
C = readtable(clinicalCsv, 'TextType', 'string', 'VariableNamingRule', 'preserve');
C.Properties.VariableNames = matlab.lang.makeValidName(C.Properties.VariableNames);

needed = {'patient_id', 'session_number', 'duration_hms'};
assert(all(ismember(needed, C.Properties.VariableNames)), ...
    'Clinical CSV is missing one of: %s', strjoin(needed, ', '));

clin_patient = double(C.patient_id) + PATIENT_ID_OFFSET;
clin_session = double(C.session_number);
clin_dur_sec = hms_to_seconds(C.duration_hms);

ok = isfinite(clin_patient) & isfinite(clin_session) & isfinite(clin_dur_sec);
assert(any(ok), 'No usable rows in clinical CSV (check duration_hms parsing).');
clin_patient = clin_patient(ok);
clin_session = clin_session(ok);
clin_dur_sec = clin_dur_sec(ok);

clinKey = string(clin_patient) + "_" + string(clin_session);
[uKey, iu] = unique(clinKey, 'stable');
if numel(uKey) < numel(clinKey)
    warning('Clinical CSV has %d duplicate (patient,session) rows; keeping first.', ...
        numel(clinKey) - numel(uKey));
end
durMap = containers.Map(cellstr(uKey), num2cell(clin_dur_sec(iu)));

%% ===== Find probability CSVs (same two patterns as before) =====
f1 = dir(fullfile(mainDir, 'sub-Penn*', 'ses-*', 'sub-Penn*_ses-*_task-EEG_eeg_*.csv'));
f2 = dir(fullfile(mainDir, 'sub-Penn*', 'ses-*', 'sub-Penn*_ses-*_combined_*.csv'));
filesAll = [f1; f2];
assert(~isempty(filesAll), 'No CSVs found under: %s', mainDir);

filePaths = string(fullfile({filesAll.folder}, {filesAll.name}));
[~, uniqIdx] = unique(filePaths, 'stable');
files = filesAll(uniqIdx);

nFound = numel(files);
if isfinite(MAX_FILES) && nFound > MAX_FILES
    files = files(1:MAX_FILES);
    fprintf('TEST RUN: processing first %d of %d files (set MAX_FILES = Inf for all).\n', ...
        numel(files), nFound);
else
    fprintf('Processing all %d files.\n', nFound);
end
nToDo = numel(files);

%% ===== Schema =====
varNames = {'EEG_Name', 'Patient', 'Session', ...
    'EDF_Duration_sec', 'Recorded_Duration_sec', 'Deadtime_sec', ...
    'n_spikes_deadtime', 'n_spikes_recording', 'n_spikes_total', ...
    'rate_deadtime_per_min', 'rate_recording_per_min', 'Flag'};
varTypes = [{'string'}, repmat({'double'}, 1, 10), {'string'}];

Summary = table('Size', [0 numel(varNames)], ...
    'VariableTypes', varTypes, 'VariableNames', varNames);

%% ===== Process each file =====
tStart = tic;
for k = 1:nToDo
    fpath = fullfile(files(k).folder, files(k).name);
    fname = files(k).name;

    % --- Parse Patient/Session from path ---
    tok = regexp(fpath, 'sub-Penn(\d+).+?ses-(\d+)', 'tokens', 'once');
    if isempty(tok)
        tok = regexp(fname, 'sub-Penn(\d+)_ses-(\d+)', 'tokens', 'once');
    end
    if isempty(tok)
        warning('Skipping %s (could not parse Patient/Session).', fname);
        continue
    end
    patientNum = str2double(tok{1});
    sessionNum = str2double(tok{2});

    pairKey = char(string(patientNum) + "_" + string(sessionNum));

    % --- Read probability trace ---
    try
        T = readtable(fpath, 'TextType', 'string');
    catch ME
        warning('Failed to read %s. Skipping. (%s)', fname, ME.message);
        continue
    end
    if ~ismember('SN2', T.Properties.VariableNames)
        warning('Skipping %s: column SN2 not found.', fname);
        continue
    end

    SN2 = double(T.SN2(:));
    SN2(~isfinite(SN2)) = -inf;      % NaN treated as below threshold
    N = numel(SN2);
    edf_dur = N / FS;

    % --- Detections (indices, not just count) ---
    detIdx = detect_spikes_wholefile(SN2, THRESH, WIN_SAMP);
    nTotal = numel(detIdx);

    % --- Look up recorded duration ---
    flag = "";
    if isKey(durMap, pairKey)
        rec_dur = durMap(pairKey);
    else
        rec_dur = NaN;
        flag = "no_clinical_match";
    end

    dead_sec = edf_dur - rec_dur;

    if isnan(dead_sec)
        nDeadSamp = NaN;
        nDead = NaN; nRec = NaN;
    else
        if dead_sec < -NEG_DEAD_TOL
            flag = strtrim(flag + " negative_deadtime");
        end
        dead_sec  = max(dead_sec, 0);
        dead_sec  = min(dead_sec, edf_dur);
        nDeadSamp = round(dead_sec * FS);

        nDead = sum(detIdx <= nDeadSamp);
        nRec  = nTotal - nDead;
        assert(nDead + nRec == nTotal, ...
            'Dead/recording split does not sum to total for %s', fname);
    end

    rec_sec_eff = max(edf_dur - dead_sec, 0);
    rate_dead = 60 * nDead / max(dead_sec, eps);
    rate_rec  = 60 * nRec  / max(rec_sec_eff, eps);

    Summary = [Summary; {string(fname), patientNum, sessionNum, ...
        edf_dur, rec_dur, dead_sec, nDead, nRec, nTotal, ...
        rate_dead, rate_rec, flag}]; %#ok<AGROW>

    if mod(k, PRINT_EVERY) == 0 || k == nToDo
        elapsed = toc(tStart);
        perFile = elapsed / k;
        remain  = perFile * (nToDo - k);
        fprintf('  [%d/%d] %.0f%% | elapsed %s | ~%s remaining | %.2f s/file\n', ...
            k, nToDo, 100*k/nToDo, fmt_dur(elapsed), fmt_dur(remain), perFile);
    end
end

%% ===== Sanity checks =====
assert(height(Summary) > 0, 'No files were successfully processed.');
matched = ~ismissing(Summary.Recorded_Duration_sec) & isfinite(Summary.Recorded_Duration_sec);
assert(all(Summary.n_spikes_deadtime(matched) + Summary.n_spikes_recording(matched) ...
    == Summary.n_spikes_total(matched)), 'Spike counts do not sum to total.');
fprintf('Matched clinical durations for %d / %d files.\n', sum(matched), height(Summary));
fprintf('Median dead time: %.1f min (IQR %.1f-%.1f)\n', ...
    median(Summary.Deadtime_sec(matched))/60, ...
    prctile(Summary.Deadtime_sec(matched), 25)/60, ...
    prctile(Summary.Deadtime_sec(matched), 75)/60);
fprintf('Median spike rate: dead %.2f/min vs recording %.2f/min\n', ...
    median(Summary.rate_deadtime_per_min(matched)), ...
    median(Summary.rate_recording_per_min(matched)));

%% ===== Save =====
Summary = sortrows(Summary, {'Patient', 'Session', 'EEG_Name'});

% Don't let a test run clobber the full-run output
if isfinite(MAX_FILES)
    [p, n, e] = fileparts(outCsv);
    outCsv = fullfile(p, sprintf('%s_first%d%s', n, MAX_FILES, e));
end
writetable(Summary, outCsv);
fprintf('Saved to: %s  (rows=%d)\n', outCsv, height(Summary));

%% ===== Local functions =====
function detIdx = detect_spikes_wholefile(SN2, THRESHOLD, WIN_SAMP)
% Returns the sample index of each clustered detection, using the identical
% 1-second skip rule as the original counting routine. numel(detIdx) equals
% the count the original code would have produced.

    Nloc   = numel(SN2);
    detIdx = zeros(Nloc, 1);
    nDet   = 0;
    i = 1;
    while i <= Nloc
        if SN2(i) > THRESHOLD
            nDet = nDet + 1;
            detIdx(nDet) = i;
            i = i + WIN_SAMP;      % skip ahead ~1 second
        else
            i = i + 1;
        end
    end
    detIdx = detIdx(1:nDet);
end

function s = fmt_dur(secs)
% Compact h:mm:ss / m:ss formatting for progress messages.
    if ~isfinite(secs)
        s = '--';
        return
    end
    h = floor(secs/3600);
    m = floor(mod(secs, 3600)/60);
    sec = floor(mod(secs, 60));
    if h > 0
        s = sprintf('%dh%02dm', h, m);
    elseif m > 0
        s = sprintf('%dm%02ds', m, sec);
    else
        s = sprintf('%ds', sec);
    end
end

function secs = hms_to_seconds(x)
% Parse " 0:36:51" style strings (or existing duration/numeric) into seconds.

    if isduration(x)
        secs = seconds(x);
        return
    end
    if isnumeric(x)
        secs = double(x);
        return
    end

    s = strtrim(string(x));
    secs = nan(numel(s), 1);
    for j = 1:numel(s)
        if strlength(s(j)) == 0 || ismissing(s(j))
            continue
        end
        parts = str2double(split(s(j), ":"));
        if any(isnan(parts))
            continue
        end
        switch numel(parts)
            case 3
                secs(j) = parts(1)*3600 + parts(2)*60 + parts(3);
            case 2
                secs(j) = parts(1)*60 + parts(2);
            case 1
                secs(j) = parts(1);
        end
    end
end
