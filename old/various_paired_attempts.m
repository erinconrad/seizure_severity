
function pairs = run_paired_spikerate_by_hasSz_NEW(SpikeSummaryTable, ReportTable, PatientTypingAll, ...
    badTypes, NESD_LABEL, EPS_PER_MIN)
% run_paired_spikerate_by_hasSz
%   Uses ALREADY-LOADED and ALREADY-OUTPATIENT-FILTERED:
%       - SpikeSummaryTable (with SpikeRate_perHour / SpikeRate_perMin)
%       - ReportTable       (with patient_id, session_number, visit_dates_deid, visit_hasSz, etc.)
%       - PatientTypingAll  (from build_patient_metrics_from_report)
%
%     For each patient:
%       1) Look across ALL visits and keep the patient only if they have at
%          least one visit with HasSz==0 and at least one visit with HasSz==1.
%       2) For that patient, define two EEG sets:
%            - HasSz==1 set: all EEGs within [ -preGapDays, +postGapDays ]
%              of ANY visit with HasSz==1 (union over visits).
%            - HasSz==0 set: same but for HasSz==0 visits.
%       3) Compute per-patient mean spike rates (spikes/hour) over these
%          two EEG sets → paired values for that patient.
%
%   Stats: Wilcoxon sign-rank on RAW spikes/hour.
%   Plot:  log10(spikes/min) with jittered points + paired lines.
%


    % ---------------- CONFIG LOCAL TO THIS ANALYSIS ----------------
    pairsOut    = '../output/eeg_visit_pairs_pairedmeans.csv';
    saveAudit   = true;
    preGapDays  = 365;   % EEG may be up to this many days BEFORE the visit
    postGapDays = 30;    % EEG may be up to this many days AFTER  the visit

    % EPS in spikes/hour and spikes/min
    EPS_PLOT_perMin = EPS_PER_MIN;          % already in per-minute units from main script
    EPS_PLOT_perHr  = EPS_PLOT_perMin * 60; % equivalent in spikes/hour

    % ---------------- COPY TABLES LOCALLY ----------------
    S = SpikeSummaryTable;
    R = ReportTable;

    % ---------------- PATIENT-LEVEL INCLUSION: "ALL EPILEPSY" ----------------
    % Use PatientTypingAll (already built in main script)
    et = string(PatientTypingAll.EpilepsyType);
    et_norm = lower(strtrim(et));
    isEmpty = ismissing(et) | strlength(strtrim(et))==0;
    isNESD  = et_norm == lower(strtrim(NESD_LABEL));
    isBad   = isNESD | ismember(et_norm, badTypes) | isEmpty;

    validP = double(PatientTypingAll.Patient(~isBad));

    S = S(ismember(S.Patient,    validP), :);
    R = R(ismember(R.patient_id, validP), :);

    % ---------------- EEG DATE FROM start_time_deid ----------------
    rawStart = string(R.start_time_deid);
    okStart  = ~ismissing(rawStart) & strlength(strtrim(rawStart))>0;
    R.EEG_Date = NaT(height(R),1);
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart)), 'InputFormat','yyyy-MM-d HH:mm:ss');

    % ---------------- PARSE visit_dates_deid + visit_hasSz ----------------
    R.VisitDates = cell(height(R),1);
    R.HasSzVec   = cell(height(R),1);
    for i = 1:height(R)
        % visit_dates_deid
        vstr = strtrim(string(R.visit_dates_deid(i)));
        if strlength(vstr)>0 && vstr~="[]"
            vcell = jsondecode(vstr);
            vdt   = datetime(vcell, 'InputFormat','yyyy-MM-dd');
            R.VisitDates{i} = vdt(:);
        else
            R.VisitDates{i} = NaT(0,1);
        end

        % visit_hasSz
        hstr = strtrim(string(R.visit_hasSz(i)));
        if strlength(hstr)>0 && hstr~="[]"
            hv = jsondecode(char(hstr));
            has_sz = double(hv(:));
            has_sz(has_sz==2) = nan;
            R.HasSzVec{i} = has_sz;
        else
            R.HasSzVec{i} = nan(0,1);
        end
    end

    % ---------------- JOIN S + R ON (Patient, Session) ----------------
    rightVars = {'EEG_Date','VisitDates','HasSzVec','epilepsy_type'};
    if ~ismember('acquired_on',          R.Properties.VariableNames), R.acquired_on = strings(height(R),1); end
    if ~ismember('report_PATIENT_CLASS', R.Properties.VariableNames), R.report_PATIENT_CLASS = strings(height(R),1); end
    rightVars{end+1} = 'acquired_on';
    rightVars{end+1} = 'report_PATIENT_CLASS';

    JR = innerjoin(S, R, ...
        'LeftKeys',{'Patient','Session'}, ...
        'RightKeys',{'patient_id','session_number'}, ...
        'RightVariables', rightVars);

    % keep rows with valid EEG date and at least one visit date
    JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) & any(~isnat(v)), JR.VisitDates), :);

    if isempty(JR)
        warning('run_paired_spikerate_by_hasSz: No EEG rows with valid visit arrays; skipping.');
        return;
    end

    %   * Aggregate all visits per patient
    %   * Require at least one HasSz==0 and one HasSz==1 visit
    %   * For each patient, define EEG sets that are within the gap of ANY
    %     visit of each type (union of EEGs per condition)
    %   * Then compute per-patient mean spike rates for HasSz==0 and HasSz==1

    % Preallocate containers
    allPatients = unique(JR.Patient);
    nPatients   = numel(allPatients);

    mu0_raw = nan(nPatients,1);  % mean spikes/hour for HasSz==0 set
    mu1_raw = nan(nPatients,1);  % mean spikes/hour for HasSz==1 set
    n0      = zeros(nPatients,1); % # EEGs in HasSz==0 set
    n1      = zeros(nPatients,1); % # EEGs in HasSz==1 set
    pid_vec = allPatients;        % keep patient IDs aligned

    % Audit table: one row per (patient, EEG, condition)
    pairs = table('Size',[0 9], ...
        'VariableTypes', {'double','double','string','datetime','double','datetime','double','double','string'}, ...
        'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour', ...
                          'Visit_Date','HasSz','GapDays_abs','EpilepsyType'});

    for ip = 1:nPatients
        pid = allPatients(ip);

        rowsP = JR(JR.Patient == pid, :);
        if isempty(rowsP), continue; end

        % ----- Collect ALL visits for this patient across rows -----
        vdates_all = NaT(0,1);
        has_all    = [];
        for k = 1:height(rowsP)
            vdates_all = [vdates_all; rowsP.VisitDates{k}(:)]; %#ok<AGROW>
            has_all    = [has_all;    rowsP.HasSzVec{k}(:)];   %#ok<AGROW>
        end

        % Keep only visits with HasSz in {0,1} and valid dates
        validVisitMask = isfinite(has_all) & (has_all==0 | has_all==1) & ~isnat(vdates_all);
        vdates_all = vdates_all(validVisitMask);
        has_all    = has_all(validVisitMask);

        if isempty(vdates_all)
            continue;
        end

        % Require at least one HasSz==0 visit and one HasSz==1 visit
        if ~any(has_all==0) || ~any(has_all==1)
            % Skip this patient entirely
            continue;
        end

        % Split visits by HasSz type
        v0 = vdates_all(has_all==0);
        v1 = vdates_all(has_all==1);

        % EEGs for this patient
        eDates   = rowsP.EEG_Date;
        eRates   = rowsP.SpikeRate_perHour;
        eSess    = rowsP.Session;
        eNames   = rowsP.EEG_Name;
        eEpiType = string(rowsP.epilepsy_type);

        nEEG = height(rowsP);

        isNear0   = false(nEEG,1);
        bestGap0  = nan(nEEG,1);
        bestV0    = NaT(nEEG,1);

        isNear1   = false(nEEG,1);
        bestGap1  = nan(nEEG,1);
        bestV1    = NaT(nEEG,1);

        % For each EEG, see if it is within the allowed gap of ANY visit of
        % each type (union across visits)
        for e = 1:nEEG
            % --- HasSz == 0 visits ---
            if ~isempty(v0)
                d0 = days(eDates(e) - v0);
                elig0 = (d0 >= -preGapDays) & (d0 <= postGapDays);
                if any(elig0)
                    isNear0(e) = true;
                    % representative visit date: nearest HasSz==0 visit
                    d0_elig = d0(elig0);
                    idx0_rel = find(elig0);
                    [~, ix0] = min(abs(d0_elig));
                    idx0 = idx0_rel(ix0);
                    bestGap0(e) = abs(d0(idx0));
                    bestV0(e)   = v0(idx0);
                end
            end

            % --- HasSz == 1 visits ---
            if ~isempty(v1)
                d1 = days(eDates(e) - v1);
                elig1 = (d1 >= -preGapDays) & (d1 <= postGapDays);
                if any(elig1)
                    isNear1(e) = true;
                    % representative visit date: nearest HasSz==1 visit
                    d1_elig = d1(elig1);
                    idx1_rel = find(elig1);
                    [~, ix1] = min(abs(d1_elig));
                    idx1 = idx1_rel(ix1);
                    bestGap1(e) = abs(d1(idx1));
                    bestV1(e)   = v1(idx1);
                end
            end
        end


        % EEG sets (UNION over visits) for each condition
        eegIdx0 = find(isNear0);
        eegIdx1 = find(isNear1);


        % Require at least one EEG in each set; otherwise skip the patient
        if isempty(eegIdx0) || isempty(eegIdx1)
            continue;
        end

        % Per-patient mean spike rates (spikes/hour) for each condition
        mu0_raw(ip) = mean(eRates(eegIdx0), 'omitnan');
        mu1_raw(ip) = mean(eRates(eegIdx1), 'omitnan');
        n0(ip)      = numel(eegIdx0);
        n1(ip)      = numel(eegIdx1);

        % ---- Audit rows: one row per (patient, EEG, condition) ----
        % HasSz == 0 EEGs
        for ii = eegIdx0'
            pairs = [pairs; {pid, eSess(ii), eNames(ii), eDates(ii), eRates(ii), ...
                             bestV0(ii), 0, bestGap0(ii), eEpiType(ii)}]; %#ok<AGROW>
        end

        % HasSz == 1 EEGs
        for ii = eegIdx1'
            pairs = [pairs; {pid, eSess(ii), eNames(ii), eDates(ii), eRates(ii), ...
                             bestV1(ii), 1, bestGap1(ii), eEpiType(ii)}]; %#ok<AGROW>
        end
    end

    % ---------------- PER-PATIENT FILTERING (BOTH STATES PRESENT) ----------------
    % At this point:
    %   mu0_raw, mu1_raw, n0, n1 are NaN/0 for patients that were skipped.

    isMu0NaN = isnan(mu0_raw);
    isN0zero = (n0 == 0);
    % Consistency check: should match except for patients we never visited
    assert(all(isMu0NaN == isN0zero), ...
        'Mismatch: Some patients have mu0_raw = NaN but n0 > 0, or vice versa.');

    isMu1NaN = isnan(mu1_raw);
    isN1zero = (n1 == 0);
    assert(all(isMu1NaN == isN1zero), ...
        'Mismatch: Some patients have mu1_raw = NaN but n1 > 0, or vice versa.');

    % Only keep patients with BOTH sets non-empty and finite means
    keep = isfinite(mu0_raw) & isfinite(mu1_raw) & (n0>0) & (n1>0);

    yNoSz_raw     = mu0_raw(keep);          % HasSz==0 mean spike rate (spikes/hour)
    ySz_raw       = mu1_raw(keep);          % HasSz==1 mean spike rate (spikes/hour)
    patients_kept = pid_vec(keep); %#ok<NASGU>

    Npairs = numel(ySz_raw);

    if Npairs == 0
        warning('run_paired_spikerate_by_hasSz: No patients with EEG coverage near both HasSz==0 and HasSz==1 visits.');
        return;
    end

    % Wilcoxon sign-rank on RAW spikes/hour
    diffs_raw   = ySz_raw - yNoSz_raw;
    p_signrank  = signrank(ySz_raw, yNoSz_raw, 'method','approx');
    medDeltaRaw = median(diffs_raw, 'omitnan');
    med0        = median(yNoSz_raw, 'omitnan');
    med1        = median(ySz_raw,   'omitnan');

    % ---------------- PLOT: log10(spikes/min) WITH JITTER ----------------
    % Convert raw means from spikes/hour → spikes/min
    yNoSz_perMin = yNoSz_raw / 60;
    ySz_perMin   = ySz_raw   / 60;

    % log10(spikes/min) with EPS guard (per-minute)
    logNoSz = log10(max(yNoSz_perMin, EPS_PLOT_perMin));
    logSz   = log10(max(ySz_perMin,   EPS_PLOT_perMin));

    jitterWidth = 0.08;   % 0.05–0.12 looks nice

    % Jittered x-positions for points (but NOT for lines)
    x0_jit = -jitterWidth/2 + jitterWidth*rand(Npairs,1);      % around 0
    x1_jit =  1 - jitterWidth/2 + jitterWidth*rand(Npairs,1);  % around 1

    % "zero spikes" line in log10(spikes/min)
    Y_ZERO = log10(EPS_PLOT_perMin);

    % Choose y-limits: a bit below Y_ZERO, a bit above max data
    yDataMax = max([logNoSz; logSz]);
    Y_LIMS   = [Y_ZERO - 0.4, yDataMax + 0.4];

    figure('Color','w','Position',[200 200 750 580]);
    hold on; grid on; box off; set(gca,'FontSize',16);

    % Paired lines (no jitter)
    for i = 1:Npairs
        plot([0 1], [logNoSz(i) logSz(i)], '-', 'Color',[0.7 0.7 0.7]);
    end

    % Jittered scatter points
    scatter(x0_jit, logNoSz, 50, 'filled', 'MarkerFaceAlpha',0.8);
    scatter(x1_jit, logSz,   50, 'filled', 'MarkerFaceAlpha',0.8);

    % Horizontal "zero spike rate" line
    yline(Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    % Axes / labels
    xlim([-0.25 1.25]);
    ylim(Y_LIMS);
    set(gca,'XTick',[0 1], ...
            'XTickLabel',{'Visits without seizures','Visits with seizures'});
    ylabel('log_{10} spike rate (spikes/min)');
    title(sprintf('All epilepsy patients'), ...
        'FontSize',18,'FontWeight','bold');

    subtitle(sprintf(['n=%d  median HasSz=0=%.3f/hr, HasSz=1=%.3f/hr, ' ...
                      'median Δ=%.3f/hr,  p=%.3g'], ...
        Npairs, med0, med1, medDeltaRaw, p_signrank));

    % --------- Significance bar with p-value above the groups ----------
    % Put it a bit below the top of the y-axis
    yBar = Y_LIMS(2) - 0.07 * range(Y_LIMS);
    tick = 0.03 * range(Y_LIMS);

    % Horizontal bar from x=0 to x=1
    plot([0 0 1 1], [yBar - tick, yBar, yBar, yBar - tick], ...
        'k-', 'LineWidth', 1.3);

    % p-value text centered between the two groups
    text(0.5, yBar + 0.02*range(Y_LIMS), ...
         sprintf('p = %.3g', p_signrank), ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', ...
         'FontSize',16, 'FontWeight','bold');

    % ---------------- CONSOLE SUMMARY + AUDIT ----------------
    fprintf('\nPaired sign-rank on per-patient mean spike rates (HasSz=1 minus HasSz=0):\n');
    fprintf('  n patients with both states (and EEGs near both) = %d\n', Npairs);
    fprintf('  Median HasSz=0 mean spike rate      = %.3f spikes/hour\n', med0);
    fprintf('  Median HasSz=1 mean spike rate      = %.3f spikes/hour\n', med1);
    fprintf('  Median difference (HasSz=1 - HasSz=0) = %.3f spikes/hour\n', medDeltaRaw);
    fprintf('  Wilcoxon sign-rank p = %.3g\n', p_signrank);

    if saveAudit
        if ~exist(fileparts(pairsOut),'dir'), mkdir(fileparts(pairsOut)); end
        writetable(pairs, pairsOut);
        fprintf('Saved EEG–visit unions audit to: %s\n', pairsOut);
    end
end




function pairs = run_paired_spikerate_by_hasSz_OLD(SpikeSummaryTable, ReportTable, PatientTypingAll, ...
    badTypes, NESD_LABEL, EPS_PER_MIN)
% run_paired_spikerate_by_hasSz
%   Uses ALREADY-LOADED and ALREADY-OUTPATIENT-FILTERED:
%       - SpikeSummaryTable (with SpikeRate_perHour / SpikeRate_perMin)
%       - ReportTable       (with patient_id, session_number, visit_dates_deid, visit_hasSz, etc.)
%       - PatientTypingAll  (from build_patient_metrics_from_report)
%
%   For each patient:
%       compares mean spike rate across visits with HasSz==0 vs HasSz==1
%       (per-patient paired analysis), for ALL valid epilepsy patients
%       (NESD + badTypes excluded).
%
%   Stats: Wilcoxon sign-rank on RAW spikes/hour.
%   Plot:  log10(spikes/min) with jittered points + paired lines.
%
%   Cohort filters respected:
%       - only_amb (0/all, 1/amb only, 2/routine only)
%       - outpatient filter has already been applied upstream.

    % ---------------- CONFIG LOCAL TO THIS ANALYSIS ----------------
    pairsOut    = '../output/eeg_visit_pairs_pairedmeans.csv';
    saveAudit   = true;
    preGapDays  = 365;   % EEG may be up to this many days BEFORE the visit
    postGapDays = 30;    % EEG may be up to this many days AFTER  the visit

    % EPS in spikes/hour and spikes/min
    EPS_PLOT_perMin = EPS_PER_MIN;          % already in per-minute units from main script
    EPS_PLOT_perHr  = EPS_PLOT_perMin * 60; % equivalent in spikes/hour


    % ---------------- COPY TABLES LOCALLY ----------------
    S = SpikeSummaryTable;
    R = ReportTable;

    % ---------------- PATIENT-LEVEL INCLUSION: "ALL EPILEPSY" ----------------
    % Use PatientTypingAll (already built in main script)
    et = string(PatientTypingAll.EpilepsyType);
    et_norm = lower(strtrim(et));
    isEmpty = ismissing(et) | strlength(strtrim(et))==0;
    isNESD  = et_norm == lower(strtrim(NESD_LABEL));
    isBad   = isNESD | ismember(et_norm, badTypes) | isEmpty;

    validP = double(PatientTypingAll.Patient(~isBad));

    S = S(ismember(S.Patient,     validP), :);
    R = R(ismember(R.patient_id,  validP), :);


    % ---------------- EEG DATE FROM start_time_deid ----------------
    rawStart = string(R.start_time_deid);
    okStart  = ~ismissing(rawStart) & strlength(strtrim(rawStart))>0;
    R.EEG_Date = NaT(height(R),1);
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart)), 'InputFormat','yyyy-MM-d HH:mm:ss');
    

    % ---------------- PARSE visit_dates_deid + visit_hasSz ----------------
    R.VisitDates = cell(height(R),1);
    R.HasSzVec   = cell(height(R),1);
    for i = 1:height(R)
        % visit_dates_deid
        vstr = strtrim(string(R.visit_dates_deid(i)));
        if strlength(vstr)>0 && vstr~="[]"
            vcell = jsondecode(vstr);
            vdt   = datetime(vcell, 'InputFormat','yyyy-MM-dd');
            R.VisitDates{i} = vdt(:);
        else
            R.VisitDates{i} = NaT(0,1);
        end

        % visit_hasSz
        hstr = strtrim(string(R.visit_hasSz(i)));
        if strlength(hstr)>0 && hstr~="[]"
            hv = jsondecode(char(hstr));
            has_sz = double(hv(:));
            has_sz(has_sz==2) = nan;
            R.HasSzVec{i} = has_sz;
        else
            R.HasSzVec{i} = nan(0,1);
        end
    end

    % ---------------- JOIN S + R ON (Patient, Session) ----------------
    rightVars = {'EEG_Date','VisitDates','HasSzVec','epilepsy_type'};
    if ~ismember('acquired_on',          R.Properties.VariableNames), R.acquired_on = strings(height(R),1); end
    if ~ismember('report_PATIENT_CLASS', R.Properties.VariableNames), R.report_PATIENT_CLASS = strings(height(R),1); end
    rightVars{end+1} = 'acquired_on';
    rightVars{end+1} = 'report_PATIENT_CLASS';

    JR = innerjoin(S, R, ...
        'LeftKeys',{'Patient','Session'}, ...
        'RightKeys',{'patient_id','session_number'}, ...
        'RightVariables', rightVars);

    % keep rows with valid EEG date and at least one visit date
    JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) & any(~isnat(v)), JR.VisitDates), :);

    if isempty(JR)
        warning('run_paired_spikerate_by_hasSz: No EEG rows with valid visit arrays; skipping.');
        return;
    end

    % ---------------- CONSTRUCT EEG–VISIT PAIRS (NEAREST VISIT IN ASYMMETRIC WINDOW) ----------------
    pairs = table('Size',[0 9], ...
        'VariableTypes', {'double','double','string','datetime','double','datetime','double','double','string'}, ...
        'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour', ...
                          'Visit_Date','HasSz','GapDays_abs','EpilepsyType'});

    % Loop over EEGs
    for i = 1:height(JR)
        eegDate = JR.EEG_Date(i); % get eeg date
        vdates  = JR.VisitDates{i}; % get all visit dates
        hasV    = JR.HasSzVec{i};
        pid     = JR.Patient(i);

        if isempty(vdates) || isempty(hasV), continue; end
        assert(numel(hasV)==numel(vdates))
        nAlign = min(numel(vdates), numel(hasV));
        if nAlign==0, continue; end
        vdates = vdates(1:nAlign);
        hasV   = hasV(1:nAlign);

        % compute the difference in days between EEG and visit date
        dSigned = days(eegDate - vdates);
    
        % base eligibility: within asymmetric time window
        inWindow = (dSigned >= -preGapDays) & (dSigned <= postGapDays);
    
        % require HasSz == 0 or 1 for eligibility
        hasValidHasSz = isfinite(hasV) & (hasV == 0 | hasV == 1);
    
        elig = inWindow & hasValidHasSz;
    
        % if no visit in the window with HasSz==0/1, skip this EEG
        if ~any(elig)
            continue
        end
        
    
        % among eligible visits, choose the nearest
        dSigned_elig = dSigned(elig);
        [~, idxRel]  = min(abs(dSigned_elig));
        idxVec       = find(elig);
        j            = idxVec(idxRel);
    
        thisHas = hasV(j);  % guaranteed 0 or 1 here
        

        pairs = [pairs; {pid, JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perHour(i), ...
                         vdates(j), double(thisHas), abs(dSigned(j)), string(JR.epilepsy_type(i))}]; %#ok<AGROW>
    end

    fprintf('Matched %d EEG–visit pairs under [%d days pre, %d days post].\n', ...
        height(pairs), preGapDays, postGapDays);

    orig_pairs = pairs;

    % Require ≥2 pairs per patient (at least two EEG–visit pairs)
    gc = groupcounts(pairs,'Patient');
    keepP2 = gc.Patient(gc.GroupCount>=2);
    pairs  = pairs(ismember(pairs.Patient, keepP2), :);


    % ---------------- PER-PATIENT MEANS: HasSz==0 vs HasSz==1 ----------------
    % Group by patient
    [G, pid] = findgroups(double(pairs.Patient));

    % patient level means for visits without seizures
    mu0_raw = splitapply(@(y,hs) mean(y(hs==0), 'omitnan'), pairs.SpikeRate_perHour, pairs.HasSz, G); % average spike rate for each patient for visits without szs
    mu1_raw = splitapply(@(y,hs) mean(y(hs==1), 'omitnan'), pairs.SpikeRate_perHour, pairs.HasSz, G);

    % number of visits of each type (has sz==1 and has sz == 0)
    n0      = splitapply(@(hs) sum(hs==0),       pairs.HasSz, G); 
    n1      = splitapply(@(hs) sum(hs==1),       pairs.HasSz, G);

    if 0
        table(mu0_raw,mu1_raw,n0,n1)
    end
    % Logical masks
    isMu0NaN = isnan(mu0_raw);
    isN0zero = (n0 == 0);
    
    % Assert that NaNs in mu0_raw correspond exactly to n0 == 0
    assert(isequal(isMu0NaN, isN0zero), ...
        'Mismatch: Some patients have mu0_raw = NaN but n0 > 0, or vice versa.');

    % Logical masks
    isMu1NaN = isnan(mu1_raw);
    isN1zero = (n1 == 0);
    
    % Assert consistency
    assert(isequal(isMu1NaN, isN1zero), ...
        'Mismatch: Some patients have mu1_raw = NaN but n1 > 0, or vice versa.');


    % only keep patients if there are visits with seizures AND visits
    % without seizures
    keep = isfinite(mu0_raw) & isfinite(mu1_raw) & (n0>0) & (n1>0);
    yNoSz_raw     = mu0_raw(keep);       % HasSz==0 mean spike rate (spikes/hour)
    ySz_raw       = mu1_raw(keep);       % HasSz==1 mean spike rate (spikes/hour)
    patients_kept = pid(keep); %#ok<NASGU>

    Npairs = numel(ySz_raw);


    % Wilcoxon sign-rank on RAW spikes/hour
    diffs_raw   = ySz_raw - yNoSz_raw;
    p_signrank  = signrank(ySz_raw, yNoSz_raw, 'method','approx');
    medDeltaRaw = median(diffs_raw, 'omitnan');
    med0        = median(yNoSz_raw, 'omitnan');
    med1        = median(ySz_raw,   'omitnan');

      % ---------------- PLOT: log10(spikes/min) WITH JITTER ----------------
    % Convert raw means from spikes/hour → spikes/min
    yNoSz_perMin = yNoSz_raw / 60;
    ySz_perMin   = ySz_raw   / 60;

    % log10(spikes/min) with EPS guard (per-minute)
    logNoSz = log10(max(yNoSz_perMin, EPS_PLOT_perMin));
    logSz   = log10(max(ySz_perMin,   EPS_PLOT_perMin));

    jitterWidth = 0.08;   % 0.05–0.12 looks nice

    % Jittered x-positions for points (but NOT for lines)
    x0_jit = -jitterWidth/2 + jitterWidth*rand(Npairs,1);      % around 0
    x1_jit =  1 - jitterWidth/2 + jitterWidth*rand(Npairs,1);  % around 1

    % "zero spikes" line in log10(spikes/min)
    Y_ZERO = log10(EPS_PLOT_perMin);

    % Choose y-limits: a bit below Y_ZERO, a bit above max data
    yDataMax = max([logNoSz; logSz]);
    Y_LIMS   = [Y_ZERO - 0.4, yDataMax + 0.4];

    figure('Color','w','Position',[200 200 750 580]);
    hold on; grid on; box off; set(gca,'FontSize',16);

    % Paired lines (no jitter)
    for i = 1:Npairs
        plot([0 1], [logNoSz(i) logSz(i)], '-', 'Color',[0.7 0.7 0.7]);
    end

    % Jittered scatter points
    scatter(x0_jit, logNoSz, 50, 'filled', 'MarkerFaceAlpha',0.8);
    scatter(x1_jit, logSz,   50, 'filled', 'MarkerFaceAlpha',0.8);

    % Horizontal "zero spike rate" line
    yline(Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    % Axes / labels
    xlim([-0.25 1.25]);
    ylim(Y_LIMS);
    set(gca,'XTick',[0 1], ...
            'XTickLabel',{'Visits without seizures','Visits with seizures'});
    ylabel('log_{10} spike rate (spikes/min)');
    title(sprintf('All epilepsy patients'), ...
        'FontSize',18,'FontWeight','bold');

    subtitle(sprintf(['n=%d  median HasSz=0=%.3f/hr, HasSz=1=%.3f/hr, ' ...
                      'median Δ=%.3f/hr,  p=%.3g'], ...
        Npairs, med0, med1, medDeltaRaw, p_signrank));

    % --------- Significance bar with p-value above the groups ----------
    % Put it a bit below the top of the y-axis
    yBar = Y_LIMS(2) - 0.07 * range(Y_LIMS);
    tick = 0.03 * range(Y_LIMS);

    % Horizontal bar from x=0 to x=1
    plot([0 0 1 1], [yBar - tick, yBar, yBar, yBar - tick], ...
        'k-', 'LineWidth', 1.3);

    % p-value text centered between the two groups
    text(0.5, yBar + 0.02*range(Y_LIMS), ...
         sprintf('p = %.3g', p_signrank), ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', ...
         'FontSize',16, 'FontWeight','bold');


    % ---------------- CONSOLE SUMMARY + AUDIT ----------------
    fprintf('\nPaired sign-rank on per-patient mean spike rates (HasSz=1 minus HasSz=0):\n');
    fprintf('  n patients with both states = %d\n', Npairs);
    fprintf('  Median HasSz=0 mean spike rate      = %.3f spikes/hour\n', med0);
    fprintf('  Median HasSz=1 mean spike rate      = %.3f spikes/hour\n', med1);
    fprintf('  Median difference (HasSz=1 - HasSz=0) = %.3f spikes/hour\n', medDeltaRaw);
    fprintf('  Wilcoxon sign-rank p = %.3g\n', p_signrank);

    if saveAudit
        if ~exist(fileparts(pairsOut),'dir'), mkdir(fileparts(pairsOut)); end
        writetable(pairs, pairsOut);
        fprintf('Saved EEG–visit pairs audit to: %s\n', pairsOut);
    end
end




function run_paired_spikerate_by_hasSz(SpikeSummaryTable, ReportTable, PatientTypingAll, ...
    badTypes, NESD_LABEL, EPS_PER_MIN)
% run_paired_spikerate_by_hasSz
%   Uses ALREADY-LOADED and ALREADY-OUTPATIENT-FILTERED:
%       - SpikeSummaryTable (with SpikeRate_perHour / SpikeRate_perMin)
%       - ReportTable       (with patient_id, session_number, visit_dates_deid, visit_hasSz, etc.)
%       - PatientTypingAll  (from build_patient_metrics_from_report)
%
%   For each patient:
%       compares mean spike rate across visits with HasSz==0 vs HasSz==1
%       (per-patient paired analysis), for ALL valid epilepsy patients
%       (NESD + badTypes excluded).
%
%   Stats:
%       - Wilcoxon sign-rank on RAW spikes/hour.
%       - Effect size = Hodges–Lehmann (HL) estimator of the paired
%         difference in spikes/hour (HasSz==1 minus HasSz==0).
%
%   Plot:
%       - Left panel: paired plot of log10(spikes/min) with jitter
%                    + p-value + HL Δ in title/subtitle.
%       - Right panel: sensitivity of HL vs. allowable EEG–visit gap
%                     (two lines, significance marked by * where p<0.05).
%
%   Cohort filters respected:
%       - only_amb (0/all, 1/amb only, 2/routine only)
%       - outpatient filter has already been applied upstream.

    % ---------------- CONFIG LOCAL TO THIS ANALYSIS ----------------
    pairsOut    = '../output/eeg_visit_pairs_pairedmeans.csv';
    saveAudit   = true;
    preGapDays  = 365;   % BASELINE: EEG may be up to this many days BEFORE the visit
    postGapDays = 30;    % BASELINE: EEG may be up to this many days AFTER  the visit

    % EPS in spikes/hour and spikes/min
    EPS_PLOT_perMin = EPS_PER_MIN;          % already in per-minute units from main script


    % ---------------- COPY TABLES LOCALLY ----------------
    S = SpikeSummaryTable;
    R = ReportTable;

    % Sanity: need these columns
    needS = {'Patient','Session','Duration_sec','SpikeRate_perHour'};
    if ~all(ismember(needS, S.Properties.VariableNames))
        error('SpikeSummaryTable is missing required columns: %s', ...
            strjoin(setdiff(needS, S.Properties.VariableNames), ', '));
    end

    needR = {'patient_id','session_number','start_time_deid','visit_dates_deid', ...
             'visit_hasSz','epilepsy_type','epilepsy_specific'};
    if ~all(ismember(needR, R.Properties.VariableNames))
        error('ReportTable is missing required columns: %s', ...
            strjoin(setdiff(needR, R.Properties.VariableNames), ', '));
    end

    % Ensure numeric IDs
    if ~isnumeric(S.Patient),        S.Patient        = double(str2double(string(S.Patient)));        end
    if ~isnumeric(S.Session),        S.Session        = double(str2double(string(S.Session)));        end
    if ~isnumeric(R.patient_id),     R.patient_id     = double(str2double(string(R.patient_id)));     end
    if ~isnumeric(R.session_number), R.session_number = double(str2double(string(R.session_number))); end

    % ---------------- PATIENT-LEVEL INCLUSION: "ALL EPILEPSY" ----------------
    % Use PatientTypingAll (already built in main script)
    et = string(PatientTypingAll.EpilepsyType);
    et_norm = lower(strtrim(et));
    isEmpty = ismissing(et) | strlength(strtrim(et))==0;
    isNESD  = et_norm == lower(strtrim(NESD_LABEL));
    isBad   = isNESD | ismember(et_norm, badTypes) | isEmpty;

    validP = double(PatientTypingAll.Patient(~isBad));

    S = S(ismember(S.Patient,     validP), :);
    R = R(ismember(R.patient_id,  validP), :);

    % ---------------- EEG DATE FROM start_time_deid ----------------
    rawStart = string(R.start_time_deid);
    okStart  = ~ismissing(rawStart) & strlength(strtrim(rawStart))>0;
    R.EEG_Date = NaT(height(R),1);
    R.EEG_Date(okStart) = datetime(strtrim(rawStart(okStart)), 'InputFormat','yyyy-MM-d HH:mm:ss');

    % ---------------- PARSE visit_dates_deid + visit_hasSz ----------------
    R.VisitDates = cell(height(R),1);
    R.HasSzVec   = cell(height(R),1);
    for i = 1:height(R)
        % visit_dates_deid
        vstr = strtrim(string(R.visit_dates_deid(i)));
        if strlength(vstr)>0 && vstr~="[]"
            vcell = jsondecode(vstr);
            vdt   = datetime(vcell, 'InputFormat','yyyy-MM-dd');
            R.VisitDates{i} = vdt(:);
        else
            R.VisitDates{i} = NaT(0,1);
        end

        % visit_hasSz
        hstr = strtrim(string(R.visit_hasSz(i)));
        if strlength(hstr)>0 && hstr~="[]"
            hv = jsondecode(char(hstr));
            has_sz = double(hv(:));
            has_sz(has_sz==2) = nan;
            R.HasSzVec{i} = has_sz;
        else
            R.HasSzVec{i} = nan(0,1);
        end
    end

    % ---------------- JOIN S + R ON (Patient, Session) FOR BASELINE ----------------
    rightVars = {'EEG_Date','VisitDates','HasSzVec','epilepsy_type','acquired_on','report_PATIENT_CLASS'};
    JR = innerjoin(S, R, ...
        'LeftKeys',{'Patient','Session'}, ...
        'RightKeys',{'patient_id','session_number'}, ...
        'RightVariables', rightVars);

    % keep rows with valid EEG date and at least one visit date
    JR = JR(~isnat(JR.EEG_Date) & cellfun(@(v)~isempty(v) & any(~isnat(v)), JR.VisitDates), :);

    if isempty(JR)
        warning('run_paired_spikerate_by_hasSz: No EEG rows with valid visit arrays; skipping.');
        return;
    end

    % ---------------- BASELINE EEG–VISIT PAIRS (NEAREST VISIT IN ASYMMETRIC WINDOW) ----------------
    pairs = table('Size',[0 9], ...
        'VariableTypes', {'double','double','string','datetime','double','datetime','double','double','string'}, ...
        'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour', ...
                          'Visit_Date','HasSz','GapDays_abs','EpilepsyType'});

    for i = 1:height(JR)
        eegDate = JR.EEG_Date(i);
        vdates  = JR.VisitDates{i};
        hasV    = JR.HasSzVec{i};

        if isempty(vdates) || isempty(hasV), continue; end
        assert(numel(hasV)==numel(vdates))
        nAlign = min(numel(vdates), numel(hasV));
        if nAlign==0, continue; end
        vdates = vdates(1:nAlign);
        hasV   = hasV(1:nAlign);

        % compute the difference in days between eeg and visit date
        dSigned = days(eegDate - vdates);

        % determine which are within allowable gap
        elig = (dSigned >= -preGapDays) & (dSigned <= postGapDays);
        if ~any(elig), continue; end

        % Find the closest
        [~, idxRel] = min(abs(dSigned(elig)));
        idxVec = find(elig);
        j = idxVec(idxRel);

        thisHas = hasV(j);
        if ~(isfinite(thisHas) && (thisHas==0 || thisHas==1))
            continue
        end

        pairs = [pairs; {JR.Patient(i), JR.Session(i), JR.EEG_Name(i), eegDate, JR.SpikeRate_perHour(i), ...
                         vdates(j), double(thisHas), abs(dSigned(j)), string(JR.epilepsy_type(i))}]; %#ok<AGROW>
    end

    fprintf('Matched %d EEG–visit pairs under [%d days pre, %d days post].\n', ...
        height(pairs), preGapDays, postGapDays);

    % Require ≥2 pairs per patient (at least two EEG–visit pairs)
    gc = groupcounts(pairs,'Patient');
    keepP2 = gc.Patient(gc.GroupCount>=2);
    pairs  = pairs(ismember(pairs.Patient, keepP2), :);

    % ---------------- PER-PATIENT MEANS: HasSz==0 vs HasSz==1 ----------------
    % Group by patient
    [G, pid] = findgroups(double(pairs.Patient));

    % patient level means for visits without seizures
    mu0_raw = splitapply(@(y,hs) mean(y(hs==0), 'omitnan'), pairs.SpikeRate_perHour, pairs.HasSz, G); % HasSz==0
    mu1_raw = splitapply(@(y,hs) mean(y(hs==1), 'omitnan'), pairs.SpikeRate_perHour, pairs.HasSz, G); % HasSz==1

    % number of visits of each type (has sz==1 and has sz == 0)
    n0      = splitapply(@(hs) sum(hs==0), pairs.HasSz, G); 
    n1      = splitapply(@(hs) sum(hs==1), pairs.HasSz, G);

    % Consistency checks
    isMu0NaN = isnan(mu0_raw);
    isN0zero = (n0 == 0);
    assert(isequal(isMu0NaN, isN0zero), ...
        'Mismatch: Some patients have mu0_raw = NaN but n0 > 0, or vice versa.');

    isMu1NaN = isnan(mu1_raw);
    isN1zero = (n1 == 0);
    assert(isequal(isMu1NaN, isN1zero), ...
        'Mismatch: Some patients have mu1_raw = NaN but n1 > 0, or vice versa.');

    % only keep patients if there are visits with seizures AND visits without seizures
    keep = isfinite(mu0_raw) & isfinite(mu1_raw) & (n0>0) & (n1>0);
    yNoSz_raw     = mu0_raw(keep);       % HasSz==0 mean spike rate (spikes/hour)
    ySz_raw       = mu1_raw(keep);       % HasSz==1 mean spike rate (spikes/hour)
    patients_kept = pid(keep); %#ok<NASGU>

    Npairs = numel(ySz_raw);

    % Wilcoxon sign-rank on RAW spikes/hour
    diffs_raw   = ySz_raw - yNoSz_raw;
    if Npairs >= 1
        p_signrank  = signrank(ySz_raw, yNoSz_raw, 'method','approx');
    else
        p_signrank = NaN;
    end

    % Effect size = Hodges–Lehmann estimator (paired difference)
    HL_raw      = median(diffs_raw, 'omitnan');        % spikes/hour
    med0        = median(yNoSz_raw, 'omitnan');
    med1        = median(ySz_raw,   'omitnan');

    % ---------------- PLOT: log10(spikes/min) WITH JITTER + SENSITIVITY PANEL ----------------
    % Convert raw means from spikes/hour → spikes/min
    yNoSz_perMin = yNoSz_raw / 60;
    ySz_perMin   = ySz_raw   / 60;

    % log10(spikes/min) with EPS guard (per-minute)
    logNoSz = log10(max(yNoSz_perMin, EPS_PLOT_perMin));
    logSz   = log10(max(ySz_perMin,   EPS_PLOT_perMin));

    jitterWidth = 0.08;   % 0.05–0.12 looks nice

    % Jittered x-positions for points (but NOT for lines)
    x0_jit = -jitterWidth/2 + jitterWidth*rand(Npairs,1);      % around 0
    x1_jit =  1 - jitterWidth/2 + jitterWidth*rand(Npairs,1);  % around 1

    % "zero spikes" line in log10(spikes/min)
    Y_ZERO = log10(EPS_PLOT_perMin);

    % Choose y-limits: a bit below Y_ZERO, a bit above max data
    yDataMax = max([logNoSz; logSz]);
    Y_LIMS   = [Y_ZERO - 0.4, yDataMax + 0.4];

    % ---- 2-panel figure: paired plot + sensitivity ----
    fPair = figure('Color','w','Position',[200 200 1200 580]);
    tl = tiledlayout(fPair,1,2,'Padding','compact','TileSpacing','compact');

    % ===== LEFT PANEL: PAIRED PLOT =====
    ax1 = nexttile(tl,1);
    hold(ax1,'on'); grid(ax1,'on'); box(ax1,'off'); set(ax1,'FontSize',16);

    % Paired lines (no jitter)
    for i = 1:Npairs
        plot(ax1, [0 1], [logNoSz(i) logSz(i)], '-', 'Color',[0.7 0.7 0.7]);
    end

    % Jittered scatter points
    scatter(ax1, x0_jit, logNoSz, 50, 'filled', 'MarkerFaceAlpha',0.8);
    scatter(ax1, x1_jit, logSz,   50, 'filled', 'MarkerFaceAlpha',0.8);

    % Horizontal "zero spike rate" line
    yline(ax1, Y_ZERO, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    % Axes / labels
    xlim(ax1, [-0.25 1.25]);
    ylim(ax1, Y_LIMS);
    set(ax1,'XTick',[0 1], ...
            'XTickLabel',{'Visits without seizures','Visits with seizures'});
    ylabel(ax1, 'log_{10} spike rate (spikes/min)');
    title(ax1, sprintf('All epilepsy patients'), ...
        'FontSize',18,'FontWeight','bold');

    % Subtitle with HL (spikes/hour) and p
    subtitle(ax1, sprintf(['n=%d  median HasSz=0=%.3f/hr, HasSz=1=%.3f/hr, ' ...
                           'HL \\Delta=%.3f/hr,  p=%.3g'], ...
        Npairs, med0, med1, HL_raw, p_signrank));

    % Significance bar with p-value above the groups
    yBar = Y_LIMS(2) - 0.07 * range(Y_LIMS);
    tick = 0.03 * range(Y_LIMS);

    plot(ax1, [0 0 1 1], [yBar - tick, yBar, yBar, yBar - tick], ...
        'k-', 'LineWidth', 1.3);

    text(ax1, 0.5, yBar + 0.02*range(Y_LIMS), ...
         sprintf('p = %.3g', p_signrank), ...
         'HorizontalAlignment','center', ...
         'VerticalAlignment','bottom', ...
         'FontSize',16, 'FontWeight','bold');

    % ===== RIGHT PANEL: SENSITIVITY OF HL TO GAP WINDOWS =====
    ax2 = nexttile(tl,2);
    hold(ax2,'on'); grid(ax2,'on'); box(ax2,'off'); set(ax2,'FontSize',16);

    % 1) Fix preGap = 365, vary postGap
    pre_fixed   = 365;
    %post_values = [0 1 10 15 30 45 60 90 180 365];
    post_values = 0:10:600;

    HL_post   = nan(size(post_values));
    p_post    = nan(size(post_values));
    N_post    = nan(size(post_values));

    for k = 1:numel(post_values)
        [HL_post(k), p_post(k), N_post(k)] = gap_effect(pre_fixed, post_values(k));
    end

    % 2) Fix postGap = 30, vary preGap
    post_fixed = 30;
    %pre_values = [730 365 180 60 30 15];  % as requested
    pre_values = [600:-10:0];

    HL_pre   = nan(size(pre_values));
    p_pre    = nan(size(pre_values));
    N_pre    = nan(size(pre_values));

    for k = 1:numel(pre_values)
        [HL_pre(k), p_pre(k), N_pre(k)] = gap_effect(pre_values(k), post_fixed);
    end

    % Plot lines
    h1=plot(ax2, post_values, HL_post, '-o', 'LineWidth',2, 'MarkerSize',7, ...
        'DisplayName','vary post-gap (pre=365d)');
    h2=plot(ax2, pre_values,  HL_pre,  '-s', 'LineWidth',2, 'MarkerSize',7, ...
        'DisplayName','vary pre-gap (post=30d)');

    % Reference line at zero effect
    yline(ax2, 0, ':', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    % Axis labels
    xlabel(ax2, 'Gap (days)');
    ylabel(ax2, 'Hodges–Lehmann \Delta spike rate (spikes/hour)');
    title(ax2, 'Sensitivity of HL to EEG–visit gap','FontSize',18,'FontWeight','bold');

    % Nice limits
    allHL = [HL_post(:); HL_pre(:)];
    if all(isnan(allHL))
        yMin = -1; yMax = 1;
    else
        yMin = min(allHL(~isnan(allHL))) - 0.2*range(allHL(~isnan(allHL)));
        yMax = max(allHL(~isnan(allHL))) + 0.2*range(allHL(~isnan(allHL)));
        if yMin == yMax
            yMin = yMin - 0.5;
            yMax = yMax + 0.5;
        end
    end
    ylim(ax2, [yMin yMax]);

    % Asterisks over significant points (p<0.05)
    yr = diff(ylim(ax2));
        % Color of each line
    col_post = get(h1, 'Color');   % color of post-gap curve
    col_pre  = get(h2, 'Color');   % color of pre-gap curve
    
    % Asterisks for post-gap line
    for k = 1:numel(post_values)
        if ~isnan(p_post(k)) && p_post(k) < 0.05 && ~isnan(HL_post(k))
            text(ax2, post_values(k), HL_post(k) + 0.04*yr, '*', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',18, 'FontWeight','bold', ...
                'Color', col_post);
        end
    end
    
    % Asterisks for pre-gap line
    for k = 1:numel(pre_values)
        if ~isnan(p_pre(k)) && p_pre(k) < 0.05 && ~isnan(HL_pre(k))
            text(ax2, pre_values(k), HL_pre(k) + 0.04*yr, '*', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',18, 'FontWeight','bold', ...
                'Color', col_pre);
        end
    end


    legend(ax2,'Location','best');

    % ---------------- CONSOLE SUMMARY + AUDIT ----------------
    fprintf('\nPaired sign-rank on per-patient mean spike rates (HasSz=1 minus HasSz=0):\n');
    fprintf('  n patients with both states = %d\n', Npairs);
    fprintf('  Median HasSz=0 mean spike rate      = %.3f spikes/hour\n', med0);
    fprintf('  Median HasSz=1 mean spike rate      = %.3f spikes/hour\n', med1);
    fprintf('  Hodges–Lehmann (HL) difference      = %.3f spikes/hour\n', HL_raw);
    fprintf('  Wilcoxon sign-rank p                = %.3g\n', p_signrank);

    if saveAudit
        if ~exist(fileparts(pairsOut),'dir'), mkdir(fileparts(pairsOut)); end
        writetable(pairs, pairsOut);
        fprintf('Saved EEG–visit pairs audit to: %s\n', pairsOut);
    end

        % ======= NESTED HELPER: HL & p FOR ARBITRARY GAP SETTINGS =======
    function [HL_gap, p_gap, N_gap] = gap_effect(preD, postD)
        % Computes Hodges–Lehmann (HL) and Wilcoxon p-value for the same
        % S/R cohort, but with a different asymmetric gap window:
        %   -preD days <= (EEG_Date - VisitDate) <= postD
        %
        % Uses the already-built JR table (same S/R cohort) and repeats the
        % pairing + per-patient averaging logic.

        % Initialize output
        HL_gap = NaN;
        p_gap  = NaN;
        N_gap  = 0;

        if isempty(JR)
            return;
        end

        % ---- Build pairs for this gap window ----
        pairs_tmp = table('Size',[0 9], ...
            'VariableTypes', {'double','double','string','datetime','double','datetime','double','double','string'}, ...
            'VariableNames', {'Patient','Session','EEG_Name','EEG_Date','SpikeRate_perHour', ...
                              'Visit_Date','HasSz','GapDays_abs','EpilepsyType'});

        for ii = 1:height(JR)
            eegDate = JR.EEG_Date(ii);
            vdates  = JR.VisitDates{ii};
            hasV    = JR.HasSzVec{ii};

            if isempty(vdates) || isempty(hasV), continue; end
            assert(numel(hasV)==numel(vdates))
            nAlign = min(numel(vdates), numel(hasV));
            if nAlign==0, continue; end
            vdates = vdates(1:nAlign);
            hasV   = hasV(1:nAlign);

            dSigned = days(eegDate - vdates);
            elig    = (dSigned >= -preD) & (dSigned <= postD);
            if ~any(elig), continue; end

            [~, idxRel] = min(abs(dSigned(elig)));
            idxVec = find(elig);
            j      = idxVec(idxRel);

            thisHas = hasV(j);
            if ~(isfinite(thisHas) && (thisHas==0 || thisHas==1))
                continue;
            end

            pairs_tmp = [pairs_tmp; {JR.Patient(ii), JR.Session(ii), JR.EEG_Name(ii), eegDate, JR.SpikeRate_perHour(ii), ...
                                     vdates(j), double(thisHas), abs(dSigned(j)), string(JR.epilepsy_type(ii))}]; %#ok<AGROW>
        end

        if isempty(pairs_tmp)
            return;
        end

        % Require ≥2 pairs per patient
        gc_tmp   = groupcounts(pairs_tmp,'Patient');
        keepP2_t = gc_tmp.Patient(gc_tmp.GroupCount>=2);
        pairs_tmp = pairs_tmp(ismember(pairs_tmp.Patient, keepP2_t), :);

        if isempty(pairs_tmp)
            return;
        end

        % ---- Per-patient mean spike rates for HasSz==0 vs HasSz==1 ----
        [G_t, pid_t] = findgroups(double(pairs_tmp.Patient));

        mu0_t = splitapply(@(y,hs) mean(y(hs==0), 'omitnan'), ...
                           pairs_tmp.SpikeRate_perHour, pairs_tmp.HasSz, G_t);
        mu1_t = splitapply(@(y,hs) mean(y(hs==1), 'omitnan'), ...
                           pairs_tmp.SpikeRate_perHour, pairs_tmp.HasSz, G_t);

        n0_t  = splitapply(@(hs) sum(hs==0), pairs_tmp.HasSz, G_t);
        n1_t  = splitapply(@(hs) sum(hs==1), pairs_tmp.HasSz, G_t);

        % keep patients with at least one 0-visit and one 1-visit
        keep_t = isfinite(mu0_t) & isfinite(mu1_t) & (n0_t>0) & (n1_t>0);

        if ~any(keep_t)
            return;
        end

        y0_t = mu0_t(keep_t);
        y1_t = mu1_t(keep_t);

        diffs_t = y1_t - y0_t;
        N_gap   = numel(diffs_t);

        % HL = median of paired differences
        HL_gap  = median(diffs_t, 'omitnan');

        % Wilcoxon sign-rank p-value (if enough pairs)
        if N_gap >= 1
            p_gap = signrank(y1_t, y0_t, 'method','approx');
        else
            p_gap = NaN;
        end
    end % gap_effect

end % run_paired_spikerate_by_hasSz

