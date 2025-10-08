%% severity vs spikes — add avg_sz_freq to T without join

redcap_table = '../data/Routineeegpec-Deidreport_DATA_LABELS_2025-10-07_1611.csv';

% Load with exact headers and string parsing
T = readtable(redcap_table, 'TextType','string', 'VariableNamingRule','preserve');

% Ensure patient_id is numeric (as your head(T) shows). If it's not, convert:
% T.patient_id = double(T.patient_id);

% --- Helper: parse one sz_freqs cell into numeric vector (ignores nulls) ---
function v = parseSzFreqCell(s)
    v = [];
    if ~exist('s','var') || ismissing(s) || strlength(s)==0
        return
    end
    s = strtrim(s);
    try
        x = jsondecode(s);   % e.g., [0.7, null, 0.8, ...]
        if isnumeric(x)
            v = double(x(:));
        elseif iscell(x)
            % convert each element; allow numeric or string
            v = nan(numel(x),1);
            for ii = 1:numel(x)
                xi = x{ii};
                if isempty(xi)
                    v(ii) = NaN;
                elseif isnumeric(xi)
                    v(ii) = double(xi);
                else
                    % strings like "null" or "0.724..."
                    xi = string(xi);
                    if xi == "null" || xi == ""
                        v(ii) = NaN;
                    else
                        v(ii) = str2double(xi);
                    end
                end
            end
        else
            v = [];
        end
    catch
        % Fallback: strip brackets, split commas, drop 'null'
        s2 = replace(s, ["[", "]"], "");
        parts = strtrim(split(s2, ","));
        parts(parts=="null" | parts=="") = [];
        if ~isempty(parts)
            v = str2double(parts);
        end
    end
    v = v(~isnan(v));  % keep only valid numbers
end

% --- Expand to one value per row ---
valid_pid = ~ismissing(T.patient_id);
allPid  = [];
allFreq = [];

for i = find(valid_pid).'
    vec = parseSzFreqCell(T.sz_freqs(i));
    if ~isempty(vec)
        allPid  = [allPid;  repmat(T.patient_id(i), numel(vec), 1)];
        allFreq = [allFreq; vec(:)];
    end
end

% Handle case where no numeric values were found
if isempty(allFreq)
    T.avg_sz_freq = NaN(height(T),1);
    disp('No valid numeric entries in sz_freqs; avg_sz_freq set to NaN.');
    % writetable(T, 'with_avg.csv');  % optional
    return
end

TF = table(allPid, allFreq, 'VariableNames', {'patient_id','sz_freq'});

% --- Per-patient mean ---
G = groupsummary(TF, "patient_id", "mean", "sz_freq");
G.Properties.VariableNames{'mean_sz_freq'} = 'avg_sz_freq';

% --- Map back WITHOUT join (avoids missing-key restriction) ---
avg_col = nan(height(T),1);
[tf, loc] = ismember(T.patient_id, G.patient_id);
avg_col(tf) = G.avg_sz_freq(loc(tf));
T.avg_sz_freq = avg_col;

% Quick peek
disp(T(1:8, {'Record ID','patient_id','avg_sz_freq','report_SPORADIC_EPILEPTIFORM_DISCHARGES'}));

%% Bar + error bars (mean ± SEM) for avg_sz_freq by SPORADIC_EPILEPTIFORM_DISCHARGES

% Ensure the group variable is string, normalized
g = string(T.report_SPORADIC_EPILEPTIFORM_DISCHARGES);
g = lower(strtrim(g));

% Keep rows with a label and a numeric avg
keep = ~ismissing(g) & ~isnan(T.avg_sz_freq);

% Define groups (adjust the strings here if your labels differ)
present_vals = T.avg_sz_freq( keep & g == "present" );
absent_vals  = T.avg_sz_freq( keep & g == "absent"  );

% Basic stats
mu  = [ mean(present_vals,'omitnan'),  mean(absent_vals,'omitnan') ];
sd  = [  std(present_vals,'omitnan'),   std(absent_vals,'omitnan') ];
n   = [ sum(~isnan(present_vals)),      sum(~isnan(absent_vals))   ];
sem = sd ./ sqrt(max(n,1));  % guard against n=0

% Plot
figure; 
b = bar(1:2, mu); hold on
e = errorbar(1:2, mu, sem, 'k', 'LineStyle','none', 'LineWidth', 1.5);
set(gca, 'XTick', 1:2, 'XTickLabel', {'present','absent'});
ylabel('Average seizure frequency');
title('Avg seizure frequency by SPORADIC EPILEPTIFORM DISCHARGES (mean ± SEM)');

% Annotate with n
yl = ylim;
ytext = mu + sem + 0.05*range(yl);
for i = 1:2
    text(i, ytext(i), sprintf('n=%d', n(i)), ...
         'HorizontalAlignment','center', 'VerticalAlignment','bottom');
end
box off; hold off

% Optional: Welch’s t-test (unequal variances)
if n(1) > 1 && n(2) > 1
    [~, p, ~, stats] = ttest2(present_vals, absent_vals, 'Vartype','unequal');
    fprintf('Welch t-test: t = %.3f, df = %.1f, p = %.3g\n', stats.tstat, stats.df, p);
end

