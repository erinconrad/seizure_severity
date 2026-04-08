function compare_report_csvs_for_outpatient_discrepancies(reportCsv1, reportCsv2, outCsv)
% Compare two report CSVs row-by-row at the EEG key level
% and identify why outpatient classification differs.
%
% Usage:
% compare_report_csvs_for_outpatient_discrepancies( ...
%     '../data/clinical_data_deidentified_A.csv', ...
%     '../data/clinical_data_deidentified_B.csv', ...
%     '../output/report_compare_outpatient_discrepancies.csv');

T1 = readtable(reportCsv1, 'TextType','string', 'VariableNamingRule','preserve');
T2 = readtable(reportCsv2, 'TextType','string', 'VariableNamingRule','preserve');

req = ["patient_id","session_number","acquired_on","report_PATIENT_CLASS","jay_in_or_out"];
require_cols_local(T1, req, "T1");
require_cols_local(T2, req, "T2");

% Keep only columns needed for comparison
T1 = T1(:, req);
T2 = T2(:, req);

% Rename non-key columns so they can coexist after join
T1 = renamevars(T1, ["acquired_on","report_PATIENT_CLASS","jay_in_or_out"], ...
                   ["acquired_on_1","report_PATIENT_CLASS_1","jay_in_or_out_1"]);
T2 = renamevars(T2, ["acquired_on","report_PATIENT_CLASS","jay_in_or_out"], ...
                   ["acquired_on_2","report_PATIENT_CLASS_2","jay_in_or_out_2"]);

% Full outer join on EEG key
J = outerjoin(T1, T2, ...
    'Keys', {'patient_id','session_number'}, ...
    'MergeKeys', true, ...
    'Type', 'full');

% Presence flags
J.in_csv1 = ~ismissing(J.acquired_on_1) | ~ismissing(J.report_PATIENT_CLASS_1) | ~ismissing(J.jay_in_or_out_1);
J.in_csv2 = ~ismissing(J.acquired_on_2) | ~ismissing(J.report_PATIENT_CLASS_2) | ~ismissing(J.jay_in_or_out_2);

% Standardize strings
acq1   = lower(strtrim(fill_missing_string(J.acquired_on_1)));
acq2   = lower(strtrim(fill_missing_string(J.acquired_on_2)));
class1 = lower(strtrim(fill_missing_string(J.report_PATIENT_CLASS_1)));
class2 = lower(strtrim(fill_missing_string(J.report_PATIENT_CLASS_2)));
jay1   = lower(strtrim(fill_missing_string(J.jay_in_or_out_1)));
jay2   = lower(strtrim(fill_missing_string(J.jay_in_or_out_2)));

% Recompute outpatient logic exactly like your pipeline
J.isOutpt_site_1  = contains(acq1,"spe") | contains(acq1,"radnor");
J.isOutpt_site_2  = contains(acq2,"spe") | contains(acq2,"radnor");

J.isOutpt_class_1 = (class1 == "outpatient");
J.isOutpt_class_2 = (class2 == "outpatient");

J.isOutpt_jay_1   = (jay1 == "out");
J.isOutpt_jay_2   = (jay2 == "out");

J.isOutpt_any_1 = J.isOutpt_site_1 | J.isOutpt_class_1 | J.isOutpt_jay_1;
J.isOutpt_any_2 = J.isOutpt_site_2 | J.isOutpt_class_2 | J.isOutpt_jay_2;

% Field-level differences
J.diff_acquired_on = acq1 ~= acq2;
J.diff_patient_class = class1 ~= class2;
J.diff_jay_in_or_out = jay1 ~= jay2;
J.diff_outpt_call = J.isOutpt_any_1 ~= J.isOutpt_any_2;

% Human-readable discrepancy reason
reason = strings(height(J),1);

for i = 1:height(J)
    if ~J.in_csv1(i)
        reason(i) = "missing from csv1";
    elseif ~J.in_csv2(i)
        reason(i) = "missing from csv2";
    elseif J.diff_outpt_call(i)
        parts = strings(0,1);
        if J.diff_acquired_on(i)
            parts(end+1) = "acquired_on differs";
        end
        if J.diff_patient_class(i)
            parts(end+1) = "report_PATIENT_CLASS differs";
        end
        if J.diff_jay_in_or_out(i)
            parts(end+1) = "jay_in_or_out differs";
        end
        if isempty(parts)
            parts = "outpatient classification differs for unclear reason";
        end
        reason(i) = strjoin(parts, "; ");
    elseif J.diff_acquired_on(i) || J.diff_patient_class(i) || J.diff_jay_in_or_out(i)
        parts = strings(0,1);
        if J.diff_acquired_on(i)
            parts(end+1) = "acquired_on differs";
        end
        if J.diff_patient_class(i)
            parts(end+1) = "report_PATIENT_CLASS differs";
        end
        if J.diff_jay_in_or_out(i)
            parts(end+1) = "jay_in_or_out differs";
        end
        reason(i) = "fields differ but outpatient classification same: " + strjoin(parts, "; ");
    else
        reason(i) = "";
    end
end

J.discrepancy_reason = reason;

% Keep only rows that differ in some important way
keep = ~J.in_csv1 | ~J.in_csv2 | ...
       J.diff_outpt_call | ...
       J.diff_acquired_on | J.diff_patient_class | J.diff_jay_in_or_out;

D = J(keep, :);

% Sort so the most important rows are first
D = sortrows(D, {'diff_outpt_call','in_csv1','in_csv2'}, {'descend','ascend','ascend'});

% Print summary
fprintf('\n=== Summary comparing report CSVs ===\n');
fprintf('Rows/keys in csv1 only: %d\n', nnz(~J.in_csv2 & J.in_csv1));
fprintf('Rows/keys in csv2 only: %d\n', nnz(~J.in_csv1 & J.in_csv2));
fprintf('Shared keys with different outpatient call: %d\n', nnz(J.in_csv1 & J.in_csv2 & J.diff_outpt_call));
fprintf('Shared keys with same outpatient call but different fields: %d\n', ...
    nnz(J.in_csv1 & J.in_csv2 & ~J.diff_outpt_call & ...
        (J.diff_acquired_on | J.diff_patient_class | J.diff_jay_in_or_out)));

% Show the rows that actually change outpatient inclusion
D_out = D(D.diff_outpt_call | ~D.in_csv1 | ~D.in_csv2, :);
disp(D_out(:, {'patient_id','session_number', ...
    'acquired_on_1','acquired_on_2', ...
    'report_PATIENT_CLASS_1','report_PATIENT_CLASS_2', ...
    'jay_in_or_out_1','jay_in_or_out_2', ...
    'isOutpt_any_1','isOutpt_any_2', ...
    'discrepancy_reason'}));

if nargin >= 3 && strlength(string(outCsv)) > 0
    if ~exist(fileparts(outCsv), 'dir')
        mkdir(fileparts(outCsv));
    end
    writetable(D, outCsv);
    fprintf('Wrote discrepancy table: %s\n', outCsv);
end

end

function s = fill_missing_string(s)
s = string(s);
s(ismissing(s)) = "";
end

function require_cols_local(T, cols, name)
missing = setdiff(string(cols), string(T.Properties.VariableNames));
if ~isempty(missing)
    error('%s is missing required columns: %s', name, strjoin(missing, ", "));
end
end