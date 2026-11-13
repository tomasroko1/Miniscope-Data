function result = caliali_result_file(p)
% CALIALI_RESULT_FILE  The final CNMF-E result of a day (p from caliali_day_paths): the newest
% checkpoint in the extraction folder, else the copy exported to <NAS day>/CaliAli_results.
L = dir(fullfile(p.extraction, 'frames_*', 'LOGS_*', '*.mat'));
L = L(~strcmp({L.name}, 'intermediate_results.mat'));
if ~isempty(L)
    [~, k] = max([L.datenum]);
    result = fullfile(L(k).folder, L(k).name);
else
    result = p.result;
end
assert(isfile(result), 'No extraction result found for %s; run caliali_extract first.', p.day);
end
