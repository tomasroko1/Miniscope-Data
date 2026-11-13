function results = caliali_run(animal_dir, step)
% CALIALI_RUN  Run one unattended step on every day of an animal, from the MATLAB command window.
%   caliali_run('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26', 'mc')       % downsample + motion correction
%   caliali_run('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26', 'align')    % after caliali_crop_animal
%   caliali_run('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26', 'extract')  % after caliali_set_seeds per day
%   caliali_run('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26', 'summary')  % rebuild the summaries of
%       the extracted days with the current caliali_settings.m (e.g. after changing smin); always
%       redone, no extraction
% Each day is skipped if the step is already done or the day is not ready for it. An error on one
% day is caught, logged with its full report, and the run continues with the next day. Everything
% printed is also written to a log file in the animal's work folder. At the end (and in the
% returned table) you get one line per day: ok / skipped / not ready / FAILED.
%
% Full workflow per animal:
%   1. caliali_run(animal, 'mc')        unattended (can run over a weekend)
%   2. caliali_crop_animal(animal)      interactive, once per animal
%   3. caliali_run(animal, 'align')     unattended; alignment QC saved as images
%   4. caliali_set_seeds(day)           interactive, once per day
%   5. caliali_run(animal, 'extract')   unattended; also writes the summary per day
% caliali_status(animal) shows at any time which steps are done.
step = validatestring(step, {'mc', 'align', 'extract', 'summary'});
caliali_setup();
days = caliali_animal_days(animal_dir);
assert(~isempty(days), 'No day folders (with */My_V4_Miniscope) in %s', animal_dir);

p1 = caliali_day_paths(days{1});
log_file = fullfile(fileparts(p1.basedir), ...
    sprintf('caliali_run_%s_%s.log', step, datestr(now, 'yyyymmdd_HHMMSS')));
diary(log_file);
stop_log = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('caliali_run(%s) on %d days of %s\nLog: %s\n', step, numel(days), animal_dir, log_file);

results = cell(numel(days), 3);
for k = 1:numel(days)
    p  = caliali_day_paths(days{k});
    st = caliali_day_status(p);
    fprintf('\n===== %s  (%d/%d, %s) =====\n', p.day, k, numel(days), datestr(now, 'ddd HH:MM'));
    [ready, why] = check_day(step, st);
    if ~ready
        fprintf('%s\n', why);
        results(k, :) = {p.day, why, 0};
        continue
    end
    t = tic;
    try
        switch step
            case 'mc',      caliali_mc(days{k});
            case 'align',   caliali_align(days{k});
            case 'extract', caliali_extract(days{k}); caliali_summary(days{k});
            case 'summary', caliali_summary(days{k});
        end
        results(k, :) = {p.day, 'ok', toc(t) / 60};
    catch err
        fprintf(2, '\nFAILED %s:\n%s\n', p.day, getReport(err, 'extended', 'hyperlinks', 'off'));
        results(k, :) = {p.day, ['FAILED: ' err.message], toc(t) / 60};
    end
end

results = cell2table(results, 'VariableNames', {'day', 'result', 'minutes'});
fprintf('\n===== caliali_run(%s) finished %s =====\n', step, datestr(now, 'ddd HH:MM'));
disp(results);
fprintf('Log: %s\n', log_file);
end

function [ready, why] = check_day(step, st)
ready = false;
if st.locked
    why = 'skipped: locked (someone else is working on this day)';  return
end
switch step
    case 'mc'
        if st.mc,          why = 'skipped: already motion-corrected';                   return, end
    case 'align'
        if st.aligned,     why = 'skipped: already aligned';                            return, end
        if ~st.mc,         why = 'not ready: motion correction missing';                return, end
        if ~st.cropped,    why = 'not ready: not cropped (caliali_crop_animal)';        return, end
    case 'summary'
        if ~st.extracted,  why = 'not ready: not extracted';                            return, end
    case 'extract'
        if st.extracted,   why = 'skipped: already extracted (with the current seeds)';  return, end
        if ~st.aligned,    why = 'not ready: not aligned';                              return, end
        if ~st.seeds,      why = 'not ready: seeds not set (caliali_set_seeds)';        return, end
end
ready = true;
why = '';
end
