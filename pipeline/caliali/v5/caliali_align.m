function CaliAli_options = caliali_align(basedir)
% CALIALI_ALIGN  Unattended step 3 for one day: detrending, projections, inter-session alignment.
%   caliali_align('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19')
% Needs the crop to be done first (caliali_crop_animal). Instead of opening figures, it saves the
% alignment QC to <NAS day>/CaliAli_results/: CaliAli_alignment_<tag>.txt (BV score, metrics),
% _scores.png (correlation per session) and _projections.png (sessions x alignment stages).
% Skipped if the aligned file is complete (rerunning alignment would erase the stored seeds).
% Days with a single session go through the same CaliAli function: it detrends the session and
% skips the inter-session steps, so the day gets the same files (_det, _Aligned) as any other day
% and is ready for extraction and for a later cross-day registration. Needs the catpad_centered fix.
caliali_setup();
p = caliali_day_paths(basedir);
assert(isfile(p.crop_done), ['Day %s is not cropped yet: run caliali_crop_animal first.\n' ...
    'To align without a manual crop, create the marker:  fclose(fopen(''%s'', ''w''))'], p.day, p.crop_done);
lock = caliali_lock(p); %#ok<NASGU>

if is_aligned(p)
    fprintf('Alignment already done: %s\n', p.aligned);
    CaliAli_options = CaliAli_load(p.aligned, 'CaliAli_options');
else
    mc = p.align_input;               % _mc_crop.mat (or _mc.mat for days cropped in place earlier)
    assert(all(cellfun(@caliali_file_complete, mc)), 'Day %s: input files missing or incomplete.', p.day);
    % CaliAli takes any existing _det.mat as finished, also an empty one pre-allocated by a run
    % that was stopped: delete those first so they are redone
    det = strrep(mc, '.mat', '_det.mat');
    for k = 1:numel(det)
        if isfile(det{k}) && ~caliali_file_complete(det{k})
            fprintf('Unfinished file from an interrupted run, deleted to be redone: %s\n', det{k});
            delete(det{k});
        end
    end
    CaliAli_options = CaliAli_parameters(CaliAli_my_parameters());
    CaliAli_options.inter_session_alignment.input_files = mc;
    CaliAli_options = CaliAli_align_sessions(CaliAli_options);
    assert(strcmp(CaliAli_options.inter_session_alignment.out_aligned_sessions, p.aligned), ...
        'Aligned file is not where caliali_day_paths expects it: %s', ...
        CaliAli_options.inter_session_alignment.out_aligned_sessions);
end
save(p.options, 'CaliAli_options');
try
    save_alignment_qc(p, CaliAli_options);
catch err
    warning('Alignment QC figures could not be saved: %s', err.message);
end
if numel(p.sessions) == 1
    fprintf('Single-session day %s: no inter-session alignment, movie prepared for extraction.\n', p.day);
else
    fprintf('Alignment done: %s (BV score %s)\n', p.day, num2str(CaliAli_options.inter_session_alignment.BV_score, 3));
end
end

function tf = is_aligned(p)
tf = false;
if isfile(p.aligned)
    try, tf = logical(CaliAli_load(p.aligned, 'alignment_completed')); catch, tf = false; end
end
end

function save_alignment_qc(p, o)
a = o.inter_session_alignment;
if ~isfolder(p.export), mkdir(p.export); fileattrib(p.export, '+w', 'g'); end
base = fullfile(p.export, ['CaliAli_alignment_' p.tag]);
fid = fopen([base '.txt'], 'w');
if numel(p.sessions) == 1
    fprintf(fid, 'Single session: no inter-session alignment (BV score and metrics not applicable).\n');
else
    fprintf(fid, 'BV score: %s\n\n%s', num2str(a.BV_score, 5), evalc('disp(a.alignment_metrics)'));
end
fclose(fid);

vis = get(groot, 'DefaultFigureVisible');                 % draw everything invisibly
set(groot, 'DefaultFigureVisible', 'off');
restore = onCleanup(@() set(groot, 'DefaultFigureVisible', vis)); %#ok<NASGU>

if numel(p.sessions) > 1                % scores compare sessions: meaningless for one session
    plot_alignment_scores(o);
    set(gcf, 'InvertHardcopy', 'off');  % save the dark-mode colours as drawn (else white labels on white)
    print(gcf, [base '_scores.png'], '-dpng', '-r110');
    close(gcf);
end

P = a.P;                                                   % table: one column per alignment stage
n_stage = size(P, 2);
n_sess  = size(P.(1)(1, :).(5){1, 1}, 4);
fig = figure('Position', [0 0 320 * n_stage, 280 * n_sess]);
for i = 1:n_sess
    for j = 1:n_stage
        ax = subplot(n_sess, n_stage, (i - 1) * n_stage + j);
        image(ax, mat2gray(P.(j)(1, :).(5){1, 1}(:, :, :, i)));
        axis(ax, 'image'); set(ax, 'XTick', [], 'YTick', []);
        if i == 1, title(ax, P.Properties.VariableNames{j}, 'Interpreter', 'none'); end
        if j == 1, ylabel(ax, sprintf('session %d', i)); end
    end
end
print(fig, [base '_projections.png'], '-dpng', '-r100');
close(fig);
end
