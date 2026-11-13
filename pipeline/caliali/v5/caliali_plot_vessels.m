function caliali_plot_vessels(basedir, which_session, save_png, style)
% CALIALI_PLOT_VESSELS  Neurons and blood vessels as CaliAli saw them, with the extracted neurons.
%   caliali_plot_vessels(day)              % CaliAli's projection image, averaged over the sessions
%   caliali_plot_vessels(p)                % or the struct from caliali_day_paths(day, false)
%   caliali_plot_vessels(day, 2)           % session 2 only (e.g. to compare sessions)
%   caliali_plot_vessels(day, 'mean', true)   % also save CaliAli_vessels_<animal>_<day>.png
%   caliali_plot_vessels(day, 'mean', false, 'vessels')   % grey vessel image only, instead
%
% Left:  CaliAli's own fused image of the final alignment stage: the neuron projection and the
%        blood-vessel image in two colour channels, exactly the image stored during processing
%        (the one in the alignment QC, CaliAli_alignment_*_projections.png).
% Right: the same image with the contour of every extracted neuron (50% of peak) in white.
% The two panels are linked: zoom or pan in one and the other follows.
%
% Reads only results on the NAS (CaliAli_results: CaliAli_summary_*.mat, CaliAli_options_*.mat),
% so it works on any computer that mounts the NAS at the same path; nothing is created.
% Needs caliali_extract and caliali_summary to have run for the day.
if nargin < 2 || isempty(which_session), which_session = 'mean'; end
if nargin < 3 || isempty(save_png), save_png = false; end
if nargin < 4 || isempty(style), style = 'fused'; end
style = validatestring(style, {'fused', 'vessels'});
if isstruct(basedir), p = basedir; else, p = caliali_day_paths(basedir, false); end
f_sum = fullfile(p.export, ['CaliAli_summary_' p.tag '.mat']);
f_opt = fullfile(p.export, ['CaliAli_options_' p.tag '.mat']);
assert(isfile(f_sum) && isfile(f_opt), ...
    'Missing CaliAli_summary / CaliAli_options in %s: run caliali_extract and caliali_summary first.', p.export);
sm = load(f_sum, 'summary');  sm = sm.summary;
o  = load(f_opt, 'CaliAli_options');  P = o.CaliAli_options.inter_session_alignment.P;

% images of the LAST alignment stage (what the aligned movie and the extraction are based on):
%   P.(stage)(1,:).(5) = fused neuron + vessel image, rows x cols x 3 (RGB) x sessions
%   P.(stage)(1,:).(2) = blood-vessel image, rows x cols x sessions   (both as in CaliAli's BV_gray2RGB)
stage = P.(size(P, 2));
if strcmp(style, 'fused')
    A = double(stage(1, :).(5){1, 1});
    n_ses = size(A, 4);
else
    A = double(stage(1, :).(2){1, 1});
    n_ses = size(A, 3);
end
if ischar(which_session) || isstring(which_session)
    assert(strcmpi(which_session, 'mean'), 'Second argument: a session number or ''mean''.');
    if strcmp(style, 'fused')
        img = mat2gray(mean(A, 4));
    else
        img = zeros(size(A, 1), size(A, 2));
        for k = 1:n_ses, img = img + mat2gray(A(:, :, k)); end   % each session scaled on its own
        img = mat2gray(img / n_ses);
    end
    what = sprintf('mean of %d session(s)', n_ses);
else
    assert(which_session >= 1 && which_session <= n_ses, 'The day has %d session(s).', n_ses);
    if strcmp(style, 'fused'), img = mat2gray(A(:, :, :, which_session));
    else, img = mat2gray(A(:, :, which_session)); end
    what = sprintf('session %d: %s', which_session, sm.sessions{which_session});
end
if ~isequal([size(img, 1), size(img, 2)], sm.dims_hw)
    warning(['The projection image is %d x %d but the neurons are in a %d x %d frame; the overlay ' ...
             'may be misplaced.'], size(img, 1), size(img, 2), sm.dims_hw(1), sm.dims_hw(2));
end

fig = figure('Name', sprintf('Neurons and vessels: %s / %s', p.animal, p.day), 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [0.05 0.1 0.9 0.75], 'Color', 'w');
ax1 = subplot(1, 2, 1);
show(ax1, img, style);
if strcmp(style, 'fused')
    title(sprintf('CaliAli projection image: neurons + vessels (%s)', what), 'Interpreter', 'none');
else
    title(sprintf('Blood vessels (%s)', what), 'Interpreter', 'none');
end
ax2 = subplot(1, 2, 2);
show(ax2, img, style);  hold(ax2, 'on');
for k = 1:sm.n_neurons
    c = sm.contours_xy{k};
    if ~isempty(c), plot(ax2, c(:, 1), c(:, 2), 'Color', [1 1 1], 'LineWidth', 0.8); end
end
title(ax2, sprintf('%d extracted neurons (contours at %g of peak)', sm.n_neurons, sm.contour_level));
linkaxes([ax1, ax2]);
sgtitle(sprintf('%s / %s', p.animal, p.day), 'Interpreter', 'none');

if save_png
    out = fullfile(p.export, ['CaliAli_vessels_' p.tag '.png']);
    set(fig, 'InvertHardcopy', 'off');
    print(fig, out, '-dpng', '-r150');
    fprintf('Saved %s\n', out);
end
end

function show(ax, img, style)
if strcmp(style, 'fused')
    image(ax, img);
else
    imagesc(ax, img);  colormap(ax, gray);
end
axis(ax, 'image');  axis(ax, 'off');
end
