function summary = caliali_summary(basedir, min_events, contour_level, neuron_in, smin)
% CALIALI_SUMMARY  Contours, neuron counts and activity per session for one extracted day.
%   summary = caliali_summary('/mnt/NAS/.../R004_M4_26/T1_SD_CNO_2026_06_19')
%   summary = caliali_summary(basedir, 3, 0.5)
%     min_events, contour_level, smin : default to caliali_settings.m. If smin differs from CaliAli's
%                     (-5), the saved result's C and S are deconvolved again with it first; the
%                     result file itself is not changed, summary and traces carry the new C and S
%     neuron_in     : used by caliali_redeconv (a neuron already deconvolved with smin)
% Activity at stripe frames (caliali_stripes_day) is set to NaN for every neuron, and those
% frames are left out of the event counts, rates and SNR.
% Writes to <NAS day>/CaliAli_results/:
%   CaliAli_summary_<animal>_<day>.mat  contours, centroids, sizes, per-session activity (MATLAB struct)
%     sizes: diameter_px = diameter of a disc with the area inside the contour;
%            diameter_4sigma_px = 4 x the footprint's spatial SD (pixels below 10% of peak ignored);
%            compare with 3 x gSig (CNMF-E: gSig is about 1/3 of the neuron diameter)
%   CaliAli_summary_<animal>_<day>.csv  one row per neuron
%   CaliAli_summary_<animal>_<day>.png  overview figure
%   CaliAli_traces_<animal>_<day>.mat   C, C_raw, S per session with NaN at stripe frames,
%                                       plus the raw video frame of every trace frame
st = caliali_settings();                               % pipeline settings (caliali_settings.m)
if nargin < 2 || isempty(min_events),    min_events = st.min_events;       end
if nargin < 3 || isempty(contour_level), contour_level = st.contour_level; end
if nargin < 4, neuron_in = []; end
if nargin < 5 || isempty(smin), smin = st.smin; end
caliali_setup();
p = caliali_day_paths(basedir);

% ---- load the final result (newest CNMF-E checkpoint; else the exported copy) -----------------
result = caliali_result_file(p);
fprintf('Result: %s\n', result);
if isempty(neuron_in)
    R = load(result, 'neuron');
    neuron = R.neuron;
    if smin ~= st.smin_caliali        % the saved result is deconvolved with CaliAli's smin
        fprintf('Deconvolving again with smin = %g (CaliAli used %g) ...\n', smin, st.smin_caliali);
        neuron = postprocessDeconvolvedTraces(neuron, 'foopsi', 'ar2', smin);
    end
else
    neuron = neuron_in;
end

opt = neuron.CaliAli_options;
d1 = neuron.options.d1;  d2 = neuron.options.d2;
K  = size(neuron.A, 2);
F  = double(opt.inter_session_alignment.F(:)');                 % frames per session
assert(sum(F) == size(neuron.C, 2), 'Session frame counts do not add up to the trace length');
fps = opt.downsampling.sf / opt.downsampling.temporal_ds;       % frame rate of the traces
edges = [0 cumsum(F)];
S_n = numel(F);
sess_names = cellfun(@session_name, p.nas_sessions, 'UniformOutput', false);   % e.g. A_AV1

% ---- stripe frames: NaN for all neurons (as mask_frame_values in the MiniAn pipeline) ----------
stripe = false(1, sum(F));
raw_frame = cell(1, S_n);  n_stripe = zeros(1, S_n);  stripe_masking = true;
for s = 1:S_n
    if isfile(p.stripes{s})
        info = load(p.stripes{s});
        assert(numel(info.trace_raw_frame) == F(s), 'Stripe info of session %d does not match its frames', s);
        stripe(edges(s) + info.trace_stripe_frames) = true;
        raw_frame{s} = info.trace_raw_frame;
        n_stripe(s)  = numel(info.trace_stripe_frames);
    else
        warning('No stripe-frame info for %s (processed before stripe handling): NOT masked.', p.sessions{s});
        stripe_masking = false;
    end
end
C = full(neuron.C);  C_raw = full(neuron.C_raw);  S = full(neuron.S);
C(:, stripe) = NaN;  C_raw(:, stripe) = NaN;  S(:, stripe) = NaN;

% ---- per neuron and session: events, rate, SNR, active (stripe frames excluded) ---------------
events = zeros(K, S_n);  rate = zeros(K, S_n);  snr = zeros(K, S_n);
for s = 1:S_n
    idx = edges(s)+1:edges(s+1);
    ok  = idx(~stripe(idx));                                            % frames that count
    on  = S(:, ok) > 0;
    events(:, s) = sum(diff([zeros(K, 1) on], 1, 2) == 1, 2);          % event onsets
    rate(:, s)   = events(:, s) / (numel(ok) / fps / 60);               % events per minute
    noise        = std(C_raw(:, ok) - C(:, ok), 0, 2);
    snr(:, s)    = max(C(:, ok), [], 2) ./ max(noise, eps);
end
active = events >= min_events;
n_active = sum(active, 2);

% ---- contours and centroids (own implementation; CNMF-E's get_contours stops at the first
%      empty contour and leaves all later neurons without one) ------------------------------------
contours = cell(K, 1);  centroid = nan(K, 2);  area = zeros(K, 1);
diam50 = nan(K, 1);  sigma = nan(K, 1);
[xx, yy] = meshgrid(1:d2, 1:d1);
for k = 1:K
    a = reshape(full(neuron.A(:, k)), d1, d2);
    if max(a(:)) <= 0, continue; end
    centroid(k, :) = [sum(xx(:) .* a(:)), sum(yy(:) .* a(:))] / sum(a(:));   % [x y], weighted
    w = a .* (a >= 0.1 * max(a(:)));                                     % ignore faint tails
    sigma(k) = sqrt(sum(w(:) .* ((xx(:) - centroid(k, 1)).^2 + (yy(:) - centroid(k, 2)).^2)) ...
                    / sum(w(:)) / 2);                                    % spatial SD per axis (px)
    bw = a >= contour_level * max(a(:));
    bw = bwareafilt(bw, 1);                                              % largest connected piece
    area(k) = nnz(bw);
    diam50(k) = 2 * sqrt(area(k) / pi);                                  % diameter of a disc of that area
    B = bwboundaries(bw, 'noholes');
    if ~isempty(B)
        contours{k} = [B{1}(:, 2), B{1}(:, 1)];                          % [x y] pixel coordinates
    end
end

% ---- save ------------------------------------------------------------------------------------
summary = struct('animal', p.animal, 'day', p.day, 'sessions', {sess_names}, ...
    'frames_per_session', F, 'fps', fps, 'dims_hw', [d1 d2], 'n_neurons', K, ...
    'contours_xy', {contours}, 'contour_level', contour_level, 'centroid_xy', centroid, ...
    'area_px', area, 'diameter_px', diam50, 'footprint_sigma_px', sigma, 'diameter_4sigma_px', 4 * sigma, ...
    'gSig', neuron.options.gSig, 'events', events, 'events_per_min', rate, 'snr', snr, ...
    'min_events', min_events, 'active', active, 'n_sessions_active', n_active, ...
    'stripe_masking', stripe_masking, 'stripe_frames_per_session', n_stripe, ...
    'correlation_image', opt.inter_session_alignment.Cn, 'result_file', result, 'smin', smin);
if ~isfolder(p.export), mkdir(p.export); fileattrib(p.export, '+w', 'g'); end
base = fullfile(p.export, ['CaliAli_summary_' p.tag]);
save([base '.mat'], 'summary');

% traces per session, NaN at stripe frames, with the raw video frame of every trace frame
traces = struct('session', sess_names(:), 'C', [], 'C_raw', [], 'S', [], 'raw_frame', [], 'stripe_frames', []);
for s = 1:S_n
    idx = edges(s)+1:edges(s+1);
    traces(s).C             = single(C(:, idx));            % neurons x frames, denoised
    traces(s).C_raw         = single(C_raw(:, idx));        % neurons x frames, before denoising
    traces(s).S             = single(S(:, idx));            % neurons x frames, deconvolved events
    traces(s).raw_frame     = raw_frame{s};                 % 0-based frame in the session's raw AVIs
    traces(s).stripe_frames = find(stripe(idx));            % 1-based trace frames set to NaN
end
traces_note = ['Rows = neurons (same order as the summary), columns = frames at ' num2str(fps) ...
    ' fps. NaN = stripe frame (replaced before processing). raw_frame maps each column to the ' ...
    'raw video frame (0-based, all AVIs of the session concatenated, as in the MiniAn pipeline).'];
save(fullfile(p.export, ['CaliAli_traces_' p.tag '.mat']), 'traces', 'traces_note', 'fps', '-v7.3');

T = table((1:K)', centroid(:, 1), centroid(:, 2), area, round(diam50, 2), round(4 * sigma, 2), n_active, ...
          'VariableNames', {'neuron', 'centroid_x', 'centroid_y', 'area_px', 'diameter_px', ...
                            'diameter_4sigma_px', 'n_sessions_active'});
for s = 1:S_n
    tag = matlab.lang.makeValidName(sess_names{s});
    T.([tag '_events'])      = events(:, s);
    T.([tag '_events_min'])  = round(rate(:, s), 3);
    T.([tag '_snr'])         = round(snr(:, s), 2);
    T.([tag '_active'])      = active(:, s);
end
writetable(T, [base '.csv']);

% ---- overview figure -------------------------------------------------------------------------
fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 1600 1000]);
cols = lines(S_n + 1);
ax = subplot(2, 3, [1 4]);
imagesc(ax, opt.inter_session_alignment.Cn); axis(ax, 'image'); colormap(ax, gray); hold(ax, 'on');
for k = 1:K
    if ~isempty(contours{k})
        plot(ax, contours{k}(:, 1), contours{k}(:, 2), 'Color', cols(n_active(k) + 1, :), 'LineWidth', 1);
    end
end
h = gobjects(1, S_n + 1);
for n = 0:S_n, h(n + 1) = plot(ax, NaN, NaN, 'Color', cols(n + 1, :), 'LineWidth', 2); end
legend(ax, h, arrayfun(@(n) sprintf('active in %d', n), 0:S_n, 'UniformOutput', false), 'Location', 'southoutside');
title(ax, sprintf('%d neurons (contours at %g of peak)', K, contour_level), 'Interpreter', 'none');
set(ax, 'XTick', [], 'YTick', []);

ax = subplot(2, 3, 2);
bar(ax, sum(active, 1)); set(ax, 'XTickLabel', sess_names, 'TickLabelInterpreter', 'none');
ylabel(ax, 'active neurons'); title(ax, sprintf('Active per session (>= %d events)', min_events));

ax = subplot(2, 3, 3);
bar(ax, 0:S_n, histcounts(n_active, -0.5:1:S_n + 0.5));
xlabel(ax, 'number of sessions active'); ylabel(ax, 'neurons'); title(ax, 'Activity across sessions');

ax = subplot(2, 3, [5 6]);
[~, order] = sortrows([-n_active, -sum(rate, 2)]);
Cn_ = C(order, :) ./ max(max(C(order, :), [], 2), eps);
imagesc(ax, (0:size(C, 2) - 1) / fps / 60, 1:K, Cn_); colormap(ax, hot); caxis(ax, [0 1]);
hold(ax, 'on');
for s = 2:S_n, xline(ax, edges(s) / fps / 60, 'c', 'LineWidth', 1.5); end
xlabel(ax, 'time (min, sessions concatenated)'); ylabel(ax, 'neuron (sorted by sessions active)');
title(ax, 'Denoised activity (each neuron scaled to its maximum)');
print(fig, [base '.png'], '-dpng', '-r110');
close(fig);

% ---- console summary ---------------------------------------------------------------------------
fprintf('\n%s / %s: %d neurons (active = at least %d events in a session)\n', p.animal, p.day, K, min_events);
fprintf('  deconvolution threshold smin = %g (events of at least %g noise units)\n', smin, abs(smin));
fprintf(['  neuron size: median diameter %.1f px at the contour (%g of peak), %.1f px as 4 x footprint SD;\n' ...
         '  gSig %.1f corresponds to a diameter of about %.1f px (3 x gSig)\n'], median(diam50, 'omitnan'), ...
         contour_level, median(4 * sigma, 'omitnan'), neuron.options.gSig, 3 * neuron.options.gSig);
for s = 1:S_n
    fprintf('  %-12s %4d active, median %.2f events/min among active, %d stripe frames masked\n', ...
        sess_names{s}, sum(active(:, s)), median(rate(active(:, s), s)), n_stripe(s));
end
fprintf('  active in all %d sessions: %d | in 2 or more: %d | in exactly 1: %d | in none: %d\n', S_n, ...
    sum(n_active == S_n), sum(n_active >= 2), sum(n_active == 1), sum(n_active == 0));
fprintf('Saved %s.mat / .csv / .png and CaliAli_traces_%s.mat\n', base, p.tag);
end

function n = session_name(miniscope_dir)
[~, n] = fileparts(fileparts(miniscope_dir));   % .../<session>/My_V4_Miniscope -> <session>
end
