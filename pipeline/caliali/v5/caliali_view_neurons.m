function caliali_view_neurons(basedir, start)
% CALIALI_VIEW_NEURONS  Browse the extracted neurons of one day, one at a time.
%   caliali_view_neurons(day)          % day folder, e.g. from caliali_animal_days(animal)
%   caliali_view_neurons(p)            % or the struct from caliali_day_paths(day)
%   caliali_view_neurons(day, 25)      % start at neuron 25
% Needs the outputs of caliali_summary (CaliAli_summary_*.mat and CaliAli_traces_*.mat).
%
% Left:   correlation image with all contours (grey), the selected neuron (yellow) and its
%         nearest neighbour (orange); a zoom around them; the neuron's numbers per session.
% Right:  its activity in every session (same y-scale in all sessions): C_raw (grey),
%         denoised C (colour), deconvolved events (ticks). Gaps = stripe frames (NaN).
% Keys:   right/down arrow = next, left/up arrow = previous, Home/End = first/last.
% Order:  neuron number, sessions active, SNR, events, size, nearest pairs (closest first:
%         each neuron followed by its nearest neighbour; close + correlated = possibly one cell
%         split in two), or spatial (snake path through the field, neighbours consecutive).
% Click:  on the overview or the zoom to select the neuron closest to the click (if nothing
%         happens, switch off the zoom/pan tools of the axes toolbar).
if nargin < 2 || isempty(start), start = 1; end
% reads only the results on the NAS: no CaliAli needed, nothing created (works on any computer
% that mounts the NAS at the same path)
if isstruct(basedir), p = basedir; else, p = caliali_day_paths(basedir, false); end
f_sum = fullfile(p.export, ['CaliAli_summary_' p.tag '.mat']);
f_tr  = fullfile(p.export, ['CaliAli_traces_'  p.tag '.mat']);
assert(isfile(f_sum) && isfile(f_tr), 'Run caliali_summary for %s first.', p.day);
sm = load(f_sum, 'summary');  sm = sm.summary;
assert(isfield(sm, 'diameter_px'), ['The summary of %s predates the size measures; ' ...
    'run caliali_summary for this day again (it does not redo the extraction).'], p.day);
T  = load(f_tr, 'traces', 'fps');  tr = T.traces;  fps = T.fps;

K = sm.n_neurons;  S_n = numel(tr);
Cn = sm.correlation_image;
cols = lines(max(S_n, 1));
ORANGE = [0.95 0.55 0.1];
ZOOM = 30;                                         % half-width of the zoom window (px)
order = (1:K)';  pos = min(max(1, start), K);

% nearest neighbour of every neuron (centroid distance) and trace correlation helper
xy = sm.centroid_xy;
Dm = hypot(xy(:, 1) - xy(:, 1)', xy(:, 2) - xy(:, 2)');
Dm(1:K + 1:end) = inf;
[nn_dist, nn] = min(Dm, [], 2);
C_all = cat(2, tr.C);                              % neurons x all frames of all sessions
ok_frames = ~any(isnan(C_all), 1);                 % stripe frames (NaN) left out

% ---- layout --------------------------------------------------------------------------------
fig = figure('Name', sprintf('Neurons: %s / %s', p.animal, p.day), 'NumberTitle', 'off', ...
    'Units', 'normalized', 'Position', [0.04 0.06 0.92 0.84], 'Color', 'w', 'KeyPressFcn', @on_key);
axFull = axes('Parent', fig, 'Position', [0.02 0.45 0.31 0.50]);
axZoom = axes('Parent', fig, 'Position', [0.02 0.08 0.15 0.30]);
info   = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.18 0.06 0.19 0.34], ...
    'HorizontalAlignment', 'left', 'FontName', 'FixedWidth', 'FontSize', 9, 'BackgroundColor', 'w');
axT = gobjects(1, S_n);
h_tr = 0.86 / S_n;
for s = 1:S_n
    axT(s) = axes('Parent', fig, 'Position', [0.40, 0.97 - s * h_tr, 0.58, h_tr - 0.045]);
end
uicontrol(fig, 'Style', 'pushbutton', 'String', '< Previous', 'Units', 'normalized', ...
    'Position', [0.40 0.01 0.07 0.04], 'Callback', @(src, ~) step(src, -1), 'KeyPressFcn', @on_key);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Next >', 'Units', 'normalized', ...
    'Position', [0.48 0.01 0.07 0.04], 'Callback', @(src, ~) step(src, +1), 'KeyPressFcn', @on_key);
uicontrol(fig, 'Style', 'text', 'String', 'Go to neuron', 'Units', 'normalized', ...
    'Position', [0.57 0.012 0.06 0.03], 'BackgroundColor', 'w', 'HorizontalAlignment', 'right');
h_goto = uicontrol(fig, 'Style', 'edit', 'String', '', 'Units', 'normalized', ...
    'Position', [0.635 0.01 0.05 0.04], 'Callback', @on_goto);
uicontrol(fig, 'Style', 'text', 'String', 'Order by', 'Units', 'normalized', ...
    'Position', [0.70 0.012 0.05 0.03], 'BackgroundColor', 'w', 'HorizontalAlignment', 'right');
uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [0.755 0.01 0.17 0.04], ...
    'String', {'neuron number', 'sessions active', 'SNR (best session)', 'events (all sessions)', ...
               'size (4 x SD, largest first)', 'nearest pairs (closest first)', ...
               'spatial (bands, top-left to bottom-right)'}, ...
    'Callback', @on_order);

% ---- static overview: correlation image with every contour ---------------------------------
imagesc(axFull, Cn); axis(axFull, 'image'); colormap(axFull, gray); hold(axFull, 'on');
for k = 1:K
    c = sm.contours_xy{k};
    if ~isempty(c), plot(axFull, c(:, 1), c(:, 2), 'Color', [0.55 0.65 0.8], 'LineWidth', 0.5, 'HitTest', 'off'); end
end
h_nn  = plot(axFull, NaN, NaN, 'Color', ORANGE, 'LineWidth', 1.5, 'HitTest', 'off');
h_sel = plot(axFull, NaN, NaN, 'y', 'LineWidth', 2, 'HitTest', 'off');
h_box = plot(axFull, NaN, NaN, 'c-', 'LineWidth', 1, 'HitTest', 'off');
set(axFull, 'XTick', [], 'YTick', []);
title(axFull, sprintf('%s: %d neurons (click on the image to select a neuron)', p.day, K), ...
    'Interpreter', 'none', 'FontSize', 9);
% clicks: on the image, and on the axes (clicks on a contour line are passed to the axes, not to
% the image below it); MATLAB's built-in axes interactions would also intercept clicks
set(findobj(axFull, 'Type', 'image'), 'ButtonDownFcn', @on_click);
set([axFull, axZoom], 'ButtonDownFcn', @on_click);
disableDefaultInteractivity(axFull);
disableDefaultInteractivity(axZoom);

go(pos);

% ---- drawing --------------------------------------------------------------------------------
    function go(new_pos)
        pos = min(max(1, new_pos), K);
        k = order(pos);  j = nn(k);
        c = sm.contours_xy{k};  cj = sm.contours_xy{j};
        cx = xy(k, 1);  cy = xy(k, 2);
        set_line(h_sel, c);  set_line(h_nn, cj);
        x0 = cx - ZOOM;  x1 = cx + ZOOM;  y0 = cy - ZOOM;  y1 = cy + ZOOM;
        set(h_box, 'XData', [x0 x1 x1 x0 x0], 'YData', [y0 y0 y1 y1 y0]);

        % zoom: neighbours grey, nearest orange, selected yellow
        cla(axZoom);
        him = imagesc(axZoom, Cn); axis(axZoom, 'image'); colormap(axZoom, gray); hold(axZoom, 'on');
        set(him, 'ButtonDownFcn', @on_click);           % clicking in the zoom selects a neuron too
        near = find(abs(xy(:, 1) - cx) < 2 * ZOOM & abs(xy(:, 2) - cy) < 2 * ZOOM)';
        for i = near
            ci = sm.contours_xy{i};
            if isempty(ci) || i == k || i == j, continue, end
            plot(axZoom, ci(:, 1), ci(:, 2), 'Color', [0.55 0.65 0.8], 'LineWidth', 0.8, 'HitTest', 'off');
        end
        if ~isempty(cj), plot(axZoom, cj(:, 1), cj(:, 2), 'Color', ORANGE, 'LineWidth', 1.5, 'HitTest', 'off'); end
        if ~isempty(c),  plot(axZoom, c(:, 1),  c(:, 2),  'y', 'LineWidth', 2, 'HitTest', 'off'); end
        xlim(axZoom, [x0 x1]);  ylim(axZoom, [y0 y1]);
        set(axZoom, 'XTick', [], 'YTick', []);
        title(axZoom, sprintf('neuron %d (yellow), nearest %d (orange)', k, j), 'FontSize', 9);

        % trace correlation with the nearest neighbour, over all sessions
        r = corrcoef(double(C_all(k, ok_frames)), double(C_all(j, ok_frames)));
        r = r(1, 2);

        % traces, same y-scale in every session
        lo = inf;  hi = -inf;
        for s = 1:S_n
            lo = min(lo, min(tr(s).C_raw(k, :), [], 'omitnan'));
            hi = max(hi, max(tr(s).C_raw(k, :), [], 'omitnan'));
        end
        if ~isfinite(lo) || ~isfinite(hi) || hi <= lo, lo = 0; hi = 1; end
        pad = 0.08 * (hi - lo);
        for s = 1:S_n
            ax = axT(s);  cla(ax);  hold(ax, 'on');
            t = (0:size(tr(s).C, 2) - 1) / fps / 60;
            plot(ax, t, tr(s).C_raw(k, :), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
            plot(ax, t, tr(s).C(k, :), 'Color', cols(s, :), 'LineWidth', 1);
            ev = find(tr(s).S(k, :) > 0);
            if ~isempty(ev)
                plot(ax, [t(ev); t(ev)], repmat([lo - pad; lo - 0.4 * pad], 1, numel(ev)), ...
                    'Color', [0.1 0.1 0.1], 'LineWidth', 0.5);
            end
            xlim(ax, [0, max(t(end), eps)]);  ylim(ax, [lo - pad, hi + pad]);
            if s == S_n, xlabel(ax, 'time in session (min)'); else, set(ax, 'XTickLabel', []); end
            ylabel(ax, 'C');
            if sm.active(k, s), act = 'ACTIVE'; else, act = 'not active'; end
            title(ax, sprintf('%s   |   %d events, %.2f /min, SNR %.1f   |   %s   |   %d stripe frames', ...
                sm.sessions{s}, sm.events(k, s), sm.events_per_min(k, s), sm.snr(k, s), act, ...
                numel(tr(s).stripe_frames)), 'Interpreter', 'none', 'FontSize', 9, 'FontWeight', 'normal');
            box(ax, 'on');
        end

        % numbers
        lines_ = {sprintf('Neuron %d   (%d of %d)', k, pos, K), '', ...
                  sprintf('centroid x, y : %.0f, %.0f', cx, cy), ...
                  sprintf('area          : %d px', sm.area_px(k)), ...
                  sprintf('diameter      : %.1f px (contour)', sm.diameter_px(k)), ...
                  sprintf('4 x SD        : %.1f px', sm.diameter_4sigma_px(k)), ...
                  sprintf('3 x gSig      : %.1f px', 3 * sm.gSig), '', ...
                  sprintf('nearest       : %d at %.1f px', j, nn_dist(k)), ...
                  sprintf('trace corr    : %.2f', r), '', ...
                  sprintf('active in     : %d of %d sessions', sm.n_sessions_active(k), S_n), ...
                  sprintf('(active = >= %d events)', sm.min_events), '', ...
                  sprintf('%-11s %5s %6s %5s', 'session', 'ev', 'ev/min', 'SNR')};
        for s = 1:S_n
            name = sm.sessions{s};  if numel(name) > 11, name = name(1:11); end
            lines_{end + 1} = sprintf('%-11s %5d %6.2f %5.1f', name, sm.events(k, s), ...
                sm.events_per_min(k, s), sm.snr(k, s)); %#ok<AGROW>
        end
        set(info, 'String', lines_);
        set(h_goto, 'String', num2str(k));
    end

    function set_line(h, c)
        if isempty(c), set(h, 'XData', NaN, 'YData', NaN);
        else, set(h, 'XData', c(:, 1), 'YData', c(:, 2)); end
    end

% ---- controls -------------------------------------------------------------------------------
    function release_focus(src)
        % take the keyboard focus away from a control, so the arrow keys reach the figure again
        set(src, 'Enable', 'off');  drawnow;  set(src, 'Enable', 'on');
        figure(fig);
    end

    function step(src, d)
        go(pos + d);
        release_focus(src);
    end

    function on_key(~, evt)
        switch evt.Key
            case {'rightarrow', 'downarrow'}, go(pos + 1);
            case {'leftarrow', 'uparrow'},    go(pos - 1);
            case 'home',                      go(1);
            case 'end',                       go(K);
        end
    end

    function on_goto(src, ~)
        k = round(str2double(src.String));
        if isfinite(k) && k >= 1 && k <= K
            go(find(order == k, 1));
        end
        release_focus(src);
    end

    function on_order(src, ~)
        k_now = order(pos);
        switch src.Value
            case 1, order = (1:K)';
            case 2, [~, order] = sortrows([-sm.n_sessions_active, -max(sm.snr, [], 2)]);
            case 3, [~, order] = sort(max(sm.snr, [], 2), 'descend');
            case 4, [~, order] = sort(sum(sm.events, 2), 'descend');
            case 5
                d = sm.diameter_4sigma_px;  d(isnan(d)) = -inf;
                [~, order] = sort(d, 'descend');
            case 6, order = nearest_pairs();
            case 7, order = spatial_order();
        end
        order = order(:);
        go(find(order == k_now, 1));                   % stay on the same neuron
        release_focus(src);
    end

    function o = nearest_pairs()
        % closest pairs first: each neuron followed by its nearest neighbour (if not shown yet)
        [~, idx] = sort(nn_dist);
        seen = false(K, 1);  o = zeros(K, 1);  n = 0;
        for k = idx'
            for i = [k, nn(k)]
                if ~seen(i)
                    n = n + 1;  o(n) = i;  seen(i) = true;
                end
            end
        end
    end

    function o = spatial_order()
        % snake path through the field: horizontal bands about two cell diameters high, top to
        % bottom, alternating left->right and right->left, so consecutive neurons are neighbours
        band_h = 2 * 3 * sm.gSig;                      % 2 x expected neuron diameter (px)
        b = floor(xy(:, 2) / band_h);
        x = xy(:, 1);
        x(mod(b, 2) == 1) = -x(mod(b, 2) == 1);        % odd bands run right to left
        [~, o] = sortrows([b, x]);
    end

    function on_click(~, evt)
        pt = evt.IntersectionPoint(1:2);               % jump to the neuron closest to the click
        [~, k] = min(hypot(xy(:, 1) - pt(1), xy(:, 2) - pt(2)));
        go(find(order == k, 1));
    end
end
