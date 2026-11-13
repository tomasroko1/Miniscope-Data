function caliali_set_seeds(basedir)
% CALIALI_SET_SEEDS  Step 2 of 3 (server, no App Designer): interactive seed settings.
%   caliali_set_seeds('/mnt/NAS/.../R004_M4_26/T1_SD_CNO_2026_06_19')
% Normal MATLAB figure (works over x2go) with the same controls as CaliAli's app:
%   min_corr, min_pnr, gSig  -> type a value + Enter, the seeds update immediately
%   Draw mask / Whole field  -> seed mask on the correlation image
%   Save & close             -> writes them into the day's _Aligned.mat (read by CaliAli_cnmfe)
% Seeds are computed exactly as CaliAli does (get_seeds / Check_initialization_parameters).
caliali_setup();
p = caliali_day_paths(basedir);
assert(isfile(p.aligned), 'Not aligned yet; run caliali_run(animal_dir, ''align'') first.');
lock = caliali_lock(p); %#ok<NASGU>

o   = CaliAli_load(p.aligned, 'CaliAli_options');
Cn  = o.inter_session_alignment.Cn;
PNR = o.inter_session_alignment.PNR;
if isempty(o.cnmf.seed_mask), o.cnmf.seed_mask = true(size(Cn)); end
o.cnmf.seed_mask = logical(o.cnmf.seed_mask);
vmax = []; vmax_gSig = NaN; do_save = false;

fig = figure('Name', ['Seeds: ' p.animal ' / ' p.day], 'NumberTitle', 'off', ...
             'Units', 'normalized', 'Position', [0.05 0.15 0.9 0.65]);
ax = gobjects(1, 3);
for k = 1:3
    ax(k) = axes('Parent', fig, 'Position', [0.03 + (k-1)*0.325, 0.22, 0.3, 0.72]);
end
h.corr = add_edit('min_corr', o.cnmf.min_corr, 0.03);
h.pnr  = add_edit('min_pnr',  o.cnmf.min_pnr,  0.18);
h.gsig = add_edit('gSig',     o.cnmf.gSig,     0.33);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Update', 'Units', 'normalized', ...
          'Position', [0.48 0.06 0.06 0.06], 'Callback', @update);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Draw mask', 'Units', 'normalized', ...
          'Position', [0.55 0.06 0.08 0.06], 'Callback', @draw_mask);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Whole field', 'Units', 'normalized', ...
          'Position', [0.64 0.06 0.08 0.06], 'Callback', @clear_mask);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Save & close', 'Units', 'normalized', ...
          'Position', [0.76 0.06 0.09 0.06], 'Callback', @save_close);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', 'Units', 'normalized', ...
          'Position', [0.86 0.06 0.08 0.06], 'Callback', @(~,~) delete(fig));
h.info = uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'FontSize', 11, ...
                   'Position', [0.03 0.13 0.94 0.04], 'HorizontalAlignment', 'left');
fig.CloseRequestFcn = @(~,~) delete(fig);

update();
% Wait until the window is deleted (Save & close, Cancel or closing it). Not uiwait:
% drawpolygon's wait() uses uiwait/uiresume on the same figure, and finishing the
% polygon would also end an outer uiwait and close the window.
waitfor(fig);

if do_save
    if o.cnmf.gSig ~= o.inter_session_alignment.gSig
        warning('gSig for extraction (%.1f) differs from the one used for the projections (%.1f).', ...
                o.cnmf.gSig, o.inter_session_alignment.gSig);
    end
    CaliAli_save(p.aligned, 'CaliAli_options', o);
    fclose(fopen(p.seeds_done, 'w'));                 % marker for caliali_run / caliali_status
    s = CaliAli_load(p.aligned, 'CaliAli_options');
    fprintf('Saved: min_corr = %g, min_pnr = %g, gSig = %g, mask pixels = %d of %d\n', ...
        s.cnmf.min_corr, s.cnmf.min_pnr, s.cnmf.gSig, nnz(s.cnmf.seed_mask), numel(s.cnmf.seed_mask));
    if isfolder(p.extraction)
        fprintf(2, 'Old CNMF-E checkpoints exist; delete before extracting:\n  %s\n', p.extraction);
    end
    fprintf('Next: caliali_run(animal_dir, ''extract'') or caliali_extract(''%s'')\n', p.nas);
else
    fprintf('Nothing saved.\n');
end

    function e = add_edit(label, value, x)
        uicontrol(fig, 'Style', 'text', 'String', label, 'Units', 'normalized', ...
                  'Position', [x 0.065 0.06 0.04], 'HorizontalAlignment', 'right');
        e = uicontrol(fig, 'Style', 'edit', 'String', num2str(value), 'Units', 'normalized', ...
                      'Position', [x + 0.065 0.06 0.07 0.06], 'Callback', @update);
    end

    function update(~, ~)
        try
            do_update();
        catch err
            h.info.String = ['Update failed: ' err.message];
            fprintf(2, 'caliali_set_seeds update failed:\n%s\n', getReport(err));
        end
        drawnow;
    end

    function do_update()
        v = str2double({h.corr.String, h.pnr.String, h.gsig.String});
        if any(isnan(v)) || v(3) <= 0
            h.info.String = 'Invalid value: numbers only, gSig > 0';
            return
        end
        o.cnmf.min_corr = v(1); o.cnmf.min_pnr = v(2); o.cnmf.gSig = v(3);
        if o.cnmf.gSig ~= vmax_gSig                   % local maxima depend only on gSig
            vmax = CaliAli_get_local_maxima(o);
            vmax_gSig = o.cnmf.gSig;
        end
        % same as CaliAli's preview app (estimate_Corr_PNR): each panel shows its own threshold
        cand    = (vmax == Cn .* PNR) & o.cnmf.seed_mask;          % local maxima inside the mask
        cn_ind  = cand & Cn  >= o.cnmf.min_corr;                   % passes min_corr
        pnr_ind = cand & PNR >= o.cnmf.min_pnr;                    % passes min_pnr
        seed    = cn_ind & pnr_ind;                                % seeds CNMF-E will use
        imgs  = {Cn, PNR, Cn .* PNR};
        dots  = {cn_ind, pnr_ind, seed};
        cols  = {'r.', 'r.', 'm.'};
        clims = {[0 1], [0 prctile(PNR(:), 99.5)], [0 prctile(Cn(:) .* PNR(:), 99.5)]};   % saturate only the top 0.5%
        names = {sprintf('Correlation >= %g: %d', o.cnmf.min_corr, nnz(cn_ind)), ...
                 sprintf('PNR >= %g: %d', o.cnmf.min_pnr, nnz(pnr_ind)), ...
                 sprintf('Both (seeds): %d', nnz(seed))};
        for kk = 1:3
            imagesc(ax(kk), imgs{kk}); axis(ax(kk), 'image'); colormap(ax(kk), parula);
            if all(isfinite(clims{kk})) && clims{kk}(2) > 0, caxis(ax(kk), clims{kk}); end
            colorbar(ax(kk));
            hold(ax(kk), 'on');
            if ~all(o.cnmf.seed_mask(:))
                contour(ax(kk), double(o.cnmf.seed_mask), [0.5 0.5], 'w', 'LineWidth', 1);
            end
            [y, x] = find(dots{kk});
            plot(ax(kk), x, y, cols{kk}, 'MarkerSize', 10);
            hold(ax(kk), 'off');
            title(ax(kk), names{kk}); set(ax(kk), 'XTick', [], 'YTick', []);
        end
        h.info.String = sprintf('%d seeds  (min_corr %g, min_pnr %g, gSig %g; mask %d%% of field)', ...
            nnz(seed), o.cnmf.min_corr, o.cnmf.min_pnr, o.cnmf.gSig, round(100 * mean(o.cnmf.seed_mask(:))));
        if mod(max(1, round(o.cnmf.gSig)), 2) == 0
            h.info.String = [h.info.String sprintf(['   WARNING: round(gSig) = %d is even, CaliAli''s ' ...
                'local-maximum window is then off-centre and seeds become blobs'], round(o.cnmf.gSig))];
        end
        fprintf('%s\n', h.info.String);
    end

    function draw_mask(~, ~)
        h.info.String = 'Click points around the usable tissue on the LEFT image; double-click inside to finish';
        roi = drawpolygon(ax(1));
        wait(roi);
        o.cnmf.seed_mask = createMask(roi, size(Cn, 1), size(Cn, 2));
        delete(roi);
        update();
    end

    function clear_mask(~, ~)
        o.cnmf.seed_mask = true(size(Cn));
        update();
    end

    function save_close(~, ~)
        do_save = true;
        delete(fig);
    end
end
