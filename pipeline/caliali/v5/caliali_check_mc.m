function T = caliali_check_mc(animal_dir)
% CALIALI_CHECK_MC  Verify, for every session of every day, that the steps before alignment are
% complete and were done in the right order (stripe frames replaced BEFORE motion correction).
%   T = caliali_check_mc('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26')
% Columns:
%   con_frames      frames in the session's _con.mat (downsampled movie)
%   stripe_source   'MiniAn CSV' (same list as the MiniAn pipeline) / 'detected' / 'MISSING'
%   raw_stripes     stripe frames in the raw video; kept = among the downsampled frames (replaced)
%   stripes_current stripe record belongs to this very _con.mat (date and size match)
%   mc_after        the session's movie (_con_mc_crop.mat once cropped, else _con_mc.mat) was
%                   written after the stripe replacement
%   mc_frames_ok    that movie has as many frames as _con.mat
%   mc_filled       that movie was written to the end (not an empty shell of a stopped run)
%   cropped / det   crop marker of the day / detrended file present
% A session is OK when the record is current, motion correction came after it, and frames match.
caliali_setup();
days = caliali_animal_days(animal_dir);
rows = {};
for d = 1:numel(days)
    p = caliali_day_paths(days{d});
    for k = 1:numel(p.sessions)
        [~, sess] = fileparts(fileparts(p.sessions{k}));
        con = p.con{k};
        if isfile(p.crop{k}), mc = p.crop{k}; else, mc = p.mc{k}; end   % current movie of the session
        det = strrep(p.align_input{k}, '.mat', '_det.mat');
        con_frames = NaN; if isfile(con), con_frames = nframes(con); end
        src = 'MISSING'; raw_n = NaN; kept = NaN; current = false; mc_after = false;
        if isfile(p.stripes{k})
            s = load(p.stripes{k});
            if startsWith(s.source, 'MiniAn'), src = 'MiniAn CSV'; else, src = 'detected'; end
            if s.rejected, src = [src ' (rejected)']; end
            raw_n = numel(s.raw_stripe_frames);  kept = numel(s.trace_stripe_frames);
            if isfile(con)
                dc = dir(con);
                current = abs(s.con_datenum - dc.datenum) < 1e-9 && s.con_bytes == dc.bytes;
            end
            if isfile(mc)
                dm = dir(mc);
                mc_after = dm.datenum >= s.con_datenum;
            end
        end
        mc_ok = isfile(mc) && isfile(con) && nframes(mc) == con_frames;
        mc_filled = caliali_file_complete(mc);                % not an empty pre-allocated shell
        ok = current && mc_after && mc_ok && mc_filled;
        rows(end+1, :) = {p.day, sess, con_frames, src, raw_n, kept, current, mc_after, mc_ok, mc_filled, ...
                          isfile(p.crop_done), isfile(det), ok}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'day', 'session', 'con_frames', 'stripe_source', ...
    'raw_stripes', 'kept_replaced', 'stripes_current', 'mc_after', 'mc_frames_ok', 'mc_filled', 'cropped', 'det', 'OK'});
disp(T);
fprintf('%d of %d sessions OK\n', sum(T.OK), height(T));
end

function n = nframes(f)
d = get_data_dimension(f);
n = d(3);
end
