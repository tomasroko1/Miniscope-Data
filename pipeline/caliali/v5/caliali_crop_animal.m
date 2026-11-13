function caliali_crop_animal(animal_dir)
% CALIALI_CROP_ANIMAL  One crop for every session of every day of one animal.
%   caliali_crop_animal('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26')
% Run after caliali_run(animal_dir, 'mc'). Shows the middle frames of all sessions of all days
% (centred in a common frame, as CaliAli_crop does) and applies the ONE rectangle you draw to all
% of them. Keep it generous (exclude only what is outside the lens, stay inside every session's
% valid area): for cross-day registration, tissue cropped away on any day is lost for all days.
%
% Nothing is changed in place: every session's cropped movie is written to a NEW file,
%   My_V4_Miniscope_con_mc_crop.mat   (next to My_V4_Miniscope_con_mc.mat)
% and checked (frame count, last frame written). Only when ALL sessions of the run are cropped are
% the _mc.mat files deleted and the days marked as cropped. If the run stops before that, every
% _mc.mat is still intact; the next run discards the unfinished crop files and starts over.
% Days already marked as cropped are left out.
caliali_setup();
days = caliali_animal_days(animal_dir);
mc = {};  crop = {};  todo = {};  locks = {};
for d = 1:numel(days)
    p = caliali_day_paths(days{d});
    if isfile(p.crop_done)
        fprintf('Already cropped, skipped: %s\n', p.day);
        continue
    end
    leftovers = [strrep(p.mc, '.mat', '_tmp.mat'), strrep(p.mc, '.mat', '_bak.mat')];
    if any(cellfun(@isfile, leftovers))
        error(['Day %s has leftovers of an interrupted motion correction (_mc_tmp/_mc_bak); ' ...
               'run caliali_run(animal_dir, ''mc'') first, it redoes those sessions.'], p.day);
    end
    if ~all(cellfun(@caliali_file_complete, p.mc))
        error('Day %s is not (completely) motion-corrected; run caliali_run(animal_dir, ''mc'') first.', p.day);
    end
    det = [strrep(p.mc, '.mat', '_det.mat'), strrep(p.crop, '.mat', '_det.mat')];
    if isfile(p.aligned) || any(cellfun(@isfile, det))
        error('Day %s already has detrending/alignment files; cropping now would make them stale.', p.day);
    end
    locks{end+1} = caliali_lock(p); %#ok<AGROW>   released when this function ends
    mc = [mc, p.mc];  crop = [crop, p.crop];  todo{end+1} = p; %#ok<AGROW>
end
if isempty(mc)
    fprintf('Nothing to crop.\n');
    return
end

% crop files of an earlier run that did not finish (its days were never marked): start over
for k = 1:numel(crop)
    for f = {crop{k}, strrep(crop{k}, '.mat', '_tmp.mat')}
        if isfile(f{1})
            fprintf('Unfinished crop file from an earlier run, deleted: %s\n', f{1});
            delete(f{1});
        end
    end
end

% preview: middle frame of every session, centred in a common frame (as catpad_centered)
n = numel(mc);
sz = zeros(n, 3);
for k = 1:n
    sz(k, :) = size(matfile(mc{k}), 'Y');
end
D = max(sz(:, 1:2), [], 1);
off = floor((D - sz(:, 1:2)) / 2);
P = nan([D, n]);
fprintf('Reading the middle frame of %d sessions ...\n', n);
for k = 1:n
    m = matfile(mc{k});
    P(off(k, 1) + (1:sz(k, 1)), off(k, 2) + (1:sz(k, 2)), k) = double(m.Y(:, :, round(sz(k, 3) / 2)));
end
fig = figure('Name', sprintf('Crop: %d sessions from %d days', n, numel(todo)));
imagesc(mean(P, 3, 'omitnan')); axis image; colormap gray;
title('Draw the area to keep (generous), then double-click inside it');
impixelinfo;
roi = drawrectangle;
wait(roi);
pr = roi.Position;                                           % [xmin ymin width height]
close(fig);
R = max(1, ceil(pr(2))):min(D(1), floor(pr(2) + pr(4)));     % rows, common frame
C = max(1, ceil(pr(1))):min(D(2), floor(pr(1) + pr(3)));     % cols, common frame
assert(~isempty(R) && ~isempty(C), 'Empty crop rectangle.');
fprintf('Crop rectangle (common frame %d x %d): rows %d:%d, cols %d:%d\n', D(1), D(2), R(1), R(end), C(1), C(end));

% crop every session into its new file
boxes = zeros(n, 4);
for k = 1:n
    r = R - off(k, 1);  r = r(r >= 1 & r <= sz(k, 1));
    c = C - off(k, 2);  c = c(c >= 1 & c <= sz(k, 2));
    assert(~isempty(r) && ~isempty(c), 'The rectangle does not overlap %s', mc{k});
    boxes(k, :) = [r(1), r(end), c(1), c(end)];
    crop_to_new_file(mc{k}, crop{k}, r, c, sz(k, 3));
    fprintf('  %2d/%d  %s  rows %d:%d, cols %d:%d\n', k, n, crop{k}, r(1), r(end), c(1), c(end));
end

% all cropped and checked: now delete the _mc files and mark the days (with a record of the crop)
for k = 1:n
    delete(mc{k});
end
for d = 1:numel(todo)
    p = todo{d};
    fid = fopen(p.crop_done, 'w');
    fprintf(fid, 'cropped %s by %s\ncommon frame %d x %d, rectangle rows %d:%d, cols %d:%d\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM'), getenv('USER'), D(1), D(2), R(1), R(end), C(1), C(end));
    for k = find(ismember(crop, p.crop))
        fprintf(fid, '%s: rows %d:%d, cols %d:%d of %d x %d\n', crop{k}, boxes(k, :), sz(k, 1), sz(k, 2));
    end
    fclose(fid);
end
fprintf('Done: %d sessions cropped, %d days marked. Next: caliali_run(''%s'', ''align'')\n', ...
    n, numel(todo), animal_dir);
end

function crop_to_new_file(src, dst, rows, cols, n_frames)
% Cropped copy of src in dst (via dst_tmp.mat, renamed when complete); src is not modified.
tmp = strrep(dst, '.mat', '_tmp.mat');
if isfile(tmp), delete(tmp); end                       % never reuse a stale temporary file
s = load(src, 'CaliAli_options');
m = matfile(src);
mt = matfile(tmp, 'Writable', true);
mt.Y(numel(rows), numel(cols), n_frames) = cast(0, class(m.Y(1, 1, 1)));
for t1 = 1:2000:n_frames
    t2 = min(t1 + 1999, n_frames);
    mt.Y(:, :, t1:t2) = m.Y(rows(1):rows(end), cols(1):cols(end), t1:t2);
end
M0 = s.CaliAli_options.motion_correction.Mask;         % stored mask, consistent with the new Y
if isequal(size(M0), [size(m, 'Y', 1), size(m, 'Y', 2)])
    s.CaliAli_options.motion_correction.Mask = M0(rows, cols);
else
    s.CaliAli_options.motion_correction.Mask = true(numel(rows), numel(cols));
end
save(tmp, '-append', '-struct', 's');
szt = size(matfile(tmp), 'Y');
assert(isequal(szt, [numel(rows), numel(cols), n_frames]) && caliali_file_complete(tmp), ...
    'Cropped file %s is incomplete', tmp);
movefile(tmp, dst);
end
