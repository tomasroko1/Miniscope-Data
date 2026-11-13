function caliali_stripes_day(p, temporal_ds)
% CALIALI_STRIPES_DAY  Stripe frames for every session of one day: find, replace, remember.
% Called by caliali_mc after downsampling and BEFORE motion correction. Same logic as
% run_minian_batch_stripes_v4.py:
%   1. stripe frames (raw frame numbers, 0-based global over the session's AVIs, as in Python):
%      read from MiniAn's  <NAS session>/My_V4_Miniscope/stripe_frames_detected.csv  if present
%      (so both pipelines use the identical list), otherwise detected here with the same
%      first-column periodogram detector (threshold 20, peaks above 0.03 cycles/pixel) and saved
%      in the same CSV format in the local work folder; lists covering > 50% of frames are rejected
%   2. mapped to frames of the session's _con.mat (CaliAli keeps every temporal_ds-th frame of
%      each AVI, so a stripe frame is either dropped or kept unchanged)
%   3. each block of stripe frames in _con.mat is replaced by the same number of good frames,
%      previous block preferred, otherwise following (as in Python); the movie keeps its length
%   4. everything is saved in <session>/My_V4_Miniscope_stripes.mat, used at the end to set the
%      activity of all neurons to NaN at those frames (caliali_summary)
% Safe to rerun: a session whose _con.mat already has its stripes replaced is skipped.
THRESHOLD = 20;          % MINIAN_STRIPE_DETECTION_THRESHOLD
MAX_FRACTION = 0.5;      % MINIAN_MAX_AUTO_STRIPE_FRACTION
CSV_NAME = 'stripe_frames_detected.csv';

for k = 1:numel(p.sessions)
    con = p.con{k};
    info_file = p.stripes{k};
    if is_current(info_file, con)
        fprintf('Stripe frames already replaced: %s\n', con);
        continue
    end
    if isfile(strrep(con, '.mat', '_mc.mat')) || isfile(strrep(con, '.mat', '_mc_crop.mat'))
        error(['%s was motion-corrected before stripe replacement. Clean this day with ' ...
               'caliali_clean_day(day, ''all'') and run it again.'], con);
    end

    % raw AVIs in the same order as MiniAn (numeric names, natural order) and their frame counts
    avi = {dir(fullfile(p.sessions{k}, '*.avi')).name};
    avi = natsortfiles(avi(~cellfun(@isempty, regexp(avi, '^[0-9]+\.avi$'))));
    n_raw = zeros(1, numel(avi));
    for a = 1:numel(avi)
        n_raw(a) = VideoReader(fullfile(p.sessions{k}, avi{a})).NumFrames;
    end
    n_total = sum(n_raw);

    % 1. stripe frame list (raw, 0-based global)
    nas_csv = fullfile(p.nas_sessions{k}, CSV_NAME);
    if isfile(nas_csv)
        raw_bad = read_stripe_csv(nas_csv);
        source = sprintf('MiniAn CSV %s', nas_csv);
    else
        [raw_bad, scores] = detect_session(p.sessions{k}, avi, THRESHOLD);
        write_stripe_csv(fullfile(p.sessions{k}, CSV_NAME), raw_bad, scores, n_raw, avi, THRESHOLD);
        source = sprintf('detected by caliali_stripes_day (threshold %g)', THRESHOLD);
    end
    raw_bad = unique(raw_bad(raw_bad >= 0 & raw_bad < n_total));
    rejected = numel(raw_bad) > floor(MAX_FRACTION * n_total);
    if rejected
        warning('%s: %d of %d frames listed as stripe frames (> %g%%); list rejected, no replacement.', ...
            p.sessions{k}, numel(raw_bad), n_total, 100 * MAX_FRACTION);
        raw_bad = [];
    end

    % 2. raw frame of every _con.mat frame, and the stripe frames among them
    starts = [0 cumsum(n_raw(1:end-1))];
    trace_raw_frame = cell2mat(arrayfun(@(a) starts(a) + (0:temporal_ds:n_raw(a)-1), ...
                                        1:numel(avi), 'UniformOutput', false));
    n_con = get_data_dimension(con) * [0; 0; 1];
    assert(numel(trace_raw_frame) == n_con, ...
        '%s: %d frames expected from the AVIs (temporal_ds %d) but %d found; cannot map stripe frames.', ...
        con, numel(trace_raw_frame), temporal_ds, n_con);
    con_bad = find(ismember(trace_raw_frame, raw_bad));                % 1-based frames of _con.mat

    % 3. replace in place
    map = replace_blocks(con, con_bad);

    % 4. remember
    d = dir(con);
    info = struct('session', p.sessions{k}, 'source', source, 'rejected', rejected, ...
        'raw_stripe_frames', raw_bad(:)', 'raw_frames_dropped_by_ds', setdiff(raw_bad, trace_raw_frame), ...
        'trace_stripe_frames', con_bad(:)', 'trace_raw_frame', trace_raw_frame, ...
        'replacement', map, 'replacement_columns', {{'trace_frame', 'source_trace_frame', 'method'}}, ...
        'replacement_methods', {{'previous', 'following', 'previous_nonlocal', 'following_nonlocal', ...
                                 'previous_sparse', 'following_sparse', 'padded'}}, ...
        'avi_files', {avi}, 'avi_frames', n_raw, 'temporal_ds', temporal_ds, ...
        'con_datenum', d.datenum, 'con_bytes', d.bytes);
    save(info_file, '-struct', 'info');
    fprintf('%s: %d stripe frames (%s); %d kept by downsampling and replaced, %d dropped by it\n', ...
        p.sessions{k}, numel(raw_bad), source, numel(con_bad), numel(info.raw_frames_dropped_by_ds));
end
end

% ------------------------------------------------------------------------------------------------
function tf = is_current(info_file, con)
tf = false;
if isfile(info_file) && isfile(con)
    s = load(info_file, 'con_datenum', 'con_bytes');
    d = dir(con);
    tf = abs(s.con_datenum - d.datenum) < 1e-9 && s.con_bytes == d.bytes;   % this very _con.mat
end
end

function frames = read_stripe_csv(csv_file)
% Same rules as read_stripe_frames_csv in Python: column global_frame (or frame / bad_frame);
% rows explicitly marked false in mask/use/keep/include are ignored.
T = readtable(csv_file, 'TextType', 'string');
names = T.Properties.VariableNames;
col = intersect({'global_frame', 'frame', 'bad_frame'}, names, 'stable');
if isempty(col) || height(T) == 0
    frames = [];
    return
end
keep = true(height(T), 1);
for flag = {'mask', 'use', 'keep', 'include'}
    if ismember(flag{1}, names)
        v = lower(strtrim(string(T.(flag{1}))));
        keep = keep & ~(v ~= "" & ~ismember(v, ["1", "true", "yes", "y"]));
    end
end
frames = double(T.(col{1})(keep));
frames = unique(round(frames(isfinite(frames))))';
end

function [bad, scores] = detect_session(folder, avi, threshold)
% First-column periodogram detector, as stripe_peak_score_column in Python.
scores = [];
for a = 1:numel(avi)
    v = VideoReader(fullfile(folder, avi{a}));
    n = v.NumFrames;
    s = zeros(1, n);
    for t1 = 1:500:n
        t2 = min(n, t1 + 499);
        fr = read(v, [t1 t2]);                                % H x W x C x frames
        col = double(squeeze(fr(:, 1, 1, :)));                % first column, first channel
        for j = 1:size(col, 2)
            s(t1 + j - 1) = stripe_score(col(:, j), threshold);
        end
    end
    scores = [scores, s]; %#ok<AGROW>
end
bad = find(scores >= threshold) - 1;                          % 0-based global, as in Python
end

function score = stripe_score(x, threshold)
% scipy.signal.periodogram(x, 1) (mean removed, one-sided density), sqrt, first fifth of the
% bins, local maxima >= threshold above 0.03 cycles/pixel; score = largest such peak, else 0.
score = 0;
N = numel(x);
if N < 4 || ~any(isfinite(x)), return, end
x = x - mean(x);
X = fft(x);
P = abs(X(1:floor(N/2) + 1)).^2 / N;
if mod(N, 2) == 0, P(2:end-1) = 2 * P(2:end-1); else, P(2:end) = 2 * P(2:end); end
L = floor(numel(P) / 5);
if L < 3, return, end
spec = sqrt(P(1:L));
freq = (0:L-1)' / N;
i = (2:L-1)';
pk = i(spec(i) > spec(i-1) & spec(i) > spec(i+1));
pk = pk(spec(pk) >= threshold & freq(pk) > 0.03);
if ~isempty(pk), score = max(spec(pk)); end
end

function write_stripe_csv(csv_file, bad, scores, n_raw, avi, threshold)
% Same columns as save_detected_stripe_frames_csv in Python.
starts = [0 cumsum(n_raw(1:end-1))];
fid = fopen(csv_file, 'w');
fprintf(fid, 'detection_label,global_frame,video_index,video_file,frame_in_video,stripe_peak_score,stripe_detection_threshold\n');
for g = bad(:)'
    a = find(starts <= g, 1, 'last');
    fprintf(fid, 'auto_periodogram_first_column,%d,%d,%s,%d,%.6g,%g\n', ...
        g, a - 1, avi{a}, g - starts(a), scores(g + 1), threshold);
end
fclose(fid);
end

function map = replace_blocks(con, bad)
% replace_stripe_frames_with_neighbor_chunks from Python, on the frames of one _con.mat
map = zeros(0, 3);                                            % [frame, source_frame, method_id]
if isempty(bad), return, end
m = matfile(con, 'Writable', true);
T = size(m, 'Y', 3);
good = true(1, T);  good(bad) = false;
if ~any(good)
    warning('%s: every frame is a stripe frame; nothing replaced.', con);
    return
end
d = diff(bad(:)');  e = [find(d > 1), numel(bad)];  s = [1, find(d > 1) + 1];
for b = 1:numel(s)
    first = bad(s(b));  last = bad(e(b));  len = last - first + 1;
    [src, method] = good_chunk(good, first, last, len);
    for i = 1:len
        m.Y(:, :, first + i - 1) = m.Y(:, :, src(i));
    end
    map = [map; (first:last)', src(:), repmat(method, len, 1)]; %#ok<AGROW>
end
end

function [src, method] = good_chunk(good, first, last, len)
% find_contiguous_good_chunk from Python (1-based). method: 1 previous, 2 following,
% 3 previous_nonlocal, 4 following_nonlocal, 5 previous_sparse, 6 following_sparse, 7 padded
T = numel(good);
if first - len >= 1 && all(good(first-len:first-1)), src = first-len:first-1; method = 1; return, end
if last + len <= T && all(good(last+1:last+len)),    src = last+1:last+len;   method = 2; return, end
for c = first - len - 1:-1:1
    if all(good(c:c+len-1)), src = c:c+len-1; method = 3; return, end
end
for c = last + 2:T - len + 1
    if all(good(c:c+len-1)), src = c:c+len-1; method = 4; return, end
end
gp = find(good);  prev = gp(gp < first);  next = gp(gp > last);
if numel(prev) >= len, src = prev(end-len+1:end); method = 5; return, end
if numel(next) >= len, src = next(1:len);         method = 6; return, end
if ~isempty(prev), base = prev; else, base = next; end
src = [base, repmat(base(end), 1, len - numel(base))];  method = 7;
end
