function tf = caliali_file_complete(f)
% CALIALI_FILE_COMPLETE  True if a CaliAli movie file (_mc.mat, _det.mat, ...) exists AND was
% actually written to the end.
% CaliAli pre-allocates its output files (full size, all zeros) before processing and, on a
% rerun, skips any file that exists. A run that is stopped (crash, kill, Ctrl+C) can therefore
% leave an empty file that every later run takes as finished. Same criterion as CaliAli's own
% remove_corrupted_output: the last frame must not be entirely zero.
tf = false;
if ~isfile(f), return, end
try
    m = matfile(f);
    sz = size(m, 'Y');
    if numel(sz) < 3 || sz(3) < 1, return, end
    last = m.Y(:, :, sz(3));
    tf = any(last(:));
catch
    tf = false;
end
end
