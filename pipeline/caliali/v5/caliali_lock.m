function unlock = caliali_lock(p)
% CALIALI_LOCK  Stop two people (or machines) from working on the same day at once.
%   lock = caliali_lock(p);   % released automatically when the calling function ends or fails
if isfile(p.lock)
    error(['Day %s is already being processed: %s\n' ...
           'If that run is dead (e.g. MATLAB crashed), delete:\n  %s'], p.day, strtrim(fileread(p.lock)), p.lock);
end
[~, host] = system('hostname');
fid = fopen(p.lock, 'w');
fprintf(fid, '%s on %s since %s', getenv('USER'), strtrim(host), datestr(now, 'yyyy-mm-dd HH:MM'));
fclose(fid);
unlock = onCleanup(@() delete(p.lock));
end
