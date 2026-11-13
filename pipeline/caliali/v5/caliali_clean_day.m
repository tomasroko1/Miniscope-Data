function caliali_clean_day(basedir, mode, where)
% CALIALI_CLEAN_DAY  List the CaliAli files of one day and delete them after confirmation.
%   caliali_clean_day(basedir)                   % 'all': everything, to restart from the AVIs
%   caliali_clean_day(basedir, 'intermediate')   % free disk: keep _det files, aligned file, results
%   caliali_clean_day(basedir, 'all', 'nas')     % old files written on the NAS before the local setup
% Only CaliAli's own file types are touched; raw AVIs, MiniAn outputs and other files never are.
%
% File types (per session folder <day>/<session>/):
%   My_V4_Miniscope/N_ds.mat            per-AVI downsampled pieces (also in other camera folders)
%   My_V4_Miniscope_con.mat             concatenated session              (downsampling)
%   My_V4_Miniscope_stripes.mat         stripe frames + replacement map   (kept by 'intermediate')
%   My_V4_Miniscope_con_mc.mat          motion corrected                  (motion correction)
%   My_V4_Miniscope_con_mc_crop.mat     cropped (replaces _con_mc.mat)    (caliali_crop_animal)
%   My_V4_Miniscope_con_mc_det.mat      detrended + projections           (alignment, part 1)
%   *_con*_tmp.mat / *_con*_bak.mat     leftovers of an interrupted crop
%   My_V4_Miniscope_con_mc_Aligned.mat  all sessions aligned (last session only)
%   ..._Aligned_source_extraction/      CNMF-E: data_*.mat (a full data copy) + LOGS_* (results)
% Day folder: CaliAli_options_<animal>_<day>.mat (snapshot), CaliAli_crop_done_<...>.txt (marker)
if nargin < 2, mode = 'all'; end
if nargin < 3, where = 'work'; end
mode  = validatestring(mode,  {'all', 'intermediate'});
where = validatestring(where, {'work', 'nas'});
p = caliali_day_paths(basedir);
if strcmp(where, 'nas')                         % legacy layout: everything next to the NAS AVIs
    p.sessions  = p.nas_sessions;
    p.options   = fullfile(p.nas, ['CaliAli_options_' p.tag '.mat']);
    p.crop_done = fullfile(p.nas, ['CaliAli_crop_done_' p.tag '.txt']);
end
if isfile(p.lock)
    error('Day %s is locked (%s); not cleaning while someone works on it.', p.day, strtrim(fileread(p.lock)));
end
keep_all = strcmp(mode, 'all');

files = {};  folders = {};
sess_dirs = unique(cellfun(@fileparts, p.sessions, 'UniformOutput', false));
for si = 1:numel(sess_dirs)
    sdir = sess_dirs{si};
    files = [files, list(fullfile(sdir, '*', '*_ds.mat'))];                 %#ok<AGROW> per-AVI pieces
    con = list(fullfile(sdir, '*_con*.mat'));
    if ~keep_all                                                            % intermediate: keep these
        con = con(~endsWith(con, {'_det.mat', '_Aligned.mat'}));
    end
    files = [files, con];                                                   %#ok<AGROW>
    if keep_all
        files = [files, list(fullfile(sdir, '*_stripes.mat'))];             %#ok<AGROW> stripe info
    end
    ext = dir(fullfile(sdir, '*_source_extraction'));
    ext = ext([ext.isdir]);
    for e = 1:numel(ext)
        edir = fullfile(ext(e).folder, ext(e).name);
        if keep_all
            folders{end+1} = edir;                                          %#ok<AGROW> incl. results
        else
            files = [files, list(fullfile(edir, 'data_*.mat'))];            %#ok<AGROW> only the data copy
        end
    end
end
if keep_all
    files = [files, list(p.options), list(p.crop_done), list(p.seeds_done)];
end
files = files(~endsWith(lower(files), '.avi'));                             % paranoia: never raw videos

if isempty(files) && isempty(folders)
    fprintf('Nothing to clean for %s (%s).\n', p.day, mode);
    return
end
fprintf('\nCaliAli files for %s / %s  (mode: %s)\n', p.animal, p.day, mode);
total = 0;
for k = 1:numel(files)
    b = dir(files{k}).bytes;  total = total + b;
    fprintf('  %9.1f MB  %s\n', b / 1e6, files{k});
end
for k = 1:numel(folders)
    d = dir(fullfile(folders{k}, '**', '*'));
    b = sum([d(~[d.isdir]).bytes]);  total = total + b;
    fprintf('  %9.1f MB  %s%s  (folder)\n', b / 1e6, folders{k}, filesep);
end
fprintf('Total: %d files, %d folders, %.2f GB\n', numel(files), numel(folders), total / 1e9);

if ~strcmpi(strtrim(input('Delete all of the above? [y/N] ', 's')), 'y')
    fprintf('Nothing deleted.\n');
    return
end
for k = 1:numel(files),   delete(files{k});        end
for k = 1:numel(folders), rmdir(folders{k}, 's');  end
fprintf('Deleted %d files and %d folders (%.2f GB).\n', numel(files), numel(folders), total / 1e9);
end

function f = list(pattern)
d = dir(pattern);
d = d(~[d.isdir]);
f = cell(1, numel(d));
for k = 1:numel(d)
    f{k} = fullfile(d(k).folder, d(k).name);
end
end
