function p = caliali_day_paths(basedir, make_mirror)
% CALIALI_DAY_PATHS  Every path CaliAli needs for one day, derived from the NAS day folder.
%   p = caliali_day_paths('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19')
% All processing happens in a LOCAL work folder that mirrors the day on the NAS:
%   <work root>/<animal>/<day>/<session>/My_V4_Miniscope/N.avi  -> symbolic links to the NAS AVIs
% CaliAli writes each output next to its input, so every intermediate lands locally, while the
% raw videos are only read from the NAS (never copied). Final results are exported to
%   <NAS day>/CaliAli_results/
% Work root: environment variable CALIALI_WORK, default /mnt/data1/CaliAli_work.
%   p = caliali_day_paths(basedir, false)   only computes the paths and creates nothing (e.g. for
%                                           reading results on a computer without the work folder)

if nargin < 2 || isempty(make_mirror), make_mirror = true; end
nas = char(regexprep(basedir, '[\\/]+$', ''));
assert(isfolder(nas), 'Day folder not found: %s', nas);
[animal_dir, day] = fileparts(nas);
[~, animal] = fileparts(animal_dir);
work_root = getenv('CALIALI_WORK');
if isempty(work_root), work_root = '/mnt/data1/CaliAli_work'; end
work = fullfile(work_root, animal, day);

p.nas    = nas;
p.animal = animal;
p.day    = day;
p.tag    = regexprep([animal '_' day], '\W', '_');
p.basedir = work;                                   % where everything is processed
p.export  = fullfile(nas, 'CaliAli_results');       % where final results are copied

% sessions on the NAS (sorted = chronological) and their local mirrors
p.nas_sessions = unique({dir(fullfile(nas, '*', 'My_V4_Miniscope')).folder});
assert(~isempty(p.nas_sessions), 'No */My_V4_Miniscope folders in %s', nas);
p.sessions = cell(size(p.nas_sessions));
for k = 1:numel(p.nas_sessions)
    [sess_dir, cam] = fileparts(p.nas_sessions{k});
    [~, sess] = fileparts(sess_dir);
    p.sessions{k} = fullfile(work, sess, cam);
    if make_mirror
        mirror_avis(p.nas_sessions{k}, p.sessions{k});
    end
end

% our own bookkeeping, one set per day
p.options   = fullfile(work, ['CaliAli_options_' p.tag '.mat']);
p.crop_done = fullfile(work, ['CaliAli_crop_done_' p.tag '.txt']);

% one file name per processing stage (inside the local work folder), per session:
%   _con.mat      downsampled, stripe frames replaced        (caliali_mc)
%   _con_mc.mat   motion-corrected                           (caliali_mc; deleted once cropped)
%   _con_mc_crop.mat  cropped                                (caliali_crop_animal)
%   <input>_det.mat, <last input>_Aligned.mat                (caliali_align, CaliAli's naming)
p.con     = cellfun(@(s) fullfile(fileparts(s), 'My_V4_Miniscope_con.mat'), p.sessions, 'UniformOutput', false);
p.stripes = cellfun(@(c) strrep(c, '_con.mat', '_stripes.mat'), p.con, 'UniformOutput', false);
p.mc      = cellfun(@(c) strrep(c, '_con.mat', '_con_mc.mat'), p.con, 'UniformOutput', false);
p.crop    = cellfun(@(c) strrep(c, '_con.mat', '_con_mc_crop.mat'), p.con, 'UniformOutput', false);
% alignment input: the cropped files; for days cropped in place by the earlier pipeline (or
% marked as cropped without cropping) the _mc files, so their existing results keep their names
p.cropped_in_place = isfile(p.crop_done) && ~all(cellfun(@isfile, p.crop));
if p.cropped_in_place, p.align_input = p.mc; else, p.align_input = p.crop; end
p.aligned    = [p.align_input{end}(1:end-4) '_Aligned.mat'];   % CaliAli: next to the last input
p.extraction = [p.aligned(1:end-4) '_source_extraction'];      % CNMF-E
p.seeds_done = fullfile(work, ['CaliAli_seeds_done_' p.tag '.txt']);
p.result    = fullfile(p.export, ['CaliAli_result_' p.tag '.mat']);
p.lock      = fullfile(work, ['CaliAli_RUNNING_' p.tag '.lock']);
end

function mirror_avis(src, dst)
% Local folder with symbolic links to the raw AVIs (created once, safe to repeat).
if ~isfolder(dst)
    mkdir(dst);
    fileattrib(dst, '+w', 'g');                     % other users on the server can work here too
end
avi = dir(fullfile(src, '*.avi'));
for k = 1:numel(avi)
    target = fullfile(src, avi(k).name);
    link   = fullfile(dst, avi(k).name);
    if ~isfile(link)
        [status, msg] = system(sprintf('ln -s "%s" "%s"', target, link));
        assert(status == 0, 'Could not create link %s: %s', link, msg);
    end
end
end
