function caliali_mc(basedir)
% CALIALI_MC  Unattended step 1 for one day: downsampling, stripe frames, motion correction.
%   caliali_mc('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26/T1_SD_CNO_2026_06_19')
% Normally called for all days by caliali_run(animal_dir, 'mc'). Rerunning skips finished work.
caliali_setup();
p = caliali_day_paths(basedir);
lock = caliali_lock(p); %#ok<NASGU>
if all(cellfun(@caliali_file_complete, p.crop))
    fprintf('Day %s is already motion-corrected and cropped.\n', p.day);
    return
end
CaliAli_options = CaliAli_parameters(CaliAli_my_parameters());
fprintf('Day %s / %s: %d sessions, working in %s\n', p.animal, p.day, numel(p.sessions), p.basedir);

% downsampling -> one _con.mat per session (skipped once verified and cleaned up)
CaliAli_options.downsampling.input_files = p.sessions;
if all(cellfun(@isfile, p.con)) && ~any_ds_pieces(p.sessions)
    fprintf('Downsampling already done.\n');
    CaliAli_options.downsampling.output_files = p.con;
else
    CaliAli_options.downsampling.output_files = {};         % CaliAli appends, so start empty
    CaliAli_options = CaliAli_downsample_batch(CaliAli_options);
    delete_ds_pieces(p.sessions, p.con);
end

% stripe frames: find (or reuse MiniAn's list), replace in each _con.mat, remember for masking
caliali_stripes_day(p, CaliAli_options.downsampling.temporal_ds);

% motion correction (CaliAli skips finished files; includes its automatic border crop).
% CaliAli takes any existing _mc.mat as finished, also an empty one pre-allocated by a run that
% was stopped, or one whose automatic border crop was interrupted (it leaves _mc_tmp.mat or
% _mc_bak.mat behind): delete those sessions' files first so they are redone.
for k = 1:numel(p.mc)
    tmp = strrep(p.mc{k}, '.mat', '_tmp.mat');  bak = strrep(p.mc{k}, '.mat', '_bak.mat');
    if isfile(tmp) || isfile(bak)
        fprintf('Interrupted motion correction (border crop) in %s: redone\n', fileparts(p.mc{k}));
        for f = {p.mc{k}, tmp, bak}
            if isfile(f{1}), delete(f{1}); end
        end
    end
end
delete_unfinished(p.mc);
CaliAli_options.motion_correction.input_files = CaliAli_options.downsampling.output_files;
CaliAli_options = CaliAli_motion_correction(CaliAli_options); %#ok<NASGU>
fprintf('Motion correction done: %s\n', p.day);
end

function parts = ds_pieces(s)
avi   = {dir(fullfile(s, '*.avi')).name};
parts = fullfile(s, strrep(avi, '.avi', '_ds.mat'));         % only pieces made from these AVIs
parts = parts(cellfun(@isfile, parts));
end

function tf = any_ds_pieces(sessions)
tf = any(cellfun(@(s) ~isempty(ds_pieces(s)), sessions));
end

function delete_ds_pieces(sessions, con)
nfr = @(f) get_data_dimension(f) * [0; 0; 1];
for k = 1:numel(sessions)
    parts = ds_pieces(sessions{k});
    n_avi = numel(dir(fullfile(sessions{k}, '*.avi')));
    if isempty(parts), continue; end
    if isfile(con{k}) && numel(parts) == n_avi && nfr(con{k}) == sum(cellfun(nfr, parts))
        cellfun(@delete, parts);
        fprintf('Deleted %d _ds pieces in %s\n', numel(parts), sessions{k});
    else
        warning('Kept _ds pieces in %s (missing or incomplete %s)', sessions{k}, con{k});
    end
end
end

function delete_unfinished(files)
for k = 1:numel(files)
    if isfile(files{k}) && ~caliali_file_complete(files{k})
        fprintf('Unfinished file from an interrupted run, deleted to be redone: %s\n', files{k});
        delete(files{k});
    end
end
end
