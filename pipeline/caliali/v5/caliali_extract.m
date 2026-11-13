function File_path = caliali_extract(basedir)
% CALIALI_EXTRACT  Step 3 of 3 (server): CNMF-E on the day's aligned file.
%   File_path = caliali_extract('/mnt/NAS/.../R004_M4_26/T1_SD_CNO_2026_06_19')
caliali_setup();
p = caliali_day_paths(basedir);
assert(isfile(p.aligned), 'Not aligned yet; run caliali_run(animal_dir, ''align'') first.');
lock = caliali_lock(p); %#ok<NASGU>
opt = CaliAli_load(p.aligned, 'CaliAli_options');
if isempty(opt.cnmf.seed_mask), m = 'none (whole field)'; else, m = sprintf('%d pixels', nnz(opt.cnmf.seed_mask)); end
fprintf('Using thresholds stored in the file: min_corr = %g, min_pnr = %g, seed mask: %s\n', ...
    opt.cnmf.min_corr, opt.cnmf.min_pnr, m);
% CNMF-E copies the video into <extraction folder>/data_*.mat and REUSES an existing copy
% without checking it. A copy left by an earlier (possibly interrupted) run would be used
% instead of the current aligned file, so always remove it; it is rebuilt from p.aligned.
old_data = dir(fullfile(p.extraction, 'data_*.mat'));
for k = 1:numel(old_data)
    fprintf('Removing old CNMF-E data copy: %s\n', fullfile(old_data(k).folder, old_data(k).name));
    delete(fullfile(old_data(k).folder, old_data(k).name));
end
File_path = CaliAli_cnmfe(p.aligned);
CaliAli_options = CaliAli_load(p.aligned, 'CaliAli_options');
save(p.options, 'CaliAli_options');

% export the final result and the options snapshot to the NAS day folder
result = File_path; if iscell(result), result = result{1}; end
if ~isfolder(p.export), mkdir(p.export); fileattrib(p.export, '+w', 'g'); end
copyfile(result,    fullfile(p.export, ['CaliAli_result_' p.tag '.mat']));
copyfile(p.options, fullfile(p.export, ['CaliAli_options_' p.tag '.mat']));
fprintf('Results exported to %s\n', p.export);
end
