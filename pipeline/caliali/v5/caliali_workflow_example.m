%% CaliAli workflow for one animal: example script
% Run the sections one by one (Ctrl+Enter in the editor). Sections marked UNATTENDED can be left
% running (e.g. over a weekend); sections marked INTERACTIVE need you at the screen for a few minutes.
% Everything is driven by the animal folder on the NAS; all heavy files are written locally in
% /mnt/data1/CaliAli_work/<animal>/<day>/, only results go back to <NAS day>/CaliAli_results/.

%% 0. Setup (once per MATLAB session)
addpath('/mnt/NAS/Miniscopes/CaliAli/v5');                       % folder with the caliali_*.m files
animal = '/mnt/NAS/Miniscopes/Reg_CA1/GDi/R005_M5_26';        % the animal: all its day folders
animal = '/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26';
caliali_status(animal)                                         % table: which steps each day has done
% Parameters (spatial_ds, temporal_ds, gSig, ...) come from CaliAli_my_parameters.m and must stay
% the same for all days of an animal, otherwise days cannot be registered against each other later.

%% 1. UNATTENDED: downsampling, stripe frames, motion correction (all days)
% Per day and session: downsamples the AVIs into one My_V4_Miniscope_con.mat (deleting the per-AVI
% pieces once verified), finds the stripe frames (MiniAn's stripe_frames_detected.csv if present,
% otherwise the same detector), replaces them with neighbouring good frames, saves the list for the
% NaN masking at the end, and motion-corrects (_con_mc.mat). Finished days are skipped; a failing
% day is logged and the run continues. The table at the end shows ok / skipped / FAILED per day.
r = caliali_run(animal, 'mc');

%% 2. INTERACTIVE: one crop for the whole animal
% Shows the middle frames of all sessions of all days on top of each other. Draw ONE generous
% rectangle (exclude only what is outside the lens; stay inside every session's valid area) and
% double-click it. It is applied to every session of every day not cropped yet.
caliali_crop_animal(animal)

%% 3. UNATTENDED: detrending + inter-session alignment (all days)
% Per day: detrends each session (_det.mat, also the input for a later cross-day registration),
% aligns the sessions (blood vessels + neurons, non-rigid) and writes the aligned movie
% (_Aligned.mat, in the last session's folder). Single-session days get the same files without
% inter-session alignment. The alignment QC is saved, not shown:
%   <NAS day>/CaliAli_results/CaliAli_alignment_<animal>_<day>.txt              BV score, metrics
%   <NAS day>/CaliAli_results/CaliAli_alignment_<animal>_<day>_scores.png       correlation per session
%   <NAS day>/CaliAli_results/CaliAli_alignment_<animal>_<day>_projections.png  sessions x stages
r = caliali_run(animal, 'align');
% Check the QC images of every day before continuing: vessels and neurons of all sessions should
% overlay without visible shifts in the last column of the projections image.

%% 4. INTERACTIVE: seed thresholds and mask, one day at a time
% For each day: a window with the correlation, PNR and correlation x PNR images. Type min_corr /
% min_pnr (keep gSig as is), click Update, draw a mask around the usable tissue if needed, then
% "Save & close". The settings go into that day's _Aligned.mat, where extraction reads them.
days = caliali_animal_days(animal);
caliali_set_seeds(days{6})
% caliali_set_seeds(days{2})
% ... one call per day; caliali_status(animal) shows which days still need seeds.

%% 5. UNATTENDED: extraction + summary (all days with seeds)
% Per day: CNMF-E on the aligned movie, result exported to <NAS day>/CaliAli_results/
% (CaliAli_result_*.mat), then the summary:
%   CaliAli_summary_<animal>_<day>.mat/.csv/.png  contours, neuron counts, activity per session
%   CaliAli_traces_<animal>_<day>.mat             C, C_raw, S per session, NaN at stripe frames,
%                                                 with the raw video frame of every trace frame
% Note: CaliAli clears the variables of the base workspace at each extraction.
r = caliali_run(animal, 'extract');

%% 6. Look at one day's results
day = days{1};
p = caliali_day_paths(day, false);      % false: only compute the paths (works on any computer)
load(fullfile(p.export, ['CaliAli_summary_' p.tag '.mat']), 'summary');
load(fullfile(p.export, ['CaliAli_traces_'  p.tag '.mat']), 'traces');
fprintf('%d neurons; active per session: %s\n', summary.n_neurons, mat2str(sum(summary.active, 1)));
imshow(imread(fullfile(p.export, ['CaliAli_summary_' p.tag '.png'])))       % overview figure

% Browse the neurons: footprint in the tissue, activity in every session, numbers per session.
% Arrow keys move between neurons; the menu at the bottom changes the order (e.g. by SNR).
caliali_view_neurons(p)

% CaliAli's neuron + vessel projection image with the extracted neuron contours; reads results only
caliali_plot_vessels(p)

% For your own analyses: traces(s).C / C_raw / S are neurons x frames (NaN = stripe frame), and
% traces(s).raw_frame gives the raw video frame of every column (e.g. to line up with behaviour).

%% 6b. Choosing the deconvolution threshold smin (events smaller than abs(smin) noise units are dropped)
% Compare values on one day (nothing is saved), look at the result in caliali_view_neurons, then
% make the chosen value the pipeline's: edit s.smin in caliali_settings.m and rebuild the summaries
% of the extracted days. New extractions use it automatically.
% caliali_redeconv(day, [-5 -4 -3 -2])
% caliali_redeconv(day, -3)                 % try one value on one day
% edit caliali_settings                     % set s.smin = -3
% caliali_run(animal, 'summary')

%% 7. Free disk space once a day is checked
% % Keeps the _det files (for cross-day registration), the aligned movie, stripe info and results;
% % deletes everything else that can be regenerated. Lists the files and asks before deleting.
% caliali_clean_day(day, 'intermediate')
% % To redo a day completely from the AVIs instead:  caliali_clean_day(day, 'all')
