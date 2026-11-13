# CaliAli day pipeline

## Install (once)
1. CaliAli at `/mnt/data1/CaliAli/CaliAli-main`, with the patches listed at the end of this file.
2. This folder (all `caliali_*.m` files) on the MATLAB path, e.g. `/mnt/NAS/Miniscopes/CaliAli/`,
   added in `startup.m`. Your `CaliAli_my_parameters.m` too.
3. Work root: `mkdir -p /mnt/data1/CaliAli_work && chmod 2775 /mnt/data1/CaliAli_work`

## Files
| File | Called by you? | Role |
|---|---|---|
| caliali_run.m | yes | runs an unattended step ('mc', 'align', 'extract') on all days of an animal |
| caliali_crop_animal.m | yes | one crop for all days of an animal (interactive); writes _mc_crop.mat, then deletes _mc.mat |
| caliali_set_seeds.m | yes | thresholds + seed mask for one day (interactive) |
| caliali_status.m | yes | which steps are done for every day |
| caliali_check_mc.m | yes | verify stripe replacement and motion correction of every session |
| caliali_view_neurons.m | yes | browse neurons: footprint in the tissue, activity in every session, arrow keys |
| caliali_plot_vessels.m | yes | CaliAli's neuron + vessel projection image with the extracted neuron contours |
| caliali_clean_day.m | yes | list and delete a day's CaliAli files |
| caliali_summary.m | (also by caliali_run) | contours, activity per session, masked traces |
| caliali_redeconv.m | yes | redo only the deconvolution (C, S) with another smin; scan mode compares values |
| caliali_mc.m | by caliali_run | downsampling, stripe frames, motion correction (one day) |
| caliali_stripes_day.m | by caliali_mc | find / replace / remember stripe frames |
| caliali_align.m | by caliali_run | detrending + alignment (one day; single-session days included), QC images |
| caliali_extract.m | by caliali_run | CNMF-E (one day), exports the result |
| caliali_day_paths.m | helper | all paths of a day + local mirror with links to the AVIs |
| caliali_day_status.m | helper | done/not done per step, from the files |
| caliali_animal_days.m | helper | day folders of an animal |
| caliali_lock.m | helper | one person per day at a time |
| caliali_settings.m | edit | pipeline settings that are not CaliAli parameters: smin, min_events, contour_level |
| caliali_result_file.m | helper | the final CNMF-E result file of a day |
| caliali_file_complete.m | helper | is a CaliAli movie file written to the end (not an empty pre-allocated shell)? |
| caliali_setup.m | helper | adds CaliAli to the path |
| caliali_workflow_example.m | example | the whole workflow, commented |

## Workflow per animal
1. `caliali_run(animal, 'mc')`
2. `caliali_crop_animal(animal)`
3. `caliali_run(animal, 'align')`
4. `caliali_set_seeds(day)` for each day
5. `caliali_run(animal, 'extract')`

## Patches applied in place to CaliAli (CaliAli-main)
Each patch block in the CaliAli code is marked `PATCH Pn (caliali pipeline)`, with the original
lines kept as comments where something was replaced.
| Patch | File | Why |
|---|---|---|
| P1 | Motion_Correction/Vessel_MC/CaliAli_motion_correction.m | rerun with all batches done failed ("Unrecognized ... 'out'") |
| P2 | Motion_Correction/apply_crop_on_disk.m | already-cropped files failed on rerun; stored mask kept consistent |
| P3 | Other_codes/Utilities/CaliAli_crop.m | files from options; drawrectangle figure instead of the crop app; mask kept consistent |
| P4 | Align_sessions/detrend_batch_and_calculate_projections.m | an error deleted every finished _det file; now only broken ones |
| P5 | Other_codes/Utilities/remove_corrupted_output.m | return -> continue: check every file |
| P6 | Other_codes/Utilities/plot_darkmode.m | find the tiledlayout among the figure's children |
| P7 | CNMF-e/ca_source_extraction/endoscope/get_patch_data.m | Source() -> Source (warning) |
| P8 | Other_codes/Utilities/plot_P.m | clf(fig) before each frame |
| P9 | Align_sessions/CaliAli_align_sessions.m | single-session days: prepared for extraction without inter-session alignment |
Earlier single-session attempts also edited `Other_codes/catpad_centered.m` (line 77) and
`Align_sessions/match_video_size.m` (if numel(Mask) > 1). With P9 they are no longer reached for
single-session days; they are harmless for multi-session days and can stay.

Note: the animal crop no longer uses CaliAli_crop (patch P3 is unused since this version); it is
done by caliali_crop_animal itself and writes a new file per session instead of changing _mc.mat.
