function st = caliali_day_status(p)
% CALIALI_DAY_STATUS  Which steps are done for one day (p from caliali_day_paths), from its files.
% motion-corrected: _mc.mat written to the end without leftovers of an interrupted border crop
% (_mc_tmp/_mc_bak), or already cropped into a complete _mc_crop.mat
leftover = cellfun(@isfile, strrep(p.mc, '.mat', '_tmp.mat')) | cellfun(@isfile, strrep(p.mc, '.mat', '_bak.mat'));
mc_ok    = cellfun(@caliali_file_complete, p.mc) & ~leftover;
st.mc    = all(mc_ok | cellfun(@caliali_file_complete, p.crop));
st.cropped   = isfile(p.crop_done);
st.aligned   = false;
if isfile(p.aligned)
    try, st.aligned = logical(CaliAli_load(p.aligned, 'alignment_completed')); catch, end
end
st.seeds     = isfile(p.seeds_done);
st.extracted = isfile(p.result);
if st.extracted && isfile(p.seeds_done)          % seeds set again after the extraction: redo it
    st.extracted = dir(p.result).datenum >= dir(p.seeds_done).datenum;
end
st.summary   = isfile(fullfile(p.export, ['CaliAli_summary_' p.tag '.png']));
st.locked    = isfile(p.lock);
end
