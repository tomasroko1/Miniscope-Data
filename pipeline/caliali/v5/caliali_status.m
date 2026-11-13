function T = caliali_status(animal_dir)
% CALIALI_STATUS  Table of which steps are done for every day of one animal.
%   caliali_status('/mnt/NAS/Miniscopes/Reg_CA1/GDi/R004_M4_26')
caliali_setup();
days = caliali_animal_days(animal_dir);
rows = cell(numel(days), 1);
for k = 1:numel(days)
    p  = caliali_day_paths(days{k});
    st = caliali_day_status(p);
    rows{k} = {p.day, st.mc, st.cropped, st.aligned, st.seeds, st.extracted, st.summary, st.locked};
end
T = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'day', 'mc', 'cropped', 'aligned', 'seeds', 'extracted', 'summary', 'locked'});
disp(T);
end
