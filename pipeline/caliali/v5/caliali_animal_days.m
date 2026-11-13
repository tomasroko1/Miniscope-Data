function days = caliali_animal_days(animal_dir)
% CALIALI_ANIMAL_DAYS  Day folders of one animal (folders containing */My_V4_Miniscope), sorted.
animal_dir = char(regexprep(animal_dir, '[\\/]+$', ''));
d = dir(animal_dir);
d = d([d.isdir] & ~startsWith({d.name}, {'.', '@', '#'}));
days = {};
for k = 1:numel(d)
    day = fullfile(animal_dir, d(k).name);
    if ~isempty(dir(fullfile(day, '*', 'My_V4_Miniscope')))
        days{end+1} = day; %#ok<AGROW>
    end
end
days = sort(days);
end
