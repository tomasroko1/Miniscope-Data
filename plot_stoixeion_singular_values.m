function plot_stoixeion_singular_values(resultsDir)
% Recreate singular_values.png from an existing Stoixeion singular_values.csv.

if nargin < 1 || isempty(resultsDir)
    resultsDir = fullfile(pwd, 'results', 'stoixeion', 'phasewise');
end
csvPath = fullfile(resultsDir, 'singular_values.csv');
if exist(csvPath, 'file') ~= 2
    error('No existe %s', csvPath);
end
T = readtable(csvPath);
required = {'animal','day_name','session','phase','singular_rank', ...
    'singular_value','selected_by_stoixeion'};
if ~all(ismember(required, T.Properties.VariableNames))
    error('singular_values.csv no tiene las columnas esperadas.');
end

keys = cell(height(T), 1);
for row = 1:height(T)
    keys{row} = sprintf('%s|%s|%s|%s', textAt(T.animal, row), ...
        textAt(T.day_name, row), textAt(T.session, row), textAt(T.phase, row));
end
[~, firstRows] = unique(keys, 'stable');
for group = firstRows(:)'
    mask = strcmp(keys, keys{group});
    ranks = double(T.singular_rank(mask));
    values = double(T.singular_value(mask));
    selected = selectedFlags(T.selected_by_stoixeion(mask));
    [ranks, order] = sort(ranks);
    values = values(order);
    selected = selected(order);
    if isempty(ranks), continue; end

    animal = textAt(T.animal, group);
    dayName = textAt(T.day_name, group);
    session = textAt(T.session, group);
    phase = textAt(T.phase, group);
    [~, sessionStem] = fileparts(session);
    phaseDir = fullfile(resultsDir, animal, dayName, sessionStem, phase);
    if exist(phaseDir, 'dir') ~= 7, mkdir(phaseDir); end

    fig = figure('Visible', 'off', 'Color', 'w');
    semilogx(ranks, values, 'k-', 'LineWidth', 1.5);
    hold on;
    if any(selected)
        semilogx(ranks(selected), values(selected), 'ro', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 6);
        legend({'Valores singulares','Rangos seleccionados'}, 'Location', 'northeast');
    end
    xlim([1 max(10, max(ranks))]);
    ylim([0 max(1, max(values) * 1.05)]);
    xlabel('Rango singular (escala log)');
    ylabel('Valor singular');
    title(sprintf('%s | %s | SVD de similitud', animal, phase));
    box on;
    saveas(fig, fullfile(phaseDir, 'singular_values.png'));
    close(fig);
    fprintf('Figura: %s\n', fullfile(phaseDir, 'singular_values.png'));
end
end

function value = textAt(column, row)
if iscell(column)
    value = column{row};
elseif ischar(column)
    value = deblank(column(row, :));
else
    value = char(column(row));
end
end

function flags = selectedFlags(column)
if isnumeric(column) || islogical(column)
    flags = column ~= 0;
else
    flags = false(numel(column), 1);
    for row = 1:numel(column)
        value = textAt(column, row);
        flags(row) = strcmpi(value, 'true') || strcmp(value, '1');
    end
end
flags = flags(:);
end
