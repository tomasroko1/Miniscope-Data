function run_stoixeion_crossday(dataDir, mappingCsv, outDir, animal, stoixeionDir)
% Cross-day ensemble fits using the MATLAB Stoixeion implementation.
% Requires mappings.csv from run_cross_registration_by_animal.py.

if nargin < 1 || isempty(dataDir), dataDir = '/mnt/NAS/Miniscopes/Reg_CA1/DataBase'; end
if nargin < 2 || isempty(mappingCsv), error('Pass the cross-registration mappings.csv path.'); end
if nargin < 3 || isempty(outDir), outDir = fullfile(pwd, 'results', 'stoixeion', 'cross_day'); end
if nargin < 4 || isempty(animal), animal = 'R005'; end
if nargin < 5 || isempty(stoixeionDir), stoixeionDir = getenv('STOIXEION_DIR'); end
if ~isempty(stoixeionDir) && exist(fullfile(stoixeionDir, 'Stoixeion.m'), 'file') == 2
    addpath(genpath(stoixeionDir));
end
animal = char(animal);
animalDir = fullfile(dataDir, animal);
if exist(animalDir, 'dir') ~= 7, error('No existe la carpeta: %s', animalDir); end
if exist(mappingCsv, 'file') ~= 2, error('No existe mappings.csv: %s', mappingCsv); end
if exist('Stoixeion', 'file') ~= 2, error('Stoixeion.m no esta en el MATLAB path.'); end
if nargout('Stoixeion') < 2, error('Aplica stoixeion_exports.patch antes de correr.'); end
if exist(outDir, 'dir') ~= 7, mkdir(outDir); end

map = readCrossDayMapping(mappingCsv);
oldVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
visibilityCleanup = onCleanup(@() set(groot, 'DefaultFigureVisible', oldVisibility)); %#ok<NASGU>
files = dir(fullfile(animalDir, '*_merged.mat'));
days = struct('file', {}, 'name', {});
for k = 1:numel(files)
    matPath = fullfile(animalDir, files(k).name);
    try
        day = loadDay(matPath, map);
        days(end + 1).file = matPath; %#ok<AGROW>
        days(end).name = day.name;
        fprintf('[DIA] %s: %s\n', files(k).name, day.name);
    catch ME
        fprintf('[OMITIDA] %s: %s\n', files(k).name, ME.message);
    end
end
if numel(days) < 2, error('Se necesitan dos .mat mergeados compatibles con el mapping.'); end

summaryPath = fullfile(outDir, 'cross_day_pairs.csv');
fid = fopen(summaryPath, 'w');
if fid < 0, error('No se pudo crear %s', summaryPath); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'animal,source_file,source_phase,target_file,target_phase,n_shared_cells,status,n_frames,n_significant_vectors,n_ensembles,output_dir\n');

for ia = 1:(numel(days)-1)
    for ib = (ia+1):numel(days)
        source = loadDay(days(ia).file, map);
        target = loadDay(days(ib).file, map);
        pairs = [(1:4)' (1:4)'];
        if strcmp(source.phaseNames{3}, 'TEST') && strcmp(target.phaseNames{4}, 'OF2')
            pairs(end + 1, :) = [3 4]; %#ok<AGROW>
        end
        for q = 1:size(pairs, 1)
            pa = pairs(q, 1); pb = pairs(q, 2);
            a = source.phases{pa}; b = target.phases{pb};
            common = intersect(a.globalIds, b.globalIds);
            if numel(common) < 20
                writeSummary(fid, animal, source.file, source.phaseNames{pa}, target.file, ...
                    target.phaseNames{pb}, numel(common), 'skipped_too_few_shared_cells', 0, 0, 0, '');
                continue;
            end
            [~, iaCol] = ismember(common, a.globalIds);
            [~, ibCol] = ismember(common, b.globalIds);
            spikesA = binarizeS(a.S(:, a.columns(iaCol)), 3);
            spikesB = binarizeS(b.S(:, b.columns(ibCol)), 3);
            spikes = [spikesA spikesB];
            ranges = [1 size(spikesA, 2); size(spikesA, 2)+1 size(spikes, 2)];
            pairName = sprintf('%s_%s_to_%s_%s', fileStem(source.file), ...
                source.phaseNames{pa}, fileStem(target.file), target.phaseNames{pb});
            pairName = regexprep(pairName, '[^A-Za-z0-9_.-]', '_');
            pairOut = fullfile(outDir, pairName);
            if exist(pairOut, 'dir') ~= 7, mkdir(pairOut); end

            coords = [(1:numel(common))' zeros(numel(common), 1)];
            close all force;
            [pools, diag] = Stoixeion(spikes, coords, []);
            if size(diag.ensemble_vectors, 1) ~= numel(diag.significant_frames)
                error('Vectores y marcos significativos no alinean.');
            end
            saveNativeFigures(pairOut, 8);
            saveTimeline(diag, ranges, source.name, source.phaseNames{pa}, ...
                target.name, target.phaseNames{pb}, pairOut, numel(common));
            saveCoreTable(pools, diag, common, pairOut);
            close all force;

            nFactors = size(diag.ensemble_vectors, 2);
            writeSummary(fid, animal, source.file, source.phaseNames{pa}, target.file, ...
                target.phaseNames{pb}, numel(common), 'ok', size(spikes, 2), ...
                numel(diag.significant_frames), nFactors, pairOut);
            fprintf('[OK] %s %s -> %s %s: %d celulas compartidas, %d ensembles.\n', ...
                source.name, source.phaseNames{pa}, target.name, target.phaseNames{pb}, numel(common), nFactors);
        end
        clear source target;
    end
end
fprintf('Resumen cross-day: %s\n', summaryPath);
end

function map = readCrossDayMapping(path)
fid = fopen(path, 'r');
if fid < 0, error('No se pudo abrir %s', path); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
header = parseCsvLine(fgetl(fid));
idCol = find(strcmp(header, 'global_cell_id'), 1);
if isempty(idCol), error('mappings.csv no contiene global_cell_id.'); end
qc = {'n_sessions','max_centroid_distance_pixels','mean_centroid_distance_pixels', ...
    'mean_centroid_height','mean_centroid_width'};
map.columns = [];
map.keys = {};
for c = 1:numel(header)
    if c ~= idCol && ~any(strcmp(header{c}, qc))
        map.columns(end + 1) = c; %#ok<AGROW>
        map.keys{end + 1} = sessionKey(header{c}); %#ok<AGROW>
    end
end
map.ids = [];
map.units = [];
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    row = parseCsvLine(line);
    if numel(row) < numel(header), row(end + 1:numel(header)) = {''}; end
    map.ids(end + 1, 1) = str2double(row{idCol}); %#ok<AGROW>
    values = nan(1, numel(map.columns));
    for c = 1:numel(map.columns), values(c) = str2double(row{map.columns(c)}); end
    map.units(end + 1, :) = values; %#ok<AGROW>
end
for c = 1:size(map.units, 2)
    present = map.units(isfinite(map.units(:, c)), c);
    if numel(unique(present)) ~= numel(present)
        error('unit_id repetido en una columna de session de mappings.csv.');
    end
end
end

function fields = parseCsvLine(line)
fields = {};
value = '';
quoted = false;
i = 1;
while i <= numel(line)
    ch = line(i);
    if ch == '"'
        if quoted && i < numel(line) && line(i+1) == '"'
            value(end+1) = '"'; %#ok<AGROW>
            i = i + 1;
        else
            quoted = ~quoted;
        end
    elseif ch == ',' && ~quoted
        fields{end+1} = value; %#ok<AGROW>
        value = '';
    else
        value(end+1) = ch; %#ok<AGROW>
    end
    i = i + 1;
end
fields{end+1} = value;
end

function day = loadDay(path, map)
d = load(path, 'act', 'sess');
if ~isfield(d, 'act') || ~isfield(d, 'sess'), error('falta act o sess'); end
a = d.act; s = d.sess;
if ~isfield(a, 'mapping') || ~isfield(a, 'S') || ~isfield(a, 't')
    error('falta act.mapping, act.S o act.t');
end
phaseS = cellPhases(a.S); phaseT = cellPhases(a.t);
if numel(phaseS) ~= 4 || numel(phaseT) ~= 4, error('se requieren cuatro fases'); end
if size(a.mapping, 2) ~= 4, error('act.mapping no tiene cuatro columnas'); end
paths = textPhases(s.sess_paths);
if numel(paths) ~= 4, error('sess.sess_paths no tiene cuatro rutas'); end
day.name = fileStem(path);
day.phaseNames = {'OF1','SAMPLE','TEST','OF2'};
if isfield(s, 'day_name')
    day.name = textValue(s.day_name);
    if strcmpi(day.name, 'HabL'), day.phaseNames = {'OF1','OF2','OF3','OF4'}; end
end
day.file = path; day.phases = cell(1, 4);
for p = 1:4
    S = double(phaseS{p}); t = double(phaseT{p}(:));
    unitIds = sort(double(a.mapping(isfinite(a.mapping(:,p)), p)));
    if size(S,1) ~= numel(t) || size(S,2) ~= numel(unitIds)
        error('S/t/mapping no alinean en fase %d', p);
    end
    if numel(unique(unitIds)) ~= numel(unitIds), error('unit_id duplicado en fase %d', p); end
    col = find(strcmp(map.keys, sessionKey(paths{p})), 1);
    if isempty(col), error('adquisicion %s no aparece en mappings.csv', paths{p}); end
    globalIds = nan(1, numel(unitIds));
    for k = 1:numel(unitIds)
        r = find(map.units(:,col) == unitIds(k), 1);
        if ~isempty(r), globalIds(k) = map.ids(r); end
    end
    valid = isfinite(globalIds);
    if numel(unique(globalIds(valid))) ~= sum(valid), error('global_cell_id duplicado en %s', paths{p}); end
    day.phases{p}.S = S;
    day.phases{p}.columns = find(valid);
    day.phases{p}.globalIds = globalIds(valid);
end
end

function c = cellPhases(x)
if iscell(x), c = reshape(x, 1, []); else, c = {x}; end
end

function c = textPhases(x)
if iscell(x), values = reshape(x, 1, []); else, values = {x}; end
c = cell(1, numel(values));
for k = 1:numel(values), c{k} = textValue(values{k}); end
end

function value = textValue(x)
if iscell(x) && ~isempty(x), x = x{1}; end
if isstring(x), x = char(x); end
if ischar(x), value = strtrim(x); else, value = ''; end
end

function key = sessionKey(path)
path = strrep(textValue(path), '\', '/');
parts = regexp(path, '/+', 'split');
parts = parts(~cellfun('isempty', parts));
if numel(parts) < 2, error('ruta de adquisicion inesperada: %s', path); end
key = [parts{end-1} '/' parts{end}];
end

function spikes = binarizeS(S, thresholdSD)
spikes = false(size(S,2), size(S,1));
for c = 1:size(S,2)
    signal = S(:,c);
    finite = signal(isfinite(signal));
    if numel(finite) < 2, continue; end
    threshold = thresholdSD * std(finite, 0);
    signal(~isfinite(signal)) = 0;
    spikes(c,:) = signal > threshold;
end
end

function saveNativeFigures(folder, maxFigure)
figs = findall(groot, 'Type', 'figure');
for k = 1:numel(figs)
    n = get(figs(k), 'Number');
    if n >= 1 && n <= maxFigure, saveas(figs(k), fullfile(folder, sprintf('stoixeion_%02d.png', n))); end
end
end

function saveTimeline(diag, ranges, dayA, phaseA, dayB, phaseB, folder, nCells)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 600]);
hold on;
nFactors = size(diag.ensemble_vectors, 2);
if nFactors == 0, nFactors = 1; end
colors = [0.82 0.90 0.98; 0.98 0.88 0.78];
for p = 1:2
    x1 = ranges(p,1)-0.5; x2 = ranges(p,2)+0.5;
    patch([x1 x2 x2 x1], [0.5 0.5 nFactors+0.5 nFactors+0.5], colors(p,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.7);
end
for f = 1:size(diag.ensemble_vectors,2)
    v = find(diag.ensemble_vectors(:,f) > 0);
    frames = diag.significant_frames(v);
    for k = 1:numel(frames)
        line([frames(k) frames(k)], [f-0.34 f+0.34], 'Color', [0.1 0.1 0.1], 'LineWidth', 0.7);
    end
end
boundary = ranges(1,2)+0.5;
line([boundary boundary], [0.5 nFactors+0.5], 'Color', 'k', 'LineWidth', 1.3);
set(gca, 'YDir', 'reverse', 'YTick', 1:size(diag.ensemble_vectors,2), ...
    'YTickLabel', ensembleLabels(size(diag.ensemble_vectors,2)));
ylim([0.2 nFactors+0.5]); xlim([0.5 ranges(2,2)+0.5]);
set(gca, 'XTick', [ranges(1,1) boundary ranges(2,2)], ...
    'XTickLabel', {'1', num2str(boundary), num2str(ranges(2,2))});
xlabel('Frame concatenado (limite vertical entre sesiones)'); ylabel('Ensemble');
title(sprintf('Stoixeion cross-day | %s %s -> %s %s | %d celulas', dayA, phaseA, dayB, phaseB, nCells));
text(mean(ranges(1,:)), 0.35, [dayA ' | ' phaseA], 'HorizontalAlignment', 'center');
text(mean(ranges(2,:)), 0.35, [dayB ' | ' phaseB], 'HorizontalAlignment', 'center');
set(gca, 'Position', [0.07 0.12 0.90 0.72]);
saveas(fig, fullfile(folder, 'cross_day_ensemble_timeline.png'));
close(fig);
end

function labels = ensembleLabels(n)
labels = cell(1,n);
for k = 1:n, labels{k} = sprintf('E%d',k); end
end

function saveCoreTable(pools, diag, ids, folder)
fid = fopen(fullfile(folder, 'core_neurons.csv'), 'w');
fprintf(fid, 'factor,rank,global_cell_id\n');
for f = 1:size(diag.ensemble_vectors,2)
    cells = round(pools(:,3,f));
    cells = unique(cells(isfinite(cells) & cells >= 1 & cells <= numel(ids)));
    for k = 1:numel(cells), fprintf(fid, '%d,%d,%d\n', f, k, ids(cells(k))); end
end
fclose(fid);
end

function writeSummary(fid, animal, fileA, phaseA, fileB, phaseB, nCells, status, nFrames, nVectors, nFactors, outputDir)
fprintf(fid, '%s,%s,%s,%s,%s,%d,%s,%d,%d,%d,%s\n', animal, fileStem(fileA), ...
    phaseA, fileStem(fileB), phaseB, nCells, status, nFrames, nVectors, nFactors, outputDir);
end

function stem = fileStem(path)
[~, stem] = fileparts(path);
end


