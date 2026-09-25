function run_stoixeion_miniscope(selection, dataDir, stoixeionDir, outDir)
% Run Stoixeion on CaImAn activity, independently for each phase.
% Selections: pilot, sd, habl, primary, or an animal_session stem.

rootDir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(selection), selection = 'pilot'; end
selection = char(selection);
if nargin < 2 || isempty(dataDir), dataDir = getenv('MINISCOPE_DATA_DIR'); end
if isempty(dataDir), dataDir = fullfile(rootDir, 'data'); end
if nargin < 3 || isempty(stoixeionDir), stoixeionDir = getenv('STOIXEION_DIR'); end
if nargin < 4 || isempty(outDir)
    outDir = fullfile(rootDir, 'results', 'stoixeion', char(selection));
end

if isempty(stoixeionDir) || exist(fullfile(stoixeionDir, 'Stoixeion.m'), 'file') ~= 2
    error('Set STOIXEION_DIR to the folder containing Stoixeion.m.');
end
if exist('pdist2', 'file') == 0 || exist('perfcurve', 'file') == 0
    error('Stoixeion requires Statistics and Machine Learning Toolbox functions pdist2 and perfcurve.');
end
addpath(genpath(stoixeionDir));
if exist('Stoixeion', 'file') ~= 2
    error('Stoixeion.m was not found on the MATLAB path.');
end

oldVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
restoreVisibility = onCleanup(@() set(groot, 'DefaultFigureVisible', oldVisibility)); %#ok<NASGU>
if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
writeProtocol(outDir, selection, dataDir, stoixeionDir);

summaryRows = {};
coreRows = {};
overlapRows = {};
animals = dir(fullfile(dataDir, 'R*'));
animals = animals([animals.isdir]);

for ai = 1:numel(animals)
    animal = animals(ai).name;
    if isempty(regexp(animal, '^R\d+$', 'once')), continue; end
    files = dir(fullfile(dataDir, animal, '*.mat'));
    for fi = 1:numel(files)
        sessionFile = files(fi).name;
        [~, sessionStem] = fileparts(sessionFile);
        sessionID = [animal '_' sessionStem];
        loaded = load(fullfile(dataDir, animal, sessionFile), 'act', 'sess');
        if ~isfield(loaded, 'act') || ~isfield(loaded, 'sess'), continue; end
        dayName = fieldText(loaded.sess, 'day_name');
        if ~selectedSession(selection, animal, dayName, sessionID), continue; end

        activity = loaded.act;
        if ~isfield(activity, 'S') || ~isfield(activity, 'C') || ~isfield(activity, 't')
            fprintf('[OMITIDA] %s: falta C, S o t.\n', sessionID);
            continue;
        end
        phaseS = toPhases(activity.S);
        phaseC = toPhases(activity.C);
        phaseT = toPhases(activity.t);
        nPhases = min([numel(phaseS), numel(phaseC), numel(phaseT)]);
        if nPhases ~= 4
            fprintf('[OMITIDA] %s: se esperaban 4 fases y hay %d.\n', sessionID, nPhases);
            continue;
        end

        mappingOK = false;
        phaseRows = cell(1, nPhases);
        if isfield(activity, 'mapping')
            map = double(activity.mapping);
            if size(map, 2) >= nPhases
                mappingOK = true;
                for pi = 1:nPhases
                    localIDs = map(:, pi);
                    rows = find(isfinite(localIDs));
                    if numel(rows) ~= size(phaseS{pi}, 2) || ...
                            numel(unique(localIDs(rows))) ~= numel(rows)
                        mappingOK = false;
                        break;
                    end
                    [~, order] = sort(localIDs(rows));
                    phaseRows{pi} = rows(order);
                end
            end
        end

        if mappingOK
            commonRows = phaseRows{1};
            for pi = 2:nPhases
                commonRows = commonRows(ismember(commonRows, phaseRows{pi}));
            end
            if numel(commonRows) < 30
                fprintf('[OMITIDA] %s: solo %d neuronas comunes entre fases.\n', sessionID, numel(commonRows));
                continue;
            end
            phaseColumns = cell(1, nPhases);
            for pi = 1:nPhases
                [~, phaseColumns{pi}] = ismember(commonRows, phaseRows{pi});
            end
            nCellsUsed = numel(commonRows);
        else
            commonRows = [];
            phaseColumns = cell(1, nPhases);
            for pi = 1:nPhases
                phaseColumns{pi} = 1:size(phaseS{pi}, 2);
            end
            nCellsUsed = NaN;
            fprintf('[AVISO] %s: mapping no verificable; núcleos quedarán locales por fase.\n', sessionID);
        end

        if strcmp(dayName, 'HabL')
            phaseNames = {'OF1', 'OF2', 'OF3', 'OF4'};
        else
            phaseNames = {'OF1', 'SAMPLE', 'TEST', 'OF2'};
        end

        phaseCoreSets = cell(1, nPhases);
        phaseSuccess = false(1, nPhases);
        dayOut = fullfile(outDir, animal, dayName, sessionStem);
        for pi = 1:nPhases
            phaseOut = fullfile(dayOut, phaseNames{pi});
            if exist(phaseOut, 'dir') ~= 7, mkdir(phaseOut); end
            try
                C = double(phaseC{pi});
                S = double(phaseS{pi});
                t = double(phaseT{pi}(:));
                cols = phaseColumns{pi};
                C = C(:, cols);
                S = S(:, cols);
                if size(C, 1) ~= size(S, 1) || size(S, 1) ~= numel(t)
                    error('C, S y t no coinciden en frames.');
                end
                [spikes, activeFraction, fallbackCells] = binarizePhase(C, S, t, 3);
                nPhaseCells = size(spikes, 1);
                coords = [(1:nPhaseCells)' zeros(nPhaseCells, 1)];

                close all force;
                pools = Stoixeion(spikes, coords, []);
                saveStoixeionFigures(phaseOut);
                close all force;

                if isempty(pools)
                    nFactors = 0;
                else
                    nFactors = size(pools, 3);
                end
                phaseCoreSets{pi} = cell(1, nFactors);
                coreCount = 0;
                for factor = 1:nFactors
                    poolIndex = round(pools(:, 3, factor));
                    localCells = unique(poolIndex(isfinite(poolIndex) & poolIndex >= 1 & poolIndex <= nPhaseCells));
                    if mappingOK
                        globalIDs = commonRows(localCells) - 1;
                    else
                        globalIDs = [];
                    end
                    phaseCoreSets{pi}{factor} = globalIDs;
                    coreCount = coreCount + numel(localCells);
                    for ci = 1:numel(localCells)
                        if mappingOK
                            globalID = commonRows(localCells(ci)) - 1;
                        else
                            globalID = NaN;
                        end
                        coreRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                            factor, localCells(ci), globalID}; %#ok<AGROW>
                    end
                end
                summaryRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                    nPhaseCells, nCellsUsed, activeFraction, fallbackCells, mappingOK, ...
                    nFactors, coreCount, 'ok'}; %#ok<AGROW>
                phaseSuccess(pi) = true;
                fprintf('[OK] %s %s: %d células, %d factores.\n', ...
                    sessionID, phaseNames{pi}, nPhaseCells, nFactors);
            catch ME
                close all force;
                summaryRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                    numel(phaseColumns{pi}), nCellsUsed, NaN, NaN, mappingOK, NaN, NaN, ME.message}; %#ok<AGROW>
                fprintf('[FALLO] %s %s: %s\n', sessionID, phaseNames{pi}, ME.message);
            end
            writeOutputs(outDir, summaryRows, coreRows, overlapRows);
        end

        if mappingOK && all(phaseSuccess)
            newOverlap = comparePhaseCores(phaseCoreSets, phaseNames, commonRows - 1, ...
                animal, dayName, sessionFile, 999);
            if ~isempty(newOverlap)
                if isempty(overlapRows)
                    overlapRows = newOverlap;
                else
                    overlapRows = [overlapRows; newOverlap]; %#ok<AGROW>
                end
            end
            writeOutputs(outDir, summaryRows, coreRows, overlapRows);
        end
    end
end
writeOutputs(outDir, summaryRows, coreRows, overlapRows);
fprintf('Resultados: %s\n', outDir);
end

function yes = selectedSession(selection, animal, dayName, sessionID)
switch lower(selection)
    case 'pilot'
        yes = strcmp(animal, 'R005') && strcmp(dayName, 'T1_SD_VEH');
    case 'sd'
        yes = ~isempty(strfind(dayName, '_SD_')); %#ok<STREMP>
    case 'habl'
        yes = strcmp(dayName, 'HabL');
    case 'primary'
        yes = strcmp(dayName, 'HabL') || ~isempty(strfind(dayName, '_SD_')); %#ok<STREMP>
    otherwise
        yes = strcmp(sessionID, selection);
end
end

function text = fieldText(s, fieldName)
text = '';
if ~isfield(s, fieldName), return; end
value = s.(fieldName);
if iscell(value) && ~isempty(value), value = value{1}; end
if isstring(value), value = char(value); end
if ischar(value), text = strtrim(value); end
end

function phases = toPhases(value)
if iscell(value)
    phases = reshape(value, 1, []);
else
    phases = {value};
end
end

function [spikes, activeFraction, fallbackCells] = binarizePhase(C, S, t, thresholdSD)
dt = diff(t);
if isempty(dt) || any(~isfinite(dt) | dt <= 0)
    error('Los tiempos deben crecer de forma estricta.');
end
dC = diff(C, 1, 1) ./ dt;
spikes = false(size(C, 2), size(C, 1));
fallbackCells = 0;
for cellIndex = 1:size(C, 2)
    derivative = dC(:, cellIndex);
    quiet = S(2:end, cellIndex) <= 0 & isfinite(derivative);
    noise = derivative(quiet);
    if numel(noise) < 100
        noise = derivative(isfinite(derivative));
        fallbackCells = fallbackCells + 1;
    end
    if numel(noise) < 2, continue; end
    cutoff = mean(noise) + thresholdSD * std(noise, 0);
    spikes(cellIndex, 2:end) = isfinite(derivative) & derivative > cutoff;
end
activeFraction = mean(spikes(:));
end

function saveStoixeionFigures(folder)
figures = findall(groot, 'Type', 'figure');
for k = 1:numel(figures)
    number = get(figures(k), 'Number');
    if number < 1 || number > 9, continue; end
    saveas(figures(k), fullfile(folder, sprintf('stoixeion_%02d.png', number)));
end
end

function rows = comparePhaseCores(coreSets, phaseNames, populationIDs, animal, dayName, sessionFile, nShuffles)
rows = {};
rng(20260925, 'twister');
for p1 = 1:numel(coreSets)
    for p2 = p1 + 1:numel(coreSets)
        for e1 = 1:numel(coreSets{p1})
            A = unique(coreSets{p1}{e1});
            for e2 = 1:numel(coreSets{p2})
                B = unique(coreSets{p2}{e2});
                if isempty(A) || isempty(B), continue; end
                nIntersection = numel(intersect(A, B));
                nUnion = numel(union(A, B));
                observed = nIntersection / max(nUnion, 1);
                nullValues = zeros(nShuffles, 1);
                for sh = 1:nShuffles
                    shuffledB = populationIDs(randperm(numel(populationIDs), min(numel(B), numel(populationIDs))));
                    nullValues(sh) = numel(intersect(A, shuffledB)) / ...
                        max(numel(union(A, shuffledB)), 1);
                end
                pValue = (1 + sum(nullValues >= observed)) / (nShuffles + 1);
                rows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{p1}, e1, ...
                    phaseNames{p2}, e2, nIntersection, nUnion, observed, pValue}; %#ok<AGROW>
            end
        end
    end
end
end

function writeOutputs(outDir, summaryRows, coreRows, overlapRows)
if ~isempty(summaryRows)
    T = cell2table(summaryRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'n_cells','n_common_cells','event_fraction','fallback_cells','mapping_valid', ...
        'n_ensembles','n_core_memberships','status'});
    writetable(T, fullfile(outDir, 'phase_summary.csv'));
end
if ~isempty(coreRows)
    T = cell2table(coreRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'ensemble','local_cell_index','mapped_cell_id'});
    writetable(T, fullfile(outDir, 'core_members.csv'));
end
if ~isempty(overlapRows)
    T = cell2table(overlapRows, 'VariableNames', {'animal','day_name','session', ...
        'phase_a','ensemble_a','phase_b','ensemble_b','n_overlap','n_union', ...
        'jaccard','p_permutation_unadjusted'});
    writetable(T, fullfile(outDir, 'core_overlap.csv'));
end
end

function writeProtocol(outDir, selection, dataDir, stoixeionDir)
fid = fopen(fullfile(outDir, 'protocol.txt'), 'w');
fprintf(fid, 'Selection: %s\nData: %s\nStoixeion: %s\n', selection, dataDir, stoixeionDir);
fprintf(fid, 'Each phase is analyzed independently on cells mapped across all four phases when available.\n');
fprintf(fid, 'Raster: positive dC/dt > quiescent-frame mean + 3 SD; quiet frames are S <= 0.\n');
fprintf(fid, 'Stoixeion defaults are retained: automatic pks/scut, TF-IDF, hcut=0.28, SVD and core selection.\n');
fprintf(fid, 'Core overlap is a separate exploratory Jaccard test with 999 identity permutations; p values are unadjusted.\n');
fprintf(fid, 'Figure 10 is omitted because source cell centroids are not in the merged MAT files.\n');
fclose(fid);
end
