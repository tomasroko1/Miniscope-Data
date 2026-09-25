function run_stoixeion_miniscope(selection, dataDir, stoixeionDir, outDir)
% Run Stoixeion phasewise and on a four-phase concatenated population raster.
% Selections: pilot, global, phasewise, sd, xss, habl, all4, primary, or a session stem.

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
if nargout('Stoixeion') < 2
    error('Apply stoixeion_exports.patch to the cloned Stoixeion folder, then rerun.');
end

oldVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
restoreVisibility = onCleanup(@() set(groot, 'DefaultFigureVisible', oldVisibility)); %#ok<NASGU>
if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
writeProtocol(outDir, selection, dataDir, stoixeionDir);

summaryRows = {};
coreRows = {};
ensembleRows = {};
singularRows = {};
ensembleActivityRows = {};
overlapRows = {};
thresholdRows = {};
coreShuffleRows = {};
globalFactorPhaseRows = {};
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
        if ~isfield(activity, 'S') || ~isfield(activity, 't')
            fprintf('[OMITIDA] %s: falta S o t.\n', sessionID);
            continue;
        end
        phaseS = toPhases(activity.S);
        phaseT = toPhases(activity.t);
        nPhases = min(numel(phaseS), numel(phaseT));
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
            fprintf('[AVISO] %s: mapping no verificable; los cores quedan locales por fase.\n', sessionID);
        end

        if strcmpi(selection, 'global') && ~mappingOK
            fprintf('[OMITIDA GLOBAL] %s: hace falta mapping valido en las cuatro fases.\n', sessionID);
            continue;
        end

        if strcmpi(dayName, 'HabL')
            phaseNames = {'OF1', 'OF2', 'OF3', 'OF4'};
        else
            phaseNames = {'OF1', 'SAMPLE', 'TEST', 'OF2'};
        end

        phaseCoreSets = cell(1, nPhases);
        phaseSuccess = false(1, nPhases);
        dayOut = fullfile(outDir, animal, dayName, sessionStem);
        for pi = 1:nPhases
            if strcmpi(selection, 'global')
                phaseSuccess(pi) = true;
                continue;
            end
            phaseOut = fullfile(dayOut, phaseNames{pi});
            if exist(phaseOut, 'dir') ~= 7, mkdir(phaseOut); end
            try
                S = double(phaseS{pi});
                t = double(phaseT{pi}(:));
                cols = phaseColumns{pi};
                S = S(:, cols);
                if size(S, 1) ~= numel(t)
                    error('S y t no coinciden en frames.');
                end
                [spikes, activeFraction, thresholds] = binarizePhase(S, 3);
                nPhaseCells = size(spikes, 1);
                coords = [(1:nPhaseCells)' zeros(nPhaseCells, 1)];

                close all force;
                [pools, diagnostics] = Stoixeion(spikes, coords, []);
                saveStoixeionFigures(phaseOut);
                close all force;
                nFactors = size(diagnostics.ensemble_vectors, 2);
                newShuffleRows = plotCoreVsShuffled(spikes, pools, diagnostics, t, phaseOut, ...
                    animal, dayName, sessionFile, phaseNames{pi}, 199);
                if ~isempty(newShuffleRows)
                    if isempty(coreShuffleRows)
                        coreShuffleRows = newShuffleRows;
                    else
                        coreShuffleRows = [coreShuffleRows; newShuffleRows]; %#ok<AGROW>
                    end
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
                    allCells = find(diagnostics.all_ensemble_cells(:, factor) > 0);
                    for ci = 1:numel(allCells)
                        cellIndex = allCells(ci);
                        if mappingOK
                            globalID = commonRows(cellIndex) - 1;
                        else
                            globalID = NaN;
                        end
                        ensembleRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                            factor, cellIndex, globalID, ismember(cellIndex, localCells)}; %#ok<AGROW>
                    end
                end
                for rank = 1:numel(diagnostics.singular_values)
                    singularRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                        rank, diagnostics.singular_values(rank), ...
                        ismember(rank, diagnostics.selected_singular_ranks)}; %#ok<AGROW>
                end
                for vec = 1:size(diagnostics.ensemble_vectors, 1)
                    activeFactors = find(diagnostics.ensemble_vectors(vec, :) > 0);
                    for factor = activeFactors
                        frame = diagnostics.significant_frames(vec);
                        ensembleActivityRows(end + 1, :) = {animal, dayName, sessionFile, ...
                            phaseNames{pi}, frame, t(frame), factor}; %#ok<AGROW>
                    end
                end
                for ci = 1:nPhaseCells
                    if mappingOK
                        globalID = commonRows(ci) - 1;
                    else
                        globalID = NaN;
                    end
                    thresholdRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                        ci, globalID, thresholds(ci)}; %#ok<AGROW>
                end
                validThresholds = thresholds(isfinite(thresholds));
                if isempty(validThresholds), medianThreshold = NaN;
                else, medianThreshold = median(validThresholds); end
                summaryRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                    nPhaseCells, nCellsUsed, activeFraction, medianThreshold, mappingOK, ...
                    diagnostics.pks, numel(diagnostics.significant_frames), diagnostics.scut, ...
                    diagnostics.hcut, nFactors, coreCount, 'ok'}; %#ok<AGROW>
                phaseSuccess(pi) = true;
                fprintf('[OK] %s %s: %d celulas, %d factores.\n', ...
                    sessionID, phaseNames{pi}, nPhaseCells, nFactors);
            catch ME
                close all force;
                summaryRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                    numel(phaseColumns{pi}), nCellsUsed, NaN, NaN, mappingOK, ...
                    NaN, NaN, NaN, NaN, NaN, NaN, ME.message}; %#ok<AGROW>
                fprintf('[FALLO] %s %s: %s\n', sessionID, phaseNames{pi}, ME.message);
            end
            writeOutputs(outDir, summaryRows, coreRows, ensembleRows, singularRows, ...
                ensembleActivityRows, overlapRows, thresholdRows, coreShuffleRows, globalFactorPhaseRows);
        end

        if mappingOK && all(phaseSuccess)
            if ~strcmpi(selection, 'global')
                newOverlap = comparePhaseCores(phaseCoreSets, phaseNames, commonRows - 1, ...
                    animal, dayName, sessionFile, 999);
                if ~isempty(newOverlap)
                    if isempty(overlapRows)
                        overlapRows = newOverlap;
                    else
                        overlapRows = [overlapRows; newOverlap]; %#ok<AGROW>
                    end
                end
            end

            if ~strcmpi(selection, 'phasewise')
                try
                    globalResult = runGlobalConcatenated(phaseS, phaseT, phaseNames, phaseColumns, ...
                        commonRows, animal, dayName, sessionFile, nCellsUsed, dayOut);
                    if ~isempty(globalResult.summaryRows)
                        summaryRows = [summaryRows; globalResult.summaryRows]; %#ok<AGROW>
                    end
                    if ~isempty(globalResult.coreRows)
                        coreRows = [coreRows; globalResult.coreRows]; %#ok<AGROW>
                    end
                    if ~isempty(globalResult.ensembleRows)
                        ensembleRows = [ensembleRows; globalResult.ensembleRows]; %#ok<AGROW>
                    end
                    if ~isempty(globalResult.singularRows)
                        singularRows = [singularRows; globalResult.singularRows]; %#ok<AGROW>
                    end
                    if ~isempty(globalResult.ensembleActivityRows)
                        ensembleActivityRows = [ensembleActivityRows; globalResult.ensembleActivityRows]; %#ok<AGROW>
                    end
                    if ~isempty(globalResult.thresholdRows)
                        thresholdRows = [thresholdRows; globalResult.thresholdRows]; %#ok<AGROW>
                    end
                    if ~isempty(globalResult.coreShuffleRows)
                        coreShuffleRows = [coreShuffleRows; globalResult.coreShuffleRows]; %#ok<AGROW>
                    end
                    if ~isempty(globalResult.globalFactorPhaseRows)
                        globalFactorPhaseRows = [globalFactorPhaseRows; globalResult.globalFactorPhaseRows]; %#ok<AGROW>
                    end
                catch ME
                    summaryRows(end + 1, :) = {animal, dayName, sessionFile, 'GLOBAL_CONCAT', ...
                        nCellsUsed, nCellsUsed, NaN, NaN, mappingOK, NaN, NaN, NaN, NaN, NaN, NaN, ME.message}; %#ok<AGROW>
                    fprintf('[FALLO GLOBAL] %s: %s\n', sessionID, ME.message);
                end
            end
            writeOutputs(outDir, summaryRows, coreRows, ensembleRows, singularRows, ...
                ensembleActivityRows, overlapRows, thresholdRows, coreShuffleRows, globalFactorPhaseRows);
        end
    end
end
writeOutputs(outDir, summaryRows, coreRows, ensembleRows, singularRows, ...
    ensembleActivityRows, overlapRows, thresholdRows, coreShuffleRows, globalFactorPhaseRows);
fprintf('Resultados: %s\n', outDir);
end

function yes = selectedSession(selection, animal, dayName, sessionID)
switch lower(selection)
    case 'pilot'
        yes = strcmp(animal, 'R005') && strcmp(dayName, 'T1_SD_VEH');
    case 'sd'
        yes = ~isempty(strfind(lower(dayName), '_sd_')); %#ok<STREMP>
    case 'xss'
        yes = ~isempty(strfind(lower(dayName), '_xss_')); %#ok<STREMP>
    case 'habl'
        yes = strcmpi(dayName, 'HabL');
    case 'all4'
        yes = true;
    case 'global'
        lowerDay = lower(dayName);
        yes = strcmpi(dayName, 'HabL') || ~isempty(strfind(lowerDay, '_sd_')) || ...
            ~isempty(strfind(lowerDay, '_xss_')); %#ok<STREMP>
    case 'phasewise'
        lowerDay = lower(dayName);
        yes = strcmpi(dayName, 'HabL') || ~isempty(strfind(lowerDay, '_sd_')) || ...
            ~isempty(strfind(lowerDay, '_xss_')); %#ok<STREMP>
    case 'primary'
        lowerDay = lower(dayName);
        yes = strcmpi(dayName, 'HabL') || ~isempty(strfind(lowerDay, '_sd_')) || ...
            ~isempty(strfind(lowerDay, '_xss_')); %#ok<STREMP>
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

function [spikes, activeFraction, thresholds] = binarizePhase(S, thresholdSD)
nCells = size(S, 2);
nFrames = size(S, 1);
spikes = false(nCells, nFrames);
thresholds = nan(1, nCells);
for cellIndex = 1:nCells
    signal = S(:, cellIndex);
    finiteSignal = signal(isfinite(signal));
    if numel(finiteSignal) < 2, continue; end
    thresholds(cellIndex) = thresholdSD * std(finiteSignal, 0);
    signal(~isfinite(signal)) = 0;
    spikes(cellIndex, :) = signal > thresholds(cellIndex);
end
activeFraction = mean(spikes(:));
end

function saveStoixeionFigures(folder, maxFigure)
if nargin < 2, maxFigure = 9; end
figures = findall(groot, 'Type', 'figure');
for k = 1:numel(figures)
    number = get(figures(k), 'Number');
    if number < 1 || number > maxFigure, continue; end
    saveas(figures(k), fullfile(folder, sprintf('stoixeion_%02d.png', number)));
end
end

function result = runGlobalConcatenated(phaseS, phaseT, phaseNames, phaseColumns, ...
    commonRows, animal, dayName, sessionFile, nCellsUsed, dayOut)
result.summaryRows = {};
result.coreRows = {};
result.ensembleRows = {};
result.singularRows = {};
result.ensembleActivityRows = {};
result.thresholdRows = {};
result.coreShuffleRows = {};
result.globalFactorPhaseRows = {};

nPhases = numel(phaseNames);
phaseRaw = cell(1, nPhases);
phaseTimes = cell(1, nPhases);
frameRanges = zeros(nPhases, 2);
phaseDurations = zeros(1, nPhases);
cursor = 1;
for pi = 1:nPhases
    S = double(phaseS{pi});
    S = S(:, phaseColumns{pi});
    t = double(phaseT{pi}(:));
    if size(S, 1) ~= numel(t)
        error('S y t no coinciden en la fase %s.', phaseNames{pi});
    end
    phaseRaw{pi} = S;
    phaseTimes{pi} = t;
    frameRanges(pi, :) = [cursor, cursor + size(S, 1) - 1];
    cursor = frameRanges(pi, 2) + 1;
    dt = diff(t);
    dt = dt(isfinite(dt) & dt > 0);
    if isempty(dt), medianDt = 0.05; else, medianDt = median(dt); end
    finiteTimes = t(isfinite(t));
    if isempty(finiteTimes)
        phaseDurations(pi) = size(S, 1) * medianDt;
    else
        phaseDurations(pi) = max(finiteTimes) - min(finiteTimes) + medianDt;
    end
end

allS = vertcat(phaseRaw{:});
[spikes, activeFraction, thresholds] = binarizePhase(allS, 3);
nCells = size(spikes, 1);
coords = [(1:nCells)' zeros(nCells, 1)];
globalOut = fullfile(dayOut, 'GLOBAL_CONCATENATED');
if exist(globalOut, 'dir') ~= 7, mkdir(globalOut); end

close all force;
[pools, diagnostics] = Stoixeion(spikes, coords, []);
if size(diagnostics.ensemble_vectors, 1) ~= numel(diagnostics.significant_frames)
    error('Stoixeion devolvio vectores y marcos significativos desalineados.');
end
saveStoixeionFigures(globalOut, 8);
close all force;

nFactors = size(diagnostics.ensemble_vectors, 2);
coreCount = 0;
for factor = 1:nFactors
    poolIndex = round(pools(:, 3, factor));
    localCells = unique(poolIndex(isfinite(poolIndex) & poolIndex >= 1 & poolIndex <= nCells));
    coreCount = coreCount + numel(localCells);
    for ci = 1:numel(localCells)
        localCell = localCells(ci);
        mappedID = commonRows(localCell) - 1;
        result.coreRows(end + 1, :) = {animal, dayName, sessionFile, ...
            'GLOBAL_CONCAT', factor, localCell, mappedID}; %#ok<AGROW>
    end
    allCells = find(diagnostics.all_ensemble_cells(:, factor) > 0);
    for ci = 1:numel(allCells)
        localCell = allCells(ci);
        mappedID = commonRows(localCell) - 1;
        result.ensembleRows(end + 1, :) = {animal, dayName, sessionFile, ...
            'GLOBAL_CONCAT', factor, localCell, mappedID, ismember(localCell, localCells)}; %#ok<AGROW>
    end
end

for rank = 1:numel(diagnostics.singular_values)
    result.singularRows(end + 1, :) = {animal, dayName, sessionFile, ...
        'GLOBAL_CONCAT', rank, diagnostics.singular_values(rank), ...
        ismember(rank, diagnostics.selected_singular_ranks)}; %#ok<AGROW>
end
validThresholds = thresholds(isfinite(thresholds));
if isempty(validThresholds), medianThreshold = NaN;
else, medianThreshold = median(validThresholds); end
result.summaryRows(1, :) = {animal, dayName, sessionFile, 'GLOBAL_CONCAT', ...
    nCells, nCellsUsed, activeFraction, medianThreshold, true, diagnostics.pks, ...
    numel(diagnostics.significant_frames), diagnostics.scut, diagnostics.hcut, ...
    nFactors, coreCount, 'ok'};
for ci = 1:nCells
    result.thresholdRows(end + 1, :) = {animal, dayName, sessionFile, ...
        'GLOBAL_CONCAT', ci, commonRows(ci) - 1, thresholds(ci)}; %#ok<AGROW>
end

rateMatrix = zeros(nFactors, nPhases);
for factor = 1:nFactors
    factorFrames = diagnostics.significant_frames(diagnostics.ensemble_vectors(:, factor) > 0);
    for pi = 1:nPhases
        firstFrame = frameRanges(pi, 1);
        lastFrame = frameRanges(pi, 2);
        inPhase = factorFrames >= firstFrame & factorFrames <= lastFrame;
        nFactorVectors = sum(inPhase);
        nPhaseVectors = sum(diagnostics.significant_frames >= firstFrame & ...
            diagnostics.significant_frames <= lastFrame);
        if nPhaseVectors == 0, fraction = NaN;
        else, fraction = nFactorVectors / nPhaseVectors; end
        duration = phaseDurations(pi);
        if duration > 0, rate = nFactorVectors / duration * 60;
        else, rate = NaN; end
        rateMatrix(factor, pi) = rate;
        result.globalFactorPhaseRows(end + 1, :) = {animal, dayName, sessionFile, ...
            factor, phaseNames{pi}, size(phaseRaw{pi}, 1), duration, nPhaseVectors, ...
            nFactorVectors, rate, fraction}; %#ok<AGROW>
    end
end
saveGlobalFactorPhaseFigure(rateMatrix, phaseNames, globalOut, animal, dayName);

for vec = 1:size(diagnostics.ensemble_vectors, 1)
    globalFrame = diagnostics.significant_frames(vec);
    pi = find(globalFrame >= frameRanges(:, 1) & globalFrame <= frameRanges(:, 2), 1);
    if isempty(pi), continue; end
    localFrame = globalFrame - frameRanges(pi, 1) + 1;
    activeFactors = find(diagnostics.ensemble_vectors(vec, :) > 0);
    for factor = activeFactors
        result.ensembleActivityRows(end + 1, :) = {animal, dayName, sessionFile, ...
            ['GLOBAL_' phaseNames{pi}], globalFrame, phaseTimes{pi}(localFrame), factor}; %#ok<AGROW>
    end
end

for pi = 1:nPhases
    firstFrame = frameRanges(pi, 1);
    lastFrame = frameRanges(pi, 2);
    vectorMask = diagnostics.significant_frames >= firstFrame & ...
        diagnostics.significant_frames <= lastFrame;
    phaseDiagnostics = diagnostics;
    phaseDiagnostics.significant_frames = diagnostics.significant_frames(vectorMask) - firstFrame + 1;
    phaseDiagnostics.ensemble_vectors = diagnostics.ensemble_vectors(vectorMask, :);
    phaseSpikes = spikes(:, firstFrame:lastFrame);
    phaseOut = fullfile(globalOut, phaseNames{pi});
    if exist(phaseOut, 'dir') ~= 7, mkdir(phaseOut); end
    phaseShuffleRows = plotCoreVsShuffled(phaseSpikes, pools, phaseDiagnostics, ...
        phaseTimes{pi}, phaseOut, animal, dayName, sessionFile, ...
        ['GLOBAL_' phaseNames{pi}], 199);
    if ~isempty(phaseShuffleRows)
        result.coreShuffleRows = [result.coreShuffleRows; phaseShuffleRows]; %#ok<AGROW>
    end
end
end

function saveGlobalFactorPhaseFigure(rateMatrix, phaseNames, folder, animal, dayName)
fig = figure('Visible', 'off', 'Color', 'w');
if isempty(rateMatrix)
    axis off;
    text(0.5, 0.5, 'Stoixeion no selecciono factores globales', ...
        'HorizontalAlignment', 'center');
else
    imagesc(rateMatrix);
    set(gca, 'YDir', 'normal', 'XTick', 1:numel(phaseNames), 'XTickLabel', phaseNames, ...
        'YTick', 1:size(rateMatrix, 1));
    xlabel('Fase');
    ylabel('Ensemble global');
    title(sprintf('%s | %s | vectores significativos por minuto', animal, dayName));
    colorbar;
end
saveas(fig, fullfile(folder, 'global_factor_phase_activity.png'));
close(fig);
end

function rows = plotCoreVsShuffled(spikes, pools, diagnostics, t, folder, animal, dayName, sessionFile, phase, nShuffles)
rows = {};
if isempty(pools), return; end
nCells = size(spikes, 1);
nFrames = size(spikes, 2);
if numel(t) ~= nFrames || nCells < 2, return; end
dt = median(diff(t));
if ~isfinite(dt) || dt <= 0, dt = 1; end
binFrames = max(1, round(0.5 / dt));
nBins = floor(nFrames / binFrames);
if nBins < 1, return; end
nKeep = nBins * binFrames;
binTimes = mean(reshape(t(1:nKeep), binFrames, nBins), 1);

for factor = 1:size(diagnostics.ensemble_vectors, 2)
    poolIndex = round(pools(:, 3, factor));
    coreCells = unique(poolIndex(isfinite(poolIndex) & poolIndex >= 1 & poolIndex <= nCells));
    nCore = numel(coreCells);
    if nCore == 0, continue; end
    syncThreshold = max(1, ceil(nCore / 2));
    ensembleFrames = diagnostics.significant_frames( ...
        diagnostics.ensemble_vectors(:, factor) > 0);
    ensembleFrames = ensembleFrames(ensembleFrames >= 1 & ensembleFrames <= nFrames);
    ensembleOn = false(1, nFrames);
    ensembleOn(ensembleFrames) = true;
    ensembleBins = mean(reshape(double(ensembleOn(1:nKeep)), binFrames, nBins), 1);
    coreCounts = sum(spikes(coreCells, 1:nKeep), 1);
    coreBins = mean(reshape(coreCounts, binFrames, nBins), 1);
    coreSyncFraction = mean(coreBins >= syncThreshold);
    ensembleCoreCorr = NaN;
    if std(ensembleBins) > 0 && std(coreBins) > 0
        ensembleCoreCorr = corr(ensembleBins(:), coreBins(:));
    end
    identityBins = zeros(nShuffles, nBins);
    shiftBins = zeros(nShuffles, nBins);
    identitySync = zeros(nShuffles, 1);
    shiftSync = zeros(nShuffles, 1);

    seed = 20260925 + sum(double(animal)) + sum(double(dayName)) + factor;
    rng(seed, 'twister');
    for sh = 1:nShuffles
        randomCells = randperm(nCells, nCore);
        randomCounts = sum(spikes(randomCells, 1:nKeep), 1);
        identityBins(sh, :) = mean(reshape(randomCounts, binFrames, nBins), 1);
        identitySync(sh) = mean(identityBins(sh, :) >= syncThreshold);

        shifted = spikes(coreCells, 1:nKeep);
        for ci = 1:nCore
            shiftAmount = randi(max(nKeep - 1, 1));
            shifted(ci, :) = circshift(shifted(ci, :), [0 shiftAmount]);
        end
        shiftedCounts = sum(shifted, 1);
        shiftBins(sh, :) = mean(reshape(shiftedCounts, binFrames, nBins), 1);
        shiftSync(sh) = mean(shiftBins(sh, :) >= syncThreshold);
    end

    idP = (1 + sum(identitySync >= coreSyncFraction)) / (nShuffles + 1);
    shiftP = (1 + sum(shiftSync >= coreSyncFraction)) / (nShuffles + 1);
    rows(end + 1, :) = {animal, dayName, sessionFile, phase, factor, nCore, ...
        ensembleCoreCorr, coreSyncFraction, mean(identitySync), prctile(identitySync, 95), idP, ...
        mean(shiftSync), prctile(shiftSync, 95), shiftP}; %#ok<AGROW>

    idLo = prctile(identityBins, 5, 1);
    idHi = prctile(identityBins, 95, 1);
    shiftLo = prctile(shiftBins, 5, 1);
    shiftHi = prctile(shiftBins, 95, 1);
    fig = figure('Visible', 'off', 'Color', 'w');
    subplot(3, 1, 1); hold on;
    stairs(binTimes, ensembleBins, 'Color', [0.15 0.55 0.25], 'LineWidth', 1);
    ylim([-0.05 1.05]);
    ylabel('Ensemble');
    title(sprintf('%s ensemble %d: SVD activity and core checks', phase, factor));
    subplot(3, 1, 2); hold on;
    fill([binTimes fliplr(binTimes)], [idLo fliplr(idHi)], [0.82 0.82 0.82], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(binTimes, mean(identityBins, 1), 'Color', [0.35 0.35 0.35], 'LineWidth', 1);
    plot(binTimes, coreBins, 'Color', [0.85 0.20 0.16], 'LineWidth', 1.1);
    ylabel('Core cells active / 0.5 s');
    title('Core IDs vs random neuron IDs at the same times');
    legend('Random-ID 5-95%', 'Random-ID mean', 'Observed core', 'Location', 'best');
    subplot(3, 1, 3); hold on;
    fill([binTimes fliplr(binTimes)], [shiftLo fliplr(shiftHi)], [0.78 0.84 0.92], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55);
    plot(binTimes, mean(shiftBins, 1), 'Color', [0.25 0.42 0.68], 'LineWidth', 1);
    plot(binTimes, coreBins, 'Color', [0.85 0.20 0.16], 'LineWidth', 1.1);
    xlabel('Time (s)');
    ylabel('Core cells active / 0.5 s');
    title('Same core IDs, each neuron independently time-shifted');
    legend('Time-shift 5-95%', 'Time-shift mean', 'Observed core', 'Location', 'best');
    saveas(fig, fullfile(folder, sprintf('core_vs_shuffled_%02d.png', factor)));
    close(fig);
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

function writeOutputs(outDir, summaryRows, coreRows, ensembleRows, singularRows, ...
    ensembleActivityRows, overlapRows, thresholdRows, coreShuffleRows, globalFactorPhaseRows)
if ~isempty(summaryRows)
    T = cell2table(summaryRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'n_cells','n_common_cells','event_fraction','median_threshold_sd','mapping_valid', ...
        'pks','n_significant_vectors','scut','hcut','n_ensembles','n_core_memberships','status'});
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
if ~isempty(thresholdRows)
    T = cell2table(thresholdRows, 'VariableNames', {'animal','day_name','session', ...
        'phase','local_cell_index','mapped_cell_id','S_threshold'});
    writetable(T, fullfile(outDir, 'event_thresholds.csv'));
end
if ~isempty(ensembleRows)
    T = cell2table(ensembleRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'ensemble','local_cell_index','mapped_cell_id','is_core'});
    writetable(T, fullfile(outDir, 'ensemble_members.csv'));
end
if ~isempty(singularRows)
    T = cell2table(singularRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'singular_rank','singular_value','selected_by_stoixeion'});
    writetable(T, fullfile(outDir, 'singular_values.csv'));
end
if ~isempty(ensembleActivityRows)
    T = cell2table(ensembleActivityRows, 'VariableNames', {'animal','day_name','session', ...
        'phase','frame','time_s','ensemble'});
    writetable(T, fullfile(outDir, 'ensemble_activity.csv'));
end
if ~isempty(coreShuffleRows)
    T = cell2table(coreShuffleRows, 'VariableNames', {'animal','day_name','session', ...
        'phase','ensemble','n_core_cells','ensemble_core_corr','observed_sync_fraction', ...
        'random_ID_null_mean','random_ID_null_95','random_ID_p_unadjusted', ...
        'time_shift_null_mean','time_shift_null_95','time_shift_p_unadjusted'});
    writetable(T, fullfile(outDir, 'core_shuffle_summary.csv'));
end
if ~isempty(globalFactorPhaseRows)
    T = cell2table(globalFactorPhaseRows, 'VariableNames', {'animal','day_name','session', ...
        'ensemble','phase','n_phase_frames','phase_duration_s','n_significant_vectors', ...
        'factor_vectors','factor_vectors_per_min','fraction_significant_vectors'});
    writetable(T, fullfile(outDir, 'global_factor_phase_activity.csv'));
end
end

function writeProtocol(outDir, selection, dataDir, stoixeionDir)
fid = fopen(fullfile(outDir, 'protocol.txt'), 'w');
fprintf(fid, 'Selection: %s\nData: %s\nStoixeion: %s\n', selection, dataDir, stoixeionDir);
if strcmpi(selection, 'global')
    fprintf(fid, 'Only the four-phase concatenated analysis is run; standalone phasewise detection is skipped.\n');
elseif strcmpi(selection, 'phasewise')
    fprintf(fid, 'Only standalone phasewise detection is run; concatenated detection is skipped.\n');
else
    fprintf(fid, 'Each phase is analyzed independently on cells mapped across all four phases when available.\n');
end
if ~strcmpi(selection, 'global')
    fprintf(fid, 'Phasewise raster: CaImAn S > 3 * SD(S), calculated per cell and phase.\n');
end
if ~strcmpi(selection, 'phasewise')
    fprintf(fid, 'Global raster: the four phases are concatenated and thresholded with one pooled per-cell SD(S).\n');
end
fprintf(fid, 'Stoixeion defaults are retained: automatic pks/scut, TF-IDF, hcut=0.28, SVD and core selection.\n');
if ~strcmpi(selection, 'phasewise')
    fprintf(fid, 'Global factors are detected once across the four phases; activation counts are assigned back to the original phases.\n');
    fprintf(fid, 'The global null shuffles across the concatenated raster and does not preserve phase-specific rates; interpret it as exploratory and check phasewise results.\n');
    fprintf(fid, 'Global fitting is weighted by the number of significant vectors contributed by each phase; reported phase rates are normalized by duration.\n');
end
fprintf(fid, 'Core overlap uses Jaccard with 999 identity permutations.\n');
fprintf(fid, 'Core activity controls: random neuron IDs at the same times and independent circular shifts of the same core cells.\n');
fprintf(fid, 'All overlap and shuffle p values are exploratory and unadjusted.\n');
fprintf(fid, 'Stoixeion export patch exposes existing internal outputs; the detection steps are unchanged.\n');
if ~strcmpi(selection, 'phasewise')
    fprintf(fid, 'Figure 9 is omitted for the global call because sequence correlations can cross phase boundaries.\n');
end
fprintf(fid, 'Spatial figures are omitted because source cell centroids are not in the merged MAT files.\n');
fclose(fid);
end
