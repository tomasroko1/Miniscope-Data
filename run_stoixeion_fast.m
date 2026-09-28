function run_stoixeion_fast(selection, dataDir, stoixeionDir, outDir, exportFigures)
% RUN_STOIXEION_FAST Optimized, fast runner for Stoixeion.
%   Runs phasewise and global concatenated analysis for targeted sessions
%   (default: 'target9' for HabL and SD sessions in R004, R005, R006).
%   Uses 99 shuffles (denominator 100) and disables heavy MATLAB figure rendering.

rootDir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(selection), selection = 'target9'; end
selection = char(selection);
if nargin < 2 || isempty(dataDir), dataDir = getenv('MINISCOPE_DATA_DIR'); end
if isempty(dataDir), dataDir = fullfile(rootDir, 'data'); end
if nargin < 3 || isempty(stoixeionDir), stoixeionDir = getenv('STOIXEION_DIR'); end
if isempty(stoixeionDir), stoixeionDir = fullfile(rootDir, 'external', 'Stoixeion'); end
if nargin < 4 || isempty(outDir)
    outDir = fullfile(rootDir, 'results', 'stoixeion', ['fast_' char(selection)]);
end
if nargin < 5 || isempty(exportFigures), exportFigures = false; end

nShuffles = 99; % 99 shuffles gives clean p-value steps of 0.01

if isempty(stoixeionDir) || exist(fullfile(stoixeionDir, 'Stoixeion.m'), 'file') ~= 2
    error('Set STOIXEION_DIR to the folder containing Stoixeion.m.');
end
addpath(genpath(stoixeionDir));
if exist('Stoixeion', 'file') ~= 2
    error('Stoixeion.m was not found on the MATLAB path.');
end

oldVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
restoreVisibility = onCleanup(@() set(groot, 'DefaultFigureVisible', oldVisibility)); %#ok<NASGU>
if exist(outDir, 'dir') ~= 7, mkdir(outDir); end

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
        if ~isTargetSession(selection, animal, dayName, sessionID), continue; end

        activity = loaded.act;
        if ~isfield(activity, 'S') || ~isfield(activity, 't'), continue; end
        phaseS = toPhases(activity.S);
        phaseT = toPhases(activity.t);
        nPhases = min(numel(phaseS), numel(phaseT));
        if nPhases ~= 4, continue; end

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
                fprintf('[OMITIDA] %s: solo %d neuronas comunes.\n', sessionID, numel(commonRows));
                continue;
            end
            phaseColumns = cell(1, nPhases);
            for pi = 1:nPhases
                [~, phaseColumns{pi}] = ismember(commonRows, phaseRows{pi});
            end
            nCellsUsed = numel(commonRows);
        else
            continue; % Require valid mapping for within-day tracking
        end

        if strcmpi(dayName, 'HabL')
            phaseNames = {'OF1', 'OF2', 'OF3', 'OF4'};
        else
            phaseNames = {'OF1', 'SAMPLE', 'TEST', 'OF2'};
        end

        phaseCoreSets = cell(1, nPhases);
        phaseSuccess = false(1, nPhases);
        dayOut = fullfile(outDir, animal, dayName, sessionStem);

        % --- 1. PHASEWISE ANALYSIS ---
        fprintf('\n=== Procesando %s (%s) Phasewise ===\n', sessionID, dayName);
        for pi = 1:nPhases
            phaseOut = fullfile(dayOut, phaseNames{pi});
            if exportFigures && exist(phaseOut, 'dir') ~= 7, mkdir(phaseOut); end
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
                if exportFigures
                    saveStoixeionFigures(phaseOut);
                end
                close all force;
                nFactors = size(diagnostics.ensemble_vectors, 2);

                % Core vs Shuffled check (with 99 shuffles)
                newShuffleRows = computeCoreVsShuffled(spikes, pools, diagnostics, t, ...
                    animal, dayName, sessionFile, phaseNames{pi}, nShuffles, exportFigures, phaseOut);
                if ~isempty(newShuffleRows)
                    coreShuffleRows = [coreShuffleRows; newShuffleRows]; %#ok<AGROW>
                end

                phaseCoreSets{pi} = cell(1, nFactors);
                coreCount = 0;
                for factor = 1:nFactors
                    poolIndex = round(pools(:, 3, factor));
                    localCells = unique(poolIndex(isfinite(poolIndex) & poolIndex >= 1 & poolIndex <= nPhaseCells));
                    globalIDs = commonRows(localCells) - 1;
                    phaseCoreSets{pi}{factor} = globalIDs;
                    coreCount = coreCount + numel(localCells);
                    for ci = 1:numel(localCells)
                        coreRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                            factor, localCells(ci), globalIDs(ci)}; %#ok<AGROW>
                    end
                    allCells = find(diagnostics.all_ensemble_cells(:, factor) > 0);
                    for ci = 1:numel(allCells)
                        cellIndex = allCells(ci);
                        globalID = commonRows(cellIndex) - 1;
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
                    thresholdRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                        ci, commonRows(ci) - 1, thresholds(ci)}; %#ok<AGROW>
                end

                validThresholds = thresholds(isfinite(thresholds));
                if isempty(validThresholds), medianThreshold = NaN;
                else, medianThreshold = median(validThresholds); end
                summaryRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                    nPhaseCells, nCellsUsed, activeFraction, medianThreshold, true, ...
                    diagnostics.pks, numel(diagnostics.significant_frames), diagnostics.scut, ...
                    diagnostics.hcut, nFactors, coreCount, 'ok'}; %#ok<AGROW>
                phaseSuccess(pi) = true;
                fprintf('  [OK] %s: %d celulas, %d ensambles.\n', phaseNames{pi}, nPhaseCells, nFactors);
            catch ME
                close all force;
                summaryRows(end + 1, :) = {animal, dayName, sessionFile, phaseNames{pi}, ...
                    numel(phaseColumns{pi}), nCellsUsed, NaN, NaN, false, ...
                    NaN, NaN, NaN, NaN, NaN, NaN, ME.message}; %#ok<AGROW>
                fprintf('  [FALLO] %s: %s\n', phaseNames{pi}, ME.message);
            end
        end

        % Compare Phase Cores (with 99 shuffles)
        if all(phaseSuccess)
            newOverlap = comparePhaseCores(phaseCoreSets, phaseNames, commonRows - 1, ...
                animal, dayName, sessionFile, nShuffles);
            if ~isempty(newOverlap)
                overlapRows = [overlapRows; newOverlap]; %#ok<AGROW>
            end

            % --- 2. GLOBAL CONCATENATED ANALYSIS ---
            fprintf('=== Procesando %s (%s) Global Concatenado ===\n', sessionID, dayName);
            try
                globalResult = runGlobalConcatenated(phaseS, phaseT, phaseNames, phaseColumns, ...
                    commonRows, animal, dayName, sessionFile, nCellsUsed, dayOut, nShuffles, exportFigures);
                if ~isempty(globalResult.summaryRows), summaryRows = [summaryRows; globalResult.summaryRows]; end
                if ~isempty(globalResult.coreRows), coreRows = [coreRows; globalResult.coreRows]; end
                if ~isempty(globalResult.ensembleRows), ensembleRows = [ensembleRows; globalResult.ensembleRows]; end
                if ~isempty(globalResult.singularRows), singularRows = [singularRows; globalResult.singularRows]; end
                if ~isempty(globalResult.ensembleActivityRows), ensembleActivityRows = [ensembleActivityRows; globalResult.ensembleActivityRows]; end
                if ~isempty(globalResult.thresholdRows), thresholdRows = [thresholdRows; globalResult.thresholdRows]; end
                if ~isempty(globalResult.coreShuffleRows), coreShuffleRows = [coreShuffleRows; globalResult.coreShuffleRows]; end
                if ~isempty(globalResult.globalFactorPhaseRows), globalFactorPhaseRows = [globalFactorPhaseRows; globalResult.globalFactorPhaseRows]; end
                fprintf('  [OK] Global completado exitosamente.\n');
            catch ME
                fprintf('  [FALLO GLOBAL] %s: %s\n', sessionID, ME.message);
            end
        end

        % Save intermediate outputs incrementally
        writeOutputs(outDir, summaryRows, coreRows, ensembleRows, singularRows, ...
            ensembleActivityRows, overlapRows, thresholdRows, coreShuffleRows, globalFactorPhaseRows);
    end
end

writeOutputs(outDir, summaryRows, coreRows, ensembleRows, singularRows, ...
    ensembleActivityRows, overlapRows, thresholdRows, coreShuffleRows, globalFactorPhaseRows);
fprintf('\n>>> Corrida completada. Resultados en: %s <<<\n', outDir);
end

% --- TARGET SESSIONS SELECTOR ---
function yes = isTargetSession(selection, animal, dayName, sessionID)
switch lower(selection)
    case 'target9'
        % The 9 focal sessions: HabL, SD VEH, SD CNO for R004, R005, R006
        isHabL = strcmpi(dayName, 'HabL');
        isSD = ~isempty(strfind(lower(dayName), '_sd_')); %#ok<STREMP>
        yes = (isHabL || isSD) && ismember(animal, {'R004', 'R005', 'R006'});
    case {'r004', 'r005', 'r006'}
        isHabL = strcmpi(dayName, 'HabL');
        isSD = ~isempty(strfind(lower(dayName), '_sd_')); %#ok<STREMP>
        yes = (isHabL || isSD) && strcmpi(animal, selection);
    case 'sd'
        yes = ~isempty(strfind(lower(dayName), '_sd_')); %#ok<STREMP>
    case 'habl'
        yes = strcmpi(dayName, 'HabL');
    otherwise
        yes = strcmp(sessionID, selection) || strcmp(dayName, selection);
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

function rows = computeCoreVsShuffled(spikes, pools, diagnostics, t, animal, dayName, ...
    sessionFile, phase, nShuffles, exportFigures, folder)
rows = {};
nFactors = size(diagnostics.ensemble_vectors, 2);
if nFactors < 1, return; end

dt = diff(t);
dt = dt(isfinite(dt) & dt > 0);
if isempty(dt), medianDt = 0.05; else, medianDt = median(dt); end
binFrames = max(1, round(0.5 / medianDt));
nFrames = size(spikes, 2);
nBins = floor(nFrames / binFrames);
if nBins < 2, return; end

nKeep = nBins * binFrames;
nCells = size(spikes, 1);

for factor = 1:nFactors
    poolIndex = round(pools(:, 3, factor));
    coreCells = unique(poolIndex(isfinite(poolIndex) & poolIndex >= 1 & poolIndex <= nCells));
    nCore = numel(coreCells);
    if nCore < 1, continue; end

    syncThreshold = max(1, ceil(0.5 * nCore));
    coreSpikes = sum(spikes(coreCells, 1:nKeep), 1);
    coreBins = mean(reshape(coreSpikes, binFrames, nBins), 1);
    coreSyncFraction = mean(coreBins >= syncThreshold);

    ensembleFrames = diagnostics.significant_frames(diagnostics.ensemble_vectors(:, factor) > 0);
    ensembleSignal = zeros(1, nKeep);
    validFrames = ensembleFrames(ensembleFrames >= 1 & ensembleFrames <= nKeep);
    ensembleSignal(validFrames) = 1;
    ensembleBins = mean(reshape(ensembleSignal, binFrames, nBins), 1);

    coreCentered = coreBins - mean(coreBins);
    ensCentered = ensembleBins - mean(ensembleBins);
    denom = sqrt(sum(coreCentered .^ 2) * sum(ensCentered .^ 2));
    if denom > 0, ensembleCoreCorr = sum(coreCentered .* ensCentered) / denom;
    else, ensembleCoreCorr = NaN; end

    identitySync = zeros(nShuffles, 1);
    shiftSync = zeros(nShuffles, 1);

    seed = 20260925 + sum(double(animal)) + sum(double(dayName)) + factor;
    rng(seed, 'twister');
    for sh = 1:nShuffles
        randomCells = randperm(nCells, nCore);
        randomCounts = sum(spikes(randomCells, 1:nKeep), 1);
        idBins = mean(reshape(randomCounts, binFrames, nBins), 1);
        identitySync(sh) = mean(idBins >= syncThreshold);

        shifted = spikes(coreCells, 1:nKeep);
        for ci = 1:nCore
            shiftAmount = randi(max(nKeep - 1, 1));
            shifted(ci, :) = circshift(shifted(ci, :), [0 shiftAmount]);
        end
        shiftedCounts = sum(shifted, 1);
        shBins = mean(reshape(shiftedCounts, binFrames, nBins), 1);
        shiftSync(sh) = mean(shBins >= syncThreshold);
    end

    idP = (1 + sum(identitySync >= coreSyncFraction)) / (nShuffles + 1);
    shiftP = (1 + sum(shiftSync >= coreSyncFraction)) / (nShuffles + 1);
    rows(end + 1, :) = {animal, dayName, sessionFile, phase, factor, nCore, ...
        ensembleCoreCorr, coreSyncFraction, mean(identitySync), prctile(identitySync, 95), idP, ...
        mean(shiftSync), prctile(shiftSync, 95), shiftP}; %#ok<AGROW>
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

function result = runGlobalConcatenated(phaseS, phaseT, phaseNames, phaseColumns, ...
    commonRows, animal, dayName, sessionFile, nCellsUsed, dayOut, nShuffles, exportFigures)
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
    phaseRaw{pi} = S;
    phaseTimes{pi} = t;
    frameRanges(pi, :) = [cursor, cursor + size(S, 1) - 1];
    cursor = frameRanges(pi, 2) + 1;
    dt = diff(t);
    dt = dt(isfinite(dt) & dt > 0);
    if isempty(dt), medianDt = 0.05; else, medianDt = median(dt); end
    phaseDurations(pi) = size(S, 1) * medianDt;
end

concatS = cat(1, phaseRaw{:});
[spikes, activeFraction, thresholds] = binarizePhase(concatS, 3);
nCells = size(spikes, 1);
coords = [(1:nCells)' zeros(nCells, 1)];

globalOut = fullfile(dayOut, 'GLOBAL_CONCATENATED');
if exportFigures && exist(globalOut, 'dir') ~= 7, mkdir(globalOut); end

close all force;
[pools, diagnostics] = Stoixeion(spikes, coords, []);
close all force;

nFactors = size(diagnostics.ensemble_vectors, 2);
coreCount = 0;
for factor = 1:nFactors
    poolIndex = round(pools(:, 3, factor));
    localCells = unique(poolIndex(isfinite(poolIndex) & poolIndex >= 1 & poolIndex <= nCells));
    globalIDs = commonRows(localCells) - 1;
    coreCount = coreCount + numel(localCells);
    for ci = 1:numel(localCells)
        result.coreRows(end + 1, :) = {animal, dayName, sessionFile, ...
            'GLOBAL_CONCAT', factor, localCells(ci), globalIDs(ci)}; %#ok<AGROW>
    end
    allCells = find(diagnostics.all_ensemble_cells(:, factor) > 0);
    for ci = 1:numel(allCells)
        cellIndex = allCells(ci);
        result.ensembleRows(end + 1, :) = {animal, dayName, sessionFile, ...
            'GLOBAL_CONCAT', factor, cellIndex, commonRows(cellIndex) - 1, ...
            ismember(cellIndex, localCells)}; %#ok<AGROW>
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
        result.globalFactorPhaseRows(end + 1, :) = {animal, dayName, sessionFile, ...
            factor, phaseNames{pi}, size(phaseRaw{pi}, 1), duration, nPhaseVectors, ...
            nFactorVectors, rate, fraction}; %#ok<AGROW>
    end
end

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
    phaseShuffleRows = computeCoreVsShuffled(phaseSpikes, pools, phaseDiagnostics, ...
        phaseTimes{pi}, animal, dayName, sessionFile, ...
        ['GLOBAL_' phaseNames{pi}], nShuffles, exportFigures, phaseOut);
    if ~isempty(phaseShuffleRows)
        result.coreShuffleRows = [result.coreShuffleRows; phaseShuffleRows]; %#ok<AGROW>
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
if ~isempty(ensembleRows)
    T = cell2table(ensembleRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'ensemble','cell_index','mapped_cell_id','is_core'});
    writetable(T, fullfile(outDir, 'ensemble_members.csv'));
end
if ~isempty(singularRows)
    T = cell2table(singularRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'singular_rank','singular_value','selected_by_stoixeion'});
    writetable(T, fullfile(outDir, 'singular_values.csv'));
end
if ~isempty(ensembleActivityRows)
    T = cell2table(ensembleActivityRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'frame','time_s','ensemble'});
    writetable(T, fullfile(outDir, 'ensemble_activity.csv'));
end
if ~isempty(overlapRows)
    T = cell2table(overlapRows, 'VariableNames', {'animal','day_name','session','phase_a', ...
        'ensemble_a','phase_b','ensemble_b','n_overlap','n_union','jaccard','p_permutation_unadjusted'});
    writetable(T, fullfile(outDir, 'core_overlap.csv'));
end
if ~isempty(thresholdRows)
    T = cell2table(thresholdRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'local_cell_index','mapped_cell_id','threshold'});
    writetable(T, fullfile(outDir, 'event_thresholds.csv'));
end
if ~isempty(coreShuffleRows)
    T = cell2table(coreShuffleRows, 'VariableNames', {'animal','day_name','session','phase', ...
        'ensemble','n_core','ensemble_core_corr','observed_sync_fraction','mean_id_sync', ...
        'p95_id_sync','p_id_sync','mean_shift_sync','p95_shift_sync','p_shift_sync'});
    writetable(T, fullfile(outDir, 'core_shuffle_summary.csv'));
end
if ~isempty(globalFactorPhaseRows)
    T = cell2table(globalFactorPhaseRows, 'VariableNames', {'animal','day_name','session', ...
        'ensemble','phase','n_phase_frames','phase_duration_s','n_significant_vectors', ...
        'factor_vectors','factor_vectors_per_min','fraction_significant_vectors'});
    writetable(T, fullfile(outDir, 'global_factor_phase_activity.csv'));
end
end
