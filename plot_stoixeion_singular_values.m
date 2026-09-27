function plot_stoixeion_singular_values(varargin)
% Recreate singular_values.png from existing CSV exports, without rerunning Stoixeion.
for rootIndex = 1:nargin
    outDir = char(varargin{rootIndex});
    summaryFile = fullfile(outDir, 'phase_summary.csv');
    singularFile = fullfile(outDir, 'singular_values.csv');
    if exist(summaryFile, 'file') ~= 2 || exist(singularFile, 'file') ~= 2
        continue;
    end
    summary = readtable(summaryFile);
    singular = readtable(singularFile);
    for row = 1:height(summary)
        animal = summary.animal{row};
        dayName = summary.day_name{row};
        sessionFile = summary.session{row};
        phase = summary.phase{row};
        match = strcmp(singular.animal, animal) & ...
            strcmp(singular.day_name, dayName) & ...
            strcmp(singular.session, sessionFile) & ...
            strcmp(singular.phase, phase);
        if ~any(match), continue; end

        [ranks, order] = sort(singular.singular_rank(match));
        values = singular.singular_value(match);
        values = values(order);
        nRanks = min(numel(values), round(summary.n_cells(row) / 6));
        if nRanks < 1, continue; end
        use = ranks >= 1 & ranks <= nRanks & isfinite(values) & values > 0;
        if ~any(use), continue; end

        [~, sessionStem] = fileparts(sessionFile);
        phaseDir = fullfile(outDir, animal, dayName, sessionStem, phase);
        if exist(phaseDir, 'dir') ~= 7, mkdir(phaseDir); end
        fig = figure('Visible', 'off', 'Color', 'w');
        semilogy(ranks(use), values(use), 'k-', 'LineWidth', 0.75);
        title('C: Singular values');
        xlabel('Singular rank');
        ylabel('Singular value');
        xlim([1 max(nRanks, 2)]);
        logValues = log10(values(use));
        logRange = max(logValues) - min(logValues);
        margin = 0.05 * logRange;
        if margin == 0, margin = 0.05; end
        yLimits = 10 .^ [min(logValues) - margin, max(logValues) + margin];
        ylim(yLimits);
        decades = floor(log10(yLimits(1))):ceil(log10(yLimits(2)));
        yTicks = sort(reshape([1; 2; 3; 5] * 10 .^ decades, 1, []));
        yTicks = yTicks(yTicks >= yLimits(1) & yTicks <= yLimits(2));
        if numel(yTicks) > 8
            yTicks = 10 .^ decades;
            yTicks = yTicks(yTicks >= yLimits(1) & yTicks <= yLimits(2));
        end
        if ~isempty(yTicks)
            set(gca, 'YTick', yTicks, 'YTickLabel', ...
                arrayfun(@(value) sprintf('%g', value), yTicks, 'UniformOutput', false));
        end
        set(gca, 'Color', 'w');
        saveas(fig, fullfile(phaseDir, 'singular_values.png'));
        close(fig);
        fprintf('Saved %s %s %s %s: ranks 1:%d\n', ...
            animal, dayName, sessionStem, phase, nRanks);
    end
end
end
