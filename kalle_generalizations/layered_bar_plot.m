function [bc, be] = layered_bar_plot(oo, selectedRows, numBins, fig)
% layered_bar_plot - Creates stacked bar plot where layers represent normalized category counts.
%
% INPUTS:
%   oo           - Matrix, row 1 = x-locations, other rows = categorical values
%   selectedRows - Vector of rows to include (e.g., 2:6)
%   numBins      - Number of bins to group x-values (default: 20)
% OUTPUTS:
%   bc           - Count for ordered bins
%   be           - Lower edge of the ordered bins

if nargin < 3, numBins = 20; end

x = oo(1,:);
layers = oo(selectedRows, :);
if layers == 0
    warning("No data presented.");
    return
end

% Find all unique categories across selected rows
categories = unique(layers(:));
numCats = numel(categories);

%% If the bins are too tight, you get this error:
%   Error using bar (line 182)
%   XData values must be unique.
%   Error in layered_bar_plot (line 57)
% Workaround: add some padding.
binEdges = max(x) - min(x);
if (binEdges <= 10*eps*numBins)
    % Minimum bins
    binEdges = 0:numBins;
    binEdges = binEdges *10*eps;
    binEdges = binEdges + mean(x) - mean(binEdges);
    binCenters = binEdges(1:end-1);
else
    % Compute bin edges
    binEdges = linspace(min(x) - binEdges, max(x) + binEdges, numBins + 1);
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
end

% Initialize category frequency matrix
categoryCounts = zeros(numCats, numBins);
binCounts = zeros(1, numBins+1);

%% Count normalized category frequency per bin
for b = 1:numBins
    inBin = x >= binEdges(b) & x < binEdges(b+1);
    sampleCategories = layers(:, inBin);  % size: numLayers x samplesInBin
    binCounts(b) = nnz(inBin);

    for c = 1:numCats
        catVal = categories(c);
        % Count how many samples in each column contain this category
        hasCat = sum(sampleCategories == catVal, 1);  % across rows
        % Normalize by number of active layers (== length(selectedRows))
        categoryCounts(c, b) = sum(hasCat) / length(selectedRows);
    end
end

%% Plot stacked bar chart
figure(fig);
hh = bar(binCenters, categoryCounts', 'stacked');
xlabel('oo(1,:) values');
ylabel('Normalized frequency');
legend(arrayfun(@(c) sprintf('Category %d', c), categories, 'UniformOutput', false), 'Location', 'northwest');
title('Layered Categorical Bar Plot');

bc = binCounts';
be = binEdges';
end
