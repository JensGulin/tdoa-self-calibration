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

% Precompute bin edges
if size(x) < 2
    binEdges = linspace(x-1, x+1, numBins + 1);
else    
    binEdges = linspace(min(x), max(x), numBins + 1);
end
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

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
