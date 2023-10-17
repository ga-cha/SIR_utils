% null_ROI_dist2.m
%
% Gabriella Chan 29/08/23
% gabriella.chan@monash.edu
% Monash University
%
% Updated simulation of ROI distance, following S. Oldham (2022), Sci. Adv.
% https://www.science.org/doi/full/10.1126/sciadv.abm6127
%
% The approximate idea is we just run Dijkstra's through the connectivity
% matrix to estimate distances of unconnected vertices.
%

sconnLen_40 = readtable('data/sconnLen_40.csv');
sconnLen_40 = table2array(sconnLen_40);
%nodeNames = cellstr(arrayfun(@num2str, 1:size(sconnLen_40, 1), 'UniformOutput', false));
sconnLen_G = graph(sconnLen_40, "omitselfloops");
for ii = 1:height(sconnLen_G.Edges)
    for j = 1:width(sconnLen_G.Edges)
        distances = shortestpath(sconnLen_G, ii, j, 'Method', 'positive');
    end
end
A = distances;


% plot(sconnLen_G)