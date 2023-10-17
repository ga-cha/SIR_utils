function [out, medianDists_out] = interpolateFiberLengths(verts, rois, originalFiberLengthsRaw)
%% interpolateFiberLengths Estimates fiber lengths for tracts between brain regions
% Implements the method described by Ying-Qiu Zheng and colleagues in "Local
% vulnerability and global connectivity jointly shape neurodegenerative disease
% propagation." PLoS biology 17.11 (2019): e3000495
%
%
%% Input Arguments
%  verts - vertex xyz coordinates (V x 3 matrix)
%  where V is the number of vertices
%
%  rois - ROI allocations of each vertex (V x 1 vector)
%
%  originalFiberLengthsRaw - original fibre lengths between vertices (R x R matrix) 
%  where R is the number of ROIs
%
%
%% Output Arguments
%  out - original and interpolated fiber lengths between each pair of ROIs (R x R matrix)
%  where R is the number of ROIs
%  
%  medianDists_out - median distance between each pair of ROIs (R x R matrix)
%
%
%% Authors
% Mehul Gajwani, Monash University, 2023
%
%
%% ENDPUBLISH

%% Format inputs and remove 0 rois
% [verts, ~, rois] = checkVertsFacesRoisData(verts, [], rois);
% [verts, ~, rois] = trimExcludedRois(verts, [], rois, 'overrideAssertions', true);
nrois = length(unique(rois));


%% Calculate medianDists
medianDists_out = nan(nrois);
for y = 1:nrois
    medianDists_out(:,y) = splitapply(...
        @(v) median(pdist2( verts(rois==y,:),v ),'all') , verts , rois);
end
medianDists_out = medianDists_out.*~eye(size(medianDists_out)); % output square
medianDists = squareform(medianDists_out).'; % use vector for calculations


%% Format and vectorise originalFiberLengthsRaw
if ~issymmetric(originalFiberLengthsRaw)
    assert(all(tril(originalFiberLengthsRaw==0), 'all') || ...
        all(triu(originalFiberLengthsRaw==0), 'all'));
    originalFiberLengthsRaw = originalFiberLengthsRaw + originalFiberLengthsRaw.';
end
originalFiberLengths = squareform(originalFiberLengthsRaw).';

presentEdges = logical(originalFiberLengths);


%% Regress and generate output
out = originalFiberLengths;
p = polyfit(medianDists(presentEdges), originalFiberLengths(presentEdges), 1);
out(~presentEdges) = polyval(p, medianDists(~presentEdges));
out = squareform(out);


end