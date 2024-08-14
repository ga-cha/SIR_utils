% Simulate theoretical inter-ROI distance
%
% Gabriella Chan 29/08/23
% gabriella.chan@monash.edu
% Monash University
%
% Additional script to null_ROI_dist.ipynb, using MATLAB's pdist2 function
% Generates pairwise euclidean distance between voxels of ROIs
% the mean is taken into a separate array

n_rois = 100;
sconnLen = ifod_len_100;


% since the rois are of different sizes, we read ROI coordinates
% into a cell array
% initialise the array
roi_coords = cell(n_rois, 1);
% read each file into a new row

reverseStr ='';
for roi = 1:n_rois
    f = strcat('data/ROI_coordinates_100/ROI', num2str(roi), 'coords.csv');
    roi_coords{roi} = readmatrix(f);

    msg = sprintf('Reading %d/%d\n', i, n_rois);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

sc_euc = zeros(n_rois, n_rois);

reverseStr ='';
for roi1 = 1:n_rois
    for roi2 = 1:n_rois
        if roi1 == roi2
            continue
        end
        % pdist2 will return pairwise distances between each voxel,
        % i.e. a 2d array
        dist_arr = pdist2(roi_coords{roi1}, roi_coords{roi2});
        sc_euc(roi1,roi2) = mean(nonzeros(dist_arr(:)));

        msg = sprintf('Calculating distance between %d and %d\n', roi1, roi2);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end

% % reshape euclidean and empirical distances
sc_euc = reshape(sc_euc, [], 1);
sconnLen = reshape(sconnLen, [], 1);
% % Create a mask for empirically connected regions
mask = logical(sconnLen);
% % Apply the mask to both arrays
m_sc_euc = sc_euc(mask);
m_sconnLen = sconnLen(mask);

% Calculate the linear regression coefficients
coefficients = polyfit(m_sc_euc, m_sconnLen, 1);
% slope = coefficients(1); intercept = coefficients(2)

% take empirical connection lengths if possible
% otherwise extrapolate from correlation of euclidean distances
sconnLen_sim = zeros(size(sconnLen));

for i = 1:size(sconnLen,1)
    for j = 1:size(sconnLen,2)
        if (sconnLen(i,j) == 0)
            sconnLen_sim(i,j) = sc_euc(i,j) * coefficients(1) + coefficients(2);
        else 
            sconnLen_sim(i,j) = sconnLen(i,j);
        end
    end
end

% visualisation
% a scatter plot of euclidean and empirical distances
m_sconnLen_sim = sconnLen_sim(mask);
hold on
scatter(m_sc_euc, m_sconnLen); lsline;
scatter(m_sc_euc, m_sconnLen_sim); lsline;
scatter(m_sconnLen_sim, m_sconnLen); lsline;
hold off