n_rois = 41;
sconnLen(42,:) = [];
sconnLen(:,42) = [];
sconnDen(42,:) = [];
sconnDen(:,42) = [];

sconnLen_40 = zeros(41,41);
sconnDen_40 = zeros(41,41);


for roi_i = 1:n_rois
    for roi_j = 1:n_rois
        sconnLen_40(roi_i,roi_j) = sconnLen(order(roi_i),order(roi_j));
        sconnDen_40(roi_i,roi_j) = sconnDen(order(roi_i),order(roi_j));
    end
end
