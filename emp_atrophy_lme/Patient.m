classdef Patient
    properties
        id
        dataset
        site
        age
        sex
        diagnosis
        ses
        img
        patient_data
    end

    methods
        function obj = Patient(i, metadata)
            obj.id = metadata.subj_id{i};
            obj.dataset = metadata.dataset{i};
            obj.site = metadata.site_string{i};
            obj.age = metadata.age(i);
            obj.sex = metadata.sex_string{i};
            obj.diagnosis = metadata.diagnosis_string{i};
            obj.ses = metadata.ses(i);

            if ismissing(obj.ses)
                patient_dir = sprintf('/fs04/kg98/trangc/VBM/data/%s/%s/anat/', obj.dataset, obj.id);
                obj.img = fullfile(patient_dir, ['s6mwp1' obj.id '_T1w.nii']);
            else
                patient_dir = sprintf('/fs04/kg98/trangc/VBM/data/%s/%s/%s/anat/', obj.dataset, obj.id, obj.ses);
                obj.img = fullfile(patient_dir, ['s6mwp1' obj.id '_' num2str(obj.ses) '_T1w.nii']);
            end
        end

        function rois = get_patient_rois(obj, atlas, n_parcs)
            if nargin < 3
                n_parcs = 66;
            end

            % Read and process image
            gmv = double(niftiread(obj.img));

            % Filter atlas for valid ROI indices
            mask = (atlas > 0) & (atlas <= n_parcs);
            roi_indices = atlas(mask);
            gray_values = gmv(mask);

            % Group using accumarray
            rois = accumarray(roi_indices, gray_values, [], @(x) {x});
        end

        function obj = make_patient_df(obj, rois)
            % Preallocate vectors
            num_rois = length(rois);
            MGV = zeros(num_rois, 1);
            roi_ids = (1:num_rois)';

            % Compute mean gray matter values
            for r = 1:num_rois
                if ~isempty(rois{r})
                    MGV(r) = mean(rois{r});
                end
            end

            % Create table
            obj.patient_data = table( ...
                MGV, repmat({obj.id}, num_rois, 1), roi_ids, ...
                repmat({obj.diagnosis}, num_rois, 1), ...
                repmat(obj.age, num_rois, 1), ...
                repmat({obj.sex}, num_rois, 1), ...
                repmat({obj.site}, num_rois, 1), ...
                'VariableNames', {'MGV', 'subj_id', 'roi', 'diagnosis', 'age', 'sex', 'site'});
        end
    end
end
