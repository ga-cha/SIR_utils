% For some reason this takes wayyy too long

% Main script
% Load patient metadata and explicitly set 'ses' as string
opts = detectImportOptions('/fs04/kg98/trangc/VBM/data/metaVBM_SCZ.csv');
opts = setvartype(opts, 'ses', 'string');
metadata = readtable('/fs04/kg98/trangc/VBM/data/metaVBM_SCZ.csv', opts);
% metadata.diagnosis_string = double(strcmp(metadata.diagnosis_string, 'SCZ')); % 'SCZ' -> 1, 'HC' -> 0

% load atlas
s132_img = niftiread('/fs03/kg98/gchan/Atlases/Tian/Schaefer_Tian/reordered/Schaefer2018_100Parcels_7Networks_order_Tian_Subcortex_S2_MNI152NLin6Asym_1.5mm_reordered.nii.gz');
atlas = double(s132_img);

gmv = get_gmv(metadata, atlas);
lme = run_lme(gmv);
betas = get_beta(lme);
write_beta(betas);


function [gmv] = get_gmv(metadata, atlas)
    % Output dataframe for lme
    gmv = table();
    
    for i = 1:height(metadata)
        patient = Patient(i, metadata);
        rois = patient.get_patient_rois(atlas);
        patient = patient.make_patient_df(rois);
        gmv = [gmv; patient.patient_data]; %#ok<AGROW>
    end

    % gmv_all = gmv;
    gmv = gmv_all;
    gmv = gmv_all(gmv_all.site=='BrainGluSchi', :);
end

function lme = run_lme(gmv)
    gmv.roi = categorical(gmv.roi);
    gmv.sex = categorical(gmv.sex);
    gmv.site = categorical(gmv.site);
    gmv.diagnosis = double(gmv.diagnosis)-1;  % Set 0 (HC) as the baseline


    % lme = fitlme(gmv, 'MGV ~ roi + roi*diagnosis + age + age*age + sex + (1|site)');
    lme = fitlme(gmv, 'MGV ~ roi + roi*diagnosis + age + age*age + sex + (1|site)');
end

function [betas] = get_beta(lme)
    beta = fixedEffects(lme);
    
    % Get the terms corresponding to 'roi' and 'roi*diagnosis'
    terms = lme.CoefficientNames; % Get all term names
    roi_terms = contains(terms, 'roi') & ~contains(terms, 'diagnosis'); % Terms for just 'roi'
    roi_diag_terms = contains(terms, 'roi') & contains(terms, 'diagnosis'); % Terms for 'roi*diagnosis'
    diag_terms = ~contains(terms, 'roi') & contains(terms, 'diagnosis'); % Terms for just 'roi'

    % Extract beta values
    betas.roi = beta(roi_terms); 
    betas.roi_diag = beta(roi_diag_terms); 
    betas.diag = beta(diag_terms); 

    % Add roi_1 = 0
    betas.roi = [0; betas.roi]; 
    betas.roi_diag = [0; betas.roi_diag]; 
    betas.diag = [0; betas.roi_diag]; 

    betas.roi = normalize([betas.roi]); 
    betas.roi_diag = normalize([betas.roi_diag]); 
    betas.diag = normalize([betas.roi_diag]); 
end

function write_beta(betas)
    % Save the betas to a CSV file
    writematrix(betas.roi, 'betas_roi.csv');
    writematrix(betas.roi_diag, 'betas_roi_diag.csv');
    writematrix(betas.diag, 'betas_diag.csv');
    
    disp('Beta values for ROIs:');
    disp(head(betas.roi));
    
    disp('Beta values for ROI*Diagnosis interaction:');
    disp(head(betas.roi_diag));
end

function plot_lme_residuals(lme)
    % Extract residuals from the model
    resid = residuals(lme);
    fitted_values = fitted(lme);
    
    % Plot residuals with transparency
    figure('Color','white');
    scatter(1:length(resid), resid, 'filled', 'MarkerFaceAlpha', 0.005);
    xlabel('Observation'); ylabel('Residuals');

    figure('Color','white');
    scatter(fitted_values, resid, 'filled', 'MarkerFaceAlpha', 0.005);
    xlabel('Fitted Values'); ylabel('Residuals');
end



