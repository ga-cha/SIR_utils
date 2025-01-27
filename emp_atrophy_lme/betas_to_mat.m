% Load the CSV file
data = readtable('betas.csv');

% Extract the numeric data (excluding the first column)
betas = (data(:, 2:end));

% Save the numeric data to a .mat file
save('c:/Users/gabri/ga-cha/SIR_utils/emp_atrophy_lme/betas.mat', 'betas');

% Create a new table with the same column names as emp_atr and additional columns
new_table = cell2table(C, 'VariableNames', [{'risk gene', 'clearance gene'}, emp_atr.Properties.VariableNames]);

% Save the new table to a .mat file
save('c:/Users/gabri/ga-cha/SIR_utils/emp_atrophy_lme/new_table.mat', 'new_table');