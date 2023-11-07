%% Analysis of Principal Component Scores for High vs. Low Sensory Perceptions

% Clear the environment and set initial parameters
clearvars; close all; clc; 

% Define the index for the principal component of interest and the intensity level
pcIndex = 1;
intensity = 'low';

% Specify the window and sense to examine
windowToExamine = 'late';
senseToPlot = 'Discomf';  % Options: 'Sound', 'Discomf', 'Twitch', 'Pain'
%---------------------------------------------------------------------------------------------------
% Load the perception rating data for the specified intensity
dataPath = ['/Volumes/Mana_HD/GWM/GWM_setsRates_', intensity, '.mat'];
load(dataPath);
%---------------------------------------------------------------------------------------------------
% Define conditions and perceptions for analysis
conditions = {'DLPFC', 'FEF', 'PP', 'SHAM'};
perceptionsToExamine = {'Discomf', 'Twitch', 'Sound', 'Pain'};

% Find the index of the sense to plot
senseIndex = find(strcmpi(perceptionsToExamine, senseToPlot));

% Define the sensory data to examine and calculate medians
sensoryData = {Discomf, Twitch, Sound, Pain};
medians = cellfun(@(x) median(reshape(x(:, 1:3), [], 1)), sensoryData);

% Replicate medians for comparison across conditions
medianOfSenses = repmat(medians, 3, 1);  % Order: 'DLPFC', 'FEF', 'PP', 'SHAM'

% Categorize data based on medians
withSense = cell(length(conditions) - 1, length(perceptionsToExamine));
withoutSense = cell(size(withSense));

for senseIdx = 1:length(perceptionsToExamine)
    for condIdx = 1:length(conditions) - 1
        currentData = sensoryData{senseIdx};
        medianValue = medianOfSenses(condIdx, senseIdx);

        withSense{condIdx, senseIdx} = find(currentData(:, condIdx) > medianValue);

        if medianValue == 0
            withoutSense{condIdx, senseIdx} = find(currentData(:, condIdx) == medianValue);
        else
            withoutSense{condIdx, senseIdx} = find(currentData(:, condIdx) < medianValue);
        end
    end
end

% Clear temporary variables
clear condIdx

% Grouping data based on low vs high perceptions
% Extracting indices for conditions with higher-than-median sensory perception
A_DLPFC = withSense{1, senseIndex};  % DLPFC condition
A_FEF = withSense{2, senseIndex};    % FEF condition
A_PP = withSense{3, senseIndex};     % PP condition
A_Groups = {A_DLPFC, A_FEF, A_PP}; % Consolidating into a cell array
nA = sum(cellfun(@length, A_Groups)); % Total number of high perception instances

% Extracting indices for conditions with lower-than-median sensory perception
B_DLPFC = withoutSense{1, senseIndex};  % DLPFC condition
B_FEF = withoutSense{2, senseIndex};    % FEF condition
B_PP = withoutSense{3, senseIndex};     % PP condition
B_Groups = {B_DLPFC, B_FEF, B_PP};    % Consolidating into a cell array
nB = sum(cellfun(@length, B_Groups)); % Total number of low perception instances

% ---------------------------------------------------------------------------------------------------
% Analyzing Differences in Principal Components (PCs)
% ---------------------------------------------------------------------------------------------------

% Define the type of error shading for plotting
errShade = 'CI'; % Options: 'SEM' for Standard Error or 'CI' for Confidence Interval

% Load PCA data for each condition and extract the relevant principal component
pcaData = cell(1, 3); % Preallocate cell array for PCA data
conditions = {'DLPFC', 'FEF', 'PP'}; % Define the conditions

% Loop through each condition to load and store PCA data
for i = 1:length(conditions)
    filePath = ['/Volumes/Mana_HD/GWM/Analyzed/avref_', conditions{i}, '_PCA_clusters_', intensity, '_', windowToExamine];
    load(filePath, 'allSubjScores');
    pcaData{i} = allSubjScores(:, :, pcIndex);
end

% Collect PCA data for individuals with high sensory perception
dtCondWt_pca_indiv = [];
for condIdx = 1:length(A_Groups)
    % Concatenate PCA scores for individuals with high perception
    dtCondWt_pca_indiv = [dtCondWt_pca_indiv; pcaData{condIdx}(A_Groups{condIdx}, :)];
end

% Collect PCA data for individuals with low sensory perception
dtCondWOt_pca_indiv = []; 
for condIdx = 1:length(B_Groups)
    % Concatenate PCA scores for individuals with low perception
    dtCondWOt_pca_indiv = [dtCondWOt_pca_indiv; pcaData{condIdx}(B_Groups{condIdx}, :)];
end

%---------------------------------------------------------------------------------------------------
% Statistical Analysis of PCA Differences
%---------------------------------------------------------------------------------------------------

% Define time windows based on the experimental condition
timeWindows = struct(...
    'early', struct('range', 1016:1060, 'total', 1001:1060, 'offset', 16), ...
    'late', struct('range', 1060:1400, 'total', 1001:1400, 'offset', 60), ...
    'posStim', struct('range', 1016:1400, 'total', 1001:1400, 'offset', 16) ...
    );

% Select the appropriate time window based on 'windowToExamine'
selectedTimeWindow = timeWindows.(windowToExamine);

% Calculate the absolute values for PCA data
absPCAHigh = abs(dtCondWt_pca_indiv);
absPCALow = abs(dtCondWOt_pca_indiv);

% Initialize p-values array
pValues = zeros(1, size(absPCAHigh, 2));

% Perform rank-sum test for each time point
for tIdx = 1:size(absPCAHigh, 2)
    pValues(tIdx) = ranksum(absPCAHigh(:, tIdx), absPCALow(:, tIdx));
end

% Adjust p-values for multiple comparisons using False Discovery Rate (FDR)
pFDR = mafdr(pValues, 'BHFDR', true);

% Find significant time points after FDR correction
significantTimePoints = find(pFDR < 0.05) + selectedTimeWindow.offset;

% Display significant time points
disp('Significant time points after FDR correction:');
disp(significantTimePoints)

%---------------------------------------------------------------------------------------------------
% Visualization of PCA Data with Confidence Intervals or Standard Error of the Mean
%---------------------------------------------------------------------------------------------------

% Define error shading type
errShadeType = 'CI'; % Options: 'CI' for Confidence Interval, 'SEM' for Standard Error of the Mean

% Calculate Confidence Intervals or Standard Error of the Mean
switch errShadeType
    case 'CI'
        ERR_A = confidence_intervals(dtCondWt_pca_indiv, 95);
        ERR_B = confidence_intervals(dtCondWOt_pca_indiv, 95);
    case 'SEM'
        ERR_A = std(dtCondWt_pca_indiv) / sqrt(size(dtCondWt_pca_indiv, 1));
        ERR_B = std(dtCondWOt_pca_indiv) / sqrt(size(dtCondWOt_pca_indiv, 1));
end

% Prepare data for plotting
meanHighCond = mean(dtCondWt_pca_indiv, 1);
meanLowCond = mean(dtCondWOt_pca_indiv, 1);
timeBase = zeros(1, selectedTimeWindow.range(1) - 1001);

% Plot settings
plotColors = {'k', 'r'}; % Black for high condition, Purple for low condition
lineWidth = 2;
patchSaturationHigh = 0.6;
patchSaturationLow = 0.4;

% Create figure for plotting
figure('color', 'w');
hold on;

% Plot high condition with error shading
shadedErrorBar(selectedTimeWindow.total, [timeBase, meanHighCond], ...
    [timeBase, ERR_A(:, 1)'], ...
    'lineprops', {'color', plotColors{1}, 'LineWidth', lineWidth}, ...
    'patchSaturation', patchSaturationHigh);

% Plot low condition with error shading
shadedErrorBar(selectedTimeWindow.total, [timeBase, meanLowCond], ...
    [timeBase, ERR_B(:, 1)'], ...
    'lineprops', {'color', plotColors{2}, 'LineWidth', lineWidth}, ...
    'patchSaturation', patchSaturationLow);

% Axis settings based on time window
if strcmp(windowToExamine, 'early')
    xStart = 1000;
    xEnd = 1060;
    yMin = -6;
    yMax = 6;
    xTickStep = 15;
    yTickStep = 3;
elseif strcmp(windowToExamine, 'late')
    xStart = 1000;
    xEnd = 1400;
    yMin = -20;
    yMax = 20;
    xTickStep = 100;
    yTickStep = 10;
end

xlim([xStart, xEnd]);
xticks(xStart:xTickStep:xEnd);
xticklabels(0:xTickStep:(xEnd - xStart));
ylim([yMin, yMax]);
yticks(yMin:yTickStep:yMax);
set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 2.5, 'TickDir', 'out');
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');

% Finalize figure
hold off;
set(gcf, 'color', 'w');

%% -------------------------------------------------------------------------------------------------
% Analysis of Principal Component Scores for low auditory versus low or high somatosensory 
% --------------------------------------------------------------------------------------------------

% Clear the environment and set initial parameters
clearvars; close all; clc; 

% Define the index for the principal component of interest and the intensity level
pcIndex = 1;
intensity = 'low';

% Specify the window and sense to examine
windowToExamine = 'late';
senseToPlot = 'Pain';  % Options: 'Sound', 'Discomf', 'Twitch', 'Pain'
%---------------------------------------------------------------------------------------------------
% Load the perception rating data for the specified intensity
dataPath = ['/Volumes/Mana_HD/GWM/GWM_setsRates_', intensity, '.mat'];
load(dataPath);
%---------------------------------------------------------------------------------------------------

% Define conditions and perceptions for analysis
conditions = {'DLPFC', 'FEF', 'PP', 'SHAM'};
perceptionsToExamine = {'Discomf', 'Twitch', 'Sound', 'Pain'};

% Find the index of the sense to plot
senseIndex = find(strcmpi(perceptionsToExamine, senseToPlot));

% Define the sensory data to examine and calculate medians
sensoryData = {Discomf, Twitch, Sound, Pain};
medians = cellfun(@(x) median(reshape(x(:, 1:3), [], 1)), sensoryData);

% Replicate medians for comparison across conditions
medianOfSenses = repmat(medians, 3, 1);  % Order: 'DLPFC', 'FEF', 'PP', 'SHAM'

% Categorize data based on medians
withSense = cell(length(conditions) - 1, length(perceptionsToExamine));
withoutSense = cell(size(withSense));

for senseIdx = 1:length(perceptionsToExamine)
    for condIdx = 1:length(conditions) - 1
        currentData = sensoryData{senseIdx};
        medianValue = medianOfSenses(condIdx, senseIdx);

        withSense{condIdx, senseIdx} = find(currentData(:, condIdx) > medianValue);

        if medianValue == 0
            withoutSense{condIdx, senseIdx} = find(currentData(:, condIdx) == medianValue);
        else
            withoutSense{condIdx, senseIdx} = find(currentData(:, condIdx) < medianValue);
        end
    end
end

% Clear temporary variables
clear condIdx

% Grouping data based on low auditory and low or high somatosensory 
% Extracting indices for conditions with higher-than-median sensory perception
A_DLPFC = intersect(withoutSense{1,3},withSense{1,senseIndex});  % DLPFC condition
A_FEF = intersect(withoutSense{2,3},withSense{2,senseIndex});    % FEF condition
A_PP = intersect(withoutSense{3,3},withSense{3,senseIndex});   % PP condition
A_Groups = {A_DLPFC, A_FEF, A_PP}; % Consolidating into a cell array
nA = sum(cellfun(@length, A_Groups)); % Total number of high perception instances

% Extracting indices for conditions with lower-than-median sensory perception
B_DLPFC = intersect(withoutSense{1,3},withoutSense{1,senseIndex}); % DLPFC condition
B_FEF = intersect(withoutSense{2,3},withoutSense{2,senseIndex});% FEF condition
B_PP = intersect(withoutSense{3,3},withoutSense{3,senseIndex});
B_Groups = {B_DLPFC, B_FEF, B_PP};    % Consolidating into a cell array
nB = sum(cellfun(@length, B_Groups)); % Total number of low perception instances

% ---------------------------------------------------------------------------------------------------
% Analyzing Differences in Principal Components (PCs)
% ---------------------------------------------------------------------------------------------------

% Define the type of error shading for plotting
errShade = 'CI'; % Options: 'SEM' for Standard Error or 'CI' for Confidence Interval

% Load PCA data for each condition and extract the relevant principal component
pcaData = cell(1, 3); % Preallocate cell array for PCA data
conditions = {'DLPFC', 'FEF', 'PP'}; % Define the conditions

% Loop through each condition to load and store PCA data
for i = 1:length(conditions)
    filePath = ['/Volumes/Mana_HD/GWM/Analyzed/avref_', conditions{i}, '_PCA_clusters_', intensity, '_', windowToExamine];
    load(filePath, 'allSubjScores');
    pcaData{i} = allSubjScores(:, :, pcIndex);
end

% Collect PCA data for individuals with high sensory perception
dtCondWt_pca_indiv = []; 
for condIdx = 1:length(A_Groups)
    % Concatenate PCA scores for individuals with high perception
    dtCondWt_pca_indiv = [dtCondWt_pca_indiv; pcaData{condIdx}(A_Groups{condIdx}, :)];
end

% Collect PCA data for individuals with low sensory perception
dtCondWOt_pca_indiv = []; 
for condIdx = 1:length(B_Groups)
    % Concatenate PCA scores for individuals with low perception
    dtCondWOt_pca_indiv = [dtCondWOt_pca_indiv; pcaData{condIdx}(B_Groups{condIdx}, :)];
end

%---------------------------------------------------------------------------------------------------
% Statistical Analysis of PCA Differences
%---------------------------------------------------------------------------------------------------

% Define time windows based on the experimental condition
timeWindows = struct(...
    'early', struct('range', 1016:1060, 'total', 1001:1060, 'offset', 16), ...
    'late', struct('range', 1060:1400, 'total', 1001:1400, 'offset', 60), ...
    'posStim', struct('range', 1016:1400, 'total', 1001:1400, 'offset', 16) ...
    );

% Select the appropriate time window based on 'windowToExamine'
selectedTimeWindow = timeWindows.(windowToExamine);

% Calculate the absolute values for PCA data
absPCAHigh = abs(dtCondWt_pca_indiv);
absPCALow = abs(dtCondWOt_pca_indiv);

% Initialize p-values array
pValues = zeros(1, size(absPCAHigh, 2));

% Perform rank-sum test for each time point
for tIdx = 1:size(absPCAHigh, 2)
    pValues(tIdx) = ranksum(absPCAHigh(:, tIdx), absPCALow(:, tIdx));
end

% Adjust p-values for multiple comparisons using False Discovery Rate (FDR)
pFDR = mafdr(pValues, 'BHFDR', true);

% Find significant time points after FDR correction
significantTimePoints = find(pFDR < 0.05) + selectedTimeWindow.offset;

% Display significant time points
disp('Significant time points after FDR correction:');
disp(significantTimePoints)

%---------------------------------------------------------------------------------------------------
% Visualization of PCA Data with Confidence Intervals or Standard Error of the Mean
%---------------------------------------------------------------------------------------------------

% Define error shading type
errShadeType = 'CI'; % Options: 'CI' for Confidence Interval, 'SEM' for Standard Error of the Mean

% Calculate Confidence Intervals or Standard Error of the Mean
switch errShadeType
    case 'CI'
        ERR_A = confidence_intervals(dtCondWt_pca_indiv, 95);
        ERR_B = confidence_intervals(dtCondWOt_pca_indiv, 95);
    case 'SEM'
        ERR_A = std(dtCondWt_pca_indiv) / sqrt(size(dtCondWt_pca_indiv, 1));
        ERR_B = std(dtCondWOt_pca_indiv) / sqrt(size(dtCondWOt_pca_indiv, 1));
end

% Prepare data for plotting
meanHighCond = mean(dtCondWt_pca_indiv, 1);
meanLowCond = mean(dtCondWOt_pca_indiv, 1);
timeBase = zeros(1, selectedTimeWindow.range(1) - 1001);

% Plot settings
plotColors = {'k', 'r'}; % Black for high condition, Purple for low condition
lineWidth = 2;
patchSaturationHigh = 0.6;
patchSaturationLow = 0.4;

% Create figure for plotting
figure('color', 'w');
hold on;

% Plot high condition with error shading
shadedErrorBar(selectedTimeWindow.total, [timeBase, meanHighCond], ...
    [timeBase, ERR_A(:, 1)'], ...
    'lineprops', {'color', plotColors{1}, 'LineWidth', lineWidth}, ...
    'patchSaturation', patchSaturationHigh);

% Plot low condition with error shading
shadedErrorBar(selectedTimeWindow.total, [timeBase, meanLowCond], ...
    [timeBase, ERR_B(:, 1)'], ...
    'lineprops', {'color', plotColors{2}, 'LineWidth', lineWidth}, ...
    'patchSaturation', patchSaturationLow);

% Axis settings based on time window
if strcmp(windowToExamine, 'early')
    xStart = 1000;
    xEnd = 1060;
    yMin = -6;
    yMax = 6;
    xTickStep = 15;
    yTickStep = 3;
elseif strcmp(windowToExamine, 'late')
    xStart = 1000;
    xEnd = 1400;
    yMin = -20;
    yMax = 20;
    xTickStep = 100;
    yTickStep = 10;
end

xlim([xStart, xEnd]);
xticks(xStart:xTickStep:xEnd);
xticklabels(0:xTickStep:(xEnd - xStart));
ylim([yMin, yMax]);
yticks(yMin:yTickStep:yMax);
set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 2.5, 'TickDir', 'out');
xlabel('Time (ms)');
ylabel('Amplitude (\muV)');

% Finalize figure
hold off;
set(gcf, 'color', 'w');
