% -------------------------------------------------------------------------------------------------
% GMFP in low vs high sensations with the potentials from all conditions pooled in together
% -------------------------------------------------------------------------------------------------

% Clear the workspace, close all figures, and clear command window
clearvars; close all; clc;

% Define the sense to plot
senseToPlot = 'Twitch'; % Options: 'Discomf', 'Twitch', 'Sound', 'Pain'
windowToExamine = 'late';

% Define the intensity
Int = 'low';
dataPath = fullfile('/Volumes/Mana_HD/GWM', sprintf('GWM_setsRates_%s.mat', Int));
load(dataPath);

% Define conditions and perceptions to examine
conditions = {'DLPFC', 'FEF', 'PP', 'SHAM'};
perceptionsToExamine = {'Discomf', 'Twitch', 'Sound', 'Pain'};

% Find the index of the sense to plot
senseIndex = find(strcmpi(perceptionsToExamine, senseToPlot));

% Initialize cell arrays for storing indices
withSense = cell(length(conditions)-1, length(perceptionsToExamine));
withoutSense = cell(size(withSense));

% Calculate medians for each sense
% Define the sensory data to examine and calculate medians
sensoryData = {Discomf, Twitch, Sound, Pain};
medians = cellfun(@(x) median(reshape(x(:, 1:3), [], 1)), sensoryData);
% Replicate medians for comparison across conditions
medianOfSenses = repmat(medians, length(conditions)-1, 1);

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

%--------------------------------------------------------------------------------------------------- 
% Differences in GMFP
%--------------------------------------------------------------------------------------------------- 
load(['/Volumes/Mana_HD/GWM/Analyzed/avref_FieldtripTEPs_', Int, '_mfp.mat'], 'gmfp_avg', 'gmfp');% Define time windows based on the experimental condition
timeWindows = struct(...
    'early', struct('range', 1016:1060, 'total', 1000:1060, 'offset', 16), ...
    'late', struct('range', 1060:1400, 'total', 1000:1400, 'offset', 60));

selectedTimeWindow = timeWindows.(windowToExamine);

% Initialize variables to store GMFP data
dtCondWt_gmfp_indiv = [];
dtCondWOt_gmfp_indiv = [];

% Extract GMFP data for condition A
for condIdx = 1:length(A_Groups)
    for subjIdx = 1:length(A_Groups{condIdx})
        if ~isempty(gmfp{condIdx, A_Groups{condIdx}(subjIdx)})
            dtCondWt_gmfp_indiv = [dtCondWt_gmfp_indiv; gmfp{condIdx, A_Groups{condIdx}(subjIdx)}.avg(selectedTimeWindow.total)];
        end
    end
end

% Extract GMFP data for condition B
for condIdx = 1:length(B_Groups)
    for subjIdx = 1:length(B_Groups{condIdx})
        if ~isempty(gmfp{condIdx, B_Groups{condIdx}(subjIdx)})
            dtCondWOt_gmfp_indiv = [dtCondWOt_gmfp_indiv; gmfp{condIdx,  B_Groups{condIdx}(subjIdx)}.avg(selectedTimeWindow.total)];
        end
    end
end

%--------------------------------------------------------------------------------------------------- 
% GMFP statistical difference
%--------------------------------------------------------------------------------------------------- 

% Select the appropriate time window based on 'windowToExamine'


% Initialize p-values array
pValues = zeros(1, size(dtCondWt_gmfp_indiv, 2));

% Perform rank-sum test for each time point
for tIdx = 1:size(dtCondWt_gmfp_indiv, 2)
    pValues(tIdx) = ranksum(dtCondWt_gmfp_indiv(:, tIdx), dtCondWOt_gmfp_indiv(:, tIdx));
end

% Adjust p-values for multiple comparisons using False Discovery Rate (FDR)
pFDR = mafdr(pValues, 'BHFDR', true);

% Find significant time points after FDR correction
significantTimePoints = find(pFDR < 0.05);

% Display significant time points
disp('Significant time points after FDR correction:');
disp(significantTimePoints)

%---------------------------------------------------------------------------------------------------
% Visualization of GMFP Data with Confidence Intervals or Standard Error of the Mean
%---------------------------------------------------------------------------------------------------

% Define error shading type
errShadeType = 'CI'; % Options: 'CI' for Confidence Interval, 'SEM' for Standard Error of the Mean

% Calculate Confidence Intervals or Standard Error of the Mean
switch errShadeType
    case 'CI'
        ERR_A = confidence_intervals(dtCondWt_gmfp_indiv, 95);
        ERR_B = confidence_intervals(dtCondWt_gmfp_indiv, 95);
    case 'SEM'
        ERR_A = std(dtCondWt_gmfp_indiv) / sqrt(size(dtCondWt_gmfp_indiv, 1));
        ERR_B = std(dtCondWOt_gmfp_indiv) / sqrt(size(dtCondWOt_gmfp_indiv, 1));
end

% Prepare data for plotting
meanHighCond = mean(dtCondWt_gmfp_indiv, 1);
meanLowCond = mean(dtCondWOt_gmfp_indiv, 1);
timeBase = zeros(1, selectedTimeWindow.range(1) - 1001);

% Plot settings
plotColors = {'k', [0.09,0.61,0.80]}; % Black for high condition, Purple for low condition
lineWidth = 2;
patchSaturationHigh = 0.6;
patchSaturationLow = 0.4;

% Create figure for plotting
figure('color', 'w');
hold on;

% Plot high condition with error shading
shadedErrorBar(selectedTimeWindow.total, meanHighCond, ...
    ERR_A(:, 1)', ...
    'lineprops', {'color', plotColors{1}, 'LineWidth', lineWidth}, ...
    'patchSaturation', patchSaturationHigh);

% Plot low condition with error shading
shadedErrorBar(selectedTimeWindow.total, meanLowCond, ...
     ERR_B(:, 1)', ...
    'lineprops', {'color', plotColors{2}, 'LineWidth', lineWidth}, ...
    'patchSaturation', patchSaturationLow);

% Axis settings based on time window
if strcmp(windowToExamine, 'early')
    xStart = 1000;
    xEnd = 1060;
    yMin = 0.5;
    yMax = 2.5;
    xTickStep = 15;
    yTickStep = 0.5;
elseif strcmp(windowToExamine, 'late')
    xStart = 1000;
    xEnd = 1400;
    yMin = 0.5;
    yMax = 2.5;
    xTickStep = 100;
    yTickStep = 0.5;
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