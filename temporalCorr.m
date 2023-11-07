clear; close all; clc;

% ## This script computes correlations between TEPs of each condition with refCond at each electrode
% for each participant within three distinct time windows:
% - Whole Time : From baseline to 400ms
% - Early: From baseline to 60ms
% - Mid: From the P60 peak to 400ms

study = 'Adelaide';%'Monash';%'Adelaide'
RefName = 'avref';
Int = 'high';
refCond = 'FEF_TMS';
preStim = 1000;
TmsTrig = 15;
PeakT = [60	400];% Adjust according to your specified windows of interest.
% --------------------------------------------------------------------------------------------------
if strcmp(study, 'Monash') && strcmp(Int, 'high')
    pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
    load([pathOut, 'Experiment_A.mat']);
    dataToExamin = ExpA;
    condition = {'DLPFC';'FEF';'PP';'Control'};
elseif strcmp(study, 'Monash') && strcmp(Int, 'low')
    pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
    load([pathOut, 'Experiment_B.mat']);
    dataToExamin = ExpB;
    condition = {'DLPFC';'FEF';'PP';'Control'};
elseif strcmp(study, 'Adelaide')
    pathOut = '/Volumes/Mana_HD/GWM/Adelaide/Analyzed/';
    load([pathOut, 'Experiment_C.mat']);
    dataToExamin = ExpC;
    condition = {'DLPFC_TMS';'DLPFC_ES';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'M1_TMS';'M1_ES';'M1_CONTROL';'SHOULDER_TMS'};
end
% --------------------------------------------------------------------------------------------------
for conds = 1:length(condition)
    
    for idx= 1:length(ID)
        a = [];
        b = [];
        
        % Apply exeptions (remove the IDs without a specific condition)
        if ~(strcmp(ID{idx},'014')) && ~(strcmp(ID{idx},'018')) &&~(strcmp(ID{idx},'104') && (strcmp(condition{conds},'DLPFC'))) && ~(strcmp(ID{idx},'104') && strcmp(condition{conds},'SHAM')) && ~(strcmp(ID{idx},'078') && strcmp(condition{conds},'SHAM'))&& ~(strcmp(ID{idx},'117') && strcmp(condition{conds},'PP'))...
                && ~(strcmp(ID{idx},'104') && (strcmp(refCond,'DLPFC'))) && ~(strcmp(ID{idx},'104') && strcmp(refCond,'SHAM')) && ~(strcmp(ID{idx},'078') && strcmp(refCond,'SHAM')) && ~(strcmp(ID{idx},'117') && strcmp(refCond,'PP'))...
                && ~(strcmp(ID{idx},'GWM011')&&strcmp(refCond,'DLPFC_ES')) && ~(strcmp(ID{idx},'GWM011') && strcmp(condition{conds},'DLPFC_ES'))
            
            % Correlation across conditions for Each electrode each subject(at the individualized intervals based on individuals' peaks)
            a = double(dataToExamin{strcmpi(condition,refCond)}{idx});
            b = double(dataToExamin{conds}{idx});
            t1 = preStim+PeakT(1);
            t2 = preStim+PeakT(2);
            for j = 1:length(eeglabChans)
                [EachSubj_spearman_wholeTime{conds}(j,idx), EachSubj_Pval_spearman_wholeTime{conds}(j,idx)] = corr(a(j,preStim+TmsTrig:preStim+PeakT(end))',b(j,preStim+TmsTrig:preStim+PeakT(end))','type','Spearman');%excluding 15ms post trigger time
                [EachSubj_spearman_early{conds}(j,idx), EachSubj_Pval_spearman_early{conds}(j,idx)] = corr(a(j,preStim+TmsTrig:t1)',b(j,preStim+TmsTrig:t1)','type','Spearman');%excluding 15ms post trigger time
                [EachSubj_spearman_mid{conds}(j,idx), EachSubj_Pval_spearman_mid{conds}(j,idx)] = corr(a(j,t1:t2)',b(j,t1:t2)','type','Spearman');
            end
        end
    end
    
    % Fisher's r to z transformation
    zWholeTime{conds} = .5.*log((1+EachSubj_spearman_wholeTime{conds})./(1-EachSubj_spearman_wholeTime{conds}));
    zEarly{conds} = .5.*log((1+EachSubj_spearman_early{conds})./(1-EachSubj_spearman_early{conds}));
    zMid{conds} = .5.*log((1+EachSubj_spearman_mid{conds})./(1-EachSubj_spearman_mid{conds}));
    
    % Remove the excluded IDs before averaging
    exIds = find(all(zWholeTime{conds}==0,1));
    zWholeTime{conds}(:,exIds) =[];
    zEarly{conds}(:,exIds) =[];
    zMid{conds}(:,exIds) =[];
    
    %  Meansubjects from z
    meanzWholeTime{conds} = mean(zWholeTime{conds},2) ;
    meanzEarly{conds} = mean(zEarly{conds},2) ;
    meanzMid{conds} = mean(zMid{conds},2);
    
    % Transform meansubjects' z back to r for plotting
    for j = 1:length(eeglabChans)
        rmeanzWholeTime{conds}(j,1) = (exp(1)^(2.*meanzWholeTime{conds}(j))-1)/(exp(1)^(2.*meanzWholeTime{conds}(j))+1);
        rmeanzEarly{conds}(j,1) = (exp(1)^(2.*meanzEarly{conds}(j))-1)/(exp(1)^(2.*meanzEarly{conds}(j))+1);
        rmeanzMid{conds}(j,1) = (exp(1)^(2.*meanzMid{conds}(j))-1)/(exp(1)^(2.*meanzMid{conds}(j))+1);
    end
    
    % Reorder subjects and variables for the permutation test (observation x variable)
    zWholeTimeTtest{conds} = (zWholeTime{conds})';
    zEarlyTtest{conds} = (zEarly{conds})';
    zMidTtest{conds} = (zMid{conds})';
    
    % One sample permutaion test
    [pvalWholeTime{conds}] = mult_comp_perm_t1(zWholeTimeTtest{conds},10000,1);
    [pvalEarly{conds}] = mult_comp_perm_t1(zEarlyTtest{conds},10000,1);
    [pvalMid{conds}] = mult_comp_perm_t1(zMidTtest{conds},10000,1);
    
end

% save
save([pathOut,'TempCorr_With_',refCond, '_int_',Int]);