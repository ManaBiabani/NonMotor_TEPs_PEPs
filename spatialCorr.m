clear; close all; clc;
% ## calculates spatial correlations between each condition and refCond at each point of time

Int = 'high';
refCond = 'FEF_TMS';
windLength = 2000;
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
    for t = 1:windLength
        for idx = 1:length(ID)
            % Apply exeptions (remove the IDs without a specific condition)EXCLUDED IDS ARE NOT DEFINED FOR PP CONDITION
            if ~(strcmp(ID{idx},'014')) && ~(strcmp(ID{idx},'018')) && ~(strcmp(ID{idx},'104') && (strcmp(condition{conds},'DLPFC')))...
                    && ~(strcmp(ID{idx},'104') && strcmp(condition{conds},'SHAM')) && ~(strcmp(ID{idx},'078') && strcmp(condition{conds},'SHAM'))...
                    && ~(strcmp(ID{idx},'104') && (strcmp(refCond,'DLPFC'))) && ~(strcmp(ID{idx},'104') && strcmp(refCond,'SHAM'))...
                    && ~(strcmp(ID{idx},'078') && strcmp(refCond,'SHAM'))&& ~(strcmp(ID{idx},'117') && strcmp(condition{conds},'PP'))...
                    && ~(strcmp(ID{idx},'GWM011')&&strcmp(refCond,'DLPFC_ES')) && ~(strcmp(ID{idx},'GWM011') && strcmp(condition{conds},'DLPFC_ES'))
                a = [];
                b = [];
                % Correlation across conditions for Each electrode each subject(at the individualized intervals based on individuals' peaks)
                a = double(squeeze(dataToExamin{strcmp(condition,refCond)}{idx}(:,t)));
                b = double(squeeze(dataToExamin{conds}{idx}(:,t)));
                [EachSubj_spearman{conds}(idx,t), EachSubj_Pval_spearman{conds}(idx,t)] = corr(a,b,'type','Spearman');%excluding 15ms post trigger time
            end
        end

        % Fisher's r to z transformation
        fisherZ_Corr{conds}(:,t) = .5.*log((1+EachSubj_spearman{conds}(:,t))./(1-EachSubj_spearman{conds}(:,t)));
    end

    % Remove the excluded IDs before averaging
    exIds = find(all(fisherZ_Corr{conds}==0,2));
    fisherZ_Corr{conds}(exIds,:) =[];
    clearvars exIds

    % average of z scores
    meanFisherZ_Corr.(condition{conds}) = mean(fisherZ_Corr{conds},1);

    % CI of z scores
    zCI.(condition{conds}) = confidence_intervals(fisherZ_Corr{conds},95);

    % Fisher's z to r tranformation
    rFromZ_Corr.(condition{conds}) = (exp(1).^(2.*fisherZ_Corr{conds})-1)./(exp(1).^(2.*fisherZ_Corr{conds})+1);

    % Fisher's z to r tranformation for average of z scores
    rMeanFisherZ_Corr.(condition{conds})= (exp(1).^(2.*meanFisherZ_Corr.(condition{conds}))-1)./(exp(1).^(2.* meanFisherZ_Corr.(condition{conds}))+1);

    % Fisher's z to r tranformation for CI
    rCIrFromZ.(condition{conds}) = (exp(1).^(2.* zCI.(condition{conds}))-1)./(exp(1).^(2.* zCI.(condition{conds}))+1);

    % SEM of z scores
    zSEM.(condition{conds}) = std(fisherZ_Corr{conds})/sqrt(length(ID));

    % Fisher's z to r tranformation for SEM
    rSEMrFromZ.(condition{conds}) = (exp(1).^(2.* zSEM.(condition{conds}))-1)./(exp(1).^(2.* zSEM.(condition{conds}))+1);
end

% save
save([pathOut, 'SpatialCorr_With_',refCond, '_int_',Int]);