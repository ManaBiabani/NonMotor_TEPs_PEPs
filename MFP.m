%% MFP
clear; close all; clc;

% Calculates the mean field power of EEG data utilizing FieldTrip toolbox functions.
% The input data must be formatted according to FieldTrip's data structure requirements.

Int = 'high';
RefName = 'avref';
study = 'Adelaide';
% --------------------------------------------------------------------------------------------------
if strcmp(study, 'Monash') && strcmp(Int, 'high')
    pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
    load([pathOut, 'Experiment_A_fieldtrip.mat']);
    dataToExamin = ExpA;
    condition = {'DLPFC';'FEF';'PP';'Control'};
elseif strcmp(study, 'Monash') && strcmp(Int, 'low')
    pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
    load([pathOut, 'Experiment_B_fieldtrip.mat']);
    dataToExamin = ExpB;
    condition = {'DLPFC';'FEF';'PP';'Control'};
elseif strcmp(study, 'Adelaide')
    pathOut = '/Volumes/Mana_HD/GWM/Adelaide/Analyzed/';
    load([pathOut, 'Experiment_C_fieldtrip.mat']);
    dataToExamin = ExpC;
    condition = {'DLPFC_TMS';'DLPFC_ES';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'M1_TMS';'M1_ES';'M1_CONTROL';'SHOULDER_TMS'};
end
% --------------------------------------------------------------------------------------------------

% Create time-locked average
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.5 -.01];

% GMFP calculation
for conds = 1:length(condition)
    cfg = [];
    cfg.method = 'amplitude';
    gmfp_avg{conds} = ft_globalmeanfield(cfg, grandAverage.(condition{conds}) );
    % Plot GMFP
    %     figure;
    %     plot(gmfp_avg{conds}.time(900:1400), gmfp_avg{conds}.avg(900:1400),'b');
    %     xlabel('time (s)');
    %     ylabel('GMFP (uv^2)');
    %     xlim([-0.1 0.4]);
    %     ylim([0 3]);

    for idx = 1:length(ID)
        if ~(strcmp(ID{idx},'104') && (strcmp(condition{conds},'DLPFC'))) && ~(strcmp(ID{idx},'104') && strcmp(condition{conds},'SHAM')) && ~(strcmp(ID{idx},'078') && strcmp(condition{conds},'SHAM'))...
                && ~(strcmp(ID{idx},'117') && strcmp(condition{conds},'PP'))   && ~(strcmp(ID{idx},'GWM011') && (strcmp(condition{conds},'DLPFC_ES')))
            cfg = [];
            cfg.method = 'amplitude';
            gmfp{conds,idx} = ft_globalmeanfield(cfg, allData.(condition{conds}){idx});
        end
    end
end

% Load the neighbouring electrodes from the template taken from fieldtrip
neighbour = load('/neighbours.mat');
% LMFP calculation
for conds = 1:length(condition)
    Chans = [];
    if contains(condition{conds}, 'DLPFC')
        centChan = 'F3';
    elseif contains(condition{conds}, 'FEF')
        centChan = 'FC1';
    elseif contains(condition{conds}, 'PP')
        centChan = 'P3';
    elseif contains(condition{conds}, 'M1')
        centChan = 'C3';
    elseif contains(condition{conds}, 'SHAM')
        centChan = 'CZ';
    elseif contains(condition{conds}, 'SHOULDER')
        centChan = 'CZ';
    end

    Chans = [neighbour.neighbours(strcmpi({neighbour.neighbours.label},centChan)).neighblabel; centChan];
    cfg = [];
    cfg.method = 'amplitude';
    cfg.channel = Chans;

    lmfp_avg{conds} = ft_globalmeanfield(cfg, grandAverage.(condition{conds}) );
    % Plot GMFP
    %     figure;
    %     plot(lmfp_avg{conds}.time(900:1400), lmfp_avg{conds}.avg(900:1400),'b');
    %     xlabel('time (s)');
    %     ylabel('LMFP (uv^2)');
    %     xlim([-0.1 0.4]);
    %     ylim([0 3]);
    for idx = 1:length(ID)
        if ~(strcmp(ID{idx},'104') && (strcmp(condition{conds},'DLPFC'))) && ~(strcmp(ID{idx},'104') && strcmp(condition{conds},'SHAM')) && ~(strcmp(ID{idx},'078') && strcmp(condition{conds},'SHAM'))...
                && ~(strcmp(ID{idx},'117') && strcmp(condition{conds},'PP'))   && ~(strcmp(ID{idx},'GWM011') && (strcmp(condition{conds},'DLPFC_ES')))
            cfg = [];
            cfg.method = 'amplitude';
            cfg.channel = Chans;
            lmfp{conds,idx} = ft_globalmeanfield(cfg, allData.(condition{conds}){idx});
        end
    end
end

save([pathOut, Int,'_mfp.mat'],'gmfp_avg', 'gmfp', 'lmfp_avg', 'lmfp');