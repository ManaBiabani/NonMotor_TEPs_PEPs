clear; close all; clc;
%% -------------------------------------------------------------------------------------------------
% Step 1 : Epoch, Remove TMS artifact, Down sample
% --------------------------------------------------------------------------------------------------

% IDs of participants to analyse
ID = {'001';'005';'007';'008';'009';'010';'012';'013';'014';'015';'016';'017';'018';'019';'020';'021';'022';'023';'024';'024';'025';...
    '026';'027';'028';'029';'030';'031';'032'; '033'; '034'; '035'; '036'; '037'; '038'; '039'; '040'; '041';'043'; '044'; '045'; '046'; '047'; '048'; '049'; '050';...
    '051'; '053'; '054'; '055'; '056'; '057'; '058'; '059'; '061'; '062'; '063'; '064'; '065'; '067'; '068'; '069'; '071'; '072'; '073'; '074'; '075';...
    '076'; '077'; '078'; '079'; '080'; '081'; '082'; '083'; '084'; '085'; '086'; '087'; '088'; '089'; '090'; '091'; '092'; '093'; '094'; '095'; '096'; '097'; '098'; '099';...
    '101'; '102'; '103'; '104'; '105'; '106'; '107'; '108'; '109'; '110'};

% Filename identifiers
sufix = {'DLPFC';'FEF';'PP';'SHAM'};

% File path where data is stored
pathIn = '/Users/manab/Desktop/preProcessing_GWM/';
pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
addpath(genpath('/Users/manab/Desktop/Functions/eeglab2019_1'))


[ALLEEG, EEG] = eeglab;

for     idx = 1:length(ID)
    
    % Makes a subject folder
    if ~isequal(exist([pathOut,ID{idx,1}], 'dir'),7)
        mkdir(pathOut,ID{idx,1});
    end
    
        % Clear EEG
        EEG = {};
        
        % Clear ALLEEG
        ALLEEG = [];
        
        % Apply exception for the number of stimulation block
        if strcmp(ID{idx,1},'104')
            sufix = {'FEF';'PP'};
        elseif strcmp(ID{idx,1},'078')
            sufix = {'DLPFC';'FEF';'PP'};
        end
        
         for suf = 1:length(sufix)
             
        % File path where data is stored
%         filePath = [pathIn,'sub-GWM',ID{idx,1},'/TMSEEG/GWM',ID{idx,1},'_SP_',sufix{suf},'/GWM',ID{idx,1},'_SP_',sufix{suf}];
         filePath = [pathIn,'GWM',ID{idx,1},'_SP_',sufix{suf},'/GWM',ID{idx,1},'_SP_',sufix{suf}];

        
        % Load Curry files
        EEG = loadcurry( [filePath, '.dap','CurryLocations', 'False']);
        
        % Load channel locations
        EEG = pop_chanedit(EEG, 'lookup','/Users/manab/Desktop/Functions/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
        
        % Epoch around TMS pulse
        EEG = pop_epoch( EEG, {  '128'  }, [-1  1], 'newname', 'ep', 'epochinfo', 'yes');
        
        % Remove baseline
        EEG = pop_rmbase( EEG, [-500  -10]);
        
        % Remove unused channels
        EEG = pop_select( EEG,'nochannel',{'31' '32' 'Trigger'});
        
        % Save the original EEG locations for use in interpolation later
        EEG.allchan = EEG.chanlocs;
        
        % Remove TMS artifact
        EEG = pop_tesa_removedata( EEG, [-2 15] );%[-1 10]
        
        % Interpolate missing data
        EEG = pop_tesa_interpdata( EEG, 'cubic', [1 1] );
        
        % Downsample data
        EEG = pop_resample( EEG, 1000);
        
        % Label the events
        for i = 1 : EEG.trials
            EEG.event(i).type = sufix{suf};
        end
        
        % Save EEG for each block type for each participant
        [ALLEEG, EEG] = eeg_store(ALLEEG,EEG);
        
    end;
    
    % Merge all blocks for each subject
    EEG = pop_mergeset(  ALLEEG, 1:length(sufix), 0);
    
    % save data
    EEG = pop_saveset( EEG, 'filepath',[pathOut,ID{idx,1},'/'],'filename', [ID{idx,1} '_all_blocks_ds']);
    
end
%% -------------------------------------------------------------------------------------------------
% STEP 2: REMOVE BAD TRIALS, CHECK BAD ELECTRODES 
% --------------------------------------------------------------------------------------------------
clear; close all; clc;

% select the dataset
dataFrom = 'GWM';%'GWM','Adelaide','GWM_ERPs'

Int = 'high';

if strcmp(dataFrom,'GWM')
    if strcmp(Int,'low')
        ID = {'001';'005';'007';'008';'010';'012';'013';'014';'016';'018';'019';'020';'021';'022';'023';'024';'025';...
            '026';'027';'028';'029';'030';'031';'032';'033'; '034'; '037'; '038';'039';'040';'041';'043'; '044'; '045'; '046'; '047'; '048'; '050';...
            '051'; '053'; '054'; '055'; '056'; '057'; '058'; '059'; '061'; '062'; '063';'065'; '067'; '068'; '071'; '072'; '073'; '074';...
            '076';'077';'078';'079';'080';'081'; '082'; '083';'084';'085'; '086'; '087';'088'; '091'; '092'; '093'; '094'; '095'; '096'; '097'; '098'; '099';...
            '101'; '102'; '103'; '104';'105'; '106'; '107'; '108'; '109'; '110';'111';'112';'113';'114';'116';'0117'};% 117_low = 0117
    elseif strcmp(Int,'high')
        ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    end % 117 and 120 excluded
    pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
    
elseif strcmp(dataFrom,'Adelaide')
    ID = {'GWM001';'GWM002';'GWM003';'GWM004';'GWM005';'GWM006';'GWM007';'GWM008';'GWM009';'GWM010';'GWM011';'GWM012'};
    pathOut = '/Volumes/Mana_HD/GWM/Adelaide/Analyzed/';
    condition = {'DLPFC_TMS';'DLPFC_ES';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'SHOULDER_TMS'};

elseif strcmp(dataFrom,'GWM_ERPs')
    ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    pathOut = '/Volumes/Mana_HD/GWM/WM_ERPs/';
end

addpath(genpath('/Users/manab/Desktop/Functions/eeglab2021.0/'));
eeglab;

%% SEPERATE EACH BLOCK
for idx = 1:length(ID)
    
    % Load point
    for conds =  1:length(condition)
        
        % Clear EEG
        EEG = {};
        
        % Clear ALLEEG
        ALLEEG = [];
        
        if ~(strcmp(ID{idx}, 'GWM011') && strcmp(condition{conds}, 'DLPFC_ES'))
            % load the concatenated dataset
            EEG = pop_loadset('filepath',[pathOut,ID{idx,1},'/'],'filename', [ID{idx,1} '_all_blocks_ds.set']);
            
            % Extract the data from each condition
            EEG = pop_selectevent( EEG, 'type',condition{conds},'deleteevents','on','deleteepochs','on','invertepochs','off');
            
            % Save each condition seperately
            EEG = pop_saveset(EEG, 'filename', [ID{idx,1},'_', condition{conds},'_ds'],'filepath', [pathOut ID{idx,1}]);
        end
    end
end

%% REMOVE BAD TRIALS, CHECK BAD ELECTRODES
for idx = 1:length(ID)
    
    for conds = 1:length(condition)
        if ~(strcmp(ID{idx}, 'GWM011') && strcmp(condition{conds}, 'DLPFC_ES'))
            % Clear EEG
            EEG = {};
            
            % Clear ALLEEG
            ALLEEG = [];
            
            % Load point
            EEG = pop_loadset('filepath',[pathOut,ID{idx,1},'/'],'filename', [ID{idx,1},'_', condition{conds},'_ds.set']);
            
            EEG = pop_select( EEG,'nochannel',{'31' '32'});
            
            % Find bad trials and bad channels
            EEG = pop_rejkurt(EEG,1,[1:EEG.nbchan] ,5,5,2,0);
            pop_rejmenu( EEG, 1);
            R1 = input('Highlight bad trials, update marks and then press enter');
            EEG.BadTr = unique([find(EEG.reject.rejkurt==1) find(EEG.reject.rejmanual==1)]);
            
            % Reject bad trials
            EEG = pop_rejepoch(EEG,EEG.BadTr,0);
            
            % Remove bad channels
            answer = inputdlg('Enter bad channels', 'Bad channel removal', [1 50]);
            str = answer{1};
            EEG.badChan = strsplit(str);
            close all;
            EEG = pop_select( EEG,'nochannel',EEG.badChan);
            EEG = pop_saveset( EEG, 'filename', [ID{idx,1},'_', condition{conds},'_ds_reject'], 'filepath', [pathOut ID{idx,1}]);
            
        end
    end
end

%% REMOVE DECAY FROM ES CONDITION 

esCondition = {'DLPFC_ES';'FEF_ES';'PP_ES';'M1_ES'};

if strcmp(pathOut,'/Volumes/Mana_HD/GWM/Adelaide/Analyzed/')
    for idx = 1:length(ID)
        for conds = 1:length(esCondition)
            
            if ~(strcmp(ID{idx}, 'GWM011') && strcmp(esCondition{conds}, 'DLPFC_ES'))
                % Clear EEG
                EEG = {};
                
                % Clear ALLEEG
                ALLEEG = [];
                
                % Load point
                EEG = pop_loadset('filepath',[pathOut,ID{idx,1},'/'],'filename', [ID{idx,1},'_', Condition{conds},'_ds_reject.set']);
                
                % remove decay artefact  ( method used in Conde et al 2019)
                decayTime = [1005:1500];
                for j = 1:size(EEG.data, 1)
                    
                    for tr = 1:size(EEG.data, 3)
                        dt = squeeze(EEG.data(j,decayTime,tr));
                        
                        x = 1:length(dt);
                        % exponential model
                        h = fit(x', dt', 'exp2');
                        yhat = h.a*exp(h.b*x) + h.c*exp(h.d*x);
                        
                        % remove exponential from raw data
                        EEG.data(j,decayTime,tr) = dt-yhat;
                    end
                end
            end
        end
    end
    EEG = pop_saveset( EEG, 'filename', [ID{idx,1},'_', esCondition{conds},'_ds_reject'], 'filepath', [pathOut ID{idx,1}]);
end
    
%% -------------------------------------------------------------------------------------------------
% STEP 3: Remove TMS pulse artifact and run FASTICA round 1
% --------------------------------------------------------------------------------------------------
clear; close all; clc;

% select the dataset
dataFrom = 'Adelaide';%'GWM','Adelaide','GWM_ERPs'

Int = 'high';

if strcmp(dataFrom,'GWM')
    if strcmp(Int,'low')
        ID = {'001';'005';'007';'008';'010';'012';'013';'014';'016';'018';'019';'020';'021';'022';'023';'024';'025';...
            '026';'027';'028';'029';'030';'031';'032';'033'; '034'; '037'; '038';'039';'040';'041';'043'; '044'; '045'; '046'; '047'; '048'; '050';...
            '051'; '053'; '054'; '055'; '056'; '057'; '058'; '059'; '061'; '062'; '063';'065'; '067'; '068'; '071'; '072'; '073'; '074';...
            '076';'077';'078';'079';'080';'081'; '082'; '083';'084';'085'; '086'; '087';'088'; '091'; '092'; '093'; '094'; '095'; '096'; '097'; '098'; '099';...
            '101'; '102'; '103'; '104';'105'; '106'; '107'; '108'; '109'; '110';'111';'112';'113';'114';'116';'0117'};% 117_low = 0117
    elseif strcmp(Int,'high')
        ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    end % 117 and 120 excluded
    pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
    
elseif strcmp(dataFrom,'Adelaide')
    ID = {'GWM001';'GWM002';'GWM003';'GWM004';'GWM005';'GWM006';'GWM007';'GWM008';'GWM009';'GWM010';'GWM011';'GWM012'};
    pathOut = '/Volumes/Mana_HD/GWM/Adelaide/Analyzed/';
    condition = {'DLPFC_TMS';'DLPFC_ES';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'M1_TMS';'M1_ES';'M1_CONTROL';'SHOULDER_TMS'};

elseif strcmp(dataFrom,'GWM_ERPs')
    ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    pathOut = '/Volumes/Mana_HD/GWM/WM_ERPs/';
end

addpath(genpath('/Users/manab/Desktop/Functions/eeglab2021.0/'));
eeglab;

for idx = 1:length(ID)
    for conds = 1:length(condition)
        if ~(strcmp(ID{idx}, 'GWM011') && strcmp(condition{conds}, 'DLPFC_ES'))
            % Load data
            EEG = pop_loadset('filepath',[pathOut,ID{idx,1},'/'],'filename', [ID{idx,1},'_', condition{conds},'_ds_reject.set']);
            
            % Remove TMS pulse artifact
            EEG = pop_tesa_removedata( EEG, [-2 15] );
            
            % Run FastICA (round 1)
            EEG = pop_tesa_fastica( EEG,'approach', 'symm', 'g', 'tanh', 'stabilization', 'off');
            
            % Save point
            EEG = pop_saveset( EEG, 'filename', [ID{idx,1} ,'_', condition{conds},'_ds_reject_ICA1'], 'filepath', [pathOut ID{idx,1}]);
        end
    end
end
%% -------------------------------------------------------------------------------------------------
% STEP 4: Remove TMS-evoked muscle / decay 
% --------------------------------------------------------------------------------------------------
clear; close all; clc;

% select the dataset
dataFrom = 'Adelaide';%'GWM','Adelaide','GWM_ERPs'

Int = 'high';

if strcmp(dataFrom,'GWM')
    if strcmp(Int,'low')
        ID = {'001';'005';'007';'008';'010';'012';'013';'014';'016';'018';'019';'020';'021';'022';'023';'024';'025';...
            '026';'027';'028';'029';'030';'031';'032';'033'; '034'; '037'; '038';'039';'040';'041';'043'; '044'; '045'; '046'; '047'; '048'; '050';...
            '051'; '053'; '054'; '055'; '056'; '057'; '058'; '059'; '061'; '062'; '063';'065'; '067'; '068'; '071'; '072'; '073'; '074';...
            '076';'077';'078';'079';'080';'081'; '082'; '083';'084';'085'; '086'; '087';'088'; '091'; '092'; '093'; '094'; '095'; '096'; '097'; '098'; '099';...
            '101'; '102'; '103'; '104';'105'; '106'; '107'; '108'; '109'; '110';'111';'112';'113';'114';'116';'0117'};% 117_low = 0117
    elseif strcmp(Int,'high')
        ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    end % 117 and 120 excluded
    pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
    
elseif strcmp(dataFrom,'Adelaide')
    ID = {'GWM001';'GWM002';'GWM003';'GWM004';'GWM005';'GWM006';'GWM007';'GWM008';'GWM009';'GWM010';'GWM011';'GWM012'};
    pathOut = '/Volumes/Mana_HD/GWM/Adelaide/Analyzed/';
    condition = {'DLPFC_TMS';'DLPFC_ES';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'M1_TMS';'M1_ES';'M1_CONTROL';'SHOULDER_TMS'};
elseif strcmp(dataFrom,'GWM_ERPs')
    ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    pathOut = '/Volumes/Mana_HD/GWM/WM_ERPs/';
end

addpath(genpath('/Users/manab/Desktop/Functions/eeglab2021.0/'));
eeglab;


for idx = 1:length(ID)
    for conds = 1:length(condition)
        if ~(strcmp(ID{idx}, 'GWM011') && strcmp(condition{conds}, 'DLPFC_ES'))
            % Load data
            EEG = pop_loadset('filepath',[pathOut,ID{idx,1},'/'],'filename', [ID{idx,1},'_', condition{conds},'_ds_reject_ICA1.set']);
            
            % Remove
            EEG = pop_tesa_compselect( EEG,'comps',10,'figSize','small','plotTimeX',[-200 500],'plotFreqX',[1 100],...
                'tmsMuscle','on','tmsMuscleThresh',8,'tmsMuscleWin',[16 30],'tmsMuscleFeedback','off','blink','off',...
                'blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},'blinkFeedback','off','move','off','moveThresh',2,...
                'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','off','muscleThresh',0.6,'muscleFreqWin',[30 100],'muscleFeedback','off','elecNoise','off','elecNoiseThresh',4,'elecNoiseFeedback','off','compCheck','off' );
            
            % Interpolate removed data
            EEG = pop_tesa_interpdata( EEG, 'linear' );
            
            % Bandpass (1-100 Hz) and bandstop (48-52 Hz) filter data
            EEG = pop_tesa_filtbutter( EEG, 1, 100, 4, 'bandpass' );
            EEG = pop_tesa_filtbutter( EEG, 48, 52, 4, 'bandstop' );
            
            % Save point
            EEG = pop_saveset( EEG, 'filename', [ID{idx,1},'_', condition{conds}, '_ds_reject_ICA1_clean'], 'filepath', [pathOut ID{idx,1}]);
        end
    end
    
end

%% -------------------------------------------------------------------------------------------------
% STEP 5: Run FASTICA round 2
% --------------------------------------------------------------------------------------------------
clear; close all; clc;

% select the dataset
dataFrom = 'Adelaide';%'GWM','Adelaide','GWM_ERPs'

Int = 'high';

if strcmp(dataFrom,'GWM')
    if strcmp(Int,'low')
        ID = {'001';'005';'007';'008';'010';'012';'013';'014';'016';'018';'019';'020';'021';'022';'023';'024';'025';...
            '026';'027';'028';'029';'030';'031';'032';'033'; '034'; '037'; '038';'039';'040';'041';'043'; '044'; '045'; '046'; '047'; '048'; '050';...
            '051'; '053'; '054'; '055'; '056'; '057'; '058'; '059'; '061'; '062'; '063';'065'; '067'; '068'; '071'; '072'; '073'; '074';...
            '076';'077';'078';'079';'080';'081'; '082'; '083';'084';'085'; '086'; '087';'088'; '091'; '092'; '093'; '094'; '095'; '096'; '097'; '098'; '099';...
            '101'; '102'; '103'; '104';'105'; '106'; '107'; '108'; '109'; '110';'111';'112';'113';'114';'116';'0117'};% 117_low = 0117
    elseif strcmp(Int,'high')
        ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    end % 117 and 120 excluded
    pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
    
elseif strcmp(dataFrom,'Adelaide')
    ID = {'GWM001';'GWM002';'GWM003';'GWM004';'GWM005';'GWM006';'GWM007';'GWM008';'GWM009';'GWM010';'GWM011';'GWM012'};
    pathOut = '/Volumes/Mana_HD/GWM/Adelaide/Analyzed/';
    condition = {'DLPFC_TMS';'DLPFC_ES';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'M1_TMS';'M1_ES';'M1_CONTROL';'SHOULDER_TMS'};

elseif strcmp(dataFrom,'GWM_ERPs')
    ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    pathOut = '/Volumes/Mana_HD/GWM/WM_ERPs/';
end

addpath(genpath('/Users/manab/Desktop/Functions/eeglab2021.0/'));
eeglab;


for idx = 1:length(ID)
    for conds = 1:length(condition)
        if ~(strcmp(ID{idx}, 'GWM011') && strcmp(condition{conds}, 'DLPFC_ES'))
            
            % Load data
            EEG = pop_loadset('filepath',[pathOut,ID{idx,1},'/'],'filename', [ID{idx,1} ,'_', condition{conds}, '_ds_reject_ICA1_clean.set']);
            
            % Run FastICA (round 2)
            EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' );
            
            % Save point
            EEG = pop_saveset( EEG, 'filename', [ID{idx,1} ,'_', condition{conds}, '_ds_reject_ICA1_clean_ICA2'], 'filepath', [pathOut ID{idx,1}]);
        end
    end
    
end
%% -------------------------------------------------------------------------------------------------
% STEP 6: Remove all other artifacts, Interpolate, Avref
% --------------------------------------------------------------------------------------------------
clear; close all; clc;

% select the dataset
dataFrom = 'Adelaide';%'GWM','Adelaide','GWM_ERPs'

Int = 'high';

if strcmp(dataFrom,'GWM')
    if strcmp(Int,'low')
        ID = {'001';'005';'007';'008';'010';'012';'013';'014';'016';'018';'019';'020';'021';'022';'023';'024';'025';...
            '026';'027';'028';'029';'030';'031';'032';'033'; '034'; '037'; '038';'039';'040';'041';'043'; '044'; '045'; '046'; '047'; '048'; '050';...
            '051'; '053'; '054'; '055'; '056'; '057'; '058'; '059'; '061'; '062'; '063';'065'; '067'; '068'; '071'; '072'; '073'; '074';...
            '076';'077';'078';'079';'080';'081'; '082'; '083';'084';'085'; '086'; '087';'088'; '091'; '092'; '093'; '094'; '095'; '096'; '097'; '098'; '099';...
            '101'; '102'; '103'; '104';'105'; '106'; '107'; '108'; '109'; '110';'111';'112';'113';'114';'116';'0117'};% 117_low = 0117
    elseif strcmp(Int,'high')
        ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    end % 117 and 120 excluded
    pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';
    
elseif strcmp(dataFrom,'Adelaide')
    ID = {'GWM001';'GWM002';'GWM003';'GWM004';'GWM005';'GWM006';'GWM007';'GWM008';'GWM009';'GWM010';'GWM011';'GWM012'};
    pathOut = '/Volumes/Mana_HD/GWM/Adelaide/Analyzed/';
    condition = {'DLPFC_TMS';'DLPFC_ES';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'M1_TMS';'M1_ES';'M1_CONTROL';'SHOULDER_TMS'};

elseif strcmp(dataFrom,'GWM_ERPs')
    ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
    pathOut = '/Volumes/Mana_HD/GWM/WM_ERPs/';
end

addpath(genpath('/Users/manab/Desktop/Functions/eeglab2021.0/'));
eeglab;

for idx = 1:length(ID)
    for conds = 1:length(condition)
        if ~(strcmp(ID{idx}, 'GWM011') && strcmp(condition{conds}, 'DLPFC_ES'))
            
            % Load data
            EEG = pop_loadset('filepath',[pathOut,ID{idx,1},'/'],'filename', [ID{idx,1} ,'_', condition{conds}, '_ds_reject_ICA1_clean_ICA2.set']);
            
            % Run FastICA (round 2)
            EEG = pop_tesa_compselect( EEG,'comps',[],'figSize','small','plotTimeX',[-200 500],'plotFreqX',[1 100],'tmsMuscle',...
                'on','tmsMuscleThresh',8,'tmsMuscleWin',[11 30],'tmsMuscleFeedback','off','blink','on','blinkThresh',2.5,'blinkElecs',{'Fp1','Fp2'},...
                'blinkFeedback','off','move','on','moveThresh',2,'moveElecs',{'F7','F8'},'moveFeedback','off','muscle','on','muscleThresh',-0.3,...
                'muscleFreqWin',[30 100],'muscleFreqIn',[7 75],'muscleFeedback','off','elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
            
            EEG.allchan(find(strcmpi({EEG.allchan.labels},'31'))) = [];
            EEG.allchan(find(strcmpi({EEG.allchan.labels},'32'))) = [];
        
            % Interpolate missing channels
            EEG = pop_interp(EEG, EEG.allchan, 'spherical');
            
            % Save point
            EEG = pop_saveset( EEG, 'filename', [ID{idx,1} ,'_', condition{conds}, '_ds_reject_ICA1_clean_ICA2_clean'], 'filepath', [pathOut ID{idx,1}]);
            
            % Reference each condition's data to common average and save
            EEGav = pop_reref(EEG, []);
            EEGav = pop_saveset(EEGav, 'filename', [ID{idx,1},'_', condition{conds},'_avref_FINAL'],'filepath', [pathOut ID{idx,1}]);
            
        end
    end
end
%% -------------------------------------------------------------------------------------------------
% STEP 7: Averaging
% --------------------------------------------------------------------------------------------------
clear; close all; clc;


% define the re-referencing method
RefName = 'avref';

% select the dataset
study = 'Monash';%'Monash','Adelaide','GWM_ERPs'

Int = 'high';

if strcmp(study,'Monash')
if strcmp(Int,'low')
ID = {'001';'005';'007';'008';'010';'012';'013';'014';'016';'018';'019';'020';'021';'022';'023';'024';'025';...
    '026';'027';'028';'029';'030';'031';'032';'033'; '034'; '037'; '038';'039';'040';'041';'043'; '044'; '045'; '046'; '047'; '048'; '050';...
    '051'; '053'; '054'; '055'; '056'; '057'; '058'; '059'; '061'; '062'; '063';'065'; '067'; '068'; '071'; '072'; '073'; '074';...
    '076';'077';'078';'079';'080';'081'; '082'; '083';'084';'085'; '086'; '087';'088'; '091'; '092'; '093'; '094'; '095'; '096'; '097'; '098'; '099';...
    '101'; '102'; '103'; '104';'105'; '106'; '107'; '108'; '109'; '110';'111';'112';'113';'114';'116';'0117'};% 117_low = 0117
condition = {'DLPFC';'FEF';'PP';'SHAM'};
elseif strcmp(Int,'high')
ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
condition = {'DLPFC';'FEF';'PP';'SHAM'};
end % 117 and 120 excluded
pathOut = '/Volumes/Mana_HD/GWM/Analyzed/';

elseif strcmp(study,'Adelaide')
ID = {'GWM001';'GWM002';'GWM003';'GWM004';'GWM005';'GWM006';'GWM007';'GWM008';'GWM009';'GWM010';'GWM011';'GWM012'};
pathOut = '/Volumes/Mana_HD/GWM/Adelaide/Analyzed/';
condition = {'DLPFC_TMS';'DLPFC_ES';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'M1_TMS';'M1_ES';'M1_CONTROL';'SHOULDER_TMS'};

elseif strcmp(study,'GWM_ERPs')
ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
pathOut = '/Volumes/Mana_HD/GWM/WM_ERPs/';
end
% 
% addpath(genpath('/Users/manab/Desktop/Functions/eeglab2021.0/'));
% eeglab;


% Define the window of time that you want to analyse your data in the following steps
WinOfInt = [901:1400];
%--------------------------------------------------------------------------
% Calculate the mean trials for each subject/each electrode and put all the results in one big structure
for idx = 1:length(ID)
    filePath = [pathOut,ID{idx,1},'/'];
    
    % Define conditions
    if strcmp(pathOut,'/Volumes/Mana_HD/GWM/Analyzed/')
        if strcmp(ID{idx,1},'104')
            conditionB = {'FEF';'PP'};
        elseif strcmp(ID{idx,1},'078')
            conditionB = {'DLPFC';'FEF';'PP'};
        elseif strcmp(ID{idx,1},'117')
            conditionB = {'DLPFC';'FEF';'SHAM'};
        else
            conditionB = {'DLPFC';'FEF';'PP';'SHAM'};
        end
    end
    
    if strcmp(pathOut,'/Volumes/Mana_HD/GWM/Adelaide/Analyzed/')
        if strcmp(ID{idx,1},'GWM011')
            conditionB = {'DLPFC_TMS';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'M1_TMS';'M1_ES';'M1_CONTROL';'SHOULDER_TMS'};
        else
            conditionB = {'DLPFC_TMS';'DLPFC_ES';'DLPFC_CONTROL';'FEF_TMS';'FEF_ES';'FEF_CONTROL';'PP_TMS';'PP_ES';'PP_CONTROL';'M1_TMS';'M1_ES';'M1_CONTROL';'SHOULDER_TMS'};
        end
    end
    EEG = [];
    
    cd([filePath])
    for conds = 1:length(conditionB)
        
        % Put similar conditions in the same column
        condNum = find(strcmp(conditionB{conds},condition));
        
        EEG{condNum} = pop_loadset([filePath, ID{idx,1},'_',conditionB{conds}, '_',RefName,'_FINAL.set']);
        if size(EEG{condNum}.data,1) == 64
            EEG{condNum} = pop_select( EEG{condNum},'nochannel',{'31' '32'});
        end
        
        meanTrials{condNum}(idx) = {mean(EEG{condNum}.data,3)};
        mt = cell2mat(meanTrials{condNum}(idx));
        meanTrials_WinOfInt{condNum}(:,idx,:)= mt(:,WinOfInt);
        
    end
    
end

% Calculate meansubj for each condition
for conds = 1:length(condition)
    ct = [];
    for t = 1:2000
        for idx = 1:length(ID)
            if ~(strcmp(ID{idx},'104') && (strcmp(condition{conds},'DLPFC'))) && ~(strcmp(ID{idx},'104') && strcmp(condition{conds},'SHAM')) && ~(strcmp(ID{idx},'078') && strcmp(condition{conds},'SHAM'))...
                    && ~(strcmp(ID{idx},'117') && strcmp(condition{conds},'PP'))   && ~(strcmp(ID{idx},'GWM011') && (strcmp(condition{conds},'DLPFC_ES')))
                c = cell2mat(meanTrials{conds}(idx));
                ct(:,idx) = c(:,t);
            end
        end
        meanSubject{conds}(:,t) = mean(ct,2);
        sdSubject{conds}(:,t) = std(ct,0,2);
    end
    
    % exctract data within the specified window
    meanSubject_WinOfInt{conds} = meanSubject{conds}(:,WinOfInt);
    sdSubject_WinOfInt{conds} = sdSubject{conds}(:,WinOfInt);
end

% Define the channel order/names in EEGLab/your eeg data
eeglabChans = {EEG{1}.chanlocs.labels};

save([pathOut, RefName,'_TEPs_',Int,'.mat']);