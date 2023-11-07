%% =========================== E-Field Simulation =========================

clc; close all;clear all;
% Add the SimNIBS MATLAB functions to the MATLAB path

addpath(genpath('/Users/manab/Applications/SimNIBS-3.2/matlab'));
pathIn = '/Volumes/Mana_HD/GWM/scans/simnibs/Raw/';
pathOut = '/Volumes/Mana_HD/GWM/scans/simnibs/Analyzed/';

% Path for participants/coil info
pathInfo = '/Volumes/Mana_HD/GWM/';

site = {'DLPFC';'FEF';'PP'};
ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
for idx = 1:size(ID,1)
    
    for st = 1:size(site,1)
        % --------------------------------------------------------------------------------------------------
        % Start a SESSION and select a head mesh. Each session has one head mesh and output folder(one subject)
        % --------------------------------------------------------------------------------------------------
        % Initialize a session
        s = sim_struct('SESSION');
         s.map_to_fsavg = true;
        s.map_to_MNI = true;

        % Name of head mesh
        s.fnamehead = [pathIn ID{idx} '/' ID{idx} '.msh'];
        dd = [pathOut,ID{idx},'/',site{st}];
        % Makes a subject output folder
        if ~isequal(exist([pathOut,ID{idx},'/',site{st}], 'dir'),7)
            mkdir([pathOut,ID{idx},'/',site{st}]);
        else
            sprintf('Deleting the existing directory:\n%s', [pathOut,ID{idx},'/',site{st}])
            rmdir([pathOut,ID{idx},'/',site{st},'/'],'s');
            mkdir([pathOut,ID{idx},'/',site{st}]);
        end
        
        s.pathfem = [pathOut, ID{idx},'/',site{st},'/'];
        
        % --------------------------------------------------------------------------------------------------
        % Fields to be output
        % --------------------------------------------------------------------------------------------------
        % Any combination of :
        % v for voltage (in V)
        % e/E for electric field norm/vector (in V/m)
        % j/J for current density norm/vector (int A/m2).
        % e.g., 'eEjJ'; eE(default)
         s.fields = 'eE';
        
        % --------------------------------------------------------------------------------------------------
        % Setting up a TMS simulation
        % --------------------------------------------------------------------------------------------------
        
        % Initialize a list of TMS simulations
        s.poslist{1} = sim_struct('TMSLIST');
        % Select coil
        s.poslist{1}.fnamecoil = '/Users/manab/Applications/SimNIBS-3.2/simnibs_env/lib/python3.7/site-packages/simnibs/ccd-files/Deng_Brain_Stimul_2013/No29_MagVenture_C-B60_Fig8.nii.gz';
        
        % Set a position for your coil
        s.poslist{1}.pos;
        % Center the coil using the coordinates or electrode
        T = readtable([pathInfo 'participant_info_deidentified_MB.xlsx'], 'Sheet','coilCoordinates');
        targetRow = find(strcmp(T.ID, ID{idx}));
        
        centre = [T.([site{st},'x'])(targetRow) T.([site{st},'y'])(targetRow) T.([site{st},'z'])(targetRow)];% or e.g 'F3'
        coil.centre = mni2subject_coords(centre, fullfile([pathIn ID{idx}], ['m2m_', ID{idx}]));
        
        % Coil centre (Simnibs projects coil center to scalp. Therfore the coordinate is a scalp position)
        s.poslist{1}.pos(1).centre = coil.centre;
        
        % Inverse the side of z axis in brainsight
        coilAngle = [T.([site{st},'Coilxa'])(targetRow) T.([site{st},'Coilya'])(targetRow) -T.([site{st},'Coilza'])(targetRow);...
            T.([site{st},'Coilxb'])(targetRow) T.([site{st},'Coilyb'])(targetRow) -T.([site{st},'Coilzb'])(targetRow);...
            T.([site{st},'Coilxc'])(targetRow) T.([site{st},'Coilyc'])(targetRow) -T.([site{st},'Coilzc'])(targetRow)];
        
        % Set the transformation matrix (Simnibs will get the coil centre and coil to scalp distance from this matrix)
        s.poslist{1}.pos.matsimnibs = [...
            [coilAngle(1,:) coil.centre(1)]
            [coilAngle(2,:)	coil.centre(2)]
            [coilAngle(3,:)	coil.centre(3)]
            0	0	0	1];
        
        % Set Di/Dt (the unit is A/us in magventure and A/s in simnibs) 
        s.poslist{1}.pos(1).didt = [T.([site{st},'DiDt'])(targetRow)]*1e6;
        
        % just in case you need barycenters of cortex/scalp surfaces
        % m = mesh_load_gmsh4(s.fnamehead);
       
        run_simnibs(s)
    end
end
%% ===================================== E-Field Group Averaging ===================================
clc; close all;clear all;
pathIn = '/Volumes/Mana_HD/GWM/scans/simnibs/Raw/';
pathOut = '/Volumes/Mana_HD/GWM/scans/simnibs/Analyzed/';
site = {'DLPFC';'FEF';'PP'};
ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
results_folder = "fsavg_overlays";
for st = 1:size(site,1)
    for idx = 1:length(ID)
        sub = ID{idx};
        % Load normal field data
        normal_surf = ['lh.',sub,'_TMS_1-0001_No29_MagVenture_C-B60_Fig8_nii_scalar.fsavg.E.Normal'];
        m{st} = mesh_load_fsresults(char(...
            fullfile(pathOut, sub, site{st},results_folder, normal_surf)));
        % Add to cell
        normals{st,idx} = m{st}.node_data{1}.data;
        normal_surf = [];
    end
end

%% ===================================== E-Field Group Averaging ===================================
clc; close all;clear all;
pathOut = '/scans/simnibs/Analyzed/';
ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
fsavg_msh_name = '_TMS_1-0001_No29_MagVenture_C-B60_Fig8_nii_scalar_fsavg.msh';
site = {'DLPFC';'FEF';'PP'};
field_name = 'E_norm';%E_normal?

for st = 1:size(site,1)   
    for idx = 1:length(ID)
        sub = ID{idx};        
        % load mesh with results transformed to fsaverage space
        m = mesh_load_gmsh4(fullfile(pathOut, ID{idx},  site{st}, '/fsavg_overlays/',[ID{idx},fsavg_msh_name]));       
        % Save the field of each subject
        fields{st,idx} = m.node_data{get_field_idx(m, field_name, 'node')}.data;
    end
end
save([pathOut, 'groupAvg_norm']);

%% ======================================= Individual Vector E nx3  ================================
clc; close all;clear all;
pathOut = '/Volumes/Mana_HD/GWM/scans/simnibs/Analyzed/';
ID = {'118';'119';'121';'122';'123';'124';'125';'126';'127';'128';'129';'130';'131';'132';'133';'134';'136';'137';'138';'139';'140';'141';'142';'143';'145';'146';'147';'148';'149'};
fsavg_msh_name = '_TMS_1-0001_No29_MagVenture_C-B60_Fig8_nii_scalar.msh';

% takes ~30 minutes to run
for idx = 1:length(ID)

    mesh = [];
    % load mesh with results transformed to fsaverage space
    mesh = mesh_load_gmsh4(fullfile(pathOut, ID{idx},  'DLPFC/',[ID{idx},fsavg_msh_name]));
    
    % Save the field of each subject each condition
    save([pathOut, ID{idx},'_vectorE_DLPFC.mat'],'mesh');
    
end
for idx = 1:length(ID)
    
    mesh = [];
    % load mesh with results transformed to fsaverage space
    mesh = mesh_load_gmsh4(fullfile(pathOut, ID{idx},  'FEF/',[ID{idx},fsavg_msh_name]));
    % Save the field of each subject each condition
    save([pathOut, ID{idx},'_vectorE_FEF.mat'],'mesh');

end
for idx = 1:length(ID)
    
    mesh = [];
    % load mesh with results transformed to fsaverage space
    mesh = mesh_load_gmsh4(fullfile(pathOut, ID{idx},  'PP/',[ID{idx},fsavg_msh_name]));
    % Save the field of each subject each condition
    save([pathOut, ID{idx},'_vectorE_PP.mat'],'mesh');
    
end

%% =====================================  Calculate and plot averages ==============================
clc; close all;clear all;
pathOut = '/Volumes/Mana_HD/GWM/scans/simnibs/Analyzed/';
addpath(genpath('/Users/manab/Applications/SimNIBS-3.2/matlab'));
load([pathOut, 'groupAvg_norm']);

% Calculate
for st = 1:length(site)  

    st_fields = fields(st,:);
    st_fields = cell2mat(st_fields);
    avg_field = mean(st_fields, 2);
    std_field = std(st_fields, 0, 2);
    % Plot
    m.node_data = {}; %cleanup fields
    m.node_data{1}.data = avg_field; % add average field
    m.node_data{1}.name = [field_name '_avg'];
    m.node_data{2}.data = std_field; % add std field
    m.node_data{2}.name = [field_name '_std'];
    clearvars st_fields avg_field std_field
    mCond{st}= m;
    
end
st = 1;
% Plot with fields
mesh_show_surface( mCond{st}, 'field_idx', [field_name '_avg'])
mesh_show_surface( mCond{st}, 'field_idx', [field_name '_std'])% Calculate average and standard deviation of the normal at each node
