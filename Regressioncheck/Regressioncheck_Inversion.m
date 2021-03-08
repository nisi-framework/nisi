%% Regression Inversion
% This program features testcases with analytically calculated impulse responses. 
% It only checks if new or changed inversion algorithms run properly.
% 
%   run     heat flux           nameDataset             Filter
%--------------------------------------------------------------------
%   1       gauss               RegressionM1        -
%   2       gauss + noise       RegressionM2        predefined in Config
%   3       traingular          RegressionM3        -
%
% Change the inversion algoritm in config-files RegressionM1, RegressionM2,
% RegressionM3
clear all

%% Copy data for first use
copyRegData() % Subfunction to copy necessary datasets into data folder 

%% 
%Measurement Gauss
try 
  NISI('RegressionM1','M');
catch err
end
if exist('err','var') == 1
        Standard_Calibration =   [ '"' err.message '"occured in ' err.stack(1,1).name ...
        ' called by ' err.stack(2,1).name]  ;
else
        Standard_Calibration = 'No error occured';
end
fprintf('Finished regression run 1,\nPress any key to continue')
pause;

%% 
%Measurement Gauss + Noise
try 
  NISI('RegressionM2','M');
catch err
end
if exist('err','var') == 1
        D3_Calibration =   [ '"' err.message '"occured in ' err.stack(1,1).name ...
        ' called by ' err.stack(2,1).name]  ;
else
        D3_Calibration = 'No error occured';
end
fprintf('Finished regression run 2,\nPress any key to continue')
pause;

%% 
%Measurement Triangular
try 
  NISI('RegressionM3','M');
catch err
end
if exist('err','var') == 1
        NL_Calibration =   [ '"' err.message '"occured in ' err.stack(1,1).name ...
        ' called by ' err.stack(2,1).name]  ;
else
        NL_Calibration = 'No error occured';
end
fprintf('Finished regression run 3, \nPress any key to clean up data repository')
pause;

%% Clean up
load pfade.mat

fprintf('\ncheck dialog window for data clean up\n')
answer = questdlg('Delete regression data folders from local data repository?',...
    'Data Clean Up','Yes','No','Yes');

switch answer
    case 'Yes'
rmdir([pathDataFolder 'RegressionM1'], 's')
rmdir([pathDataFolder 'RegressionM2'], 's')
rmdir([pathDataFolder 'RegressionM3'], 's')
    case'No'
end

ccc;

fprintf('\n\nRegressioncheck_Inversion.m done!\n');

clear
%% Subfunctions
function copyRegData()
try 
load pfade.mat;
catch
       error(['Matlab could not find pfade.mat']); 
end
if exist('pathDataFolder','var')==1 && exist('pathNisiFolder','var')==0
        error(['pathDataFolder in pfade.mat is not defined!']);
elseif exist('pathNisiFolder','var')==0 && exist('pathDataFolder','var')==1
        error(['pathNisiFolder in pfade.mat is not defined!']);
elseif exist('pathNisiFolder','var')==0 && exist('pathDataFolder','var')==0
    error(['pfade.mat is empty']);
else
    
end


addpath(pathDataFolder);
%% regr1
% Search for NISI_regr1.mat and regr1 folder in path for NISI data
% (pathDataFolder) - if nonexistent, create dir and copy file to it from
% Regressionchek folder in NISI dir.
dir_regr1 = ([pathDataFolder 'RegressionM1\Gauss\']);
dir_regr1_mat1 = ([pathDataFolder 'RegressionM1\Gauss\RegressionM1_Gauss_C_NISI.mat']);
dir_regr1_mat2 = ([pathDataFolder 'RegressionM1\Gauss\RegressionM1_Gauss_C_Impulse_Response.mat']);
dir_regr1_mat3 = ([pathDataFolder 'RegressionM1\Gauss\RegressionM1_Gauss_C_NISI_filtered.mat']);
dir_regr1_mat4 = ([pathDataFolder 'RegressionM1\Gauss\RegressionM1_Gauss_C_NISI_Parameter.mat']);
dir_regr1_mat5 = ([pathDataFolder 'RegressionM1\Gauss\RegressionM1_Gauss_C_Solution.mat']);
if exist(dir_regr1_mat1,'file') ~= 2
    if exist(dir_regr1,'dir') ~= 7
        mkdir(dir_regr1); 
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss\NISI_Test_Gauss_C_NISI.mat'), dir_regr1_mat1)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss\NISI_Test_Gauss_C_Impulse_Response.mat'), dir_regr1_mat2)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss\NISI_Test_Gauss_C_NISI_filtered.mat'), dir_regr1_mat3)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss\NISI_Test_Gauss_C_NISI_Parameter.mat'), dir_regr1_mat4)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss\NISI_Test_Gauss_C_Solution.mat'), dir_regr1_mat5)
end
addpath(dir_regr1);

% Create Measurement Folder
dir_regr1_Meas = ([pathDataFolder 'RegressionM1\Msrt_Gauss\']);
dir_regr1_Meas_mat1 = ([pathDataFolder 'RegressionM1\Msrt_Gauss\RegressionM1_Msrt_Gauss_M_NISI.mat']);
dir_regr1_Meas_mat2 = ([pathDataFolder 'RegressionM1\Msrt_Gauss\RegressionM1_Msrt_Gauss_C_Impulse_Response.mat']);
dir_regr1_Meas_mat3 = ([pathDataFolder 'RegressionM1\Msrt_Gauss\RegressionM1_Msrt_Gauss_M_NISI_filtered.mat']);
dir_regr1_Meas_mat4 = ([pathDataFolder 'RegressionM1\Msrt_Gauss\RegressionM1_Msrt_Gauss_C_NISI_Parameter.mat']);
dir_regr1_Meas_mat5 = ([pathDataFolder 'RegressionM1\Msrt_Gauss\RegressionM1_Msrt_Gauss_M_Solution.mat']);
if exist(dir_regr1_Meas_mat1,'file') ~= 2
    if exist(dir_regr1_Meas,'dir') ~= 7
        mkdir(dir_regr1_Meas);
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss\NISI_Test_Msrt_Gauss_M_NISI.mat'), dir_regr1_Meas_mat1)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss\NISI_Test_Gauss_C_Impulse_Response.mat'), dir_regr1_Meas_mat2) 
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss\NISI_Test_Msrt_Gauss_M_NISI_filtered.mat'), dir_regr1_Meas_mat3)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss\NISI_Test_Gauss_C_NISI_Parameter.mat'), dir_regr1_Meas_mat4)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss\NISI_Test_Msrt_Gauss_M_Solution.mat'), dir_regr1_Meas_mat5)
    
end
addpath(dir_regr1_Meas);



%% regr2
% same as regr1
dir_regr2 = ([pathDataFolder 'RegressionM2\Gauss_Noise\']);
dir_regr2_mat1 = ([pathDataFolder 'RegressionM2\Gauss_Noise\RegressionM2_Gauss_Noise_C_NISI.mat']);
dir_regr2_mat2 = ([pathDataFolder 'RegressionM2\Gauss_Noise\RegressionM2_Gauss_Noise_C_Impulse_Response.mat']);
dir_regr2_mat3 = ([pathDataFolder 'RegressionM2\Gauss_Noise\RegressionM2_Gauss_Noise_C_NISI_filtered.mat']);
dir_regr2_mat4 = ([pathDataFolder 'RegressionM2\Gauss_Noise\RegressionM2_Gauss_Noise_C_NISI_Parameter.mat']);
dir_regr2_mat5 = ([pathDataFolder 'RegressionM2\Gauss_Noise\RegressionM2_Gauss_Noise_C_Solution.mat']);
if exist(dir_regr2_mat1,'file') ~= 2
    if exist(dir_regr2,'dir') ~= 7
        mkdir(dir_regr2);
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss_Noise\NISI_Test_Gauss_Noise_C_NISI.mat'), dir_regr2_mat1)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss_Noise\NISI_Test_Gauss_Noise_C_Impulse_Response.mat'), dir_regr2_mat2)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss_Noise\NISI_Test_Gauss_Noise_C_NISI_filtered.mat'), dir_regr2_mat3)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss_Noise\NISI_Test_Gauss_Noise_C_NISI_Parameter.mat'), dir_regr2_mat4)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Gauss_Noise\NISI_Test_Gauss_Noise_C_Solution.mat'), dir_regr2_mat5)

end
addpath(dir_regr2);

% Create Measurement Folder
dir_regr2_Meas = ([pathDataFolder 'RegressionM2\Msrt_Gauss_Noise\']);
dir_regr2_Meas_mat1 = ([pathDataFolder 'RegressionM2\Msrt_Gauss_Noise\RegressionM2_Msrt_Gauss_Noise_M_NISI.mat']);
dir_regr2_Meas_mat2 = ([pathDataFolder 'RegressionM2\Msrt_Gauss_Noise\RegressionM2_Msrt_Gauss_Noise_C_Impulse_Response.mat']);
dir_regr2_Meas_mat3 = ([pathDataFolder 'RegressionM2\Msrt_Gauss_Noise\RegressionM2_Msrt_Gauss_Noise_M_NISI_filtered.mat']);
dir_regr2_Meas_mat4 = ([pathDataFolder 'RegressionM2\Msrt_Gauss_Noise\RegressionM2_Msrt_Gauss_Noise_C_NISI_Parameter.mat']);
dir_regr2_Meas_mat5 = ([pathDataFolder 'RegressionM2\Msrt_Gauss_Noise\RegressionM2_Msrt_Gauss_Noise_M_Solution.mat']);
if exist(dir_regr2_Meas_mat1,'file') ~= 2
    if exist(dir_regr2_Meas,'dir') ~= 7
        mkdir(dir_regr2_Meas);
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss_Noise\NISI_Test_Msrt_Gauss_Noise_M_NISI.mat'), dir_regr2_Meas_mat1)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss_Noise\NISI_Test_Gauss_Noise_C_Impulse_Response.mat'), dir_regr2_Meas_mat2) 
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss_Noise\NISI_Test_Msrt_Gauss_Noise_M_NISI_filtered.mat'), dir_regr2_Meas_mat3)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss_Noise\NISI_Test_Gauss_Noise_C_NISI_Parameter.mat'), dir_regr2_Meas_mat4)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Gauss_Noise\NISI_Test_Msrt_Gauss_Noise_M_Solution.mat'), dir_regr2_Meas_mat5)
end
addpath(dir_regr2_Meas);


%% regr3
% same as regr1
dir_regr3 = ([pathDataFolder 'RegressionM3\Triangular\']);
dir_regr3_mat1 = ([pathDataFolder 'RegressionM3\Triangular\RegressionM3_Triangular_C_NISI.mat']);
dir_regr3_mat2 = ([pathDataFolder 'RegressionM3\Triangular\RegressionM3_Triangular_C_Impulse_Response.mat']);
dir_regr3_mat3 = ([pathDataFolder 'RegressionM3\Triangular\RegressionM3_Triangular_C_NISI_filtered.mat']);
dir_regr3_mat4 = ([pathDataFolder 'RegressionM3\Triangular\RegressionM3_Triangular_C_NISI_Parameter.mat']);
dir_regr3_mat5 = ([pathDataFolder 'RegressionM3\Triangular\RegressionM3_Triangular_C_Solution.mat']);
if exist(dir_regr3_mat1,'file') ~= 2
    if exist(dir_regr3,'dir') ~= 7
        mkdir(dir_regr3);
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Triangular\NISI_Test_Triangular_C_NISI.mat'), dir_regr3_mat1)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Triangular\NISI_Test_Triangular_C_Impulse_Response.mat'), dir_regr3_mat2)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Triangular\NISI_Test_Triangular_C_NISI_filtered.mat'), dir_regr3_mat3)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Triangular\NISI_Test_Triangular_C_NISI_Parameter.mat'), dir_regr3_mat4)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Triangular\NISI_Test_Triangular_C_Solution.mat'), dir_regr3_mat5)
end
addpath(dir_regr3);

% Create Measurement Folder - plus copy NISI_regr3_Measurement.mat as measurement case
dir_regr3_Meas = ([pathDataFolder 'RegressionM3\Msrt_Triangular\']);
dir_regr3_Meas_mat1 = ([pathDataFolder 'RegressionM3\Msrt_Triangular\RegressionM3_Msrt_Triangular_M_NISI.mat']);
dir_regr3_Meas_mat2 = ([pathDataFolder 'RegressionM3\Msrt_Triangular\RegressionM3_Msrt_Triangular_C_Impulse_Response.mat']);
dir_regr3_Meas_mat3 = ([pathDataFolder 'RegressionM3\Msrt_Triangular\RegressionM3_Msrt_Triangular_M_NISI_filtered.mat']);
dir_regr3_Meas_mat4 = ([pathDataFolder 'RegressionM3\Msrt_Triangular\RegressionM3_Msrt_Triangular_C_NISI_Parameter.mat']);
dir_regr3_Meas_mat5 = ([pathDataFolder 'RegressionM3\Msrt_Triangular\RegressionM3_Msrt_Triangular_M_Solution.mat']);
if exist(dir_regr3_Meas_mat1,'file') ~= 2
    if exist(dir_regr3_Meas,'dir') ~= 7
        mkdir(dir_regr3_Meas);        
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Triangular\NISI_Test_Msrt_Triangular_M_NISI.mat'), dir_regr3_Meas_mat1)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Triangular\NISI_Test_Triangular_C_Impulse_Response.mat'), dir_regr3_Meas_mat2) 
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Triangular\NISI_Test_Msrt_Triangular_M_NISI_filtered.mat'), dir_regr3_Meas_mat3)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Triangular\NISI_Test_Triangular_C_NISI_Parameter.mat'), dir_regr3_Meas_mat4)
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\NISI_Test\Msrt_Triangular\NISI_Test_Msrt_Triangular_M_Solution.mat'), dir_regr3_Meas_mat5)
end
addpath(dir_regr3_Meas);

end

