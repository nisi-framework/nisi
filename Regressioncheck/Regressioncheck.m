%% Regressioncheck
% This file is made to test the correct setup of the NISI code. 
% It tests a 1D NISI and 3D NISI case. 
clear all

%% Copy data for first use
FileCheck = exist('pfade.mat','file');
if FileCheck == 2
        copyRegData() % Subfunction to copy necessary datasets into data folder 
else
        error('Matlab:FileNotFound','Error: pfade.mat is either not defined or not on the Matlab search path.\nCheck NISI handbook (section 1.1 Folder Structure)',FileCheck)  
end
%% Standard Calibration Process
try 
  NISI('regr1','C')
catch err
end
if exist('err','var') == 1
        Standard_Calibration =   [ '"' err.message '"occured in ' err.stack(1,1).name ...
        ' called by ' err.stack(2,1).name]  ;
else
        Standard_Calibration = 'No error occured';
end
pause

try 
  NISI('regr1','M');
catch err
end
if exist('err','var') == 1
        Standard_Measurement_Analysis =   [ '"' err.message '"occured in ' err.stack(1,1).name ...
        ' called by ' err.stack(2,1).name]  ;
else
        Standard_Measurement_Analysis = 'No error occured';
end
pause


%% 3D Calibration Process
try 
  NISI('regr2','C');
catch err
end
if exist('err','var') == 1
        D3_Calibration =   [ '"' err.message '"occured in ' err.stack(1,1).name ...
        ' called by ' err.stack(2,1).name]  ;
else
        D3_Calibration = 'No error occured';
end
pause
try 
  NISI('regr2','M');
catch err
end
if exist('err','var') == 1
        D3_Measurement =   [ '"' err.message '"occured in ' err.stack(1,1).name ...
        ' called by ' err.stack(2,1).name]  ;
else
        D3_Measurement = 'No error occured';
end
pause


%% NL Calibration and Measurement Process

% NOT FUNCTIONAL YET - Date: 17. Jan 2018 - FH
% try 
%   NISI('regr3','C');
% catch err
% end
% if exist('err','var') == 1
%         NL_Calibration =   [ '"' err.message '"occured in ' err.stack(1,1).name ...
%         ' called by ' err.stack(2,1).name]  ;
% else
%         NL_Calibration = 'No error occured';
% end
% 
% try 
%   NISI('regr3','M');
% catch err
% end
% if exist('err','var') == 1
%         NL_Measurement =   [ '"' err.message '"occured in ' err.stack(1,1).name ...
%         ' called by ' err.stack(2,1).name]  ;
% else
%         NL_Measurement = 'No error occured';
% end
% 


%% Clean up
load pfade.mat

fprintf('\ncheck dialog window for data clean up\n')
answer = questdlg('Delete regression data folders from local data repository?',...
    'Data Clean Up','Yes','No','Yes');

switch answer
    case 'Yes'
rmdir([pathDataFolder 'regr1'], 's')
rmdir([pathDataFolder 'regr2'], 's')
    case 'No'
end

clear answer CorM 
% delete([pathDataFolder 'regr1\Measurement\regr1_Measurement_M_Impulse_Response.mat']);
% delete([pathDataFolder 'regr1\Calibration\regr1_Calibration_C_Impulse_Response.mat']);
% delete([pathDataFolder 'regr1\Calibration\regr1_Calibration_C_Impulse_Response_Raw.mat']);
% delete([pathDataFolder 'regr1\Calibration\regr1_Calibration_C_NISI_filtered.mat']);
% delete([pathDataFolder 'regr1\Measurement\regr1_Measurement_M_NISI_filtered.mat']);
% delete([pathDataFolder 'regr1\Calibration\regr1_Calibration_C_NISI_Parameter.mat']);
% delete([pathDataFolder 'regr1\Measurement\regr1_Measurement_M_Solution.mat']);
% delete([pathDataFolder 'regr1\Calibration\regr1_Calibration_C_Solution.mat']); 
% 
% delete([pathDataFolder 'regr2\Measurement\regr2_Measurement_M_Impulse_Response.mat']);
% delete([pathDataFolder 'regr2\Calibration\regr2_Calibration_C_Impulse_Response.mat']);
% delete([pathDataFolder 'regr2\Calibration\regr2_Calibration_C_Impulse_Response_Raw.mat']);
% delete([pathDataFolder 'regr2\Calibration\regr2_Calibration_C_NISI_filtered.mat']);
% delete([pathDataFolder 'regr2\Measurement\regr2_Measurement_M_NISI_filtered.mat']);
% delete([pathDataFolder 'regr2\Calibration\regr2_Calibration_C_NISI_Parameter.mat']);
% delete([pathDataFolder 'regr2\Measurement\regr2_Measurement_M_Solution.mat']);
% delete([pathDataFolder 'regr2\Calibration\regr2_Calibration_C_Solution.mat']); 

% delete([pathDataFolder 'regr3\Measurement\regr3_Measurement_M_Impulse_Response.mat']); % <-- Unvollständig
% delete([pathDataFolder 'regr3\Calibration\regr3_Calibration_C_NISI_filtered.mat']);
% delete([pathDataFolder 'regr3\Calibration\regr3_Calibration_C_NISI_Parameter.mat']);
% delete([pathDataFolder 'regr3\Measurement\regr3_Measurement_M_Solution.mat']);

ccc;

fprintf('\n\nRegressioncheck.m done!\n');



%% Subfunctions
function copyRegData()

load pfade.mat
addpath(pathDataFolder);
%% regr1
% Search for NISI_regr1.mat and regr1 folder in path for NISI data
% (pathDataFolder) - if nonexistent, create dir and copy file to it from
% Regressionchek folder in NISI dir.
dir_regr1 = ([pathDataFolder 'regr1\Calibration\']);
dir_regr1_mat = ([pathDataFolder 'regr1\Calibration\regr1_Calibration_C_NISI.mat']);
if exist(dir_regr1_mat,'file') ~= 2
    if exist(dir_regr1,'dir') ~= 7
        mkdir(dir_regr1); 
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\Testcases\regr1\regr1_Calibration_C_NISI.mat'), dir_regr1_mat)
end
addpath(dir_regr1);

% Create Measurement Folder
dir_regr1_Meas = ([pathDataFolder 'regr1\Measurement\']);
dir_regr1_Meas_mat = ([pathDataFolder 'regr1\Measurement\regr1_Measurement_M_NISI.mat']);
if exist(dir_regr1_Meas_mat,'file') ~= 2
    if exist(dir_regr1_Meas,'dir') ~= 7
        mkdir(dir_regr1_Meas);
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\Testcases\regr1\regr1_Measurement_M_NISI.mat'), dir_regr1_Meas_mat)
end
addpath(dir_regr1_Meas);



%% regr2
% same as regr1
dir_regr2 = ([pathDataFolder 'regr2\Calibration\']);
dir_regr2_mat = ([pathDataFolder 'regr2\Calibration\regr2_Calibration_C_NISI.mat']);
if exist(dir_regr2_mat,'file') ~= 2
    if exist(dir_regr2,'dir') ~= 7
        mkdir(dir_regr2);
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\Testcases\regr2\regr2_Calibration_C_NISI.mat'), dir_regr2_mat)
end
addpath(dir_regr2);

% Create Measurement Folder
dir_regr2_Meas = ([pathDataFolder 'regr2\Measurement\']);
dir_regr2_Meas_mat = ([pathDataFolder 'regr2\Measurement\regr2_Measurement_M_NISI.mat']);
if exist(dir_regr2_Meas_mat,'file') ~= 2
    if exist(dir_regr2_Meas,'dir') ~= 7
        mkdir(dir_regr2_Meas);
    end
    copyfile(fullfile(pathNisiFolder, 'Regressioncheck\Testcases\regr2\regr2_Measurement_M_NISI.mat'), dir_regr2_Meas_mat)
end
addpath(dir_regr2_Meas);


%% regr3
% same as regr1
% dir_regr3 = ([pathDataFolder 'regr3\Calibration\']);
% dir_regr3_mat = ([pathDataFolder 'regr3\Calibration\regr3_Calibration_C_NISI.mat']);
% if exist(dir_regr3_mat,'file') ~= 2
%     if exist(dir_regr3,'dir') ~= 7
%         mkdir(dir_regr3);
%     end
%     copyfile(fullfile(pathNisiFolder, 'Regressioncheck\Testcases\regr3\regr3_Calibration_C_NISI.mat'), dir_regr3_mat)
% end
% addpath(dir_regr3);
% 
% % Create Measurement Folder - plus copy NISI_regr3_Measurement.mat as measurement case
% dir_regr3_Meas = ([pathDataFolder 'regr3\Measurement\']);
% dir_regr3_Meas_mat = ([pathDataFolder 'regr3\Measurement\regr3_Measurement_M_NISI.mat']);
% if exist(dir_regr3_Meas_mat,'file') ~= 2
%     if exist(dir_regr3_Meas,'dir') ~= 7
%         mkdir(dir_regr3_Meas);        
%     end
%     copyfile(fullfile(pathNisiFolder, 'Regressioncheck\Testcases\regr3\regr3_Measurement_M_NISI.mat'), ...
%              dir_regr3_Meas_mat)
% end
% addpath(dir_regr3_Meas);

end

