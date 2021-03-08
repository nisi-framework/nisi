function[Time, Calibration_Heat_Flux, Calibration_Temperature] = E_check(nameDataset)

load('pfade.mat')


configNameBuffer = [ nameDataset '_Config' ];     %Create full config file name string from the particular handle
configLoadHandle = str2func(configNameBuffer);    %Create function handle from name string
[nameCalFolder,~] = configLoadHandle();


fname = [pathDataFolder, nameDataset '\' nameCalFolder '\'...
          nameDataset '_' nameCalFolder,'_C_NISI.mat'];
      
load(fname);
figure('Name','Imported data unfiltered','NumberTitle','off'),
if exist('Calibration_Temperature','var') == 1
	plotyy(Time,squeeze(Calibration_Temperature),Time,squeeze(Calibration_Heat_Flux));
elseif exist('Calibration_Pressure','var') == 1  
	plotyy(Time,squeeze(Calibration_Pressure),Time,squeeze(Calibration_Heat_Flux));
else
    fprintf(['\n\nWarning! Neither Calibration_Temperature nor Calibration_Pressure '...
        'exists in ' fname '.\n'] )
end
drawnow;


fname2 = [pathDataFolder, nameDataset '\' nameCalFolder '\'...
          nameDataset '_' nameCalFolder,'_C_NISI_filtered.mat'];

if exist(fname2,'file') == 2
    load(fname2);
    figure('Name','Imported data filtered','NumberTitle','off'),
    if exist('Calibration_Temperature','var') == 1
        plotyy(Time,squeeze(Calibration_Temperature),Time,squeeze(Calibration_Heat_Flux));
    elseif exist('Calibration_Pressure','var') == 1  
        plotyy(Time,squeeze(Calibration_Pressure),Time,squeeze(Calibration_Heat_Flux));
    else
        fprintf(['\n\nWarning! Neither Calibration_Temperature nor Calibration_Pressure '...
            'exists in ' fname2 '.\n'] )
    end
    drawnow;
else
    fprintf(['\n\nWarning! Filtered data file ' nameDataset '_' nameCalFolder ...
        '_C_NISI_filtered.mat does not exist and was not plotted.\n'] )
end

assignin('base','Time',Time);
if exist('Calibration_Temperature','var') == 1
    assignin('base','Calibration_Temperature',Calibration_Temperature);
elseif exist('Calibration_Pressure','var') == 1  
    assignin('base','Calibration_Pressure',Calibration_Pressure);
end    
assignin('base','Calibration_Heat_Flux',Calibration_Heat_Flux);
end



