function [tGes,Deltat] = tRawCheck(nameDataset,nameFolder)

load pfade.mat   

pathRawDataC = [pathDataFolder, nameDataset '\' nameFolder '\RAW_DATA\'];
[traceTitleCal,numberOfDatasetC] = readTraceTitleAndNumberOfDataset(pathRawDataC);

% try
[Time,~,~,~]=LeCroy_Oszilloscope_and_Laserline(nameDataset,pathDataFolder,...
                                        pathRawDataC,'C',traceTitleCal,numberOfDatasetC,...
                                        'L','T', '0','0',...
                                        0,0,nameFolder,'');  
% end

tGes = Time(end);
Deltat = Time(2) - Time(1);

end