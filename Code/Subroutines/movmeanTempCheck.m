nameDataset         = 'Stylo5mm_5micros';

load([nameDataset '_SettingsForAnalysis.mat']);

A=(50:25:250)*1000;
% A=50000:25000:250000;
% A=[50000 100000 250000 500000];
% A=[1000];

for n=1:length(A)
    movmeanTemp = A(n);
    
    
%% Calibration Data Import 
[rawTcSignal,rawLaserSignal]=LeCroy_Oszilloscope_and_Laserline(nameDataset,pathDataFolder,...
                                        pathRawData,'C',traceTitleCal,numberOfDatasetC,...
                                        C1DataType,C2DataType, C3DataType,C4DataType,...
                                        flagSave,flagPlot,exp);


%% Data Preparation (HF and Temp calculation) 
 
 % Data import
 tabname = [pathNisiFolder 'LDM500.xls'];
 Laservolt = xlsread(tabname,'A10:A20');
 LaserP = xlsread(tabname,'B10:B20');

 % Temp calculation
 Temp = rawTcSignal(1:end,2)*tcCalcFactor;
 Temp = movmean(Temp,movmeanTemp);
 Calibration_Temperature(n,1,:) = Temp(:,1);

 % Time (read in)
 Time = rawTcSignal(1:end,1);
 Time = Time - Time(1);
 
 % HF calculation
 E_L_P = pchip(Laservolt,LaserP,rawLaserSignal(:,2)); 
 HF_L = E_L_P./area * (1-r); % Net Heat flux = Laser power/area * (1-reflectivity)
 Calibration_Heat_Flux(1,1,:) = HF_L(1:end,1);
 
%  calTempDiff(n,1)=0;
%  for m = 2: length(E_L_P) % E_L_P because it is as long as Calibration_Temperature
%     calTempDiff(n,m) = abs(Calibration_Temperature(n,1,m)- Calibration_Temperature(n,1,m-1));
%  end
end


figure()
hold on
for n = 1:length(A)
% plotyy(Time,squeeze(Calibration_Temperature(n,1,:)),Time,squeeze(calTempDiff(n,:)));
plot(Time,squeeze(Calibration_Temperature(n,1,:)),'DisplayName',num2str(A(n)))
end
legend('show');
hold off
