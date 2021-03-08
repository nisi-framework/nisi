function [ Filtered_Data ] = filter_sharpPulseEdges_NISI(Input_Data,Timestep_Reduction,nPlateaus)
% Custom Filter, der ein gepulstes Lasersignal auf unterschiedlichen Spannungsplateaus
% (z.B.: 0-10%, 10-20%,usw...) mittelt. Dadurch wird der scharfe Anstieg und Abfall des Lasersignals
% nicht gedämpft und der HF-Wert (maximal) geglättet. 
% Praktisch wird das HF-Signal als perfektes Rechteckssingal interpretiert.
%
% Ursprungsversion von
% Fabian Hufgard
%
% 

% Waitbar settings
nameMWB = 'Running SPE filter on Calibration_Heat_Flux';

if isempty(nPlateaus)
    nPlateaus = 3; % Default Value: 3
    warning(['Variable nPlat_SPE was not defined in Config File and set to value 3.' ...
            ' This variable should be the 31st variable given out in the Config File.'])
end  

[m,n] = size(Input_Data);
if n==1 && m > 10
    buffer(1,:) = Input_Data(:,1);
    clear Input_Data
    Input_Data=buffer;
    warning('Data vector has been changed form (:,1) to (1,:) within shaprePulseEdges filter');
end

bufferSPE = Input_Data;
[maxVal,~] = max(bufferSPE);

for i = 1:nPlateaus     
    uG(i) = maxVal*(1/(2*nPlateaus)+ (i-1)/nPlateaus); % Untergrenze des i-ten Plateaus
end


%% Filter 
% 0-plateau
[~,idx] = find(bufferSPE(1,:)<uG(1));
bufferSPE(idx) = 0;
clear idx
multiWaitbar( nameMWB , 1/(nPlateaus+1));

% Plateaus i bis zweithöchstes
for i = 1:nPlateaus-1  
    [~,idx] = find(bufferSPE(1,:)>=uG(i) & bufferSPE(1,:)<uG(i+1));
    meanVal = mean(bufferSPE(1,idx));
    bufferSPE(idx) = meanVal;
    clear idx
    
    multiWaitbar( nameMWB , (1+i)/(nPlateaus+1));
end

% höchests Plateau
[~,idx] = find(bufferSPE(1,:)>=uG(nPlateaus));
meanVal = mean(bufferSPE(1,idx));
bufferSPE(idx) = meanVal;


%% Reduction
Filtered_Data = bufferSPE(1,1:Timestep_Reduction:end);


multiWaitbar( nameMWB , 1);
end