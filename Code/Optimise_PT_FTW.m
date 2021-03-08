function [ftw_opt] = Optimise_PT_FTW(nameDataset,nameMeasFolder,varargin)
% This function opimizes the Future Time Window (ftw). 
%
% Ursprungsversion von
% Ben Kefford
%
% Modified:
% F. Hufgard
%
%
% INPUT
%       nameDataset
%       nameMeasFolder
%       varargin (optional) - Choosable parameters:
%                                   nameCalFolder
%
% OUTPUT
%       ftw_opt - optimal future time step
%
% Data should be filtered first and then the filtering option in the Config file should
% be turned off to improve computational speed. 
% 
% Penetration Time (pt) should be chosen acording to the theoretical value (Kaviany_2001_01). 
% In order to account for uncertainty in sensor position, this theoretical value can be slightly 
% varied in this function (see definition of 'pt_list').

 
load pfade.mat 

% nameDataset = 'Maya_HkCC_5mm_K_Danzig';
% nameMeasFolder = '200525_run7_160mgs_Arb';

    % load Input Heat Flux
    nameNISIM = [pathDataFolder, nameDataset '\' nameMeasFolder '\' ...
                nameDataset '_' nameMeasFolder '_M_NISI.mat'];       
    load(nameNISIM);  

    % Load delta_t (and Time)
    nameNISIM1 = [pathDataFolder, nameDataset '\' nameMeasFolder '\' ...
                nameDataset '_' nameMeasFolder '_M_NISI_filtered.mat'];     
            try
    load(nameNISIM1);
            catch
                NISI(nameDataset,'M','nameMeasFolder',nameMeasFolder,'noPlot','noClear')
                warning('Data is not filtered yet\n-> Initial NISI run was conducted!\nThere might not be the right IR chosen for this run, but that does not matter for the rest of the code.')
                load(nameNISIM1);
            end
    
    % Calculate Time_Step_Reduction
    tsr = ceil(length(Input_Heat_Flux)/length(Time));
    delta_t = Time(2);

Input_HF(1,:) = Input_Heat_Flux(1,1,1:tsr:end); 
idx_s = find(Input_HF(1,:)>max(Input_HF)*0.03,1,'first');
idx_e = find(Input_HF(1,:)>max(Input_HF)*0.03,1,'last');
   

clearvars -except nameDataset Input_HF delta_t idx_s idx_e nameMeasFolder varargin pathDataFolder

pt_theory =  0.15; % <- 0.25=random   % Choose according to Kaviany_2001_01 --> pt = (x/3.6)^2*rho*cp/k
% pt_list = [pt_theory-delta_t pt_theory pt_theory+delta_t]; % <- minimal variation of pt <-> for the tests so far this does not yield a large variance
pt_list = pt_theory; % <- no variation of pt

for i = 1:length(pt_list)      
    pt = pt_list(i);
    B = [1 2 3]; % logspace(-3,-1,20)
    dftw = (B(end)-B(1))/length(B);
    A = [];         % initialize A
   
    while dftw > delta_t/2 && dftw > 0.01
    B = round(B./delta_t).*delta_t;         % Auf Vielfacher der Zeitschrittweite runden 
        dftw = dftw/2;
        for j = length(A)+1:length(B)
            ftw = B(j);

            NISI(nameDataset,'M','pt', pt, 'ftw', ftw,'noWB','noPlot','noClear','nameMeasFolder',nameMeasFolder, varargin{:})
            
            load([pathDataFolder nameDataset '\' nameMeasFolder '\' nameDataset '_' nameMeasFolder '_M_Solution.mat'])
            res = Msrmnt_Heat_Flux - Input_HF;  	% res = residual
            M_HF(i,j,:) = Msrmnt_Heat_Flux(1,:);
            NRMSD(i,j) = 100*sqrt((1/length(res(idx_s:idx_e))*sum(res(idx_s:idx_e).^2)))/max(Input_HF);
%             NSTD(i,j) = 100*std(res(idx_s:idx_e))/max(Input_HF);  % Give normalized Standard Deviation of residual as a percent of maximum input heat flux
            clear res
            
            if dftw > delta_t/2 && dftw > 0.01
                multiWaitbar('Optimization Progress',j/(length(B)+2)*i/length(pt_list));
            else
                multiWaitbar('Optimization Progress',j/length(B)*i/length(pt_list));
            end
        end
        [~,idx1] = min(NRMSD(i,1:length(B)));
        A = B;
        B = [B B(idx1)-dftw B(idx1)+dftw];
    end     % ftw Schleife


        if i~=1
            a = length(ftw_list(1,:));
            b = length(A);
        if a<b
            ftw_list(1:i-1,a+1:b) = NaN;
            NRMSD(1:i-1,a+1:b) = NaN;
            M_HF(1:i-1,a+1:b,:) = NaN;
        elseif a>b
            A(b+1:a) = NaN;
        end
        end

        [ftw_list(i,:),I] = sort(A);
        NRMSD(i,:) = NRMSD(i,I);
        M_HF(i,:,:) = M_HF(i,I,:);
        clear A B I 
end  % pt Schleife

%Save file
load pfade.mat
file = ['Optimize_PT_FTW_Results.mat'];
save([pathDataFolder nameDataset '\' file], 'pt_list', 'ftw_list','NRMSD');

%Plot curves
h = figure('name','Optimal Future Time Window'); hold
set(gca, 'ColorOrder', colormap(autumn(length(ftw_list(:,1))+1)));
legend show
for i = 1:length(pt_list) 
    plot(ftw_list(i,:),NRMSD(i,:),'Displayname',['Penetration Time: ' num2str(pt_list(i))])
    [val,pos] = min(NRMSD(i,:));
    plot(ftw_list(i,pos),val,'rd','HandleVisibility','off')
    ftw_opt_vonpt(i) = ftw_list(i,pos);     % find best ftw for current pt
    M_HF_opt(i,:) = M_HF(i,pos,:);
end
xlabel('Future Time Window, s')
ylabel('% NRMSD of Heat Flux')
ylim([0 max(max(NRMSD(:,end),5*min(min((NRMSD)))))]);


% find best overall ftw and its pt
[ftw_opt,idx_pt] = min(ftw_opt_vonpt);
pt_opt = pt_list(idx_pt);

ftw_opt = round(ftw_opt/delta_t)*delta_t;   % round value
% pt_opt = round(pt_opt/delta_t)*delta_t;   % round value

% calulate residual
normRes = (M_HF_opt(idx_pt,:)-Input_HF) - max(Input_HF)/4;


% Plot best Match w/o Shift
figure('name','Plot best Match w/o Shift'); hold
legend show
plot(Input_HF,'k-', 'DisplayName','Input HF')
plot(M_HF_opt(idx_pt,:),'b-', 'DisplayName',['Msrmnt HF, pt = ' num2str(pt_opt) ', ftw = ' num2str(ftw_opt)])
plot(normRes,'r-','DisplayName','Residual')
plot([1 length(normRes)],  -max(Input_HF)/4*[1 1], 'k--' ,'HandleVisibility','off')
plot([1 length(normRes)],  (max(Input_HF)*0.05 -max(Input_HF)/4)*[1 1], 'k:' ,'DisplayName','5% Deviation')
plot([1 length(normRes)], (-max(Input_HF)*0.05 -max(Input_HF)/4)*[1 1], 'k:' ,'HandleVisibility','off')
multiWaitbar('Close All');


% Write optimal parameters into Config file in fig folderC
overwriteConfig(nameDataset,'ftw',ftw_opt,'pt',pt_opt)

end

