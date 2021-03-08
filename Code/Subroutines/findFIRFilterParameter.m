function [FIR_Runs_opt,FIR_Config_opt] = findFIRFilterParameter(nameDataset,Time,rawInputSignal,...
                                     FIR_Config_Possibilities,FIR_Runs_Possibilities)
    
	flagPlot = 1;
 
    configNameBuffer = [ nameDataset '_Config' ];      %Create full config file name string from the particular handle
    configLoadHandle = str2func(configNameBuffer);    %Create function handle from name string
    [~,~,~,~,~,~,~,~                                                , ...
    Timestep_Reduction                                              , ...
    Temperature_Zero                                                , ...
    ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~  ]    ... 
    = configLoadHandle();
    lenghtFilteredSignal = ceil(length(rawInputSignal)/Timestep_Reduction);
   
   
    Delta_t = Time(2)-Time(1);
    StepsZero = double(int32(Temperature_Zero / Delta_t));
    Mean_Zero_Level  = sum(rawInputSignal(1:StepsZero)) / StepsZero;
    rawInputSignal = rawInputSignal - Mean_Zero_Level;
    rawInputSignal(1:StepsZero) = 0;
    
    
    % Variables definitions for later usage
    ExCount = 0;
    
    fTC( 1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities),1:lenghtFilteredSignal) = 0;
    dfTC(1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities),1:lenghtFilteredSignal) = 0;
    
    
    fTC_terms(       1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1,1:lenghtFilteredSignal) = 0;
    fTC_terms_mprev( 1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1,1:lenghtFilteredSignal) = 0;
    fTC_terms_mnow(  1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1,1:lenghtFilteredSignal) = 0;
    r_fTC(1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1) = 0;
    y_fTC(1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1) = 0;
    
    fTC_terms_cc(       1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1,1:lenghtFilteredSignal) = 0;
    dfTC_terms_cc(      1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1,1:lenghtFilteredSignal) = 0;
    fTC_terms_mnow_cc(  1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1,1:lenghtFilteredSignal) = 0;
    fTC_terms_mnext_cc( 1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1,1:lenghtFilteredSignal) = 0;
    dfTC_terms_mnow_cc( 1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1,1:lenghtFilteredSignal) = 0;
    dfTC_terms_mnext_cc(1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1,1:lenghtFilteredSignal) = 0;

    r_fTC_cc( 1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1) = 0;
    r_dfTC_cc(1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1) = 0;
    y_fTC_cc( 1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1) = 0;
    y_dfTC_cc(1:length(FIR_Runs_Possibilities),1:length(FIR_Config_Possibilities)-1) = 0;
    
    
    fTC_terms_rc(       1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities),1:lenghtFilteredSignal) = 0;
    dfTC_terms_rc(      1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities),1:lenghtFilteredSignal) = 0;
    fTC_terms_mnow_rc(  1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities),1:lenghtFilteredSignal) = 0;
    fTC_terms_mnext_rc( 1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities),1:lenghtFilteredSignal) = 0;
    dfTC_terms_mnow_rc( 1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities),1:lenghtFilteredSignal) = 0;
    dfTC_terms_mnext_rc(1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities),1:lenghtFilteredSignal) = 0;

    r_fTC_rc( 1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities)) = 0;
    r_dfTC_rc(1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities)) = 0;
    y_fTC_rc( 1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities)) = 0;
    y_dfTC_rc(1:length(FIR_Runs_Possibilities)-1,1:length(FIR_Config_Possibilities)) = 0;
    
    
   %% Filter Data using all filter possibilities
       
    if flagPlot
    figure('Name','Phase Plane'); % [num2str(FIR_Runs) ' FIR runs'])
    hold  
    end
    

    percstp = 1/(length(FIR_Runs_Possibilities)*length(FIR_Config_Possibilities));
    
    for ii = 1:length(FIR_Runs_Possibilities)
       perc = (ii-1)/length(FIR_Runs_Possibilities);
       for jj = 1:length(FIR_Config_Possibilities)
       multiWaitbar('Find FIR Filter Parameters',perc);
       perc = perc + percstp;
       
       FIR_Runs = FIR_Runs_Possibilities(ii);
       FIR_Config = FIR_Config_Possibilities(jj);

       [ fTC(ii,jj,:) ]  = filter_r_fir_NISI(rawInputSignal, ...
                          FIR_Config(1),FIR_Runs(1),...
                          Timestep_Reduction,1);
        for k=2:length(fTC(ii,jj,:))
        dfTC(ii,jj,k)  = fTC(ii,jj,k)-fTC(ii,jj,k-1);
        end 
        
        if flagPlot
            A(1,:)=fTC(ii,jj,:);
            B(1,:)=dfTC(ii,jj,:);
            scatter(A,B,3,'DisplayName',['FIR config: ' num2str(FIR_Config)...
                                            '; FIR runs: ' num2str(FIR_Runs)]); 
        end
        
        
        % Abbruchkriterium für rückläufige Filterqualität
        if jj>1
        for k = 1:lenghtFilteredSignal
         fTC_terms(      ii,jj,k) =  fTC(ii,jj-1,k) *  fTC(ii,jj,k);
         fTC_terms_mprev(ii,jj,k) =  fTC(ii,jj-1,k)^2.;
         fTC_terms_mnow( ii,jj,k) =  fTC(ii,jj,  k)^2.;
        end
        % cross-correlation function (Eq. 4)
        r_fTC( ii,jj) = sum( fTC_terms(ii,jj,:));
        % cross-corellation parameter (Eq. 5)
        y_fTC( ii,jj) = r_fTC( ii,jj)/(sqrt(sum( fTC_terms_mprev(ii,jj,:))*sum( fTC_terms_mnow(ii,jj,:))));
        
        if y_fTC(ii,jj) < y_fTC(ii,jj-1)
            ExCount = ExCount+1;
        else
            ExCount = 0;
        end % y_fTC_cc(ii,jj) < y_fTC_cc(ii,jj-1) 
        
        if ExCount > 4
            break
        end % ExCount > 4
        
        end %jj>1
        
       end % jj
       
    end % ii 
      legend show


      
    %% Cross correlation 
    % - as in Frankel 2014 Phase-Plane and Cross-Correlation Analysis for Estimating
    % Optimal Regularization in Inverse Heat Conduction - Eq. 4 und Eq. 5
    
    %% Cross correlation - Correlate Config 

    
    if flagPlot
    figure('Name','Cross-Correlation Parameters - correlate Config')
    hold on
    end
    
    for mi   = 1:length(FIR_Runs_Possibilities)
        FIR_Runs = FIR_Runs_Possibilities(mi);
    for mj   = 1:length(FIR_Config_Possibilities)-1
       FIR_Config = FIR_Config_Possibilities(mj);
      for k = 1:lenghtFilteredSignal
         fTC_terms_cc(      mi,mj,k) =  fTC(mi,mj,  k) *  fTC(mi,mj+1,k);
        dfTC_terms_cc(      mi,mj,k) = dfTC(mi,mj,  k) * dfTC(mi,mj+1,k);
         fTC_terms_mnow_cc( mi,mj,k) =  fTC(mi,mj,  k)^2.;
         fTC_terms_mnext_cc(mi,mj,k) =  fTC(mi,mj+1,k)^2.;
        dfTC_terms_mnow_cc( mi,mj,k) = dfTC(mi,mj,  k)^2.;
        dfTC_terms_mnext_cc(mi,mj,k) = dfTC(mi,mj+1,k)^2.;
      end
      % cross-correlation function (Eq. 4)
      r_fTC_cc( mi,mj) = sum( fTC_terms_cc(mi,mj,:));
      r_dfTC_cc(mi,mj) = sum(dfTC_terms_cc(mi,mj,:));
      % cross-corellation parameter (Eq. 5)
      y_fTC_cc( mi,mj) = r_fTC_cc( mi,mj)/(sqrt(sum( fTC_terms_mnow_cc(mi,mj,:))*sum( fTC_terms_mnext_cc(mi,mj,:))));
      y_dfTC_cc(mi,mj) = r_dfTC_cc(mi,mj)/(sqrt(sum(dfTC_terms_mnow_cc(mi,mj,:))*sum(dfTC_terms_mnext_cc(mi,mj,:))));
      
      if flagPlot
      scatter(y_fTC_cc(mi,mj),y_dfTC_cc(mi,mj),10,'DisplayName',[ 'Runs: ' num2str(FIR_Runs) ', Config: ' num2str(FIR_Config)]);
      end
    end
    [max_y_fTC(mi),max_idx(mi)] = max(y_fTC_cc(mi,:));
    y_dfTC_max_y_fTC(mi) = y_dfTC_cc(mi,max_idx(mi));
    
    
    end
    
    if flagPlot    
    
        for mi   = 1:length(FIR_Runs_Possibilities)
            FIR_Runs = FIR_Runs_Possibilities(mi);
    % % --------------------------------------------    
    %     % Kurven ohne Farbverlauf
            plot(y_fTC_cc(mi,:),y_dfTC_cc(mi,:),...
                'Color',[mi/length(FIR_Runs_Possibilities) 1-mi/length(FIR_Runs_Possibilities) 0.5],...
                'DisplayName',[ 'Runs: ' num2str(FIR_Runs)]);
    %     % Kurven mit Farbverlauf
    %     hold on
    %     for nn=1:length(y_fTC(:,1))-1
    %         p=plot(y_fTC(mi,nn:nn+1,1),y_dfTC(mi,nn:nn+1,1),'Color',[1-nn/length(y_fTC(:,1)) 0.5 nn^20/length(y_fTC(:,1))^20]);
    %     end  
    %     hold off 
    % % --------------------------------------------        

        end % mi = 1:length(FIR_Runs_Possibilities)
    

        plot(max_y_fTC,y_dfTC_max_y_fTC,'r')
        legend show
    end % flagPlot
    
    
    for i=1:length(y_dfTC_max_y_fTC(:))
        faktor(i) = max_y_fTC(i)* y_dfTC_max_y_fTC(i);
    end
    
    [~,idx_runs_final] = min(faktor);
    FIR_Runs_opt = FIR_Runs_Possibilities(idx_runs_final);
    fprintf(['\nEs wurden aus den gegebenen Möglichkeiten ' num2str(FIR_Runs_Possibilities(idx_runs_final)) ...
             ' FIR-Filterruns als optimal identifiziert!\n']); 
     
    
    FIR_Config_y_max = FIR_Config_Possibilities(max_idx(idx_runs_final));
%     fprintf(['\nBei dieser Runanzahl wurde die maximale y-Übereinstimmung für ' num2str(FIR_Config_y_max) ...
%              ' FIR-Config erzielt.\n']);     
    
   %% Cross correlation - Correlate Runs  
   
    if flagPlot     
    figure('Name','Cross-Correlation Parameters - correlate Runs')
    hold on
    end
    
    for mj   = 1:length(FIR_Config_Possibilities)
    for mi   = 1:length(FIR_Runs_Possibilities)-1
        FIR_Runs = FIR_Runs_Possibilities(mi);
        FIR_Config = FIR_Config_Possibilities(mj);
      for k = 1:lenghtFilteredSignal
         fTC_terms_rc(      mi,mj,k) =  fTC(mi,mj,  k) *  fTC(mi+1,mj,k);
        dfTC_terms_rc(      mi,mj,k) = dfTC(mi,mj,  k) * dfTC(mi+1,mj,k);
         fTC_terms_mnow_rc( mi,mj,k) =  fTC(mi,  mj,k)^2.;
         fTC_terms_mnext_rc(mi,mj,k) =  fTC(mi+1,mj,k)^2.;
        dfTC_terms_mnow_rc( mi,mj,k) = dfTC(mi,  mj,k)^2.;
        dfTC_terms_mnext_rc(mi,mj,k) = dfTC(mi+1,mj,k)^2.;
      end
      % cross-correlation function (Eq. 4)
      r_fTC_rc( mi,mj) = sum( fTC_terms_rc(mi,mj,:));
      r_dfTC_rc(mi,mj) = sum(dfTC_terms_rc(mi,mj,:));
      % cross-corellation parameter (Eq. 5)
      y_fTC_rc( mi,mj) = r_fTC_rc( mi,mj)/(sqrt(sum( fTC_terms_mnow_rc(mi,mj,:))*sum( fTC_terms_mnext_rc(mi,mj,:))));
      y_dfTC_rc(mi,mj) = r_dfTC_rc(mi,mj)/(sqrt(sum(dfTC_terms_mnow_rc(mi,mj,:))*sum(dfTC_terms_mnext_rc(mi,mj,:))));
      
      if flagPlot  
          scatter(y_fTC_rc(mi,mj),y_dfTC_rc(mi,mj),10,'DisplayName',[ 'Runs: ' num2str(FIR_Runs) ', Config: ' num2str(FIR_Config)]);
      end % flagPlot  
    end % mi   = 1:length(FIR_Runs_Possibilities)-1
    

    end % mj   = 1:length(FIR_Config_Possibilities)
    
    
    if flagPlot 
    for mj   = 1:length(FIR_Config_Possibilities)
        FIR_Config = FIR_Config_Possibilities(mj);
% % --------------------------------------------    
%     % Kurve ohne Farbverlauf
        plot(y_fTC_rc(:,mj),y_dfTC_rc(:,mj),...
            'Color',[mj/length(FIR_Config_Possibilities) 1-mj/length(FIR_Config_Possibilities) 0.5],...
            'DisplayName',[ 'Config: ' num2str(FIR_Config)]);
%     % Kurve mit Farbverlauf
%     hold on
%     for nn=1:length(y_fTC(:,1))-1
%         p=plot(y_fTC(mi,nn:nn+1,1),y_dfTC(mi,nn:nn+1,1),'Color',[1-nn/length(y_fTC(:,1)) 0.5 nn^20/length(y_fTC(:,1))^20]);
%     end  
%     hold off 
% % --------------------------------------------        
    
    end % mj = 1:length(FIR_Config_Possibilities)-1
    if idx_runs_final<length(FIR_Runs_Possibilities)
        plot(y_fTC_rc(idx_runs_final,:),y_dfTC_rc(idx_runs_final,:),'r')
    else
        idx_runs_final = idx_runs_final - 1;
        warning(['idx_runs_final wurde für ' num2str(FIR_Runs_Possibilities(idx_runs_final)) ...
                 ' FIR-Filterruns durchgeführt, da die automatische FIR_Config ' ... 
                 'Berechnung nicht mit der maximalen FIR-Filterrunanzahl funktioniert. ' ...
                 'Evtl Obergrenze der FIR_Runs erhöhen.'])
    end
    legend show
    end % flagPlot 
    
    
%     % Mittelwert, Standardabweichung und halbe Standardabweichung einzeichnen
%     for tt=1: length(y_dfTC_rc(:,1))
%         stdDev(tt) = std(y_dfTC_rc(tt,:));
%         mean_y_dfTC_rc(tt,1:length(y_fTC_rc(1,:)))=mean(y_dfTC_rc(tt,:));
%         plot(y_fTC_rc(tt,:),mean_y_dfTC_rc(tt,:),':g')
%         obereGrenze(tt,1:length(y_fTC_rc(1,:)))=mean(y_dfTC_rc(tt,:))+stdDev(tt);
%         untereGrenze(tt,1:length(y_fTC_rc(1,:)))=mean(y_dfTC_rc(tt,:))-stdDev(tt);
%         stdDevhalb(tt)=stdDev(tt)/2;
%         obereGrenzehalb(tt,1:length(y_fTC_rc(1,:)))=mean(y_dfTC_rc(tt,:))+stdDevhalb(tt);
%         untereGrenzehalb(tt,1:length(y_fTC_rc(1,:)))=mean(y_dfTC_rc(tt,:))-stdDevhalb(tt);
%         
%         plot(y_fTC_rc(tt,:),obereGrenze(tt,:),'g')
%         plot(y_fTC_rc(tt,:),untereGrenze(tt,:),'g')
%         plot(y_fTC_rc(tt,:),obereGrenzehalb(tt,:),'b')
%         plot(y_fTC_rc(tt,:),untereGrenzehalb(tt,:),'b')
%     end
        

    % Ich betrachte nur y-Werte die innerhalb der halben Standabweichung liegen, um
    % Ausreisser rauszuschmeissen:
    stdDev = std(y_dfTC_rc(idx_runs_final,:));
    stdDevhalb=stdDev/2;
    obereGrenzehalb=mean(y_dfTC_rc(idx_runs_final,:))+stdDevhalb;
    untereGrenzehalb=mean(y_dfTC_rc(idx_runs_final,:))-stdDevhalb;

    transvektor = y_dfTC_rc(idx_runs_final,:);
    idx_transvektor=find(transvektor>obereGrenzehalb & transvektor<untereGrenzehalb);
    y_fTC_rc_corrected(:)=y_fTC_rc(idx_runs_final,:);
    y_fTC_rc_corrected(idx_transvektor)=0;
    [~,idx_config_final] = max(y_fTC_rc_corrected);

    FIR_Config_opt = FIR_Config_Possibilities(idx_config_final);
    fprintf(['\nEs wurden aus den gegebenen Möglichkeiten eine FIR-Config von ' ...
                num2str(FIR_Config_Possibilities(idx_config_final)) ' als optimal identifiziert!\n']); 
   

            
multiWaitbar('CLOSEALL');           
end




%% Zusatzstuff

%     configNameBuffer = [ nameDataset '_Config' ];      %Create full config file name string from the particular handle
%     configLoadHandle = str2func(configNameBuffer);    %Create function handle from name string
%     [~,~,~,~,~,~,~,~                            , ...
%     Timestep_Reduction                          , ...
%     Temperature_Zero                            , ...
%     ~,~,~,~,~,~,~,~                             , ...
%     ~                                           , ...
%     ~                                           , ...
%     ~                                           , ...
%     ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~  ]    ... 
%     = configLoadHandle();
%
%    % Remove 50 Hz
%   Fs=1/(Time(2)-Time(1)); % Hz, sampling rate
%   wo =  50/(Fs/2);
%   
%   bw1 =   10/(Fs/2);
%   [B1,A1] = iirnotch(wo,bw1);
%     o50(1,:)=filtfilt(B1,A1,rawTcSignal(:,2)); 
%
%
%     figure
%     hold
%     plot(TC_mit50Hz,'DisplayName','TC mit50Hz')
%     plot(o50(1,1:Timestep_Reduction:end),'DisplayName','o50')
%     plot(bufferA(1,1:Timestep_Reduction:end),'DisplayName','buffer A')
%     plot(o50_double,'DisplayName','o50 double')
%     legend show
% 
%     figure;
%     deltaT         = Time(1+Timestep_Reduction)-Time(1); 
%     
%     next_FFT       = 2^nextpow2(length(TC_mit50Hz)-1);
%     Domain         = fft(TC_mit50Hz,next_FFT)/length(TC_mit50Hz);
%     Magnitude      = 2*abs(Domain(1,1:next_FFT/2+1));
%     Frequency_Axis = linspace(0,1,next_FFT/2 + 1) / (2 * deltaT);
%     semilogy(Frequency_Axis,Magnitude);
%     clear Domain Magnitude Frequency_Axis next_FFT 
%     hold
% 
%     next_FFT       = 2^nextpow2(length(o50)-1);
%     
%     Domain         = fft(o50,next_FFT)/length(o50);
%     Magnitude      = 2*abs(Domain(1,1:next_FFT/2+1));
%     Frequency_Axis = linspace(0,1,next_FFT/2 + 1) / (2 * deltaT);
%     semilogy(Frequency_Axis,Magnitude);
%     clear Domain Magnitude Frequency_Axis 
%     
%     Domain         = fft(o50_double,next_FFT)/length(o50_double);
%     Magnitude      = 2*abs(Domain(1,1:next_FFT/2+1));
%     Frequency_Axis = linspace(0,1,next_FFT/2 + 1) / (2 * deltaT);
%     semilogy(Frequency_Axis,Magnitude);
%     clear Domain Magnitude Frequency_Axis   