function [Msrmnt_Heat_Flux] = nisi_solver_sfe_linearHF(Msrmnt_Temperature,Impulse_Response,...
Delta_t,Future_Time_Window)
%% Solver using sequential function estimation (Beck's)
% f�r konstantes q_m f�r alle future_time_steps
Future_Time_Steps                         = double(int32(Future_Time_Window/Delta_t));
Datapoints                                = length(Msrmnt_Temperature(1,:));
Msrmnt_Points                             = length(Msrmnt_Temperature(:,1));


% StabCount = 0;
kk=0;
Msrmnt_Heat_Flux=zeros(Msrmnt_Points,Datapoints);
Coeff_Mat=zeros(Msrmnt_Points,Msrmnt_Points);
delta_Temperature=zeros(Msrmnt_Points,Future_Time_Steps);
sumOfPhi = zeros(Msrmnt_Points,Msrmnt_Points,Future_Time_Steps);


%------------------------------------------------------------------------
  % Creation of the coefficient matrix
%------------------------------------------------------------------------

% Summe der Impulsantworten f�r future_time_steps
for handle_Sensor = 1:Msrmnt_Points
for handle_Surface = 1:Msrmnt_Points
    for i = 1:Future_Time_Steps
        % Creation of sum [k = 1 -> i] (phi_k)   <- Summe der IR von 1 bis i (siehe Beck Gl. 4.4.34)
        sumOfPhi(handle_Sensor,handle_Surface,i) = sum(Impulse_Response(handle_Sensor,handle_Surface,2:i)); % <- IR(1) = phi_0 in Beck (Gl. 4.4.34) <- wenn man hier nicht i+1 sondern i als hintere Grenze nimmt, kommt das einer Regularisierung nach Beck gleich (siehe Beck Gl. 4.5.28 und Folgeabsatz (alpha ist hier so gro� wie das letzte Summenglied, d.h. i+1 entf�llt) - Eine Bemerkung hierzu: Es ergeben sich bessere �bereinstimmungen der Simulierten Temperatur mit der Messtemperatur f�r i statt i+1)
        
        Coeff_Mat(handle_Surface,handle_Sensor) = Coeff_Mat(handle_Surface,handle_Sensor) + ...
        sumOfPhi(handle_Sensor,handle_Surface,i)^2;    
    end
end
end


for currentTimeStep = 1:length(Msrmnt_Temperature(1,:))  - Future_Time_Steps
%------------------------------------------------------------------------
  % Determine the temperature at timestep induced by heat flux equated
  % for previous timesteps
%------------------------------------------------------------------------
           numeratorSum=zeros(Msrmnt_Points,Msrmnt_Points,Future_Time_Steps); % <- numerator in Beck Gl. 4.4.34
  for handle_Sensor = 1:Msrmnt_Points
         Induced_Temperature = 0;
         for handle_Surface = 1:Msrmnt_Points
           iA(1,:) = Msrmnt_Heat_Flux(handle_Surface,:);
           iB(1,:) = Impulse_Response(handle_Sensor,handle_Surface,:);
           ConvBuff = conv(iB,iA);
           Induced_Temperature = Induced_Temperature +ConvBuff(1:Datapoints);
         end
         
         if currentTimeStep == length(Msrmnt_Temperature(1,:))  - Future_Time_Steps
            Induced_memory(handle_Sensor,:)=Induced_Temperature;
         end

         for i = 1:Future_Time_Steps    % Hier entsteht der Z�hler in Beck Gl. 4.4.34

             delta_Temperature(handle_Sensor,1) = ...  % <- Y-Tdach aus Beck Gl. 4.4.34
                 Msrmnt_Temperature(handle_Sensor,currentTimeStep+i-1) -  ...
                 Induced_Temperature(1,currentTimeStep+i-1);

             for handle_Surface = 1:Msrmnt_Points
                 numeratorSum(handle_Sensor,handle_Surface,i)= ...
                 delta_Temperature(handle_Sensor,1)*sumOfPhi(handle_Sensor,handle_Surface,i);
             end
         end

  end
%--------------------------------------------------------------------------
  % Solving linear equation System
  %
  % delta_Temperature(s,x) = Sum( Msrmnt_Heat_Flux(k,i)*Impulse_Response(s,k,i+x)) for all s
%--------------------------------------------------------------------------

  Heat_Flux_Buffer(:,1) = ...
  Coeff_Mat\sum(sum(numeratorSum,2),3);%delta_Temperature(:,1);

% % negative HF mit 0 �berschreiben:
%  for stabilisation_counter = 1:length(Heat_Flux_Buffer(:,1))  % <- 
%       if Heat_Flux_Buffer(stabilisation_counter,1) < 0
%         Heat_Flux_Buffer(stabilisation_counter,1) = 0;        
%         StabCount = StabCount +1;
%      end
%  end

  Msrmnt_Heat_Flux(:,currentTimeStep) = Heat_Flux_Buffer(:,1);
  
  
  progress = currentTimeStep / ...          % progress increases from 0 to 1
      (length(Msrmnt_Temperature(1,:))-Future_Time_Steps);
  if progress > kk
  multiWaitbar('Inverting System',progress);
  kk=kk+0.001;                                % with kk the output rate is reduced to every 0.1%
  end
end  % for currentTimeStep = 1:length(Msrmnt_Temperature(1,:))  - Future_Time_Steps


% figure('Name','Induced to Msrmnt Temperature');hold on
% plot(Induced_memory(1,:),'b--')
% plot(Msrmnt_Temperature(1,:),'b')
% if Msrmnt_Points > 4
% plot(Induced_memory(5,:),'r--')
% plot(Msrmnt_Temperature(5,:),'r')
% legend ('I1','M1','I5','M5','I3','M3')
% end
% if Msrmnt_Points > 2
% plot(Induced_memory(3,:),'g--')
% plot(Msrmnt_Temperature(3,:),'g')
% legend ('I1','M1','I3','M3')
% end
% hold off


multiWaitbar('Inverting System',1);
end
