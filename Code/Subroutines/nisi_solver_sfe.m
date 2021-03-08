function [Msrmnt_Heat_Flux] = nisi_solver_sfe(Msrmnt_Temperature,Impulse_Response,...                                    Delta_t,Penetration_Time,Future_Time_Window)%% Solver using sequential function estimation (Beck's) Penetration_Steps                         = double(int32(Penetration_Time/Delta_t));Future_Time_Steps                         = double(int32(Future_Time_Window/Delta_t));Datapoints                                = length(Msrmnt_Temperature(1,:));Msrmnt_Points                             = length(Msrmnt_Temperature(:,1));Impulse_Response(:,:,1:Penetration_Steps) = 0;% Impulse_Response                          = Impulse_Response * Delta_t; % Korrektur zu Impulsantwort, eines Dirac Impulses von 1 W/m^2 (von 1 J/m^2)kk=0;Msrmnt_Heat_Flux=zeros(Msrmnt_Points,Datapoints);Coeff_Mat=zeros(Msrmnt_Points,Msrmnt_Points,1);delta_Temperature=zeros(Msrmnt_Points,1);  %------------------------------------------------------------------------  % Creation of the coefficient matrix  %------------------------------------------------------------------------for handle_Sensor = 1:Msrmnt_Points for handle_Surface = 1:Msrmnt_Points        Coeff_Mat(handle_Surface,handle_Sensor,1) =  ...    sum(Impulse_Response(handle_Surface,handle_Sensor,1:1 + Penetration_Steps + Future_Time_Steps));endendfor currentTimeStep = 1:length(Msrmnt_Temperature(1,:)) - Penetration_Steps - Future_Time_Steps  %------------------------------------------------------------------------  % Determine the temperature at timestep induced by heat flux equated  % for previous timesteps  %------------------------------------------------------------------------  for handle_Sensor = 1:Msrmnt_Points         Induced_Temperature = 0;         for handle_Surface = 1:Msrmnt_Points           iA(1,:) = Msrmnt_Heat_Flux(handle_Surface,:);           iB(1,:) = Impulse_Response(handle_Sensor,handle_Surface,:);           ConvBuff = conv(iB,iA);            Induced_Temperature = Induced_Temperature + ConvBuff(1:Datapoints);         end            delta_Temperature(handle_Sensor,1) = Msrmnt_Temperature( ...            handle_Sensor,currentTimeStep+Penetration_Steps+Future_Time_Steps-1) -  ...            Induced_Temperature(1,currentTimeStep+Penetration_Steps+Future_Time_Steps-1);  end  %--------------------------------------------------------------------------  % Solving linear equation System  %  % delta_Temperature(s,x) = Sum( Msrmnt_Heat_Flux(k,i)*Impulse_Response(s,k,i+x)) for all s  %--------------------------------------------------------------------------  Heat_Flux_Buffer(:,1) = ...    Coeff_Mat\delta_Temperature(:,1); for stabilisation_counter = 1:length(Heat_Flux_Buffer(:,1))      if Heat_Flux_Buffer(stabilisation_counter,1) < 0        Heat_Flux_Buffer(stabilisation_counter,1) = 0;     end end  Msrmnt_Heat_Flux(:,currentTimeStep) = Heat_Flux_Buffer(:,1);  progress = currentTimeStep / ...          % progress increases from 0 to 1      (length(Msrmnt_Temperature(1,:))-Penetration_Steps-Future_Time_Steps);  if progress > kk  multiWaitbar('Inverting System',progress);  kk=kk+0.001;                                % with kk the output rate is reduced to every 0.1%  endendmultiWaitbar('Inverting System',1);end