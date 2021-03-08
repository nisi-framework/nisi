function [Msrmnt_Heat_Flux] = nisi_solver_nl_sfe(Msrmnt_Temperature,Impulse_Response,Delta_t,A8,A9)
%% Nonlinear Solver w/o losses
%  Solver using sequential function estimation (Beck's) 
Penetration_Steps                         = double(int32(A8/Delta_t));
Future_Time_Window                        = double(int32(A9/Delta_t));
Datapoints                                = length(Msrmnt_Temperature(1,:));
Msrmnt_Points                             = length(Msrmnt_Temperature(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Impulse_Response in nonlinear is a 4-dimensional matrix:
%Impulse_Response(#1,    % Sensor
%                 #2,    % Surface
%                 #3,    % Index for Temperature at time
%                 #4)    % IR progression as equated for Temperature at #3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%



Percent_Done = 0;  
Percent_Per_Run = 100/Datapoints;

Msrmnt_Heat_Flux(1:Msrmnt_Points,1:Datapoints) = 0;

  %------------------------------------------------------------------------
  % Creation of the coefficient matrix
  %------------------------------------------------------------------------
for handle_Sensor = 1:Msrmnt_Points 
for handle_Surface = 1:Msrmnt_Points    
    Coeff_Mat(handle_Surface,handle_Sensor,1) =  ...
    sum(Impulse_Response(handle_Surface,handle_Sensor,1:1 + Penetration_Steps + Future_Time_Window));
end
end
  
  
  
for current_Timestep = 1:length(Msrmnt_Temperature(1,:)) - Penetration_Steps - Future_Time_Window
  %------------------------------------------------------------------------
  % Determine the temperature at timestep induced by heat flux equated
  % for previous timesteps
  %------------------------------------------------------------------------
  for handle_Sensor = 1:Msrmnt_Points
         Induced_Temperature = 0;
         for handle_Surface = 1:Msrmnt_Points
           iA(1,:) = Msrmnt_Heat_Flux(handle_Surface,:);
           iB(1,:) = Impulse_Response(handle_Sensor,handle_Surface,:);
           ConvBuff = conv(iB,iA); 
           Induced_Temperature = Induced_Temperature + ConvBuff(1:Datapoints);
         end
            delta_Temperature(handle_Sensor,1) = Msrmnt_Temperature( ...
            handle_Sensor,current_Timestep+Penetration_Steps+Future_Time_Window-1) -  ...
            Induced_Temperature(1,current_Timestep+Penetration_Steps+Future_Time_Window-1);
  end

  %--------------------------------------------------------------------------
  % Solving linear equation System
  %
  % delta_Temperature(s,x) = Sum( Msrmnt_Heat_Flux(k,i)*Impulse_Response(s,k,i+x)) for all s
  %--------------------------------------------------------------------------

  Heat_Flux_Buffer(:,1) = ...  
  Coeff_Mat\delta_Temperature(:,1);
  for stabilisation_counter = 1:length(Heat_Flux_Buffer(:,1))
      if Heat_Flux_Buffer(stabilisation_counter,1) < 0
         Heat_Flux_Buffer(stabilisation_counter,1) = 0;
      end
  end
  Msrmnt_Heat_Flux(:,current_Timestep) = Heat_Flux_Buffer(:,1);

  
  Percent_Done = Percent_Done + Percent_Per_Run;
  fprintf([ '\n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n' ...
            'Inversion Completed: \n' ...
             num2str(Percent_Done) ' Percent Done \n' ]) 
end

