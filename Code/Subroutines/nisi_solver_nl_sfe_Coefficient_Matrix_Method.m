function [Msrmnt_Heat_Flux] = nisi_solver_nl_sfe(Msrmnt_Temperature,Coefficient_Matrix,Delta_t,A8,A9)
%% Nonlinear solver via Coefficient Matrix Method
Penetration_Steps                         = double(int32(A8/Delta_t));
Future_Time_Window                        = double(int32(A9/Delta_t));
Datapoints                                = length(Msrmnt_Temperature(1,:));
Msrmnt_Points                             = length(Msrmnt_Temperature(:,1));

Percent_Done = 0;  
Percent_Per_Run = 100/Datapoints;
Msrmnt_Heat_Flux(1:Msrmnt_Points,1:Datapoints) = 0;

for current_Timestep = 1:length(Msrmnt_Temperature(1,:)) - Penetration_Steps - Future_Time_Window   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Coefficient Matrix Form:
    % Coefficient_Matrix  = [ 1 % Sensor
    %                         2 % Surface
    %                         3 % current timestep
    %                         4 % index of the phased impulse response
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for handle_Sensor = 1:Msrmnt_Points
        Induced_Temperature(1:1+Future_Time_Window) = 0;
        for virtual_time_steps = 1:1+Future_Time_Window
            for handle_Surface = 1:Msrmnt_Points
                Delayed_Influence = sum(Coefficient_Matrix ...
               (handle_Sensor,handle_Surface,current_Timestep, ...
                1:current_Timestep+virtual_time_steps-1));
                Induced_Temperature = Induced_Temperature + Delayed_Influence;
            end
            delta_Temperature(handle_Sensor,virtual_time_steps) = Msrmnt_Temperature( ...
            handle_Sensor,current_Timestep+Penetration_Steps+virtual_time_steps-1) -  ...
            Induced_Temperature(1,virtual_time_steps);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creation of the solver matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for handle_Sensor     = 1:Msrmnt_Points 
    for handle_Surface    = 1:Msrmnt_Points  
    for future_Time_steps = 1:Future_Time_Window+1
        Solv_Mat(handle_Surface,handle_Sensor,future_Time_steps) =  ...
        sum(Coefficient_Matrix(handle_Surface,handle_Sensor,        ...
        current_Timestep,current_Timestep:current_Timestep +        ...
        Penetration_Steps + Future_Time_Window-1));
    end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Creation of the solver matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    for future_Time_step = 1:Future_Time_Window+1
        Heat_Flux_Buffer(:,future_Time_step) = ...  
        Solv_Mat(:,:,future_Time_step)\delta_Temperature(:,future_Time_step);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % Approximating the final solution by weighted error handling of the 
    % Impulse Response
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for handle_Surface = 1:Msrmnt_Points
          Msrmnt_Heat_Flux(handle_Surface,current_Timestep) = ...
          sum(Heat_Flux_Buffer(handle_Surface,:))/(1 + Future_Time_Window );
          Coefficient_Matrix(handle_Surface,:,current_Timestep,:) = ...
          Coefficient_Matrix(handle_Surface,:,current_Timestep,:).* ...
          Msrmnt_Heat_Flux(handle_Surface,current_Timestep);
      end
      Percent_Done = Percent_Done + Percent_Per_Run;
      fprintf([ '\n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n \n' ...
                     'Inversion Completed: \n' ...
                       num2str(Percent_Done) ' Percent Done \n' ]) 
end
