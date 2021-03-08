function  [Output_Vector,Time] = filter_control(Input_Vector, ...
           Delta_t,fcHandle,Sensor, Timestep_Reduction, ... 
           Temperature_Zero,Filter_Config,Cut_Off_Frequency,FIR_Config,FIR_Runs,nPlat_SPE) 
%--------------------------------------------------------------------------
% fcHandel values for filter input data:
% 1 >>> Calibration temperature or pressure     -> Zero leveling is applied
% 2 >>> Calibration heat flux                   -> Zero leveling is applied
% 3 >>> Measurement temperature or pressure     -> Zero leveling is applied    
% 4 >>> Impulse Response
% 5 >>> Inverted heat flux
%--------------------------------------------------------------------------


switch fcHandle
    case {1,2} % Filtering calibration temperature/pressure or Calibration heat flux
        StepsZero = double(int32(Temperature_Zero / Delta_t));
        if Sensor(1) ~= 0 && Sensor(2) ~= 0
            Mean_Zero_Level  = sum(Input_Vector(Sensor(1),Sensor(2),1:StepsZero)) / StepsZero;
            Input_Vector(Sensor(1),Sensor(2),:) = Input_Vector(Sensor(1),Sensor(2),:) - Mean_Zero_Level;
            Input_Vector(Sensor(1),Sensor(2),1:StepsZero) = 0;

            filterInput(1,:) = Input_Vector(Sensor(1),Sensor(2),:);
            
            [filterOutput,Time] = chooseFilter(...
                    filterInput,Timestep_Reduction,Filter_Config,FIR_Config,FIR_Runs,...
                    fcHandle,Delta_t,Cut_Off_Frequency,nPlat_SPE);

            filterOutput(1,1:floor(StepsZero/Timestep_Reduction)) = 0;
                
            Output_Vector(1,1,:) = filterOutput(1,:); % Transfer_Vector is already reduced

        else % if Sensor(1) == 0 or Sensor(2) == 0
            for i=1:length(Input_Vector(:,1,1))
            for j=1:length(Input_Vector(1,:,1))
                Mean_Zero_Level         = sum(Input_Vector(i,j,1:StepsZero)) / StepsZero;
                Input_Vector(i,j,1:end) = Input_Vector(i,j,1:end) - Mean_Zero_Level;
                Input_Vector(i,j,1:StepsZero) = 0;

                filterInput(1,:) = Input_Vector(i,j,:);

                [filterOutput,Time] = chooseFilter(...
                        filterInput,Timestep_Reduction,Filter_Config,FIR_Config,FIR_Runs,...
                        fcHandle,Delta_t,Cut_Off_Frequency,nPlat_SPE);

                filterOutput(1,1:floor(StepsZero/Timestep_Reduction)) = 0;
                    
                Output_Vector(i,j,:) = filterOutput(1,:); % Transfer_Vector is already reduced
            end
            end

        end % if Sensor(1) == 0 or Sensor(2) == 0
 
 
    case 3 % Filtering Measurement Data
        StepsZero                        = double(int32(Temperature_Zero / Delta_t));
        if Sensor(1) ~= 0   
            Mean_Zero_Level         = sum(Input_Vector(Sensor(1),1:StepsZero)) / StepsZero;
            Input_Vector(Sensor(1),:)   = Input_Vector(Sensor(1),:) - Mean_Zero_Level;
            Input_Vector(Sensor(1),1:StepsZero) = 0;

            filterInput(1,:)    = Input_Vector(Sensor(1),:);
            
            [filterOutput,Time] = chooseFilter(...
                filterInput,Timestep_Reduction,Filter_Config,FIR_Config,FIR_Runs,...
                fcHandle,Delta_t,Cut_Off_Frequency,nPlat_SPE);

              filterOutput(1,1:floor(StepsZero/Timestep_Reduction)) = 0;
              
            Output_Vector(1,:)  = filterOutput(1,:); % Transfer_Vector is already reduced
            
        else
        
            for i=1:length(Input_Vector(:,1))
                Mean_Zero_Level         = sum(Input_Vector(i,1:StepsZero)) / StepsZero;
                Input_Vector(i,1:end)   = Input_Vector(i,1:end) - Mean_Zero_Level;
                Input_Vector(i,1:StepsZero) = 0;
                
                filterInput(1,:) 	= Input_Vector(i,:);
                
                [filterOutput,Time] = chooseFilter(...
                    filterInput,Timestep_Reduction,Filter_Config,FIR_Config,FIR_Runs,...
                    fcHandle,Delta_t,Cut_Off_Frequency,nPlat_SPE);
                
                filterOutput(1,1:floor(StepsZero/Timestep_Reduction)) = 0;
                
                Output_Vector(i,:)  = filterOutput(1,:); % Transfer_Vector is already reduced
            end
        end
 
    case 4 % Filtering Impulse Response 
        warning(['Data Reduction does not apply to the Impulse Response.\nAfter its calulation, the ' ...
                 'Impulse Response is already as long as the Temperature and Heat Flux Vector.\n'...
                 'The reduction of the Impulse Response has therefore been turned off.'], '');
        if Sensor(1) ~= 0 && Sensor(2) ~= 0
            filterInput(1,:) = Input_Vector(Sensor(1),Sensor(2),:); 
            
            [filterOutput,Time] = chooseFilter(...
                filterInput,1,Filter_Config,FIR_Config,FIR_Runs,... % Timestep_Reduction set to 1, because IR is not to be reduced!
                fcHandle,Delta_t,Cut_Off_Frequency,nPlat_SPE);
            
            filterOutput(1,1) = 0;
            Output_Vector(Sensor(1),Sensor(2),:) = filterOutput(1,:); 
        else
            for i=1:length(Input_Vector(:,1,1))
            for j=1:length(Input_Vector(1,:,1))
                filterInput(1,:)        = Input_Vector(i,j,:); 
                
                [filterOutput,Time]     = chooseFilter(...
                    filterInput,1,Filter_Config,FIR_Config,FIR_Runs,... % Timestep_Reduction set to 1, because IR is not to be reduced!
                    fcHandle,Delta_t,Cut_Off_Frequency,nPlat_SPE);
                
                filterOutput(1,1) = 0;
                Output_Vector(i,j,:)    = filterOutput(1,:);  
            end
            end
        end
                
    case 5 % Filtering Generic Input Vector <-> Inverted Heat Flux
        switch Filter_Config(fcHandle)
            case {1,2,3,4} 
                Timestep_Reduction = 1; % Only in case Filter_Config = 5 a (second) reduction is performed
        end  

        filterInput = Input_Vector;
        
        [filterOutput,Time] = chooseFilter(...
                filterInput,Timestep_Reduction,Filter_Config,FIR_Config,FIR_Runs,...
                fcHandle,Delta_t,Cut_Off_Frequency,nPlat_SPE);
            
        Output_Vector = filterOutput; 

end

 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Subfunction %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Legit Filter_Config values:
% 0 >>> no filter            | no timestep increase
% 1 >>> FIR Filter           | timestep increase for 'Handle' 1/2/3
% 2 >>> Running Gauss Filter | timestep increase for 'Handle' 1/2/3
% 3 >>> RGF + FIR            | timestep increase for 'Handle' 1/2/3
% 4 >>> FIR + RGF            | timestep increase for 'Handle' 1/2/3
% 5 >>> Postprocess Inv HF
%--------------------------------------------------------------------------

 function [filterOutput,Time] = chooseFilter(...
            filterInput,Timestep_Reduction,Filter_Config,FIR_Config,FIR_Runs,...
            fcHandle,Delta_t,Cut_Off_Frequency,nPlat_SPE)
    
    % Filter Transfer_Vector / Data Vector
    switch Filter_Config(fcHandle) 
        case 0
            filterOutput = filterInput;
            Timestep_Reduction = 1; % <- no reduction (of time vector)

        case 1                     
           	[ filterOutput ]  = filter_r_fir_NISI(filterInput, ...
                                      FIR_Config(fcHandle),FIR_Runs(fcHandle),...
                                      Timestep_Reduction,fcHandle);

        case 2               
           	[ filterOutput ]  = filter_r_gauss(filterInput, ...
                                      Cut_Off_Frequency(fcHandle),Delta_t,...
                                      Timestep_Reduction,fcHandle);

        case 3
        	[ bufferA ]   = filter_r_gauss(filterInput, ...
                                      Cut_Off_Frequency(fcHandle),Delta_t,...
                                      1,fcHandle); % 1 because no timestep reduction here, but during fir run!
        	[ filterOutput ]  = filter_r_fir_NISI(bufferA, ...
                                      FIR_Config(fcHandle),FIR_Runs(fcHandle),...
                                      Timestep_Reduction,fcHandle);

        case 4             
        	[ bufferA ]   = filter_r_fir_NISI(filterInput, ...
                                      FIR_Config(fcHandle),FIR_Runs(fcHandle),...
                                      1,fcHandle); % 1 because no timestep reduction here, but during gauss run!
          	[ filterOutput ]  = filter_r_gauss(bufferA, ...
                                      Cut_Off_Frequency(fcHandle),Delta_t,...
                                      Timestep_Reduction,fcHandle);        

        case 5
            filterOutput = filterInput(1,1:Timestep_Reduction:end); 
            
        case 7
            [ filterOutput ]  = filter_sharpPulseEdges_NISI(filterInput,Timestep_Reduction,nPlat_SPE);

    end
    
    % Reduce time vector to match Transfer_Vector
    if ismember(fcHandle,[1,2,3]) % Wenn fcHandle = 1, 2 oder 3, also Cal-TC, Cal-HF oder Meas-TC
        for k = 1:length(filterOutput(1,:))
            Time(k,1) = (k - 1) * Delta_t * Timestep_Reduction;
        end 
    else
            Time    = '';    
    end
     
 end

end