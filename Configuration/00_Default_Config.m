function [ nameCalFolder        , ...
           nameMeasFolder       , ...
           Frequency_Analysis   , ...
           Pre_Calc_Calibration , ...
           Pre_Calc_Measurement , ...
           Sensor               , ...
           Penetration_Time     , ...
           Future_Time_Window   , ...
           Timestep_Reduction_C , ...
           Temperature_Zero_C   , ...
           Auto_Param_Finder    , ...
           Max_Parameter_n_o    , ...
           Max_Parameter_d_o    , ...
           Parameter_Config_n   , ...
           Parameter_Config_d   , ...
           Manual_n_o           , ...
           Manual_d_o           , ...
           Filter_Config        , ...
           Cut_Off_Frequency    , ...
           FIR_Config           , ...
           FIR_Runs             , ...
           Plot_Configuration   , ...
           Solver_Configuration , ...
           analyzeIR            , ...
           parameterFile        , ...
           flagNL               , ...
           currentTemp          , ...
           pastTemps            , ...
           chemieIdentifikation , ...
           chemieSolver         , ...
           nPlat_SPE            , ...
           Timestep_Reduction_M , ...
           Temperature_Zero_M   , ...
           P4,P5                , ...
           P6,P7,P8,P9,P10    ] = ...
           00_Default_Config()
              %  /\ To be adjusted for each configuration


%% Path Configuration
nameCalFolder    = '...';
nameMeasFolder   = '...'; %Handle to identify different measurements in one folder
parameterFile    = '' ; %z.b. \HENDERSON_PYRAMIDE\NISI_HENDERSON_PYRAMIDE_Parameter_var_T.mat;     

%% Program Control:
Frequency_Analysis  = 0 ; %Activates Frequency Domain Visualisation                      
analyzeIR           = 0 ;  %|Suppresses inversion for faster calculation;
                          %|Additionally plots temperatur difference
                          %|between Simulated vs. Calibration Temperature

flagNL                = false ; %Calculation of nonlinear System Impulse Response            
currentTemp           = false ; %Calculation with Impulse Response of current Temperature    
pastTemps             = false ; %Calculation with Impulse Responses of past Time steps       

chemieIdentifikation  = false ; %Activates chemical terms in identification                  
chemieSolver          = false ; %Activates Chemical solver                                   

%% Program Configuration                                                                    
                                                                                              
Pre_Calc_Calibration  = 1  ; %|Perform input data calculations 
Pre_Calc_Measurement  = 1 ; %|>> Filters and reduction etc     
                                                                                              
                                                                                              
Sensor             = [ 0    % |Analyse a specific sensor / surface     
                       0 ]; % |combination                             
Penetration_Time   = 0.02;  % |Thermal penetration time in s                                 
Future_Time_Window = 0.15;  % |Width in s of the analysation window counting from          
                            % |Penetration_Time , if =0 -> Only the timestep                 
                            % |at Penetrationtime + 1*timestep is considered                  
Timestep_Reduction_C = 110  ; % |Only every Timestep_Reduction timestep is used               
Timestep_Reduction_M = 10  ; % ^same for Measurement                
Temperature_Zero_C   = 10  ; % Time in s used to equate zero level temperature  
Temperature_Zero_M   = 10  ; % ^same for Measurement    

Auto_Param_Finder  = 0    ; % Activation of the automatic parameter finder                   
Max_Parameter_n_o  = 6    ; % Amount of Parameters checked for viability                     
Max_Parameter_d_o  = 8    ; % Amount of Parameters checked for viability                     
Parameter_Config_n = 1    ; % |Configuration of the used parameters of the transfer          
Parameter_Config_d = 1    ; % |function -> 1 == D^(j/2) standard derivatives row             
                            % |For other configurations look @                               
                            % |<NISIPath>/Code/Subroutines/parameter_conf.m                  
%--------------------------------------------------------------------------                        
%Parameters for identification of the transfer function if                                         
%Auto_Param-Finder = false                                                                         
%First row defines the used derivatives                                                            
%Second row  = 0 >> Parameter defined by least squares algorithm                                   
%Second row ~= 0 >> Parameter is fixed on the given value                                          
%--------------------------------------------------------------------------                        
                                                                                                    
%Flux Terms                                                                                        
Manual_n_o =   [  -0.5  0  0.5 1 1.5                                                              
                   0    0  0   0 0  ] ;                                                            
				                                                                                     
				                                                                                     
% Temperature Terms                                                                                
Manual_d_o =   [  0  0.5  1 1.5                                                                  
                  0  0    0 0    ] ;                                                            
                                                                                                    
                                                                                                    
%% Filter Configuration                                                                           
                                                                                                    
Filter_Config = [ 4    % Calibration Temperature Data  | 0 >>> No Filter                           
                  7    % Calibration Heat Flux Data    | 1 >>> FIR Filter                          
                  4    % Measurement Temperature Data  | 2 >>> Running Gauss Filter                
                  0    % Impulse Response              | 3 >>> RGF + FIR                           
                  0 ]; % Inverted Heat Flux            | 4 >>> FIR + RGF                           
                       %                               | 5 >>> Reduction                           
                       %                               | 7 >>> Sharp Edge Filter (for NISI pulses) 
                                                                                                    
%Gauss-Filter Frequencies:                                                                         
Cut_Off_Frequency = [ 20    %Gauss-Filter cut off frequency (Hz) |Calibration T                     
                      0    %Gauss-Filter cut off frequency (Hz) |Calibration Q                     
                      50    %Gauss-Filter cut off frequency (Hz) |Measurement T                     
                      0    %Gauss-Filter cut off frequency (Hz) |Pulse Resp.                       
                      0 ]; %Gauss-Filter cut off frequency (Hz) |generic                           
                                                                                                    
%FIR-Filter:                                                                                       
%See <NISIPath>/Code/Subroutines/filter_r_fir.m for configuration types                            
FIR_Config    = [ 300    % Calibration Temperature Data                                              
                  1    % Calibration Heat Flux Data                                                
                  5    % Measurement Temperature Data                                              
                  3    % Impulse Response                                                          
                  3 ]; % Inverted Heat Flux]                                                                  
                                                                                                    
FIR_Runs      = [ 5    % Calibration Temperature Data                                             
                  10    % Calibration Heat Flux Data                                              
                  5    % Measurement Temperature Data                                              
                  1    % Impulse Response                                                          
                  5 ]; % Inverted Heat Flux]                                                                  
              
% Sharp Edge Filter
% see <NISIPath>/Code/Subroutines/filter_sharpPulseEdges_NISI.m
nPlat_SPE     = 3;     % Number of Plateaus the SPE is reducing to

%% Plot Configuration                                                                             
                       % C  M    % Column 1 of Plot_Configuration for Calib, 2 for Meas           
Plot_Configuration   = [ 1  0    % Calibration/Simulation Temperature    |Manual configuration     
                         1  0    % Calibration/Inverted Heat Flux        |of plots                 
                         0  1    % Measurement Heat Flux (Inverted)      |0 - deactivated          
                         1  1    % Impulse Response                      |1 - activated            
                         0  0    % FFT - Calibration Temperature                                   
                         0  0    % FFT - Calibration Heat Flux                                     
                         0  0    % FFT - Measurement Temperature                                   
                         0  0    % FFT - Inverted Heat Flux                                     
                         0  0 ]; % FFT - Impulse Response                                          

%% Solver Configuration                                                                           
Solver_Configuration = 1;        % 0 - phased Van Cittert                                          
                                 % 1 - sequential function estimation                              
                                                                                                    
%% Dummy variables for future use                           
P4 = '';                          
P5 = '';                          
P6 = '';                          
P7 = '';                          
P8 = '';                          
P9 = '';                          
P10= '';                          
end