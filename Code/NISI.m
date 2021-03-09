function NISI(nameDataset,CorM,varargin)
% Dies ist die Hauptroutine des NISI-Programms
%
% Ursprungsversion von
% Ulf Fuchs
%
% Modified:
% S. Loehle, F. Hufgard
% 
%
% INPUT
%       nameDataset - Welches NISI-Problem (im Data-Folder) soll bearbeitet werden
%       CorM        - Calibrierung oder Messung? - Mögliche Inputstrings: 'C' oder 'M' 
%       varargin (optional)   
%                   - Parameters which should be overwritten in the config file
%                       - Input in pairs with: >name of parameter< and >new value of parameter<
%                       - Possible paramters to overwrite:
%                           nameCalFolder           ('nameCalFolder')
%                           nameMeasFolder          ('nameMeasFolder')
%                           Penetration_Time        ('pt')
%                           Future_Time_Window      ('ft')
%                           Timestep_Reduction_C    ('tsrC')
%                           Temperature_Zero_C      ('tzC')
%                           Timestep_Reduction_M    ('tsrM')
%                           Temperature_Zero_M      ('tzM')
%                   - noPlot  - disables the plot_module
%                   - noClear - disables clearing the base workspace 
%                   - noWB    - disables closing the Waitbar at the end of the NISI code   
%
% OUTPUT (Wird in Workspace geschrieben)
%       Calibration_Temperature or Calibration_Pressure
%       Calibration_Heat_Flux
%       Msrmnt_Temperature or Msrmnt_Pressure
%       Msrmnt_Heat_Flux 
%       Impulse_Response 
%       Time
%       nameDataset
%       nameMeasurment
% 
% 
% Additional files necessary for execution: 
%       Path definition file called: pfade.mat 
%       Config file in config folder (pathNisiFolder\Configuration\)
%       CorM = 'C': Calibration data in file named: 
%                           nameDataset_NISI.mat or
%                           nameDataset_NISI_filtered.mat
%                   _NISI-file contains variables: 
%                           Calibration_Heat_Flux
%                           Calibration_Temperature or Calibration_Pressure 
%                           Time
%                   Location:
%                           pathDataFolder\nameDataset\nameCalFolder
%       CorM = 'M': Measurement data in file named: 
%                           nameDataset_nameMeasFolder_NISI.mat or
%                           nameDataset_nameMeasFolder_NISI_filtered.mat
%                   _NISI-file contains variables: 
%                           S or Msrmnt_Temperature or Msrmnt_Pressure
%                           Time
%                   Location:
%                           pathDataFolder\nameDataset\nameMeasFolder

%---------------------------------------------------------------------------------------
fprintf('NISI Copyright (C) 2021, [Stefan Loehle, Ulf Fuchs, Fabian Hufgard, Christian Duernhofer] \n')
fprintf('This program comes with ABSOLUTELY NO WARRANTY; for details open LICENSE.md. \n This is free software, and you are welcome to redistribute it \n under certain conditions; check LICENSE.md for details.') 


fprintf('\n \n \n \n')
fprintf('Program started \n')
fprintf('--------------------------------------------------------\n')
Timer_Start = cputime;              
tic                                 % take start time (to calculate program run time at the end)
%---------------------------------------------------------------------------------------

if ~sum(strcmp('noClear',varargin)) % Check for noClear option and if not set, then ...
evalin('base','clear')              % clear (Base-)Workspace
end 


%% Path Control
load('pfade.mat');                  % Contains two variables called:  
                                    % pathNisiFolder und pathDataFolder,
                                    % which represent the user specific as strings

% Set TorP as global variable
global TorP %TorP is defined in the parser




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Configuration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% nameDataset contains the string naming the configuration file. Each configuration
% file contains a function of the same name. Creating a generic function
% handle enables the generic call of the configuration function
%---------------------------------------------------------------------------------------

configNameBuffer = [ nameDataset '_Config' ];      %Create full config file name string from the particular handle
configLoadHandle = str2func(configNameBuffer);    %Create function handle from name string
fprintf('Loading configuration file... \n')

try
    [
    nameCalFolder	     , ...
    nameMeasFolder  	 , ...
    Frequency_Analysis	 , ...
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
    ~,~,~,~,~,~          , ... % <- NL parameter
    nPlat_SPE            , ... % See e.g. Leiser_Rundstab_Config.m
    Timestep_Reduction_M , ...
    Temperature_Zero_M   , ...
    ~,~,~,~,~,~,~    ] ... % <- Dummy parameter
    = configLoadHandle();

catch
    error(['Unable to call config file in expected path.\nFile name expected:\n%s\n' ...
          'Folder path expected:\n%s\n\nCheck existance of called config file.\n' ...
          'Check folder paths in your pfade.mat file.'], ...
          [configNameBuffer '.m'], [pathNisiFolder 'Configuration\']);
end

%% Overwrite Config Parameters: Process potential input variables and in case overwrite config
if nargin > 2 
% Overwrite parameters according to NISI-input
i = 1;      
while i < nargin-1
   switch varargin{i} 
       case 'nameCalFolder'                         % if 3rd, 5th, .. NISI-Input is nameCalFolder then
           nameCalFolder = varargin{(i+1)};         % overwrite nameCalFolder with the 4th, 6th, ... Input
       case 'nameMeasFolder'                         % same as above
           nameMeasFolder = varargin{(i+1)};
       case {'pt', 'PT', 'Penetration_Time'}        % same as above
           Penetration_Time = varargin{(i+1)};      
       case {'ftw', 'FTW','ft','FT','Future_Time_Window'}     % same as above
           Future_Time_Window = varargin{(i+1)};
       case {'tsrC', 'TSRC','Timestep_Reduction_C'} % same as above
           Timestep_Reduction_C = varargin{(i+1)};
       case {'tzC', 'TZC','Temperature_Zero_C'}     % same as above
           Temperature_Zero_C = varargin{(i+1)};
       case {'tsrM', 'TSRM','Timestep_Reduction_M'} % same as above
           Timestep_Reduction_M = varargin{(i+1)};
       case {'tzM', 'TZM','Temperature_Zero_M'}     % same as above
           Temperature_Zero_M = varargin{(i+1)};
       case {'noPlot', 'Plotoff', 'no_plot'}        % option to disable plot_module at the end of the code
           Plot_Configuration = Plot_Configuration*0;   % <- this will cause the plot_module to skip all plotting
           i = i-1;         % <- this is required, because noPlot is not a input pair, so i should only be increased by 1 overall
       case {'noClear', 'noWSClear', 'no_clear'}        % option to disable clearing the base workspace at the start of the code
           % No action required, but counter needs to be reduced by 1:
           i = i-1;         % <- this is required, because noClear is not a input pair, so i should only be increased by 1 overall
       case {'noWB', 'noWaitBar', 'noWaitbar'}        % option to disable closing the multiWaitbar at the start of the code
           % No action required, but counter needs to be reduced by 1:
           i = i-1;         % <- this is required, because noClear is not a input pair, so i should only be increased by 1 overall
       otherwise
           error('foo:bar',['Invalid overwrite input parameter: ''' varargin{i} '''\n\n'...
                    'Valid overwrite parameters are:\n'...
                    'nameCalFolder:                 ''nameCalFolder''\n' ...
                    'nameMeasFolder:                ''nameMeasFolder''\n' ...
                    'Penetration_Time:              ''pt''\n' ...
                    'Future_Time_Window:            ''ft''\n' ...
                    'Timestep_Reduction_C:          ''tsrC''\n' ...
                    'Temperature_Zero_C:            ''tzC''\n' ...
                    'Timestep_Reduction_M:          ''tsrM''\n' ...
                    'Temperature_Zero_M:            ''tzM''\n' ...
                    'Disable plot_module:           ''noPlot''\n' ...
                    'Disable workspace clearance:   ''noClear''\n'])
   end
   i = i+2;
end

% Overwrite the Config file in the Config folder
overwriteConfig(nameDataset, varargin{:})

clear i
end % if nargin > 2



%% Choose Timestep_Reduction for Calibration or Measurement
switch CorM 
    case 'C'
        Timestep_Reduction = Timestep_Reduction_C;
        Temperature_Zero = Temperature_Zero_C;
    case 'M'
        Timestep_Reduction = Timestep_Reduction_M;
        Temperature_Zero = Temperature_Zero_M;
end
    


%% Copy existing config-file into Cal or Meas Folder
    configNameBufferOrig  = [pathNisiFolder 'Configuration\'  nameDataset '_Config.m'];
    switch CorM
        case 'C'
            configNameBufferCopy = [pathDataFolder, nameDataset '\' nameCalFolder '\'...
                                    nameDataset '_' nameCalFolder '_C_Config_Copy.m'];
        case 'M'
            configNameBufferCopy = [pathDataFolder, nameDataset '\' nameMeasFolder '\'...
                                    nameDataset '_' nameMeasFolder '_M_Config_Copy.m'];
    end
    copyfile(configNameBufferOrig, configNameBufferCopy, 'f');

    clear configLoadHandle configNameBuffer configNameBufferCopy configNameBufferNew


%% Check for legit configuration
%---------------------------------------------------------------------------------------
 [ Legit ] = exception_handler(nameDataset,nameCalFolder,CorM,Temperature_Zero,...
             Filter_Config,FIR_Runs,Timestep_Reduction);

if Legit == 0                       % Legit = 0 means that the configuration is not legit
    return                          % This will lead to the termination of the program
end

clear Legit
                

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parser settings. Last Input:
% 'C' >> Calibration Data
% 'M' >> Measurement Data
%---------------------------------------------------------------------------------------
% Input data format:
%
% Calibration data file contains (filtered and unfiltered):
% Calibration_Heat_Flux(<Sensor_Number>,<Surface_Number>,<Datavector>)
% Calibration_Temperature(<Sensor_Number>,<Surface_Number>,<Datavector>)
% (same nomenclature for Calibration_Pressure)
% Time(<Datavector>) >>> Time(1) = 0!!!
%
% Measurement data file contains (filtered and unfiltered):
% S(<Sensor_Number>,<Datavector>)  - all surfaceses combined for each sensor
% (same nomenclature for Msrmnt_Temperature and Msrmnt_Pressure)
% Time(<Datavector>) >>> Time(1) = 0!!!
%
% Impulse_Response  file contains:
% H(<Sensor_Number>,<Surface_Number>,<Datavector>)
%---------------------------------------------------------------------------------------

fprintf('Loading data... \n')
    
switch CorM
    case 'C'
        fprintf('->Parsing calibration data \n')  
        [Calibration_Heat_Flux,Calib_TempOrPress,Time,Delta_t,~,~,Temperature_Regime] ...
            = parser(pathDataFolder,nameDataset,nameCalFolder,...
                     '',Pre_Calc_Calibration,'C');  
    case 'M'                      
        fprintf('->Parsing measurement temperature data and impulse response \n')  
        [~,Msrmnt_TempOrPress,Time,Delta_t,Impulse_Response,t_H,Temperature_Regime] ...
            = parser(pathDataFolder,nameDataset,nameCalFolder,...
                     nameMeasFolder,Pre_Calc_Measurement,'M');
end 

fprintf('Data loaded \n \n')


% Set up waitbars
%---------------------------------------------------------------------------------------
if ~sum(strcmp('noWB',varargin))
multiWaitbarControl(0,CorM,Pre_Calc_Calibration,Pre_Calc_Measurement,...
                    Auto_Param_Finder,Filter_Config,analyzeIR,'');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Create TorPsaveName for Calibration and Measurement Data
    if TorP == 'T'
        TorPsaveNameC = 'Calibration_Temperature';
        TorPsaveNameM = 'Msrmnt_Temperature';
    else %if TorP =='p'
        TorPsaveNameC = 'Calibration_Pressure';
        TorPsaveNameM = 'Msrmnt_Pressure';
    end

%---------------------------------------------------------------------------------------
% Filtering Calibration Data
%---------------------------------------------------------------------------------------

if CorM == 'C' && Pre_Calc_Calibration
   fprintf('Filtering input calibration data... \n')
   
 %--------------------------------------------------------------------------------------
 % Nonlinear calibration data
 %--------------------------------------------------------------------------------------
  T_regimes = length(size(Calibration_Heat_Flux));
  format    = size(Calibration_Heat_Flux);
  if T_regimes > 3
     calibration_buffer_Q = Calibration_Heat_Flux;
     calibration_buffer_T = Calib_TempOrPress;
     clear Calibration_Heat_Flux Calib_TempOrPress Temperature_Regime
     for i=1:format(3)
        fprintf(['->Filtering ' TorPsaveNameC ' at temperature level ' num2str(i) '\n'])
        input_Value = permute(calibration_buffer_T(:,:,i,:),[ 1 2 4 3 ]);
        Temperature_Regime(i)          = calibration_buffer_T(1,1,i,1);
       [output_Value,Time] = ...
        filter_control(input_Value,Delta_t,1,Sensor,Timestep_Reduction,Temperature_Zero,...
        Filter_Config,Cut_Off_Frequency,FIR_Config,FIR_Runs,'');
        Calib_TempOrPress(:,:,i,:) = output_Value;
         
        fprintf(['->Filtering calibration heat flux at temperature level ' num2str(i) '\n'])         
        input_Value = permute(calibration_buffer_Q(:,:,i,:),[1 2 4 3]);
       [output_Value,~]   = ...
        filter_control(input_Value,Delta_t,2,Sensor,Timestep_Reduction,Temperature_Zero,...
        Filter_Config,Cut_Off_Frequency,FIR_Config,FIR_Runs,nPlat_SPE);
        Calibration_Heat_Flux(:,:,i,:) = output_Value;

     end
   fprintf('->Saving filtered calibration data \n')
   if TorP == 'T'
       Calibration_Temperature = Calib_TempOrPress;
   else %if TorP =='p'
       Calibration_Pressure = Calib_TempOrPress;
   end
   filepath = [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                nameDataset '_' nameCalFolder '_C_NISI_filtered.mat' ];
   save(filepath,'Calibration_Heat_Flux',TorPsaveNameC,'Time','Temperature_Regime');
   
 %--------------------------------------------------------------------------------------
 % Linear calibration data
 %--------------------------------------------------------------------------------------
  else % -> if T_regimes < 4
      
    fprintf(['->Filtering ' TorPsaveNameC '\n'])  
   [Calib_TempOrPress,Time] = ...
    filter_control(Calib_TempOrPress,Delta_t,1,Sensor,Timestep_Reduction,Temperature_Zero,...
    Filter_Config,Cut_Off_Frequency,FIR_Config,FIR_Runs,'');
    
   fprintf('->Filtering calibration heat flux \n')
   [Calibration_Heat_Flux,~]   = ...
   filter_control(Calibration_Heat_Flux,Delta_t,2,Sensor,Timestep_Reduction,Temperature_Zero,...
   Filter_Config,Cut_Off_Frequency,FIR_Config,FIR_Runs,nPlat_SPE);
 
   fprintf('->Saving filtered calibration data \n') 
   if TorP == 'T'
       Calibration_Temperature = Calib_TempOrPress;
   else % i.e.: if TorP =='p'
       Calibration_Pressure = Calib_TempOrPress;
   end
   filepath = [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                nameDataset '_' nameCalFolder '_C_NISI_filtered.mat' ];
   save(filepath,'Calibration_Heat_Flux',TorPsaveNameC,'Time');
  end % T_regimes > 3
  clear filepath Steps i calibration_buffer_Q calibration_buffer_T *_Value
  fprintf('Input data filtering done\n\n')
  
  % Recalculate Delta_t for possibly reduced Temp and HF vectors 
  Delta_t    = Time(2) - Time(1);
end % CorM == 'C' && Pre_Calc_Calibration

%---------------------------------------------------------------------------------------
% Filtering Measurement Data
%---------------------------------------------------------------------------------------
if CorM == 'M' 
    if Pre_Calc_Measurement
       fprintf('Filtering input measurement data... \n') 
     fprintf('->Filtering measurement temperature \n') 
     [Msrmnt_TempOrPress,Time] = ...
      filter_control(Msrmnt_TempOrPress,Delta_t,3,Sensor,Timestep_Reduction,Temperature_Zero,...
      Filter_Config,Cut_Off_Frequency,FIR_Config,FIR_Runs,'');  


     fprintf('->>Saving filtered measurement data \n')
     if TorP == 'T'
       Msrmnt_Temperature  = Msrmnt_TempOrPress;
     else % i.e.: if TorP =='p'
       Msrmnt_Pressure = Msrmnt_TempOrPress;
     end
     filepath =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                   nameDataset '_' nameMeasFolder '_M_NISI_filtered.mat' ];
     save(filepath,TorPsaveNameM,'Time');
     clear filepath
     fprintf('Input data filtering done\n\n')
     
     % Recalculate Delta_t for possibly reduced Temp and HF vectors 
     Delta_t    = Time(2) - Time(1);
    end % Pre_Calc_Measurement
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% System Identification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch CorM % <- System Identification
case 'C' % <- System Identification
     if Temperature_Regime == 'L' % System is linear
         [Numerator_Order_Vector,Numerator_Parameter_Vector, ...
         Denominator_Order_Vector,Denominator_Parameter_Vector, ...
         Error_Numerator_Vector, Error_Denominator_Vector, ...
         Simulated_TorP_Vector, Impulse_Response ] =  ...
         identification_module(Delta_t,Auto_Param_Finder ,Max_Parameter_n_o,Max_Parameter_d_o, ...
         Parameter_Config_n,Parameter_Config_d,Manual_n_o,Manual_d_o, ...
         Calib_TempOrPress,Calibration_Heat_Flux);
  
        %-----------------------------------------------------------------------------------
        % Filter application for the impulse response if specified in Filter_Config
        %-----------------------------------------------------------------------------------
%
%
% 
        switch Filter_Config(4) % case 0 --> nothing happens
            
         case {1,2,3,4,5}       % cases 1 to 5 --> filter IR     
           % First save raw (unfiltered) IR
           pathIRC =  [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                        nameDataset '_' nameCalFolder '_C_Impulse_Response_Raw.mat' ];   
           H            = Impulse_Response;
           t_H          = Time;
           Delta_t_H    = t_H(2) - t_H(1);
           H_1J = H/Delta_t_H;
           save(pathIRC,'H','t_H','Delta_t_H','H_1J'); 
           clear pathIRC H t_H Delta_t_H
           
           % Then filter IR     
           [Impulse_Response,~] = ...    %% Access to interpolation models for the IR <--??
           filter_control(Impulse_Response,Delta_t,4,Sensor,Timestep_Reduction,Temperature_Zero,...
           Filter_Config,Cut_Off_Frequency,FIR_Config,FIR_Runs,'');
           
         case 6
           Filter_Config(5)              = 0; 
           [Impulse_Response] = ir_module(Impulse_Response,Time);
         case 66
%            Filter_Config(5)              = 0; 
           [Impulse_Response] = ir_module_quad(Impulse_Response,Time);           
        end

 
        %-----------------------------------------------------------------------------------
        % Saving Parameter Data & Impulse Response to files 
        %-----------------------------------------------------------------------------------
        pathParameter = [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                     nameDataset '_' nameCalFolder '_C_NISI_Parameter.mat' ];
        save(pathParameter,'Numerator_Order_Vector','Numerator_Parameter_Vector', ...
                      'Denominator_Order_Vector','Denominator_Parameter_Vector', ...
                      'Error_Numerator_Vector','Error_Denominator_Vector');
                  
        pathIRC =  [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                        nameDataset '_' nameCalFolder '_C_Impulse_Response.mat' ];   
        H        = Impulse_Response;
        t_H      = Time;
        Delta_t_H    = t_H(2) - t_H(1);
        H_1J = H/Delta_t_H;
        save(pathIRC,'H','t_H','Delta_t_H','H_1J'); 
        clear pathParameter pathIRC H t_H Delta_t_H H_1J
     else % Temperature_Regime == 'L'       % here NL code


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    % In contrast to the linear solver, the impulse response is calculated
    % after the measurement temperature vector is available since
    % the measurement temperature changes the impulse response
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ nl_Parameter ] = nl_identification_module(pathDataFolder,nameDataset,Time, ...
                          Auto_Param_Finder ,Max_Parameter_n_o,Max_Parameter_d_o, ...
                          Parameter_Config_n,Parameter_Config_d,Manual_n_o,Manual_d_o, ...
                          Calib_TempOrPress,Calibration_Heat_Flux,format);
      %% In nl_parameter_function_fit weitermachen (Kommentar von U. Fuchs)           
    %    [nl_Function_Parameter, nl_Derivative_Grid] = ...
    %     nl_parameter_function_fit(nl_Parameter,Temperature_Regime);               
       
     end % Temperature_Regime == 'L'

                 
    fprintf('System identification done\n\n')  

    
    
    
case 'M'     % <- System Identification
             % In the measurement case the impulse reponse needs to be given already, ...
             % no identification is performed 
     fprintf('System is already identified \n')
     Calibration_Heat_Flux        = 'Na';  %| Setting a variable to 'Na' means that it is not necessary for 
     Calib_TempOrPress      = 'Na';  %| the further program execution - yet it needs to exist because
     Simulated_TorP_Vector       = 'Na';  %| it may be part of a input for functions later on
     
    %-------------------------------------------------------------------------------
    % Impulse Response for Measurement
    %-------------------------------------------------------------------------------
    
    % Check time vectors of measurement data and Impulse Response    
    if isequal(t_H,Time)            %| Time vectors of IR and measurement data set are equal
                                    %| --> Use IR from Calibration
                                                    
        filepathIRM =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                         nameDataset '_' nameMeasFolder '_M_Impulse_Response.mat' ];
        filepathIRC =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                         nameDataset '_' nameCalFolder '_C_Impulse_Response.mat' ];
        try
            delete(filepathIRM)
        end
        
        copyfile(filepathIRC, filepathIRM, 'f');    % Create Msrmt IR from Cal IR by copying
        open(filepathIRM)
         
        
        
    else  % isequal(t_H,Time)  %| Time vectors of IR and measurement data set are NOT equal
                               %| --> Recalculate IR from Calibration Parameterset for new Time vector

                                            
        filepathIRM =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                         nameDataset '_' nameMeasFolder '_M_Impulse_Response.mat' ];        
        try
            delete(filepathIRM)
        end          

        filepathParaM =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                         nameDataset '_' nameCalFolder '_C_NISI_Parameter.mat' ];
        load(filepathParaM)
        Numerator_Order = Numerator_Order_Vector{1,1};
        Numerator_Parameter = Numerator_Parameter_Vector{1,1};        
        Denominator_Order = Denominator_Order_Vector{1,1};
        Denominator_Parameter = Denominator_Parameter_Vector{1,1};
                                 
        % IR mit altem Zeitschritt (und in der Länge angepassten Zeitvektor) und alten Parametern ausgerechnet
        Delta_t_H = t_H(2) - t_H(1);
        t_H_new = [0:Delta_t_H:Time(end)];
        Dirac_Impulse(1:length(t_H_new),1) = 0;         % Creating numerical 
        Dirac_Impulse(2,1)                                 = 1;    % Dirac-Impulse --> eigentlich [J/m^2], d.h. 1/Delta_t - das muss aber auch beim Invertieren beachtet werden (IR muss *Delta_T gesetzt werden). Die Implementierung von 1/Delta_t muss noch erfolgen.
        Dirac_Impulse_Derivatives  = derivative(Dirac_Impulse, Numerator_Order, Delta_t_H);
        [Response2] = simulation_system(Numerator_Parameter, ...
                                     Denominator_Parameter,Denominator_Order, ...
                                     Delta_t_H,Dirac_Impulse,Dirac_Impulse_Derivatives, 1);         
        
        % Skalierung um Verhältnis der Zeitschrittlängen (was dem Verhältnis der Energien der
        % Dirac-Impulse entspricht)
        Ratio_t_exact = Delta_t/Delta_t_H;
        Response3 = Response2*Ratio_t_exact;
        
        % Zeitschrittanpassung von Response3 auf neue Zeitdiskretisierung
        Response4 = pchip(t_H_new,Response3,Time);
        clear Impulse_Response t_H 
        Impulse_Response(1,1,:) = Response4(:,1);
        H = Impulse_Response;
        t_H = Time;
        Delta_t_H    = t_H(2) - t_H(1);
        H_1J = H/Delta_t_H;                    % Impuls Response for 1 J/m^2 pulse (=H/Delta_t_H) <-> das entspricht eigentliche einer Impulsantwort wie sie definiert ist.
        
            %-----------------------------------------------------------------------------------
            % Filter application for the impulse response if specified in Filter_Config
            %-----------------------------------------------------------------------------------
            switch Filter_Config(4) % case 0 --> nothing happens

             case {1,2,3,4,5}       % cases 1 to 5 --> filter IR     

               [Impulse_Response,~] = ...    %% Access to interpolation models for the IR <--??
               filter_control(Impulse_Response,Delta_t_H,4,Sensor,Timestep_Reduction,Temperature_Zero,...
               Filter_Config,Cut_Off_Frequency,FIR_Config,FIR_Runs,'');

             case 6
               Filter_Config(5)              = 0; 
               [Impulse_Response] = ir_module(Impulse_Response,Time);
             case 66
%                Filter_Config(5)              = 0; 
               [Impulse_Response] = ir_module_quad(Impulse_Response,Time); 
            end % switch Filter_Config(4)
        
    save(filepathIRM, 'H','t_H','Delta_t_H','H_1J')
    end

        clear Delta_t_H t_ges_H t_ges H t_h
     
     
end % switch CorM <- System Identification

    
    
    


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if analyzeIR && strcmp(CorM,'C')    % <- if this the case, then skip the inversion
    Msrmnt_TempOrPress(1,:) = Calib_TempOrPress(1,1,:);
    Msrmnt_Heat_Flux(1,:) = Calibration_Heat_Flux(1,1,:);
    Reconstructed_Heat_Flux=0;
else % analyzeIR  && strcmp(CorM,'C')
fprintf('Starting inversion process... \n')
    
% Impulse Response - Penetration Time
    % An dieser Stelle ist die IR entweder direkt aus der Kalibrierung übernommen, oder sie
    % wurde für einen neuen Zeitvektor neu berechnet. In beiden Fällen werden im folgenden 
    % die ersten n IR-Einträge auf 0 gesetzt (entsprechend der Penetration Steps). 
    
    Penetration_Steps = double(int32(Penetration_Time/Delta_t)); % Delta_t = Delta_t_H ?!? müsste!
    Impulse_Response(1,1,2:Penetration_Steps+1) = 0;    % +1, weil Time(1) = 0 <-> der erste Penetration Step entspricht dem 2. Eintrag im Zeitvekotr
    
    % Wenn IR-Filter 8 gewählt wurde, wird die IR nach der Penetration Time exponentiell gefittet. 
    % Dieser Filter ist noch ausführlicher zu testen
    if Filter_Config(4) == 8
               [Impulse_Response] = ir_module_exp(Impulse_Response,Time,Penetration_Steps); 
    end




 if strcmp(CorM,'C')    %| In Calibration, no Msrmnt Signal is given, but Msrmnt Signal is 
                        %| required for Inversion step
                        %| -> Create 'artificial' Measurement Signal from Calibration Signal 
        Msrmnt_TempOrPress(1:length(Calib_TempOrPress(:,1,1)), ...
                           1:length(Calib_TempOrPress(1,1,:))) = 0;

        for i = 1:length(Calib_TempOrPress(:,1,1))
        for j = 1:length(Calib_TempOrPress(1,:,1)) 
           Buffer(1,:) = Calib_TempOrPress(i,j,:); 
           Msrmnt_TempOrPress(i,:) = Msrmnt_TempOrPress(i,:) + Buffer;
        end
        end
        clear Buffer
 end % strcmp(CorM,'C')

 switch Solver_Configuration
    case 0
    [Msrmnt_Heat_Flux]   = nisi_solver_pVC(Msrmnt_TempOrPress,Impulse_Response,Delta_t,...
                                           Penetration_Time,Future_Time_Window);
    case 1
    [Msrmnt_Heat_Flux]   = nisi_solver_sfe(Msrmnt_TempOrPress,Impulse_Response,Delta_t, ...
                                           Penetration_Time,Future_Time_Window);
    %case 2
    % [Msrmnt_Heat_Flux]   = nisi_solver_nl_sfe(Msrmnt_Temperature,Coefficient_Matrix,Delta_t,Penetration_Time,Future_Time_Window);
    case 3
    [Msrmnt_Heat_Flux]   = nisi_solver_sfe_CD(Msrmnt_TempOrPress,Impulse_Response,Delta_t, ...
                                           Future_Time_Window);
    case 4
    [Msrmnt_Heat_Flux]   = nisi_solver_sfe_linearHF(Msrmnt_TempOrPress,Impulse_Response,Delta_t, ...
                                           Future_Time_Window);
    case 5
    [Msrmnt_Heat_Flux]   = nisi_solver_sfe_matrix_form(Msrmnt_TempOrPress,Impulse_Response,Delta_t, ...
                                           Future_Time_Window);         
 end

% Postprocessing Inverted Heat Flux
[Msrmnt_Heat_Flux,~] = filter_control(Msrmnt_Heat_Flux,Delta_t,5,Sensor,...
                       Timestep_Reduction,Temperature_Zero,Filter_Config,...
                       Cut_Off_Frequency,FIR_Config,FIR_Runs,nPlat_SPE); 

fprintf('Inverse heat flux solved\n \n')
end % analyzeIR && strcmp(CorM,'C')
                   
switch CorM
    case 'C'
      if TorP == 'T'
          Calibration_Temperature = Calib_TempOrPress;
          Simulated_Temperature  = Simulated_TorP_Vector;
          TorPsaveNameSim = 'Simulated_Temperature';
      else
          Calibration_Pressure = Calib_TempOrPress;
          Simulated_Pressure  = Simulated_TorP_Vector;
          TorPsaveNameSim = 'Simulated_Pressure';  
      end
      Reconstructed_Heat_Flux = Msrmnt_Heat_Flux;
      filepath =  [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                    nameDataset '_' nameCalFolder '_C_Solution.mat' ];  
      save(filepath,TorPsaveNameSim,TorPsaveNameC,...
           'Calibration_Heat_Flux','Reconstructed_Heat_Flux','Time');
    case 'M'
      if TorP == 'T'
          Msrmnt_Temperature = Msrmnt_TempOrPress;
      else
          Msrmnt_Pressure = Msrmnt_TempOrPress;
      end
      filepath =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                    nameDataset '_' nameMeasFolder '_M_Solution.mat' ];
      save(filepath,TorPsaveNameM,'Msrmnt_Heat_Flux','Time');
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Frequency Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Frequency_Analysis
    [FFT_CT,Axis_CT,FFT_CHF,Axis_CHF,FFT_MT,Axis_MT,FFT_IF,Axis_IF,FFT_IR,Axis_IR] = ...
    freq_analysis_module(Delta_t,Calib_TempOrPress,Calibration_Heat_Flux, ...
    Msrmnt_TempOrPress,Msrmnt_Heat_Flux,Impulse_Response);
else
    FFT_CT = 'Na'; Axis_CT = 'Na'; FFT_CHF = 'Na'; Axis_CHF = 'Na'; %| Setting a variable to 'Na' means that it is not necessary for 
    FFT_MT = 'Na'; Axis_MT = 'Na'; FFT_IR  = 'Na'; Axis_IR  = 'Na'; %| the further program execution - yet it needs to exist because
    FFT_IF = 'Na'; Axis_IF = 'Na';                                    %| it may be part of a input for functions later on
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_module(Frequency_Analysis,Plot_Configuration,nameDataset,nameCalFolder,nameMeasFolder,CorM,Time, ...
            Calib_TempOrPress,Simulated_TorP_Vector,Calibration_Heat_Flux,...
            Msrmnt_TempOrPress,Msrmnt_Heat_Flux,Impulse_Response, ...
            FFT_CT,Axis_CT,FFT_MT,Axis_MT,FFT_IF,Axis_IF,FFT_CHF,...
            Axis_CHF,FFT_IR,Axis_IR,analyzeIR);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Export Data to Workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Exporting data to workspace... \n')
    assignin('base','Impulse_Response', Impulse_Response); 
    assignin('base','Time', Time); 
    assignin('base','nameDataset', nameDataset);
    assignin('base','nameCalFolder', nameCalFolder);
    assignin('base','nameMeasFolder', nameMeasFolder);
    assignin('base','CorM', CorM);
    
switch CorM 
    case 'C'
        assignin('base',TorPsaveNameC, Calib_TempOrPress);
        assignin('base','Simulated_TorP_Vector', Simulated_TorP_Vector);
        assignin('base','Calibration_Heat_Flux', Calibration_Heat_Flux);
        assignin('base','Reconstructed_Heat_Flux', Reconstructed_Heat_Flux); 
    
    case 'M'
        assignin('base',TorPsaveNameM, Msrmnt_TempOrPress); 
        assignin('base','Msrmnt_Heat_Flux', Msrmnt_Heat_Flux); 
end
fprintf('Data export done \n \n')

if ~sum(strcmp('noWB',varargin))
multiWaitbar('CLOSEALL');
end

programDuration=toc;
cpuDuration = cputime - Timer_Start;
fprintf(['Program finished after %4.1f seconds (@' ...
         '%4.1f seconds CPU-Time)\n'], programDuration ,cpuDuration)
fprintf('--------------------------------------------------------\n')


end
