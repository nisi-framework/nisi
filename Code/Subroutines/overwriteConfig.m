
function overwriteConfig(nameDataset,varargin) 
% This funciton overwrites the Config file of nameDataset
%
% Origional version by
% Fabian Hufgard
%
%
% INPUT
%       nameDataset
%       varargin (optional)    
%                   - Parameters which should be overwritten in the config file
%                       - Input in pairs with: >name of parameter< and >new value of parameter<
%                       - Possible paramters to overwrite:
%                           nameCalFolder
%                           nameMeasFolder
%                           Penetration_Time
%                           Future_Time_Window
%                           Timestep_Reduction_C
%                           Temperature_Zero_C
%                           Timestep_Reduction_M
%                           Temperature_Zero_M
%
% OUTPUT
%       Config file in Config folder
%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Settings %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load pfade.mat
pathConfig = ([pathNisiFolder, 'Configuration\']);
nameConfig = [ nameDataset '_Config' ];
configLoadHandle = str2func(nameConfig);    %Create function handle from name string



% Load old configuration
    [nameCalFolder       , ...
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
    ~                    , ...
    Solver_Configuration , ...
    analyzeIR            , ...
    ~                    , ...
    ~                    , ...
    ~                    , ...
    ~                    , ...
    ~                    , ...
    ~                    , ...
    nPlat_SPE            , ...
    Timestep_Reduction_M , ...
    Temperature_Zero_M   , ...
    ~,~                  , ...
    ~,~,~,~,~          ] = ...
    configLoadHandle();


%% Process potential input variables and in case overwrite config
if ~isempty(varargin{1})
    
% Overwrite parameters according to NISI-input
i = 1;   
while i < nargin
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
       case {'FluxTerms','fluxterms'}               % change Flux terms to specified array
           Manual_n_o= varargin{(i+1)};
       case {'TempTerms','tempterms'}               % change Temperature terms to specified array
           Manual_d_o= varargin{(i+1)};
       case {'noPlot', 'Plotoff', 'no_plot'}        % option to disable plot_module at the end of the code
           % No action required, but counter needs to be reduced by 1:
           i = i-1;         % <- this is required, because noPlot is not a input pair, so i should only be increased by 1 overall
       case 'noClear'        % option to disable the clearance of the base workspace at the start of the NISI code
           % No action required, but counter needs to be reduced by 1:
           i = i-1;         % <- this is required, because noClear is not a input pair, so i should only be increased by 1 overall
       case {'noWB', 'noWaitBar', 'noWaitbar'}        % option to disable closing the multiWaitbar at the start of the NISI code
           % No action required, but counter needs to be reduced by 1:
           i = i-1;         % <- this is required, because noClear is not a input pair, so i should only be increased by 1 overall
       case {'noExport'}    % option to disable the result export to the workspace
           % No action required, but counter needs to be reduced by 1:
           i = i-1;
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
                    'Disable workspace clearance:   ''noClear''\n'...
                    'Disable data export:           ''noExport''\n'])
   end
   i = i+2;
end % while i < nargin

end % if ~isempty(varargin{1})


%-------------------------------------------------------
% Create Manual_d_o and Manual_n_o strings
%-------------------------------------------------------

strManual_n_o1='';
strManual_n_o2='';                                       % In every loop: add a value of the Manual_n_o
for i=1:length(Manual_n_o(1,:))                       % vector as a str to the str variable (strManual_n_o1)
    bufferStr1=num2str(Manual_n_o(1,i));               % eventually strManual_n_o1 is the str to be written
    strManual_n_o1 = [strManual_n_o1 bufferStr1 '  ']; % into the created config file
    bufferStr2=num2str(Manual_n_o(2,i));
    bufferStr2 = pad(bufferStr2,strlength(bufferStr1));%<--equalize both lengths
    strManual_n_o2 = [ strManual_n_o2 bufferStr2 '  '];
    clear bufferStr1 bufferStr2
end

strManual_d_o1='';
strManual_d_o2='';
for j=1:length(Manual_d_o(1,:))                        %similar to Manual_n_o
    bufferStr1=num2str(Manual_d_o(1,j));
    strManual_d_o1 = [strManual_d_o1 bufferStr1 '  '];
    bufferStr2=num2str(Manual_d_o(2,j));
    bufferStr2 = pad(bufferStr2,strlength(bufferStr1));
    strManual_d_o2 = [strManual_d_o2 bufferStr2 '  '];
    clear bufferStr1 bufferStr2
end
clear i j


%-------------------------------------------------------
% Write new file
%-------------------------------------------------------
cd(pathConfig)
fileID = fopen([nameConfig '.m'],'w');          %create new config file


fprintf(fileID,[''...
'function [ nameCalFolder        , ...\n' ...
'           nameMeasFolder       , ...\n' ...
'           Frequency_Analysis   , ...\n' ...
'           Pre_Calc_Calibration , ...\n' ...
'           Pre_Calc_Measurement , ...\n' ...
'           Sensor               , ...\n' ...
'           Penetration_Time     , ...\n' ...
'           Future_Time_Window   , ...\n' ...
'           Timestep_Reduction_C , ...\n' ...
'           Temperature_Zero_C   , ...\n' ...
'           Auto_Param_Finder    , ...\n' ...
'           Max_Parameter_n_o    , ...\n' ...
'           Max_Parameter_d_o    , ...\n' ...
'           Parameter_Config_n   , ...\n' ...
'           Parameter_Config_d   , ...\n' ...
'           Manual_n_o           , ...\n' ...
'           Manual_d_o           , ...\n' ...
'           Filter_Config        , ...\n' ...
'           Cut_Off_Frequency    , ...\n' ...
'           FIR_Config           , ...\n' ...
'           FIR_Runs             , ...\n' ...
'           Plot_Configuration   , ...\n' ...
'           Solver_Configuration , ...\n' ...
'           analyzeIR            , ...\n' ...
'           parameterFile        , ...\n' ...
'           flagNL               , ...\n' ...
'           currentTemp          , ...\n' ...
'           pastTemps            , ...\n' ...
'           chemieIdentifikation , ...\n' ...
'           chemieSolver         , ...\n' ...
'           nPlat_SPE            , ...\n' ...
'           Timestep_Reduction_M , ...\n' ...
'           Temperature_Zero_M   , ...\n' ...
'           P4,P5                , ...\n' ...
'           P6,P7,P8,P9,P10    ] = ...\n' ...
'           ' nameConfig '()\n' ...
'              %%  /\\ To be adjusted for each configuration\n' ...
'\n' ...
'\n' ...
'%%%% Path Configuration\n' ...
'nameCalFolder    = ''' nameCalFolder ''';\n' ...
'nameMeasFolder   = ''' nameMeasFolder '''; %%Handle to identify different measurements in one folder\n' ...
'parameterFile    = '''' ; %%z.b. \\HENDERSON_PYRAMIDE\\NISI_HENDERSON_PYRAMIDE_Parameter_var_T.mat; \n' ...
'\n' ...
'%%%% Program Control:\n' ...
'Frequency_Analysis  = ' num2str(Frequency_Analysis) ' ; %%Activates Frequency Domain Visualisation                      \n' ...
'analyzeIR           = ' num2str(analyzeIR) ' ;  %%|Suppresses inversion for faster calculation;\n' ...
'                          %%|Additionally plots temperatur difference\n' ...
'                          %%|between Simulated vs. Calibration Temperature\n' ...
'\n' ...
'flagNL                = false ; %%Calculation of nonlinear System Impulse Response            \n' ...
'currentTemp           = false ; %%Calculation with Impulse Response of current Temperature    \n' ...
'pastTemps             = false ; %%Calculation with Impulse Responses of past Time steps       \n' ...
'\n' ...
'chemieIdentifikation  = false ; %%Activates chemical terms in identification                  \n' ...
'chemieSolver          = false ; %%Activates Chemical solver                                   \n' ...
'\n' ...
'%%%% Program Configuration \n' ...
' \n' ...
'Pre_Calc_Calibration  = ' num2str(Pre_Calc_Calibration) '  ; %%|Perform input data calculations \n' ...
'Pre_Calc_Measurement  = ' num2str(Pre_Calc_Measurement) ' ; %%|>> Filters and reduction etc     \n' ...
' \n' ...
' \n' ...
'Sensor             = [ ' num2str(Sensor(1)) '    %% |Analyse a specific sensor / surface     \n' ...
'                       ' num2str(Sensor(2)) ' ]; %% |combination                             \n' ...
'Penetration_Time   = ' num2str(Penetration_Time) ';  %% |Thermal penetration time in s                                 \n' ...
'Future_Time_Window = ' num2str(Future_Time_Window) ';  %% |Width in s of the analysation window counting from          \n' ...
'                            %% |Penetration_Time , if =0 -> Only the timestep                 \n' ...
'                            %% |at Penetrationtime + 1*timestep is considered                 \n' ...
'Timestep_Reduction_C = ' num2str(Timestep_Reduction_C) '  ; %% |Only every Timestep_Reduction timestep is used                 \n' ...
'Timestep_Reduction_M = ' num2str(Timestep_Reduction_M) '  ; %% ^same for Measurement                   \n' ...
'Temperature_Zero_C   = ' num2str(Temperature_Zero_C) '  ; %% Time in s used to equate zero level temperature                \n' ...
'Temperature_Zero_M   = ' num2str(Temperature_Zero_M) '  ; %% ^same for Measurement                  \n' ...
' \n' ...
'Auto_Param_Finder  = ' num2str(Auto_Param_Finder) '    ; %% Activation of the automatic parameter finder                   \n' ...
'Max_Parameter_n_o  = ' num2str(Max_Parameter_n_o) '    ; %% Amount of Parameters checked for viability                     \n' ...
'Max_Parameter_d_o  = ' num2str(Max_Parameter_d_o) '    ; %% Amount of Parameters checked for viability                     \n' ...
'Parameter_Config_n = ' num2str(Parameter_Config_n) '    ; %% |Configuration of the used parameters of the transfer          \n' ...
'Parameter_Config_d = ' num2str(Parameter_Config_d) '    ; %% |function -> 1 == D^(j/2) standard derivatives row \n' ...
'                            %% |For other configurations look @                               \n' ...
'                            %% |<NISIPath>/Code/Subroutines/parameter_conf.m \n' ...
'%%-------------------------------------------------------------------------- \n' ...
'%%Parameters for identification of the transfer function if                                         \n' ...
'%%Auto_Param-Finder = false \n' ...
'%%First row defines the used derivatives \n' ...
'%%Second row  = 0 >> Parameter defined by least squares algorithm                                   \n' ...
'%%Second row ~= 0 >> Parameter is fixed on the given value                                          \n' ...
'%%-------------------------------------------------------------------------- \n' ...
' \n' ...
'%%Flux Terms \n' ...
'Manual_n_o =   [  ' strManual_n_o1 '                                                                \n' ...
'                  ' strManual_n_o2 ' ] ;                                                            \n' ...
' \n' ...
' \n' ...
'%% Temperature Terms \n' ...
'Manual_d_o =   [  ' strManual_d_o1 '                                                                \n' ...
'                  ' strManual_d_o2 ' ] ;                                                            \n' ...
' \n' ...
' \n' ...
'%%%% Filter Configuration \n' ...
' \n' ...
'Filter_Config = [ ' num2str(Filter_Config(1)) '    %% Calibration Temperature Data  | 0 >>> No Filter                           \n' ...
'                  ' num2str(Filter_Config(2)) '    %% Calibration Heat Flux Data    | 1 >>> FIR Filter                          \n' ...
'                  ' num2str(Filter_Config(3)) '    %% Measurement Temperature Data  | 2 >>> Running Gauss Filter                \n' ...
'                  ' num2str(Filter_Config(4)) '    %% Impulse Response              | 3 >>> RGF + FIR                           \n' ...
'                  ' num2str(Filter_Config(5)) ' ]; %% Inverted Heat Flux            | 4 >>> FIR + RGF \n' ...
'                       %%                               | 5 >>> Reduction                           \n' ...
'                       %%                               | 7 >>> Sharp Edge Filter (for NISI pulses) \n' ...
' \n' ...
'%%Gauss-Filter Frequencies: \n' ...
'Cut_Off_Frequency = [ ' num2str(Cut_Off_Frequency(1)) ' %%Gauss-Filter cut off frequency (Hz) |Calibration T                     \n' ...
'                      ' num2str(Cut_Off_Frequency(2)) ' %%Gauss-Filter cut off frequency (Hz) |Calibration Q                     \n' ...
'                      ' num2str(Cut_Off_Frequency(3)) ' %%Gauss-Filter cut off frequency (Hz) |Measurement T                     \n' ...
'                      ' num2str(Cut_Off_Frequency(4)) ' %%Gauss-Filter cut off frequency (Hz) |Pulse Resp.                       \n' ...
'                      ' num2str(Cut_Off_Frequency(5)) ' ]; %%Gauss-Filter cut off frequency (Hz) |generic                           \n' ...
' \n' ...
'%%FIR-Filter: \n' ...
'%%See <NISIPath>/Code/Subroutines/filter_r_fir.m for configuration types                            \n' ...
'FIR_Config    = [ ' num2str(FIR_Config(1)) '    %% Calibration Temperature Data                                              \n' ...
'                  ' num2str(FIR_Config(2)) '    %% Calibration Heat Flux Data                                                \n' ...
'                  ' num2str(FIR_Config(3)) '    %% Measurement Temperature Data                                              \n' ...
'                  ' num2str(FIR_Config(4)) '    %% Impulse Response \n' ...
'                  ' num2str(FIR_Config(5)) ' ]; %% Inverted Heat Flux] \n' ...
' \n' ...
'FIR_Runs      = [ ' num2str(FIR_Runs(1)) '    %% Calibration Temperature Data                                             \n' ...
'                  ' num2str(FIR_Runs(2)) '    %% Calibration Heat Flux Data                                              \n' ...
'                  ' num2str(FIR_Runs(3)) '    %% Measurement Temperature Data                                              \n' ...
'                  ' num2str(FIR_Runs(4)) '    %% Impulse Response \n' ...
'                  ' num2str(FIR_Runs(5)) ' ]; %% Inverted Heat Flux] \n' ...
' \n' ...
'%% Sharp Edge Filter                       \n' ...
'%% see <NISIPath>/Code/Subroutines/filter_sharpPulseEdges_NISI.m \n' ...
'nPlat_SPE     =  ' num2str(nPlat_SPE) ';     %% Number of Plateaus the SPE is reducing to  \n' ...                        \n' ...
' \n' ...
'%%%% Plot Configuration \n' ...
'                       %% C  M    %% Column 1 of Plot_Configuration for Calib, 2 for Meas           \n' ...
'Plot_Configuration   = [ 1  0    %% Calibration/Simulation Temperature    |Manual configuration     \n' ...
'                         1  0    %% Calibration/Inverted Heat Flux        |of plots                 \n' ...
'                         0  1    %% Measurement Heat Flux (Inverted)      |0 - deactivated          \n' ...
'                         1  1    %% Impulse Response                      |1 - activated            \n' ...
'                         0  0    %% FFT - Calibration Temperature                                   \n' ...
'                         0  0    %% FFT - Calibration Heat Flux                                     \n' ...
'                         0  0    %% FFT - Measurement Temperature                                   \n' ...
'                         0  0    %% FFT - Inverted Heat Flux                                     \n' ...
'                         0  0 ]; %% FFT - Impulse Response                                          \n\n' ...' \n' ...
'%%%% Solver Configuration \n' ...
'Solver_Configuration = ' num2str(Solver_Configuration) ';        %% 0 - phased Van Cittert                                          \n' ...
'                                 %% 1 - sequential function estimation                              \n' ...
'                                 %% 3 - Christian''s sequential function estimation                              \n' ...
'                                 %% 4 - Linear HF sequential function estimation                              \n' ...
' \n' ...
'%%%% Dummy variables for future use \n' ...
'P4 = '''';                          \n' ...
'P5 = '''';                          \n' ...
'P6 = '''';                          \n' ...
'P7 = '''';                          \n' ...
'P8 = '''';                          \n' ...
'P9 = '''';                          \n' ...
'P10= '''';                          \n' ... \n' ...
'end']);
% open(nameConfig);
rehash
cd([pathDataFolder nameDataset]);

end % function overwrite
