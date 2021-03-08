function [nameMWB] = multiWaitbarControl(...
    flagMWC,CorM,Pre_Calc_Calibration,Pre_Calc_Measurement,Auto_Param_Finder,...
    Filter_Config,analyzeIR,fcHandle)
% A multi waitbar displaying all sub processes (e.g. filtering, inversion) of the current NISI
% run and their current progress is set up and controlled by this function.
%
% Original version by 
% Fabian Hufgard
% 
%
% INPUT
%       flagMWC - flag for the multiWaitbarControl function
%       CorM - Calibrierung oder Messung? - Mögliche Inputstrings: 'C' oder 'M'
%   	Pre_Calc_Calibration - Pre Process (i.e. filtering) of Calibration Data
%   	Pre_Calc_Measurement - Pre Process (i.e. filtering) of Measurement Data
%   	Auto_Param_Finder - flag for activating function to find NISI parameters automatically
%   	Filter_Config - Vector containing configurations for the filter_cotrol module
%   	analyzeIR - flag for skipping HF inversion in the calibration case
%       fcHandle - Handle with which filter_control is controlled (switches between data types)
%
% OUTPUT 
%       A multi waitbar displaying all sub processes of the current NISI
%       run and their progress
%

% "load" global vairable(s)
global TorP

% Define Waitbar names
FIRCT  = 'Running FIR filter on Calibration_Temperature';
FIRCP  = 'Running FIR filter on Calibration_Pressure';
FIRCHF = 'Running FIR filter on Calibration_Heat_Flux';
FIRMT  = 'Running FIR filter on Msrmnt_Temperature';
FIRMP  = 'Running FIR filter on Msrmnt_Pressure';
FIRIR  = 'Running FIR filter on Impulse_Response';
FIRIHF = 'Running FIR filter on Inverted Heat Flux';

GaussCT  = 'Running Gauss filter on Calibration_Temperature';
GaussCP  = 'Running Gauss filter on Calibration_Pressure';
GaussCHF = 'Running Gauss filter on Calibration_Heat_Flux';
GaussMT  = 'Running Gauss filter on Msrmnt_Temperature';
GaussMP  = 'Running Gauss filter on Msrmnt_Pressure';
GaussIR  = 'Running Gauss filter on Impulse_Response';
GaussIHF = 'Running Gauss filter on Inverted Heat Flux';

SPECHF = 'Running SPE filter on Calibration_Heat_Flux';


switch flagMWC
    
    case 0  % Set up multiWaitbar at start of NISI run
        try
            multiWaitbar('CLOSEALL');
        end
                            % Hier noch eine 'busy' waitbar für das parsen
                            % der Inputdaten einbauen? Zu empfehlen falls
                            % das ab und zu länger dauert.
        switch CorM
            case 'M'
                if Pre_Calc_Measurement && Filter_Config(3)~=0
                    if TorP == 'T'
                        switch Filter_Config(3)
                            case 1
                                multiWaitbar( FIRMT, 0); 
                            case 2
                                multiWaitbar( GaussMT, 0); 
                            case 3
                                multiWaitbar( GaussMT, 0); 
                                multiWaitbar( FIRMT, 0); 
                            case 4
                                multiWaitbar( FIRMT, 0); 
                                multiWaitbar( GaussMT, 0); 
                        end
                    else %  if TorP == 'p'
                        switch Filter_Config(3)
                            case 1
                                multiWaitbar( FIRMP, 0); 
                            case 2
                                multiWaitbar( GaussMP, 0); 
                            case 3
                                multiWaitbar( GaussMP, 0); 
                                multiWaitbar( FIRMP, 0); 
                            case 4
                                multiWaitbar( FIRMP, 0); 
                                multiWaitbar( GaussMP, 0); 
                        end
                    end                       
                end
                
            case 'C'
                if Pre_Calc_Calibration && (Filter_Config(1)~=0 || Filter_Config(2)~=0)
                    if TorP == 'T'
                        switch Filter_Config(1)
                           	case 1
                                multiWaitbar( FIRCT, 0); 
                            case 2
                                multiWaitbar( GaussCT, 0); 
                            case 3
                                multiWaitbar( GaussCT, 0); 
                                multiWaitbar( FIRCT, 0); 
                            case 4
                                multiWaitbar( FIRCT, 0); 
                                multiWaitbar( GaussCT, 0); 
                        end
                    else %  if TorP == 'p'
                        switch Filter_Config(1)
                           	case 1
                                multiWaitbar( FIRCP, 0); 
                            case 2
                                multiWaitbar( GaussCP, 0); 
                            case 3
                                multiWaitbar( GaussCP, 0); 
                                multiWaitbar( FIRCP, 0); 
                            case 4
                                multiWaitbar( FIRCP, 0); 
                                multiWaitbar( GaussCP, 0); 
                        end
                    end 
                    

                    switch Filter_Config(2)
                        case 1
                            multiWaitbar( FIRCHF, 0); 
                        case 2
                            multiWaitbar( GaussCHF, 0); 
                        case 3
                            multiWaitbar( GaussCHF, 0); 
                            multiWaitbar( FIRCHF, 0); 
                        case 4
                            multiWaitbar( FIRCHF, 0); 
                            multiWaitbar( GaussCHF, 0); 
                        case 7
                            multiWaitbar( SPECHF, 0); 
                    end
                    
                end
                if Auto_Param_Finder
                    multiWaitbar('Automatic Parameter Selection');
                else
                    multiWaitbar('System Identification',0);    
                    multiWaitbar('Simulating Temperature',0);
                end
                multiWaitbar('Simulating Impulse Response',0);  
                
                switch Filter_Config(4)
                    case 1
                        multiWaitbar( FIRIR , 0);
                    case 2
                        multiWaitbar( GaussIR, 0); 
                    case 3
                        multiWaitbar( GaussIR, 0); 
                        multiWaitbar( FIRIR, 0); 
                    case 4
                        multiWaitbar( FIRIR, 0); 
                        multiWaitbar( GaussIR, 0); 
                end
        end
        

        if ~analyzeIR
            multiWaitbar('Inverting System',0);
        
            switch Filter_Config(5)
                case 1
                    multiWaitbar( FIRIHF , 0);
                case 2
                    multiWaitbar( GaussIHF, 0); 
                case 3
                    multiWaitbar( GaussIHF, 0); 
                    multiWaitbar( FIRIHF, 0); 
                case 4
                    multiWaitbar( FIRIHF, 0); 
                    multiWaitbar( GaussIHF, 0); 
            end  
        end

        
        
    case 'FIR' % Find WaitbarName for FIR Filter
        switch fcHandle
            case 1
                if TorP == 'T'
                    nameMWB = FIRCT;
                else % <-> TorP == 'p'
                    nameMWB = FIRCP;
                end
            case 2
                nameMWB = FIRCHF;
            case 3
                if TorP == 'T'
                    nameMWB = FIRMT;
                else % <-> TorP == 'p'
                    nameMWB = FIRMP;
                end
            case 4
                nameMWB = FIRIR;
            case 5 
                nameMWB = FIRIHF;
        end
        
    case 'Gauss' % Find WaitbarName for Gauss Filter
        switch fcHandle
            case 1
                if TorP == 'T'
                    nameMWB = GaussCT;
                else % <-> TorP == 'p'
                    nameMWB = GaussCP;
                end
            case 2
                nameMWB = GaussCHF;
            case 3
                if TorP == 'T'
                    nameMWB = GaussMT;
                else % <-> TorP == 'p'
                    nameMWB = GaussMP;
                end
            case 4
                nameMWB = GaussIR;
            case 5
                nameMWB = GaussIHF;
        end   
end


end
