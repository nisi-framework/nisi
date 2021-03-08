function [Heat_Flux,TempOrPress,Time,Delta_t,Impulse_Response,t_H,Temperature_Regime] = parser( ...
    pathDataFolder,nameDataset,nameCalFolder,nameMeasFolder,Pre_Calc,flagParser)

global TorP

switch flagParser
    %------------------------------------------------------------------------------
    case 'C'
    %------------------------------------------------------------------------------
    
        if Pre_Calc
            filepath =  [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                          nameDataset '_' nameCalFolder '_C_NISI.mat' ];
            load(filepath);            
        else
            filepath =  [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                          nameDataset '_' nameCalFolder '_C_NISI_filtered.mat' ];
            load(filepath);    
        end
        
        Heat_Flux   = Calibration_Heat_Flux   ;
        
        
        if exist('Calibration_Temperature','var') && exist('Calibration_Pressure','var')
            error(['Both Calibration Temperature and Calibration Pressure '...
                   'contained in %s file\n--> Choose only one.'],filepath);
        elseif ~exist('Calibration_Temperature','var') && ~exist('Calibration_Pressure','var')          
            error(['Neither Calibration Temperature nor Calibration Pressure '...
                   'contained in file\n%s\n--> Exactly one must be contained.'],filepath);
        elseif exist('Calibration_Temperature','var')
            TempOrPress = Calibration_Temperature ;
           	TorP = 'T';
        elseif exist('Calibration_Pressure','var')
            TempOrPress = Calibration_Pressure ;
            TorP = 'p';
        end
        
        if exist('Temperature_Regime','var') == 0
            Temperature_Regime = 'L';
        end
        Impulse_Response = 'TBC'; % TBC = To Be Calculated 
        t_H              = 'Na';
        Delta_t    = Time(2) - Time(1);
        
    %------------------------------------------------------------------------------
    case 'M'
    %------------------------------------------------------------------------------
    
        % Load measurement data
        %------------------------------------------------------------------------------
        if Pre_Calc
            filepath =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                            nameDataset '_' nameMeasFolder '_M_NISI.mat' ];
            load(filepath);
            
        else
            filepath =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                            nameDataset '_' nameMeasFolder '_M_NISI_filtered.mat' ];
            load(filepath);
        end
        % Calculate Delta_t of Measurement
        Delta_t    = Time(2) - Time(1);
        
        % Load Impulse Response
        %------------------------------------------------------------------------------
        filepathIRMF =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                         nameDataset '_' nameCalFolder '_C_Impulse_Response.mat' ];
        filepathIRCF =  [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                         nameDataset '_' nameCalFolder '_C_Impulse_Response.mat' ];
        try
            delete(filepathIRMF); %
        end
        
        copyfile(filepathIRCF, filepathIRMF, 'f'); % Validity check of this copy after filter block
            
        load(filepathIRMF); %
        

        Heat_Flux          = 'TBC' ; %| for NISI >> Temperature Data and a timevector,
        Impulse_Response   = H ;     %| Impulse Response Data only consists of a single vector H
        Temperature_Regime ='L';     %| in the form H(<Sensor_Number>,<Surface_Number><Datavector>)
        
        % Copy NISI parameter set to measurement folder
        %------------------------------------------------------------------------------
        filepathParaMF =  [ pathDataFolder nameDataset '\' nameMeasFolder '\' ...
                         nameDataset '_' nameCalFolder '_C_NISI_Parameter.mat' ];
        filepathParaCF =  [ pathDataFolder nameDataset '\' nameCalFolder '\' ...
                         nameDataset '_' nameCalFolder '_C_NISI_Parameter.mat' ];
                     
        copyfile(filepathParaCF, filepathParaMF, 'f'); % 
%         load(filepathParaMF); %



        if exist('S','var')     % <-- for old files where Msrmnt_Temperature was called S
            TorP = 'T';
            if size(S,3)~=1 
                for i=1:size(S,1)
                    TempOrPress(i,:) = S(i,1,:);
                end
                msg=sprintf(['Wrong dimension of input temperature vetor! Vector must be of form:\t'...
                         'S(<Sensor_Number>,<Datavector>)  '...
                         '--> A sensor does not distinguish between different surfaces.\n'...
                         'The vector has been transformed to the correct size.']);
                warning(msg);
            else
                TempOrPress = S ;
            end
            
        %| Check that exactly ONE! of Msrmnt_Temperature or Msrmnt_Pressure exists in loaded 
        %| NISI (or NISI_filtered) file
        elseif exist('Msrmnt_Temperature','var') && exist('Msrmnt_Pressure','var')
            error(['Both Measurement Temperature and Measurement Pressure '...
                   'contained in %s file\n--> Choose only one.'],filepath);
        elseif ~exist('Msrmnt_Temperature','var') && ~exist('Msrmnt_Pressure','var')          
            error(['Neither Measurement Temperature nor Measurement Pressure '...
                   'contained in file\n%s\n--> Exactly one must be contained.'],filepath);
               
        % Set TorP and correct vector format if necessary
        
        elseif exist('Msrmnt_Temperature','var')
            TorP = 'T';
            if size(Msrmnt_Temperature,3)~=1 
                for i=1:size(Msrmnt_Temperature,1)
                    TempOrPress(i,:) = Msrmnt_Temperature(i,1,:);
                end
                msg=sprintf(['Wrong dimension of input temperature vetor! Vector must be of form:\t'...
                         'Msrmnt_Temperature(<Sensor_Number>,<Datavector>)  '...
                         '--> A sensor does not distinguish between different surfaces.\n'...
                         'The vector has been transformed to the correct size.']);
                warning(msg);
            else
                TempOrPress = Msrmnt_Temperature ;
            end
            
        elseif exist('Msrmnt_Pressure','var')
            TorP = 'p';
            if size(Msrmnt_Pressure,3)~=1 
                for i=1:size(Msrmnt_Pressure,1)
                    TempOrPress(i,:) = Msrmnt_Pressure(i,1,:);
                end
                msg=sprintf(['Wrong dimension of input pressure vetor!' ' Vector must be of form:\t'...
                         'Msrmnt_Pressure(<Sensor_Number>,<Datavector>)  '...
                         '--> A sensor does not distinguish between different surfaces.\n'...
                         'The vector has been transformed to the correct size.']);
                warning(msg);
            else
                TempOrPress = Msrmnt_Pressure ;
            end
            
            
        end    
end

end
