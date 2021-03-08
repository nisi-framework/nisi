function [FFT_CT,Axis_CT,FFT_CHF,Axis_CHF,FFT_MT,Axis_MT,FFT_IF,Axis_IF,FFT_IR,Axis_IR] = ...
          freq_analysis_module(Delta_t,Calibration_Temperature,Calibration_Heat_Flux, ...
          Msrmnt_Temperature,Msrmnt_Heat_Flux,Impulse_Response)
%% Frequency Analysis
%--------------------------------------------------------------------------
% FFT_CT   == Spectral Intensitiy of the Calibration Temperature
% Axis_CT  == Frequency Axis of the Calibration Temperature
% FFT_CHF  == Spectral Intensitiy of the Calibration Heat Flux
% Axis_CHF == Frequency Axis of the Calibration Calibration Heat Flux
% FFT_MT   == Spectral Intensitiy of the Measurement Temperature
% Axis_MT  == Frequency Axis of the Measurement Temperature
% FFT_IF   == Spectral Intensitiy of the Inverse Heat Flux
% Axis_IF  == Frequency Axis of the Inverse Heat Flux
% FFT_IR   == Spectral Intensitiy of the Impulse Response
% Axis_IR  == Frequency Axis of the Impulse Response
%--------------------------------------------------------------------------

fprintf('Performing Frequency Analysis \n \n')

  if exist('Calibration_Temperature','var') == 1
      for i = 1:length(Calibration_Temperature(:,1,1))
          for j = 1:length(Calibration_Temperature(:,1,1))
              input(1,:)         = Calibration_Temperature(i,j,:);
             [output_1,output_2] = frequency_analysis(input,Delta_t);
              FFT_CT(i,j,:)      = output_1(:);
              Axis_CT(i,j,:)     = output_2(:);
          end
      end     
  else
    FFT_CT   = 'B';  Axis_CT  = 'B';    
  end
  
  clear input output_1 output_2
  if exist('Calibration_Heat_Flux','var') == 1
      for i = 1:length(Calibration_Heat_Flux(:,1,1))
          for j = 1:length(Calibration_Heat_Flux(:,1,1))
              input(1,:)         = Calibration_Heat_Flux(i,j,:);
             [output_1,output_2] = frequency_analysis(input,Delta_t);
              FFT_CHF(i,j,:)     = output_1(:);
              Axis_CHF(i,j,:)    = output_2(:);
          end
      end  
  else
    FFT_CHF  = 'B';  Axis_CHF  = 'B';    
  end
  
  clear input output_1 output_2
  if exist('Msrmnt_Temperature','var') == 1
      for i=1:length(Msrmnt_Temperature(:,1))
          input(1,:)         = Msrmnt_Temperature(i,:);
         [output_1,output_2] = frequency_analysis(input,Delta_t);
          FFT_MT(i,:)        = output_1(:);
          Axis_MT(i,:)       = output_2(:);
      end
  else
    FFT_MT   = 'B';  Axis_MT  = 'B';    
  end
 
  clear input output_1 output_2
  if exist('Msrmnt_Heat_Flux','var') == 1
      for i=1:length(Msrmnt_Heat_Flux(:,1))
          input(1,:)         = Msrmnt_Heat_Flux(i,:);
         [output_1,output_2] = frequency_analysis(input,Delta_t);
          FFT_IF(i,:)        = output_1(:);
          Axis_IF(i,:)       = output_2(:);
      end
  else
    FFT_IF   = 'B';  Axis_IF  = 'B';    
  end
  
  clear input output_1 output_2
  if exist('Impulse_Response','var') == 1
      for i = 1:length(Impulse_Response(:,1,1))
          for j = 1:length(Impulse_Response(:,1,1))
              input(1,:)         = Impulse_Response(i,j,:);
             [output_1,output_2] = frequency_analysis(input,Delta_t);
              FFT_IR(i,j,:)      = output_1(:);
              Axis_IR(i,j,:)     = output_2(:);
          end
      end  
  else
    FFT_IR   = 'B';  Axis_IR  = 'B';    
  end  
  
end