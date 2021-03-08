function [ nl_Parameter ] = nl_identification_module(N1,A0,Time,B2,B3,B4,B5,B6,B7,B8, ...
                            Calibration_Temperature,Calibration_Heat_Flux,format)
%% Control Module for nonlinear solver identification
for i = 1:format(3)
    Input_Flux = permute(Calibration_Heat_Flux(:,:,i,:),[1 2 4 3]);
    Input_Temp = permute(Calibration_Temperature(:,:,i,:),[1 2 4 3]);
    
    [Numerator_Order_Vector,Numerator_Parameter_Vector, ...
     Denominator_Order_Vector,Denominator_Parameter_Vector, ...
     Error_Numerator_Vector, Error_Denominator_Vector, ...
     Simulated_Temperature_Vector ] =  ...
     nl_ident(Time,B2,B3,B4,B5,B6,B7,B8, ...
     Input_Temp,Input_Flux);
 
     Errors(i,1) = Error_Numerator_Vector;
     Errors(i,2) = Error_Denominator_Vector;
     
     Simulated_Temperature_NL(:,:,i,:)= Simulated_Temperature_Vector(:,:,:);
     
     nl_Parameter(i,1) = Numerator_Order_Vector;
     nl_Parameter(i,2) = Numerator_Parameter_Vector;
     nl_Parameter(i,3) = Denominator_Order_Vector;
     nl_Parameter(i,4) = Denominator_Parameter_Vector;   
     clear Input_Flux Input_Temp Numerator_Order_Vector Numerator_Parameter_Vector
     clear Error_Numerator_Vector Error_Denominator_Vector Simulated_Temperature_Vector
     clear Impulse_Response
end
filepath = [ N1 A0 '/' A0 '_NISI_Parameter_nl.mat' ];
save(filepath,'Errors','Simulated_Temperature_NL','nl_Parameter');
clear Errors Simulated_Temperature_NL
