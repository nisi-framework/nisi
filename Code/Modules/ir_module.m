function [Impulse_Response_int] = ir_module(Impulse_Response,Time)
%% Module for interpolation routines for the IR
CounterIR = size(Impulse_Response);  

for i = 1:CounterIR(1)    
for j = 1:CounterIR(2)
   H(:,1) = Impulse_Response(i,j,:);
   %% Right side check for singularities
   for k = 1:length(H)
       if H(k)<0 &&  k>1
          H(1:k) = -1;
       end
   end
   for k = 1:length(H)
       if H(k)<0 && k>1
          H(k)=0;
       end
   end
   %% Left side check for singularities |< To be done if ever happens to be required 
   
   %% Fitting routine
   [Parameter,erroral] = IR_interpol(H,Time);
   [ Value Index ]         = max(H);
   IR_Buffer           = Time.^Parameter(3).*exp(-Parameter(2)./(Time.^Parameter(1)));
   IR_Buffer           = H(Index)/IR_Buffer(Index)*IR_Buffer;
   
   Impulse_Response_int(i,j,1:Index)         = IR_Buffer(1:Index);
   Impulse_Response_int(i,j,Index:length(H)) = H(Index:end,1); 
   Impulse_Response_int(i,j,1)               = 0; 
end
end

clear CounterIR IR_Buffer Value H  
end