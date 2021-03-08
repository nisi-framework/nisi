function [Impulse_Response_int] = ir_module_exp(Impulse_Response,Time,Penetration_Steps)
%% Module for interpolation routines for the IR



CounterIR = size(Impulse_Response);  

for i = 1:CounterIR(1)    
for j = 1:CounterIR(2)
    figure('name','Filtering IR'); hold
    legend show
   H(:,1) = Impulse_Response(i,j,:);
          plot(Time,H,'b','displayname','original IR')
   minH = min(H);
   
   [ ~,  Index1 ]         = max(H);
   [ Value, Index ]         = max(movmean(H,5));
   
   if Index1 < 0.9*Index
       figure('name','IR'); plot(H(1:2*Index))
       error('Please increase Timesstep Reduction\n')
   end
   

 
  % Fit new function 
  x1 = Penetration_Steps;
  x2 = 2*Penetration_Steps;     % <-- this is randomly chosen - is there a better solution?
  h = Impulse_Response(i,j,x2);
  b = log(h+1)/(x2-x1);
  
  for k = x1+1:x2
    H(k,1) = exp(b*(k-x1))-1;
  end
  
	plot(Time,H,'g--','displayname','IR after exp fit') 

   
   Impulse_Response_int(i,j,:) = H(:,1); 
   
   plot(Time,squeeze(Impulse_Response_int),'g','displayname','filtered IR')
end
end

clear CounterIR IR_Buffer Value H  
end