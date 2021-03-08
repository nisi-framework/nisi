function [Impulse_Response_filt] = ir_module_quad(Impulse_Response,Time)
%% Module for interpolation routines for the IR
CounterIR = size(Impulse_Response);  
flagPlot = 0;

for i = 1:CounterIR(1)    
for j = 1:CounterIR(2)
   H(:,1) = Impulse_Response(i,j,:);
   minH = min(H);
   
        if flagPlot
            figure('name','Filtering IR'); hold
            legend show
            plot(Time,H,'b','displayname','original IR')
        end

   
   [ ~,  Index1 ]         = max(H);
   [ Value, Index ]         = max(movmean(H,5));
   
   if Index1 < 0.9*Index
       figure('name','IR'); plot(H(1:2*Index))
       error('Please increase Timesstep Reduction\n')
   end
   
   % Flag, um den Anfangsanstieg quadratisch zu fitten
   flagQuadFit = 0;
   
   %% Set values which are <0 to 0
   for k = 1:Index1
       if H(k)<0 &&  k>1
          H(1:k) = minH;
%           plot(Time,H,'k')
       end
   end
   for k = 1:Index1
       if H(k)<0
          H(k)=0;
%           plot(Time,H,'r')
          flagQuadFit = 1; % Aktiviere quadratischen Fit des Anfangsanstiegs
       end
   end
   
       if flagPlot
        plot(Time,H,'r--','displayname','IR after 0 correction')
       end
   
   %% Überbügeln eines anfänglichen Peaks
   h = 0;
   while h == 0
   for k = 2:min(Index,Index1)
      if  H(k)-H(k-1) <0        % Wenn es zwischen dem Anfang und dem Maximum der IR negative Steigungen gibt, bedeutet das, dass hier ein anfängliches Peak vorliegt
          hact = H(k-1);
          a = find(H>hact,1);   % Finde ersten H-Wert nach dem anfänglichen Peak der wieder größer ist als der anflängliche Peak -> a ist dessen Index
          H(1:k) = 0;           % Set values up to the k-th entry to 0
          flagQuadFit = 1;
          h = 0;
          break
      else
          a = find(H>Value/10,1);
          h=1;
      end
   end
   end
   
   % Fit eines quadratischen Anfangsanstiegs

   if flagQuadFit
   % Fit der IR mit Polynom 2ten Grades bis zum Punkt a, d.h. index an dem H wieder größer als
   % der anfängliche Peak ist.
   hh = 0;
   while hh == 0 
       hhh = 0;
    Hstrichvona = mean([H(a+1)-H(a) H(a)-H(a-1)]);
%     b = (Hstrichvona*(a+1)-H(a+1))/(a+1)^2;
%     c = 2*H((a+1))/(a+1)-Hstrichvona;
%     d = -b-c;
    
    b = (Hstrichvona*(a-1)-H(a))/(a-1)^2;
    c = Hstrichvona-2*b*a;
    d = -b-c;
    
   
   for k = 1:a
       H(k) = b*(k)^2+c*(k)+d;
       if H(k)<0
           hhh = 1;
       end
   end  

   
     if hhh == 0
         hh = 1;
     else
         a = a+1;
     end
   end
            if flagPlot
                plot(Time,H,'r--','displayname','IR after Buegel correction') 
            end
   end
   
   Impulse_Response_filt(i,j,:) = H(:,1); 
   
            if flagPlot
                plot(Time,squeeze(Impulse_Response_filt),'g','displayname','filtered IR')
            end
end
end

clear CounterIR IR_Buffer Value H  
end