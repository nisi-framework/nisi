function [Parameter,erroral] = IR_interpol(H,Time)
%% Nelder-Mead simplex direct search for minimizing error
n               = 3;
Parameter       = [ 0 0 0 ]';
ParPlane        = zeros(n,n+1);
error           = zeros(1,n+1);
ParPlane(:,1)   = Parameter;    
[ Value Index ] = max(H);
error(:,1)      = IR_interp4(Parameter,H,Value,Index,Time);
shrink          = 0;

for j = 1:n
    if Parameter(j) ~= 0
       Parameter(j) = (1.05)*Parameter(j);
    else
       Parameter(j) = 0.00025;
    end
    ParPlane(:,j+1) = Parameter;
    error(1,j+1) = IR_interp4(Parameter,H,Value,Index,Time);
end

% Sort so ParPlane(1,:) has the lowest function value
[error,j] = sort(error);
ParPlane = ParPlane(:,j);

evals = 0;
while evals < 200 % Cut off criteria could use some more work 
    evals = evals+1;
   
    ParAverage    = sum(ParPlane(:,1:n), 2)/n;
    ParReflection = 2*ParAverage - ParPlane(:,end);
    Parameter(:)  = ParReflection; Error_Reflection = IR_interp4(Parameter,H,Value,Index,Time);
    
    if Error_Reflection < error(:,1)      
       ParExpansion    = 3*ParAverage - 2*ParPlane(:,end);
       Parameter(:)    = ParExpansion;
       Error_Expansion = IR_interp4(Parameter,H,Value,Index,Time);
       if Error_Expansion < Error_Reflection
          ParPlane(:,end) = ParExpansion;
          error(:,end)    = Error_Expansion;
       else
          ParPlane(:,end) = ParReflection;
          error(:,end)    = Error_Reflection;
       end   
    else
       if Error_Reflection < error(:,3)
          ParPlane(:,end) = ParReflection;
          error(:,end)    = Error_Reflection;
       else 
          if Error_Reflection < error(:,end)
             ParCenter    = (1.5)*ParAverage - 0.5*ParPlane(:,end);
             Parameter(:) = ParCenter; ErrorCenter = IR_interp4(Parameter,H,Value,Index,Time);      
             if ErrorCenter <= Error_Reflection
                ParPlane(:,end) = ParCenter;
                error(:,end)    = ErrorCenter;
                    
             else
                shrink = 1;
             end
          else
             ParCCCenter = 0.5*ParAverage + 0.5*ParPlane(:,end);
             Parameter(:) = ParCCCenter; ErrorCCCenter = IR_interp4(Parameter,H,Value,Index,Time);   
             if ErrorCCCenter < error(:,end)
                ParPlane(:,end) = ParCCCenter;
                error(:,end)    = ErrorCCCenter;
             else
                shrink = 1;
             end
          end
          if shrink > 0
             for j=2:n+1
                 ParPlane(:,j) = ParPlane(:,1) + 0.5*(ParPlane(:,j) - ParPlane(:,1));
                 Parameter(:)  = ParPlane(:,j); error(:,j) = IR_interp4(Parameter,H,Value,Index,Time);
             end
             shrink = 0;
          end
        end
    end
    [error,j] = sort(error);
    ParPlane  = ParPlane(:,j);
end   
Parameter(:) = ParPlane(:,1);
erroral      = error(:,1);


function error = IR_interp4(Parameter,H,Value,Index,Time)
a      = round(Index/2);
B      = Time.^Parameter(3).*exp(-Parameter(2)./(Time.^Parameter(1)));
B      = Value/max(B)*B;
error1 = abs((B - H'));
error  = sum((error1(a:Index)).^2);






