function    [nl_Function_Parameter, nl_Derivative_Grid] = ...
             nl_parameter_function_fit(nl_Parameter,Temperature_Regime) 
 %% Preparation of the input data to create function fit 
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % Sort parameter vectors depending on temperature level
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [ T_Regime_sorted Regime_sort_Index ] = sort(Temperature_Regime);
   
  [ Full_Numerator_Order ]   = nl_ordersize(nl_Parameter(:,1),Temperature_Regime);
  [ Full_Denominator_Order ] = nl_ordersize(nl_Parameter(:,3),Temperature_Regime);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Creating input matrix for fitting algorithm and applying fit 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i = 1:length(Full_Numerator_Order)
      for j = 1:length(Temperature_Regime)
      Mediator = nl_Parameter{Regime_sort_Index(j),1};
      Courator = nl_Parameter{Regime_sort_Index(j),2};
         if any(Full_Numerator_Order(i) == Mediator)
            Pillar(j) = Courator(find(Mediator == Full_Numerator_Order(i)));
         else
            Pillar(j) = 0;
         end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Calling fitting algorithm here
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %Continue here 
  end

end
         

         
         
function  Parameter  = nl_single_fit( func_values, func_pos )         
%% Nelder-Mead simplex 
n               = length(func_pos);
Parameter(1:n)  = 0;
ParPlane        = zeros(n,n+1);
error           = zeros(1,n+1);
ParPlane(:,1)   = Parameter;    
[ Value Index ] = max(H);
error(:,1)      = NL_interp(Parameter,H,Value,Index,Time);
shrink          = 0;

for j = 1:n
    if Parameter(j) ~= 0
       Parameter(j) = (1.05)*Parameter(j);
    else
       Parameter(j) = 0.00025;
    end
    ParPlane(:,j+1) = Parameter;
    error(1,j+1) = NL_interp(Parameter,H,Value,Index,Time);
end

% Sort so ParPlane(1,:) has the lowest function value
[error,j] = sort(error);
ParPlane = ParPlane(:,j);

evals = 0;
while evals < 200 % Cut off criteria could use some more work 
    evals = evals+1;
   
    ParAverage    = sum(ParPlane(:,1:n), 2)/n;
    ParReflection = 2*ParAverage - ParPlane(:,end);
    Parameter(:)  = ParReflection; Error_Reflection = NL_interp(Parameter,H,Value,Index,Time);
    
    if Error_Reflection < error(:,1)      
       ParExpansion    = 3*ParAverage - 2*ParPlane(:,end);
       Parameter(:)    = ParExpansion;
       Error_Expansion = NL_interp(Parameter,H,Value,Index,Time);
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
             Parameter(:) = ParCenter; ErrorCenter = NL_interp(Parameter,H,Value,Index,Time);      
             if ErrorCenter <= Error_Reflection
                ParPlane(:,end) = ParCenter;
                error(:,end)    = ErrorCenter;
                    
             else
                shrink = 1;
             end
          else
             ParCCCenter = 0.5*ParAverage + 0.5*ParPlane(:,end);
             Parameter(:) = ParCCCenter; ErrorCCCenter = NL_interp(Parameter,H,Value,Index,Time);   
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
                 Parameter(:)  = ParPlane(:,j); error(:,j) = NL_interp(Parameter,H,Value,Index,Time);
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
end


function error = NL_interp(Parameter,H,Value,Index,Time)
a      = round(Index/2);
B      = Time.^Parameter(3).*exp(-Parameter(2)./(Time.^Parameter(1)));
error1 = abs((B - H'));
error  = sum((error1(a:Index)).^2);
end





