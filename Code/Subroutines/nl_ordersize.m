function [ Full_Order ]  = nl_ordersize( P_Cell, Temperature_Regime )%% Sorting Algorithm   minBuffer = 0;   maxBuffer = 0;      for i = 1:length(Temperature_Regime)       max_v = max(P_Cell{i});       min_v = min(P_Cell{i});       if min_v < minBuffer          minBuffer = min_v;       end        if max_v > maxBuffer          maxBuffer = max_v;       end    end   Full_Order = minBuffer:0.5:maxBuffer;end