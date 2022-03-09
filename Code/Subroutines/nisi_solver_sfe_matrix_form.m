function [Msrmnt_Heat_Flux] = nisi_solver_sfe_matrix_form(Msrmnt_Temperature,Impulse_Response,...
    Delta_t,Penetration_Time,Future_Time_Window)
%% Solver using sequential function estimation in matrix form for multiple heat flux/sensor pairs(Beck's)  ,
% S.269-272 Inverse Heat Conduction / Ill-posed Problems
% für konstantes q_m für alle future_time_steps
% gleiche Anzahl Sensoren/Flächen
% Annahme: Gamma kürzt sich, da Rauschen Standardverteilt (7.2.20)
%
% Input:       Msrmnt_Temperature   - Matrix with measurment temperature,
%                                     sensor data in row
%              Impulse_Response     - Impulse response matrix (3D),
%                                     row = Surface, column = Sensor
%              Delta_t              - Time between datapoints, digit
%              Future_Time_Window   - Regularization window, digit
%
% Output:      Msrmnt_Heat_Flux     - Heat flux vector

%%Initalize vectors
t_penetration                   = double(int32(Penetration_Time/Delta_t));
r                               = double(int32(Future_Time_Window/Delta_t));   %Future_time_steps
J                               = size(Msrmnt_Temperature,1);          %Sensors
p                               = size(Msrmnt_Temperature,1);      %Areas
Data_Points                     = length(Msrmnt_Temperature);

IR_plotflag=1;              %Plot for analysis reasons
HF_plotflag=0;              %Plot for analysis reasons
Temp_plotFlag=0;
kk=0;
Msrmnt_Heat_Flux=zeros(J,Data_Points);
Time=linspace(0,length(Msrmnt_Heat_Flux)*Delta_t,length(Msrmnt_Heat_Flux));
%% Identify Penetrationtime and set it to zero
if IR_plotflag
f1=figure;hold on
colormap(f1,'hsv');
end

for Sensor=1:p
    for Surface=1:J
        sign_array=sign(Impulse_Response(Surface,Sensor,:));
        penetration=find(diff(sign_array)==2);          %find sign switch from -1 to 1
        if ~isempty(penetration)
            Impulse_Response(Surface,Sensor,1:penetration(1))=0;
        end
        if IR_plotflag
            plot(Time,squeeze(Impulse_Response(Surface,Sensor,:)))
        end
    end
end

Impulse_Response(:,:,1:t_penetration)=0;
Impulse_Response(isnan(Impulse_Response))=0;
Impulse_Response(Impulse_Response<0)=0;
%T_future=zeros(J*r,1);                  %(7.2.2a,b)
%q_future=zeros(r*p,1);                  %(7.2.3a,b)
%X=zeros(J*r,p*r);                       %(7.2.4)
a=zeros(J,p,r);                         %(7.2.5a)       Reihen-Sensoren, Spalten-Flächen

%% Matrix preparation
for FutureTime=1:r
    for Surface=1:p
        for Sensor=1:J
            a(Sensor,Surface,FutureTime)=Impulse_Response(Surface,Sensor,FutureTime); %(7.2.5a)
        end
    end
end

triangular_matrix=cell(r);
for FutureTime=1:r      %Building triangular matrix, hopefully (7.2.4)
    
    counter=1;
    while counter<=FutureTime
        triangular_matrix{FutureTime,counter}=a(:,:,FutureTime-counter+1);
        counter=counter+1;
    end
    nonsense=FutureTime+1;
    while nonsense<=r
        triangular_matrix{FutureTime,nonsense}=zeros(J,p,1);
        nonsense=nonsense+1;
    end
    A_buffer{FutureTime,1}=eye(p);         % für konstantes q (7.2.16)
end
locs=cellfun(@isempty, triangular_matrix);
X=cell2mat(triangular_matrix);
A=cell2mat(A_buffer);

Z=X*A;%(7.2.17)
T_i=zeros(J,r);
Y_i=zeros(J,r);
Z_squared=(Z'*Z);
%% Recursive sequential function estimation start
for currentTimeStep = 1:length(Msrmnt_Temperature(1,:))  - r
    %------------------------------------------------------------------------
    % Determine the temperature at timestep induced by heat flux equated
    % for previous timesteps (q_M=0)
    %------------------------------------------------------------------------
    for i=1:r
        for Sensor=1:J
            Induced_Temperature =0;
            for Surface = 1:p             
                %Calculate temperature through duhamel
                iA(1,:) = Msrmnt_Heat_Flux(Surface,:);
                iB(1,:) = Impulse_Response(Surface,Sensor,:);
                ConvBuff = conv(iB,iA);
                Induced_Temperature = Induced_Temperature +ConvBuff(1:Data_Points);
                if currentTimeStep == length(Msrmnt_Temperature(1,:))  - r
                    Simulated_Temp(Sensor,:)=Induced_Temperature;
                end
            end
            T_i(Sensor,i)= Induced_Temperature(currentTimeStep+i-1);
            Y_i(Sensor,i)= Msrmnt_Temperature(Sensor,currentTimeStep+i-1);
            
            if currentTimeStep == length(Msrmnt_Temperature(1,:))  - r
                Induced_Memory(Sensor,:)=Induced_Temperature;
            end
        end
    end
    
    
    for i=1:r
        if i==1
            T_buffer=T_i(:,1);
            Y_buffer=Y_i(:,1);
        else
            T_buffer=vertcat(T_buffer,T_i(:,i));        %(7.2.2.a,b)
            Y_buffer=vertcat(Y_buffer,Y_i(:,i));
        end
    end
    
    
    beta_estimator=Z_squared\(Z'*(Y_buffer-T_buffer));       %(7.2.19)
 %   beta_estimator=max(beta_estimator,0);   % No negative heat flux allowd
    
    Msrmnt_Heat_Flux(:,currentTimeStep)=beta_estimator;

    progress = currentTimeStep / ...          % progress increases from 0 to 1
        (length(Msrmnt_Temperature(1,:))-r);
    if progress > kk
        multiWaitbar('Inverting System',progress);
        kk=kk+0.001;                                % with kk the output rate is reduced to every 0.1%
    end
end % for currentTimeStep = 1:length(Msrmnt_Temperature(1,:))  - Future_Time_Steps


multiWaitbar('Inverting System',1);
if HF_plotflag
    if currentTimeStep==length(Msrmnt_Temperature(1,:))  - r
        figure;hold on
        for i=1:J
            plot(Msrmnt_Heat_Flux(i,:))
        end
    end
end
if Temp_plotFlag
    if currentTimeStep==length(Msrmnt_Temperature(1,:))  - r
    figure;hold on
    for i=1:size(Msrmnt_Temperature,1)
        plot(Time,Msrmnt_Temperature(i,:),'o-','MarkerIndices',int32(1:length(Msrmnt_Temperature)*0.1:length(Msrmnt_Temperature)))
        plot(Time,Induced_Memory(i,:))
    end
    legend show
    end
end
end
