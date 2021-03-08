%% Analytical solution
clear all

%--------------------------------------------------------------------------
c      = 384;  %Wärmekapazität J/(kg*K)
rho    = 8995; %Dichte kg/m^3
lambda = 394;  %Wärmeleitfähigkeit W/(m*K)
alpha  = lambda/(rho*c);
%--------------------------------------------------------------------------
tres   = 0.004; %Auflösung timesteps
tstep  = 0.004;
x      = 0.005;  %Tiefe in m
total_time=44;  % Time in s
lengt  = total_time/tstep;
hf=10000; %Wärmestrom in W/m^2
t_penetration=(x/(3.6*sqrt(lambda/(rho*c))))^2;
%% Heat Flux, triangular function | uncomment for use

% points_triangular=3000;
% straight=linspace(0,hf,points_triangular);
% straight_minus=fliplr(straight);
% test=[straight, straight_minus];
% Heat_Flux=[zeros(1,(lengt-points_triangular*2)/2),test,zeros(1,(lengt-points_triangular*2)/2)];

%% Heat Flux, gauss | uncomment for use
x_gauss=linspace(0,total_time,lengt);
center=total_time/2;                    %center of gauss
std_dev=0.3;

for i = 1:lengt
    
    Heat_Flux(i)=hf*exp(-(x_gauss(i)-center)^2/2*std_dev^2); %Gauss Funktion
    
end


%% Impulse Response (Heaviside derivative)
for i=1:length(Heat_Flux)
    Time(i) = (i)*tres;
end
h_heavy(1:length(Heat_Flux)) = 0; % impulse response array
for i = 1:tstep/tres
    [ T_1 ] = heaviside_t(x,lambda,alpha,Time(i));
    h_heavy(i) = T_1;
end
for i = tstep/tres+1:length(Heat_Flux)
    [ T_1 ]   = heaviside_t(x,lambda,alpha,Time(i)-tstep);
    [ T_2 ]   = heaviside_t(x,lambda,alpha,Time(i));
    h_heavy(i) = T_2 - T_1;
end
h_heaviside=h_heavy/tres;  %derivative of Heaviside function -> impulse response

%% IR analytic
for i = 1:length(Heat_Flux)
    Hard(i) = exp(-(x^2)/(4*alpha*Time(i)))/(sqrt(pi*lambda*rho*c*Time(i))); %impulse response
end

Temperature = conv(h_heaviside,Heat_Flux);
Temperature = Temperature(1:length(Heat_Flux));

%%Calibration Data without Noise
Calibration_Heat_Flux(1,1,:) = Heat_Flux;
Calibration_Temperature=zeros(1,1,lengt);
Calibration_Temperature(1,1,:) = Temperature;
%% Calibration Temperature with Noise
% for i=1:length(Calibration_Temperature)
% Calibration_Temperature(1,1,i) = Temperature(1,i)+max(Temperature)/30*(2*rand-1);
% end
%% Rearranging IR for NISI (start with 0)

H(1,1,:)=[0 , h_heaviside(1,1:end-1)];                         % heaviside IR
H_a(1,1,:)=[0, Hard(1,1:end-1)];    % analytical IR
Delta_t_H=Time(2)-Time(1);
t_H=[0,Time(1,1:end-1)]';

Msrmnt_Temperature(1,:)=Calibration_Temperature;
Input_Heat_Flux(:,1)=Calibration_Heat_Flux;
%% Plots
figure(1),plot(Time,squeeze(Calibration_Temperature));
xlabel('Time')
ylabel('Temperature, K')
figure(2),plot(Time,Heat_Flux);
xlabel('Time')
ylabel('Heat Flux, W/m^2')
figure(3);
hold on;
plot(Time, Hard);
plot(Time, h_heaviside,'--o','MarkerIndices',1:1000:length(h_heaviside));
legend('Analytical','Heaviside derivative')
hold off;
legend show
xlabel('Time')
ylabel('Temperature, K')

Time=Time';
clear buffer Temperature Heat_Flux buffer T_1 T_2 total_time c alpha i hf lengt rho...
    tstep tres x h_heaviside h_heavy lambda Hard

