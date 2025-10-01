%% Grandi model - Main file
clear
% close all
clc

%% Setting input parameters
% ICs & Stimulation period
% Male ICs
%load yfin_endo_0p25Hz; CL = 4000; % basic cycle length (ms)
%load yfin_endo_0p33Hz; CL = 3000; % basic cycle length (ms)
%load yfin_endo_0p5Hz; CL = 2000; % basic cycle length (ms)
load yfin_endo_1Hz; CL = 1000; % basic cycle length (ms)
%load yfin_endo_2Hz; CL = 500; % basic cycle length (ms)

% Female ICs
%load yfin_endo_female_0p25Hz; CL = 4000; % basic cycle length (ms)
%load yfin_endo_female_0p33Hz; CL = 3000; % basic cycle length (ms)
%load yfin_endo_female_0p5Hz; CL = 2000; % basic cycle length (ms)
%load yfin_endo_female_1Hz; CL = 1000; % basic cycle length (ms)
%load yfin_endo_female_2Hz; CL = 500; % basic cycle length (ms)

y0 = yfinal;
%CaMKt=y(58)
%mL=y(61)
%hL=y(60)
%hLp=y(59)
y0(58)= 0;%CaMKt
y0(59)= 1;%hLp
y0(60)= 1;%hL
y0(61)= 0;%mL

APD90s = [];

% Simulate mutation with 1, WT otherwise
mutation_flag_LQT8 = 0;
mutation_flag_LQT2 = 0;
mutation_flag_LQT3 = 0;

% Changes induced by mutation
if mutation_flag_LQT8 == 1
    mutation_change_GCaL_LQT8 = 0.05; % introduce mutation effects here!
    mutation_change_vCaL_LQT8 = 3; % introduce mutation effects here!
else
    mutation_change_GCaL_LQT8 = 0;
    mutation_change_vCaL_LQT8 = 0;
end

if mutation_flag_LQT2 ==1
    mutation_change_GKr = 0.99;
else
    mutation_change_GKr = 0;
end

if mutation_flag_LQT3 == 1
    mutation_change_GNaL = 30.0;
else
    mutation_change_GNaL= 0;
end

% Female with 1, male 0
female_flag = 0;

% Changes in female vs male
if female_flag == 1
    female_change_Gto =0;% -0.5; 
    female_change_GKr = 0;%-0.2;
    female_change_GKs = 0;%-0.2; 
    female_change_GK1 = 0;%-0.2; 
    female_change_vNCX = 0;%0.15; 
    female_change_vPMCA = 0;%0.8; 
else
	female_change_Gto = 0;
    female_change_GKr = 0;
    female_change_GKs = 0;
    female_change_GK1 = 0;
    female_change_vNCX = 0;
    female_change_vPMCA = 0;
end    
% Current multipliers
par_SA = ones(1,19);
% INa_Multiplier = par_SA(1); ICaL_Multiplier = par_SA(2);
% Itof_Multiplier = par_SA(3); Itos_Multiplier = par_SA(4);
% IKr_Multiplier = par_SA(5); IKs_Multiplier = par_SA(6);
% IKp_Multiplier = par_SA(7); IK1_Multiplier = par_SA(8);
% IClCa_Multiplier = par_SA(9); IClB_Multiplier = par_SA(10);
% INaB_Multiplier = par_SA(11); ICaB_Multiplier = par_SA(12);
% INaK_Multiplier = par_SA(13); INaCa_Multiplier = par_SA(14);
% IpmCa_Multiplier = par_SA(15); Jup_Multiplier = par_SA(16);
% Jrel_Multiplier = par_SA(17); Jleak_Multiplier = par_SA(18);

% Cell-type
cellType = 0; % EPI with 1, ENDO 0

% Parameter array for passing nondefault conditions
p = [cellType CL mutation_flag_LQT8 mutation_change_GCaL_LQT8 mutation_change_vCaL_LQT8...
    female_flag female_change_Gto female_change_GKr female_change_GKs...
    female_change_GK1 female_change_vNCX female_change_vPMCA mutation_flag_LQT2 mutation_change_GKr mutation_flag_LQT3 mutation_change_GNaL par_SA];

%% Single Run Simulation
tic
tspan = [0; 25e4]; % change simulation duration duration here (ms)
options = odeset('RelTol',1e-5,'MaxStep',1); 
[t,y] = ode15s(@Grandi_model_INaL,tspan,y0,options,p);
yfinal = y(end,:);
for i= 1:size(y,1)
    update_step_i =  Grandi_model_LQT(t(i), y(i,:), p, 'currents'); 
    ICaL(i,1) = update_step_i(1); %calcium current
    INa(i,1) = update_step_i(2);  %sodium current
    INaL(i,1) = update_step_i(3);  %sodium current
    IKr(i,1) = update_step_i(4); %rectified potassium current
    IKs(i,1)= update_step_i(5);  %slow potassium current
    Itos(i,1)= update_step_i(6);
    Itof(i,1)= update_step_i(7);
    IK1(i,1)=update_step_i(8);
    INaK(i,1)= update_step_i(9);
    INCX(i,1)= update_step_i(10);
      % plot currents = [I_Catot I_Na I_NaL I_kr I_ks I_tos I_tof I_ki I_nak I_ncx];
      
end
toc
%% Extractions
Vm= y(:,39); % column 39 is Vm in state vector
Ca_i= y(:,38); %column 41 corresponds to Ca_i
Na_i=y(:,34);
adjustedCa_i=10^6 *(Ca_i);

female_color = [23/255, 190/255, 187/255]; 
male_color = [239/255, 62/255, 54/255]; 

if female_flag == 0
    colors=male_color;
elseif female_flag == 1
    colors=female_color;
end

%% Plot
% Extract last two beats
start_index = find(t >= t(end) - 2*CL, 1);
T_last_two = t(start_index:end) - t(start_index);

lastBeatTime = t >= t(end) - CL;
APD90 = calculateAPD(t(lastBeatTime), Vm(lastBeatTime), 90);
APD90s = [APD90s, APD90]; 
disp(['APD90:', num2str(APD90)]);

Vm_last_two = Vm(start_index:end);
INaL_last_two= INaL(start_index:end);
adjustedCa_i_last_two = adjustedCa_i(start_index:end);

%% Plot Membrane Potential
figure(1), set(gcf, 'color', 'w'); hold on;
plot(T_last_two-CL, Vm_last_two, 'LineWidth', 3);%, 'LineStyle', '--'); 
xlabel('Time (ms)');
ylabel('Voltage (mV)');
xlim([-50 1000]);
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 28);
set(gca, 'LineWidth', 1.5);
hold off;
% 
% %% Plot Calcium Transient
% figure(2), set(gcf, 'color', 'w'); hold on;
% plot(T_last_two-CL, adjustedCa_i_last_two, 'Color', colors, 'LineWidth', 3);%, 'LineStyle', '--'); 
% xlabel('Time (ms)');
% ylabel('[Ca]_i (nM)');
% xlim([-50 2000]);
% set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
% set(gca, 'FontSize', 28);
% set(gca, 'LineWidth', 1.5);
% hold off;

figure(2), set(gcf, 'color', 'w'); hold on;
plot(T_last_two-CL,INaL_last_two,'LineWidth', 3);
ylabel('INaL');
xlabel('Time (ms)');
xlim([-50 1000]);
ylim([-1 0]);
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 28);
set(gca, 'LineWidth', 1.5);
hold off;

%% Saving
%save yfin_endo_0p25Hz yfinal
%save yfin_endo_0p33Hz yfinal
%save yfin_endo_0p5Hz yfinal
%save yfin_endo_1Hz yfinal
%save yfin_endo_2Hz yfinal

%save yfin_endo_female_0p25Hz yfinal
%save yfin_endo_female_0p33Hz yfinal
%save yfin_endo_female_0p5Hz yfinal
%save yfin_endo_female_1Hz yfinal
%save yfin_endo_female_2Hz yfinal
