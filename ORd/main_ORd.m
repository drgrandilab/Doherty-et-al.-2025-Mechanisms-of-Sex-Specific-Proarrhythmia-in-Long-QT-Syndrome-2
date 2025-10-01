% Copyright (c) 2011-2015 by Thomas O'Hara, Yoram Rudy, Washington University in St. Louis.
% All rights reserved.
% MATLAB Implementation of the O'Hara-Rudy dynamic (ORd) model for the
% undiseased human ventricular action potential and calcium transient
% The ORd model is described in the article "Simulation of the Undiseased
% Human Cardiac Ventricular Action Potential: Model Formulation and
% Experimental Valirudation"
% by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
% The article and supplemental materails are freely available in the
% Open Access jounal PLoS Computational Biology
% Link to Article:
% http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061
% % Email: tom.ohara@gmail.com / rudy@wustl.edu
% Web: http://rudylab.wustl.edu

clear; clc; %close all

%% Initial conditions for state variables
v=-87; nai=7; nass=nai; ki=145; kss=ki; cai=1.0e-4; cass=cai; cansr=1.2; cajsr=cansr; m=0;
hf=1; hs=1; j=1; hsp=1; jp=1; mL=0; hL=1; hLp=1; a=0; iF=1; iS=1; ap=0; iFp=1; iSp=1; d=0;
ff=1; fs=1; fcaf=1; fcas=1; jca=1; nca=0; ffp=1; fcafp=1; xrf=0; xrs=0; xs1=0; xs2=0; xk1=1;
Jrelnp=0; Jrelp=0; CaMKt=0;

%% Run model
% X0 is the vector for initial sconditions for state variables
X0=[v nai nass ki kss cai cass cansr cajsr m hf hs j hsp jp mL hL hLp a iF iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp xrf xrs xs1 xs2 xk1 Jrelnp Jrelp CaMKt]';
options = odeset('MaxStep', 1,'InitialStep',2e-2);

stim1_amp = -15;  %uA
stim1_dur = 5;  % Stimulus duration in ms
CL= 1000;     %pacing cycle length in ms
beats= 200;    %number of beats in the simulation


parameter_inputs.CL = CL;
parameter_inputs.stim1_amp = stim1_amp;
parameter_inputs.stim1_dur = stim1_dur;
parameter_inputs.INa_scale = 1;
parameter_inputs.IK1_scale = 1;
parameter_inputs.ICaL_scale = 1;
parameter_inputs.IKr_scale = 1;
parameter_inputs.INaL_scale = 1;

%% Simulate LQT Effects here
%LQT 8- L-type calcium gain of function in CACNA1C gene
parameter_inputs.mutation_flag_LQT8= 1;
 %0 mutation off 1 mutation ICaL on
parameter_inputs.mutation_change_GCaL = 1.2; % percent change
parameter_inputs.mutation_change_vCaL = 6; %constant change
%LQT 2- IKr loss of function in hERG
parameter_inputs.mutation_flag_LQT2= 0; %0 mutation off 1 mutation IKr on
parameter_inputs.mutation_change_GKr = 0.10; %0.8
%LQT 3- INaL gain of function in SCN5A gene
parameter_inputs.mutation_flag_LQT3= 0; %0 mutation off 1 mutation INa on
parameter_inputs.mutation_change_GNaL = 25.0; %+25,50,75,100

%T1: Complete time vector;  Vm1: Complete Voltage vector
% time: Last CL; v: Voltage corresponding to last CL
% Function model.m is the equation file; ODE solved using MATLAB ODE15s function

%1 male 2 female; 
gendertype= 2;
% 0 endo 1 epi
celltype = 0;

T1 = []; Vm1 = []; X_all=[];
for n= 1:beats
    [time, X] = ode15s(@model_ORd,[0 CL],X0,[],1, parameter_inputs, gendertype, celltype);

    X0 = X(size(X,1),:);
    voltage = X(:,1);
    T1 = [T1; (time(2:end) + (n-1)*CL)];
    Vm1 = [Vm1 ; voltage(2:end)]; X_all = [X_all; X((2:end),:)];
end
v=X(:,1);
X0 = X(size(X,1),:);

nai_vec= X_all(:, 2); 
nass_vec= X_all(:, 3);
ki_vec= X_all(:, 4);
kss_vec= X_all(:, 5);
cai_vec= X_all(:, 6);
cass_vec= X_all(:, 7);
cansr_vec= X_all(:, 8);
cajsr_vec= X_all(:, 9);
nca_vec= X_all(:, 31);


%[INa INaL Ito ICaL IKr IKs IK1 INaCa_i INaCa_ss INaK IKb INab ICab IpCa Jdiff JdiffNa JdiffK Jup Jleak Jtr Jrel CaMKa Istim ICaNa ICaK PhiCaL PhiCaNa PhiCaK];

INa_vec = zeros(length(T1), 1); INaL_vec = zeros(length(T1), 1); Ito_vec = zeros(length(T1), 1);
ICaL_vec = zeros(length(T1), 1); IKr_vec = zeros(length(T1), 1); IKs_vec = zeros(length(T1), 1);
IK1_vec = zeros(length(T1), 1); INaCai_vec = zeros(length(T1), 1); INaCass_vec = zeros(length(T1), 1);
INaK_vec = zeros(length(T1), 1); IKb_vec = zeros(length(T1), 1); INab_vec = zeros(length(T1), 1);
ICab_vec = zeros(length(T1), 1); IpCa_vec = zeros(length(T1), 1); Jdiff_vec = zeros(length(T1), 1);
JdiffNa_vec = zeros(length(T1), 1); JdiffK_vec = zeros(length(T1), 1); Jup_vec = zeros(length(T1), 1);
Jleak_vec = zeros(length(T1), 1); Jtr_vec = zeros(length(T1), 1); Jrel_vec = zeros(length(T1), 1);
CaMKa_vec = zeros(length(T1), 1); Istim_vec = zeros(length(T1), 1); ICaNa_vec = zeros(length(T1), 1);
ICaK_vec = zeros(length(T1), 1); PhiCaL_vec = zeros(length(T1), 1); PhiCaNa_vec = zeros(length(T1), 1);
PhiCaK_vec = zeros(length(T1), 1);

for i=1:size(T1,1)
    %IsJs=model_ORd11_Yang2012(T1(i),X_all(i,:),0, parameter_inputs, gendertype, celltype);
    %IsJs=model_ORd11_Alex2021_Yang(T1(i),X_all(i,:),0, parameter_inputs, gendertype, celltype);
    IsJs=model_ORd11(T1(i),X_all(i,:),0, parameter_inputs, gendertype, celltype);

    INa_vec(i)=IsJs(1); INaL_vec(i)=IsJs(2); Ito_vec(i)=IsJs(3); ICaL_vec(i)=IsJs(4); IKr_vec(i)=IsJs(5); 
    IKs_vec(i)=IsJs(6); IK1_vec(i)=IsJs(7); INaCai_vec(i)=IsJs(8); INaCass_vec(i)=IsJs(9); INaK_vec(i)=IsJs(10); 
    IKb_vec(i)=IsJs(11); INab_vec(i)=IsJs(12); ICab_vec(i)=IsJs(13); IpCa_vec(i)=IsJs(14); Jdiff_vec(i)=IsJs(15); 
    JdiffNa_vec(i)=IsJs(16); JdiffK_vec(i)=IsJs(17); Jup_vec(i)=IsJs(18); Jleak_vec(i)=IsJs(19); Jtr_vec(i)=IsJs(20); 
    Jrel_vec(i)=IsJs(21); CaMKa_vec(i)=IsJs(22); Istim_vec(i)=IsJs(23); ICaNa_vec(i)=IsJs(24); ICaK_vec(i)=IsJs(25);
    PhiCaL_vec(i)=IsJs(26); PhiCaNa_vec(i)=IsJs(27); PhiCaK_vec(i)=IsJs(28);

end
lastBeatTime = T1 >= T1(end) - CL;
% Extract data for the second-to-last beat
start_index = find(T1 >= (beats-6)*CL, 1);
end_index = length(T1);

T_last_beats = T1(start_index:end_index) - (beats-6)*CL;
Vm_last_beats = Vm1(start_index:end_index);
cai_last_beats = cai_vec(start_index:end_index);
APD90s = [];
for beat_num = 1:4
    beat_start_time = (beat_num-1)*CL;
    beat_end_time = beat_num*CL;
    beat_indices = find(T_last_beats >= beat_start_time & T_last_beats < beat_end_time);
    beat_time = T_last_beats(beat_indices);
    beat_voltage = Vm_last_beats(beat_indices);
    APD90 = calculateAPD(beat_time, beat_voltage, 90);
    APD90s = [APD90s, APD90];
end
disp(['APD90 values for the last 4 beats: ', num2str(APD90s)]);


female_color = [23/255, 190/255, 187/255]; 
male_color = [239/255, 62/255, 54/255]; 

if gendertype == 1
    colors=male_color;
elseif gendertype == 2
    colors=female_color;
end


%% Plot Membrane Voltage and Calcium Transient
figure(1),set(gcf,'color','w'); hold on
%plot(currents.time, currents.V,'LineWidth', 1.2, 'Color', 'b');
plot(T_last_beats, Vm_last_beats, 'Linewidth', 3, 'Color', colors);%, 'LineStyle', '--'); 
xlabel('Time (ms)');
ylabel('Voltage (mV)');
%xlim([-50 CL]);
%ylim([-100 50]);
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
set(gca, 'FontSize', 28); % Increase font size for tick labels
set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
hold off;

figure(2),set(gcf,'color','w');
hold on
%plot(currents.time, currents.cai_vec,'LineWidth', 1.2, 'Color', 'b');
plot(T_last_two-CL, cai_last_two * 1e6, 'Color', colors, 'Linewidth', 3);%, 'LineStyle', '--'); 
xlabel('Time (ms)');
ylabel('[Ca]_i (nM)');
xlim([-50 CL]);
ylim([50 200]);
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
set(gca, 'FontSize', 28); % Increase font size for tick labels
set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
hold off;

function apd = calculateAPD(time, V, threshold)
    Vmax = max(V);
    Vmin = min(V);
    Vth = Vmin + (Vmax - Vmin) * (1 - threshold / 100);
    
    aboveTh = find(V >= Vth);
    if isempty(aboveTh)
        apd = NaN; 
        return;
    end

    startTh = aboveTh(1);
    belowTh = find(V(startTh:end) < Vth, 1, 'first');
    
    if isempty(belowTh)
        apd = NaN; 
        return;
    end

    endTh = startTh + belowTh - 1;
    apd = time(endTh) - time(startTh);
    
end

