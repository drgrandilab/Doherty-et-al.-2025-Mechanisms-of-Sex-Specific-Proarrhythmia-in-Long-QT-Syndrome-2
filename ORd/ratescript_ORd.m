%% ORd Rate Script
%This script runs the model at different BCLs in order to create the rate
%dependant plots 

% Initial conditions for state variables
clear 
clc

v=-87; nai=7; nass=nai; ki=145; kss=ki; cai=1.0e-4; cass=cai; cansr=1.2; cajsr=cansr; m=0;
hf=1; hs=1; j=1; hsp=1; jp=1; mL=0; hL=1; hLp=1; a=0; iF=1; iS=1; ap=0; iFp=1; iSp=1; d=0;
ff=1; fs=1; fcaf=1; fcas=1; jca=1; nca=0; ffp=1; fcafp=1; xrf=0; xrs=0; xs1=0; xs2=0; xk1=1;
Jrelnp=0; Jrelp=0; CaMKt=0;

% Run model
% X0 is the vector for initial conditions for state variables
options = odeset('MaxStep', 1,'InitialStep',2e-2);

stim1_amp = -15;  %uA
stim1_dur = 5;  % Stimulus duration in ms
beats= 200;    %number of beats in the simulation

parameter_inputs.stim1_amp = stim1_amp;
parameter_inputs.stim1_dur = stim1_dur;
parameter_inputs.INa_scale = 1;
parameter_inputs.IK1_scale = 1;
parameter_inputs.ICaL_scale = 1;
parameter_inputs.IKr_scale = 1;
parameter_inputs.INaL_scale = 1;
%LQT 8- L-type calcium gain of function in CACNA1C gene
parameter_inputs.mutation_flag_LQT8= 0;
parameter_inputs.mutation_change_GCaL = 1.20; % percent change
parameter_inputs.mutation_change_vCaL = 6; %constant change
%LQT 2- IKr loss of function in hERG
parameter_inputs.mutation_flag_LQT2= 0;
parameter_inputs.mutation_change_GKr = 0.1; %-25,50,75,100
%LQT 3- INaL gain of function in SCN5A gene
parameter_inputs.mutation_flag_LQT3= 0;
parameter_inputs.mutation_change_GNaL = 25.0; %+25,50,75,100

cycleLengths = [500,1000,2000,3000,4000];%klinspace(200,500,21);
numCycleLengths = length(cycleLengths);
numConditions = 2;  % male endo, female endo (can alter the script to also include epi)

APD90s_beat1 = zeros(numCycleLengths, numConditions);
APD90s_beat2 = zeros(numCycleLengths, numConditions);
APD50s_beat1 = zeros(numCycleLengths, numConditions);
APD50s_beat2 = zeros(numCycleLengths, numConditions);
CaT_Amplitudes_beat1 = zeros(numCycleLengths, numConditions);
CaT_Amplitudes_beat2 = zeros(numCycleLengths, numConditions);
CaT_Amplitudes = zeros(numCycleLengths, numConditions);
CaT_Durations = zeros(numCycleLengths, numConditions);
Diastolic_Ca = zeros(numCycleLengths, numConditions);
SR_Content = zeros(numCycleLengths, numConditions);

T1 = []; Vm1 = []; X_all=[];


for iCL = 1:numCycleLengths
    for sexType = 1:2  % 1 for male, 2 for female
        CL = cycleLengths(iCL);
        parameter_inputs.CL = CL;

        X0=[v nai nass ki kss cai cass cansr cajsr m hf hs j hsp jp mL hL hLp a iF iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp xrf xrs xs1 xs2 xk1 Jrelnp Jrelp CaMKt]';

        T1 = []; Vm1 = []; X_all = []; Cai=[]; Cansr=[]; Cajsr=[];
        for n = 1:beats
            [time, X] = ode15s(@model_ORd11,[0 CL],X0,options,1, parameter_inputs, sexType, 0); % Endo cell type
            X0 = X(end,:);
            T1 = [T1; time(2:end) + (n-1)*CL];
            Vm1 = [Vm1; X(2:end,1)]; 
            Cai = [Cai; X(2:end,6)];
            Cansr = [Cansr; X(2:end,8)];
            Cajsr = [Cajsr; X(2:end,9)];

        end
        

            last_beat_indices = find(T1 >= T1(end) - CL, 1):length(Vm1);
            sec_last_beat_start_idx = find(T1 >= T1(end)-2*CL, 1);
            sec_last_beat_end_idx = find(T1 >= T1(end)-CL, 1) - 1;
            sec_last_beat_indices = sec_last_beat_start_idx:sec_last_beat_end_idx;

            APD90_beat1 = calculateAPD(T1(last_beat_indices), Vm1(last_beat_indices), 90);
            APD90_beat2 = calculateAPD(T1(sec_last_beat_indices), Vm1(sec_last_beat_indices), 90);
            APD90s_high(iCL, sexType) = max(max(APD90_beat1(:)), max(APD90_beat2(:)));
            APD90s_low(iCL, sexType) = min(min(APD90_beat1(:)), min(APD90_beat2(:)));
            APD50_beat1 = calculateAPD(T1(last_beat_indices), Vm1(last_beat_indices), 50);
            APD50_beat2 = calculateAPD(T1(sec_last_beat_indices), Vm1(sec_last_beat_indices), 50);

            CaTamp_beat1 = max(Cai(last_beat_indices)) - min(Cai(last_beat_indices));
            CaTamp_beat2 = max(Cai(sec_last_beat_indices)) - min(Cai(sec_last_beat_indices));

            sortedAPD90s = sort([APD90_beat1, APD90_beat2]);
            APD90s_beat1(iCL, sexType) = sortedAPD90s(1);
            APD90s_beat2(iCL, sexType) = sortedAPD90s(2);

            sortedAPD50s = sort([APD50_beat1, APD50_beat2]);
            APD50s_beat1(iCL, sexType) = sortedAPD50s(1);
            APD50s_beat2(iCL, sexType) = sortedAPD50s(2);

            sortedCaTAmps = sort([CaTamp_beat1, CaTamp_beat2]);
            CaT_Amplitudes_beat1(iCL, sexType) = sortedCaTAmps(1);
            CaT_Amplitudes_beat2(iCL, sexType) = sortedCaTAmps(2);

            APD90s(iCL, sexType) = calculateAPD(T1(last_beat_indices), Vm1(last_beat_indices), 90);
            APD50s(iCL, sexType) = calculateAPD(T1(last_beat_indices), Vm1(last_beat_indices), 50);

            % Calculate CaT Amplitude, Duration, Diastolic Ca2+, and SR Content
            CaT_Amplitudes(iCL, sexType) = calculateCaTAmplitude(Cai(last_beat_indices)*(1e6));
            CaT_Durations(iCL, sexType) = calculateCaTDuration(T1(last_beat_indices),Cai(last_beat_indices),50);
            Diastolic_Ca(iCL, sexType) = min(Cai(last_beat_indices)*1e6);
            SR_Content(iCL, sexType) = mean(Cansr) + mean(Cajsr);  % Total SR content
   
    end
end
v=X(:,1);
X0 = X(size(X,1),:);

female_endo_color = [23/255, 190/255, 187/255]; 
male_endo_color = [239/255, 62/255, 54/255]; 

colors = {male_endo_color, female_endo_color};
legendEntries = {'Male Endo', 'Female Endo'};

%% Plot APD90
figure(1),set(gcf,'color','w'); hold on 
xlabel('Cycle Length (ms)');
ylabel('APD_{90} (ms)');
%xlim([300 4000]);
%ylim([0 2000]);
%2 beats allows for alternan exploration
for sexType = 1:2
   plot(cycleLengths, APD90s_beat1(:, sexType), '-', 'LineWidth', 3.5, 'Color', colors{sexType});
   plot(cycleLengths, APD90s_beat2(:, sexType), '-', 'LineWidth', 3.5, 'Color', colors{sexType});
end

set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
title('APD90 vs. Cycle Length');
% hold off;

%% Plot APD50
figure(2),set(gcf,'color','w'); hold on 
xlabel('Cycle Length (ms)');
ylabel('APD_{50} (ms)');
%xlim([300 4000]);
%ylim([150 2000])
for sexType = 1:numConditions
    plot(cycleLengths, APD50s_beat1(:, sexType), 'o-', 'LineWidth', 3.5, 'Color', colors{sexType});
    plot(cycleLengths, APD50s_beat2(:, sexType), 'o-', 'LineWidth', 3.5, 'Color', colors{sexType});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
title('APD50 vs. Cycle Length');
hold off;

%% Plot CaT Amplitude
figure(3),set(gcf,'color','w'); hold on 
xlabel('Cycle Length (ms)');
ylabel('CaT Amplitude (nM)');
%xlim([300 4000]);
%ylim([0 500])
for sexType = 1:2
  plot(cycleLengths, CaT_Amplitudes_high(:, sexType), '-', 'LineWidth', 3.5, 'Color', colors{sexType});
  plot(cycleLengths, CaT_Amplitudes_low(:, sexType), '-', 'LineWidth', 3.5, 'Color', colors{sexType});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
title('CaT Amplitude vs. Cycle Length');
hold off;

%% Plot CaT Duration
figure(4),set(gcf,'color','w'); hold on 
xlabel('Cycle Length (ms)');
ylabel('CaD_{50} (ms)');
xlim([300 4000]);
%ylim([-50 2000]);
for sexType = 1:2
    plot(cycleLengths, CaT_Durations(:, sexType), 'o-', 'LineWidth', 3.5, 'Color', colors{sexType});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
title('CaT Duration 50 vs. Cycle Length');
hold off;

%% Plot Diastolic Ca2+
figure(5),set(gcf,'color','w'); hold on 
xlabel('Cycle Length (ms)');
ylabel('Diastolic [Ca^{2+}]_{cyto} (nM)');
xlim([300 4000]);
ylim([60 140]);
for sexType = 1:2
    plot(cycleLengths, Diastolic_Ca(:, sexType), 'o-', 'LineWidth', 3.5, 'Color', colors{sexType});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
title('Diastolic [Ca^{2+}]_{cyto} vs. Cycle Length');
hold off;

%% Plot SR Content
figure(6),set(gcf,'color','w'); hold on 
xlabel('Cycle Length (ms)');
ylabel('SR Content (mM)');
xlim([300 4000]);
ylim([2.2 3.0]);
for sexType = 1:2
    plot(cycleLengths, SR_Content(:, sexType), 'o-', 'LineWidth', 3.5, 'Color', colors{sexType});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
title('SR Content vs. Cycle Length');
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


function catDuration = calculateCaTDuration(time, Cai, percentage)
    % Handle if time is a cell array
    if iscell(time)
        time = cell2mat(time);
    end

    % Calculate baseline and peak calcium concentration
    Cmin = min(Cai);
    Cmax = max(Cai);

    % Calculate the threshold
    Cth = Cmin + (Cmax - Cmin) * (1 - percentage / 100);

    % Find indices where calcium exceeds the threshold
    aboveTh = find(Cai >= Cth);
    if isempty(aboveTh)
        catDuration = NaN;
        return;
    end

    % Starting index for the repolarization phase
    startTh = aboveTh(1);

    % Find the first index where calcium falls below the threshold after the peak
    belowTh = find(Cai(startTh:end) < Cth, 1, 'first');

    if isempty(belowTh)
        catDuration = NaN;  % Return NaN if Cai does not fall below the threshold
        return;
    end

    % Calculate the ending index for the repolarization phase
    endTh = startTh + belowTh - 1;

    % Check if the calculated index exceeds the length of the time array
    if endTh > length(time)
        endTh = length(time); 
    end

    % Calculate the duration from rise above to fall below the threshold
    catDuration = time(endTh) - time(startTh);
end


function catAmplitude = calculateCaTAmplitude(Cai)
    % Calculate the peak and baseline of calcium concentration
    baseline = min(Cai);  
    peak = max(Cai);      
    % Calculate the amplitude of the calcium transient
    catAmplitude = peak - baseline;
end
