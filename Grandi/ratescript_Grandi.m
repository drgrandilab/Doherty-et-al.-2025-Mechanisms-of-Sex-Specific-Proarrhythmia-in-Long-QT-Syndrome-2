%% Grandi model Rate Scripts
%This script runs the model at different BCLs in order to create the rate
%dependant plots 

% Initialization
clear
clc

% Define sex and cycle lengths
sex_types = {'male', 'female'};
cycle_lengths = [500, 1000, 2000, 3000, 4000]; % Define your desired cycle lengths
num_cycle_lengths = length(cycle_lengths);
frequency_labels = repmat({'2Hz'}, 1, 18);
%frequency_labels = {'2Hz', '1Hz', '0p5Hz', '0p33Hz', '0p25Hz'};
colors = {[239/255, 62/255, 54/255], [23/255, 190/255, 187/255]};  % Male Endo, Female Endo
plot_labels = {'Male Endo', 'Female Endo'};

% Array to store results
APD90s = zeros(num_cycle_lengths, 2);
APD50s = zeros(num_cycle_lengths, 2);
CaTAmp = zeros(num_cycle_lengths, 2);
CaD50 = zeros(num_cycle_lengths, 2);
SR_Content = zeros(num_cycle_lengths, 2);
Diastolic_Ca = zeros(num_cycle_lengths, 2);

% Main simulation loop
for iCL = 1:num_cycle_lengths
    CL = cycle_lengths(iCL);
    freq_label = frequency_labels{iCL};

    for sex_idx = 1:2
        if strcmp(sex_types{sex_idx}, 'male')
            female_flag = 0;
            load_prefix = ['yfin_endo_' freq_label]; % Construct file name with frequency label
        else
            female_flag = 1;
            load_prefix = ['yfin_endo_female_' freq_label]; % Construct file name with frequency label
        end

        % Adjust file loading based on CL
        load([load_prefix]);
        cellType = 0; % Endo

        % Initial conditions and model setup
        y0 = yfinal;
        y0(58) = 0; %CaMKt
        y0(59) = 1; %hLp
        y0(60) = 1; %hL
        y0(61) = 0; %mL
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
            mutation_change_GKr = 1.0;
        else
            mutation_change_GKr = 0;
        end

        if mutation_flag_LQT3 == 1
            mutation_change_GNaL = 20.0;
        else
            mutation_change_GNaL= 0;
        end

        % Changes in female vs male
        if female_flag == 1
            female_change_Gto = 0;%-0.5; 
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

        % Parameter array for simulation
        p = [cellType CL mutation_flag_LQT8 mutation_change_GCaL_LQT8 mutation_change_vCaL_LQT8...
        female_flag female_change_Gto female_change_GKr female_change_GKs...
        female_change_GK1 female_change_vNCX female_change_vPMCA mutation_flag_LQT2 mutation_change_GKr mutation_flag_LQT3 mutation_change_GNaL par_SA];

        tspan = [0; 200*CL];
        options = odeset('RelTol',1e-5,'MaxStep',1);
        [t, y] = ode15s(@Grandi_model_LQT, tspan, y0, options, p);
       for beat_num = 1:4
            beat_start_time = t(end) - beat_num * CL;
            beat_end_time = t(end) - (beat_num - 1) * CL;

            beat_indices = find(t >= beat_start_time & t < beat_end_time);

            if ~isempty(beat_indices)
                APD90 = calculateAPD(t(beat_indices), y(beat_indices, 39), 90);
                if isscalar(APD90)
                    APD90s(iCL, sex_idx, beat_num) = APD90;
                else
                    APD90s(iCL, sex_idx, beat_num) = NaN;
                end
            else
                APD90s(iCL, sex_idx, beat_num) = NaN;
            end
        end

        % Analyze the last two beats for CaT Amplitude
        last_beat_indices = find(t >= t(end) - CL, 1):length(y(:, 39));
        sec_last_beat_start_idx = find(t >= t(end) - 2*CL, 1);
        sec_last_beat_end_idx = find(t >= t(end) - CL, 1) - 1;
        sec_last_beat_indices = sec_last_beat_start_idx:sec_last_beat_end_idx;

        Cai_last_beat = y(last_beat_indices, 38)*1e6; 
        Cai_sec_last_beat = y(sec_last_beat_indices, 38)*1e6; 

        CaTamp_beat1 = max(Cai_last_beat) - min(Cai_last_beat);
        CaTamp_beat2 = max(Cai_sec_last_beat) - min(Cai_sec_last_beat);

        sortedCaTAmps = sort([CaTamp_beat1, CaTamp_beat2]);
        CaT_Amplitudes_beat1(iCL, sex_idx) = sortedCaTAmps(1);
        CaT_Amplitudes_beat2(iCL, sex_idx) = sortedCaTAmps(2);


        % Extract relevant data
        Vm_last_beat = y(last_beat_indices, 39);
        Ca_i_last_beat = y(last_beat_indices, 38) * 1e6; 
        Ca_sr = y(last_beat_indices, 31);
        % Calculate metrics
        APD90s(iCL, sex_idx) = calculateAPD(t(last_beat_indices), Vm_last_beat, 90);
        APD50s(iCL, sex_idx) = calculateAPD(t(last_beat_indices), Vm_last_beat, 50);
        %APD90s(iCL, sex_idx) = calculateAPD(t(second_to_last_beat_indices), Vm_second_to_last_beat, 90);
        %APD50s(iCL, sex_idx) = calculateAPD(t(second_to_last_beat_indices), Vm_second_to_last_beat, 50);
        CaTAmp(iCL, sex_idx) = max(Ca_i_last_beat) - min(Ca_i_last_beat);
        CaD50(iCL, sex_idx) = calculateCaTDuration(t(last_beat_indices), Ca_i_last_beat, 50);
        Diastolic_Ca(iCL, sex_idx) = min(Ca_i_last_beat);
        SR_Content(iCL, sex_idx) = mean(Ca_sr);

    end
end
%% Plot APD90
figure(1), set(gcf, 'color', 'w'); hold on;
xlabel('Cycle Length (ms)');
ylabel('APD90 (ms)');
title('APD90 vs. Cycle Length');
for sex_idx = 1:2
    low_APD90s = min(APD90s(:, sex_idx, :), [], 3);
    high_APD90s = max(APD90s(:, sex_idx, :), [], 3);

    plot(cycle_lengths, low_APD90s, 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
    plot(cycle_lengths, high_APD90s, 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
hold off;
 
%%  Plot APD50
figure(2),set(gcf,'color','w'); hold on 
xlabel('Cycle Length (ms)');
ylabel('APD_{50} (ms)');
xlim([300 4000]);
% APD50s(4,1)=987;
% APD50s(4,2)=1467;

for sexType = 1:2
    plot(cycle_lengths, APD50s(:, sexType), 'o-', 'LineWidth', 3.5, 'Color', colors{sexType});
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
for sex_idx = 1:2
    plot(cycle_lengths, CaT_Amplitudes_beat1(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
    plot(cycle_lengths, CaT_Amplitudes_beat2(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
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
for sexType = 1:2
    plot(cycle_lengths, CaD50(:, sexType), 'o-', 'LineWidth', 3.5, 'Color', colors{sexType});
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
for sexType = 1:2
    plot(cycle_lengths, Diastolic_Ca(:, sexType), 'o-', 'LineWidth', 3.5, 'Color', colors{sexType});
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
%ylim([0.52 0.64]);
for sexType = 1:2
    plot(cycle_lengths, SR_Content(:, sexType), 'o-', 'LineWidth', 3.5, 'Color', colors{sexType});
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
