%%  ToR-ORd Model Rate Script

%This script runs the model at different BCLs in order to create the rate
%dependant plots that are show in fig.1 2 and 4

clear;
clc;

% Define parameters
sex_types = {'male', 'female'};
cycle_lengths = [500,1000,2000,3000,4000];%linspace(200, 500, 26);
num_conditions = length(sex_types);
num_cycle_lengths = length(cycle_lengths);

%ikrMultipliers=0.90;
ikrMultipliers=1;

vCaLChanges = 6;
icalMultipliers= 1.45;

%vCaLChanges = 0;
%icalMultipliers= 1;

param = repmat(struct(), [length(ikrMultipliers), 1]);
for iParam = 1:length(ikrMultipliers)
    param(iParam).IKr_Multiplier = ikrMultipliers(iParam); 
end
param = repmat(param, [length(icalMultipliers), 1]);
for iParam = 1:length(icalMultipliers)
    param(iParam).ICaL_Multiplier = icalMultipliers(iParam);
    param(iParam).VCaL_Change = vCaLChanges;
end

% Initialize arrays for results
APD90s_beat1 = zeros(num_cycle_lengths, num_conditions);
APD90s_beat2 = zeros(num_cycle_lengths, num_conditions);
APD50s = zeros(num_cycle_lengths, num_conditions);
CaT_Amplitudes_beat1 = zeros(num_cycle_lengths, num_conditions);
CaT_Amplitudes_beat2 = zeros(num_cycle_lengths, num_conditions);
CaT_Duration50 = zeros(num_cycle_lengths, num_conditions);
Diastolic_Ca = zeros(num_cycle_lengths, num_conditions);
SR_Content = zeros(num_cycle_lengths, num_conditions);

% Loop through cycle lengths and conditions
for iCL = 1:num_cycle_lengths
    CL = cycle_lengths(iCL);

    for sex_idx = 1:num_conditions
        % Initialize parameters for simulation
        param.model = @model_Torord;
        param.sexType = sex_idx;
        param.cellType = 0; % Endo
        param.bcl = CL;

        % Initial conditions and simulation setup
        initialState = 'Torord_endo';
        X0 = getStartingState(initialState);
        beats = 400;
        ignoreFirst = beats - 2;

        [time, X] = modelRunner(X0, [], param, beats, ignoreFirst);
        currents = getCurrentsStructure(time, X, param, 0);

        beatDuration = param.bcl;
        numBeats = floor(currents.time(end) / beatDuration);

        % Analyze last 2 beats for APD90 and APD50
        APD90_values = [];
        for j = numBeats-1:numBeats
            beatStart = (j-1) * beatDuration;
            beatEnd = j * beatDuration;
            indices = find(currents.time >= beatStart & currents.time < beatEnd);

            beatTime = currents.time(indices) - beatStart;
            beatVoltage = currents.V(indices);

            APD90_values = [APD90_values, calculateAPD(beatTime, beatVoltage, 90)];
            APD50s(iCL, sex_idx) = calculateAPD(beatTime, beatVoltage, 50);
        end
        APD90s_high(iCL, sex_idx) = max(APD90_values); % Highest APD90
        APD90s_low(iCL, sex_idx) = min(APD90_values);  % Lowest APD90
        sortedAPD90s = sort(APD90_values);
        APD90s_beat1(iCL, sex_idx) = sortedAPD90s(1);
        APD90s_beat2(iCL, sex_idx) = sortedAPD90s(2);

        % Analyze last 2 beats for CaT Amplitude and Duration
        Cai_all = currents.Cai * 1e6; 
        CaT_Amplitudes = [];
        for j = numBeats-1:numBeats
            beatStart = (j-1) * beatDuration;
            beatEnd = j * beatDuration;
            indices = find(currents.time >= beatStart & currents.time < beatEnd);

            beatCai = Cai_all(indices);
            CaT_Amplitudes = [CaT_Amplitudes, max(beatCai) - min(beatCai)];

            if j == numBeats
                CaT_Duration50(iCL, sex_idx) = calculateCaTDuration(currents.time(indices), beatCai, 50);
            end
        end
        sortedCaTAmps = sort(CaT_Amplitudes);
        CaT_Amplitudes_beat1(iCL, sex_idx) = sortedCaTAmps(1);
        CaT_Amplitudes_beat2(iCL, sex_idx) = sortedCaTAmps(2);

        % Calculate Diastolic Ca and SR Content
        last_beat_indices = find(currents.time >= currents.time(end) - beatDuration);
        Diastolic_Ca(iCL, sex_idx) = min(Cai_all(last_beat_indices));
        SR_Content(iCL, sex_idx) = mean(currents.CaNSR(last_beat_indices)) + mean(currents.CaJSR(last_beat_indices));
    end
end

% Plot results
female_endo_color = [23/255, 190/255, 187/255];
male_endo_color = [239/255, 62/255, 54/255];
colors = {male_endo_color, female_endo_color};

%% APD90 Plot
figure(1), set(gcf, 'color', 'w'); hold on;
xlabel('Cycle Length (ms)');
ylabel('APD90 (ms)');
title('APD90 vs. Cycle Length');
for sex_idx = 1:num_conditions
    plot(cycle_lengths, APD90s_beat1(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
    plot(cycle_lengths, APD90s_beat2(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
hold off;

%% APD50 Plot
figure(2), set(gcf, 'color', 'w'); hold on;
xlabel('Cycle Length (ms)');
ylabel('APD50 (ms)');
title('APD50 vs. Cycle Length');
for sex_idx = 1:num_conditions
    plot(cycle_lengths, APD50s(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
hold off;

%% CaT Amplitude Plot
figure(3), set(gcf, 'color', 'w'); hold on;
xlabel('Cycle Length (ms)');
ylabel('CaT Amplitude (μM)');
title('CaT Amplitude vs. Cycle Length');
for sex_idx = 1:num_conditions
    plot(cycle_lengths, CaT_Amplitudes_beat1(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
    plot(cycle_lengths, CaT_Amplitudes_beat2(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
hold off;


%% CaT Duration Plot
figure(4), set(gcf, 'color', 'w'); hold on;
xlabel('Cycle Length (ms)');
ylabel('CaT Duration 50% (ms)');
title('CaT Duration 50% vs. Cycle Length');
for sex_idx = 1:num_conditions
    plot(cycle_lengths, CaT_Duration50(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
hold off;

%% Diastolic Ca2+ Plot
figure(5), set(gcf, 'color', 'w'); hold on;
xlabel('Cycle Length (ms)');
ylabel('Diastolic [Ca^{2+}]_{cyto} (μM)');
title('Diastolic Ca^{2+} vs. Cycle Length');
for sex_idx = 1:num_conditions
    plot(cycle_lengths, Diastolic_Ca(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
hold off;

%% SR Content Plot
figure(6), set(gcf, 'color', 'w'); hold on;
xlabel('Cycle Length (ms)');
ylabel('SR Content (mM)');
title('SR Content vs. Cycle Length');
for sex_idx = 1:num_conditions
    plot(cycle_lengths, SR_Content(:, sex_idx), 'o-', 'LineWidth', 3.5, 'Color', colors{sex_idx});
end
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
set(gca, 'FontSize', 30);
set(gca, 'LineWidth', 1.5);
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
