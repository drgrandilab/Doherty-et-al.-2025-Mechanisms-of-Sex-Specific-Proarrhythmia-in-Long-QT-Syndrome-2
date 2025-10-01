% This is a simple script which runs the control endocardial model for 100
% beats and plots membrane potential and calcium transient.

%% Setting parameters
% clear
% clc

% Param is the structure of model parameters that the user may wish to
% change compared to default simulation. The full list is given in the
% function ORdRunner, and it mainly includes cell type, current
% multipliers, extracellular ionic concentrations, or fraction of NCX and ICaL
% localisation in junctional subspace.
param.cellType = 0; %0 endo 1 epi
param.sexType = 1; %male 1 female 2
param.bcl =  1000; % basic cycle length in ms
param.model = @model_Torord;%@model_Torord; % which model is to be used - right now, use @model_Torord. In general, any model with the same format of inputs/outputs as @model_Torord may be simulated, which is useful when current formulations are changed within the model code, etc.
param.verbose = true; % printing numbers of beats simulated.

options = []; % parameters for ode15s - usually empty
beats = 200; % number of beats
ignoreFirst = beats - 4; % this many beats at the start of the simulations are ignored when extracting the structure of simulation outputs (i.e., beats - 1 keeps the last beat).

X0 = getStartingState('Torord_endo'); % starting state - can be also Torord_mid or Torord_epi for midmyocardial or epicardial cells respectively.
%X0 = getStartingState('Torord_epi');

%% Simulation and extraction of outputs
% time, X are cell arrays corresponding to stored beats (if 1 beat is
% simulated, this is 1-by-1 cell still), giving time vectors and state
% variable values at corresponding time points.
[time, X] = modelRunner(X0, options, param, beats, ignoreFirst);
% A structure of currents is computed from the state variables (see the
% function code for a list of properties extracted - also, hitting Tab
% following typing 'currents.' lists all the fields of the structure). Some
% state variables are also stored in a named way (time, V, Cai, Cass) so
% that the user can do most of necessary plotting simply via accessing the
% structure currents as shown below. 
currents = getCurrentsStructure(time, X, param, 0);

female_color = [23/255, 190/255, 187/255]; 
male_color = [239/255, 62/255, 54/255]; 
figure(1); set(gcf, 'color', 'w'); hold on;
if param.sexType == 1
    colors=male_color;
elseif param.sexType == 2
    colors=female_color;
end
beatDuration = param.bcl; % Duration of a single beat
numBeats = floor(currents.time(end) / beatDuration); % Number of full beats in the data
APD90s = []; % Initialize array for APD90 values

for i = 1:numBeats
    % Extract time and voltage for the current beat
    beatStart = (i-1) * beatDuration; % Start time of the current beat
    beatEnd = i * beatDuration; % End time of the current beat

    % Find indices corresponding to this beat
    indices = find(currents.time >= beatStart & currents.time < beatEnd);

    % Extract the time and voltage for this beat
    beatTime = currents.time(indices) - beatStart; % Shift time to start at 0
    beatVoltage = currents.V(indices);

    % Calculate APD90 for this beat
    APD90 = calculateAPD(beatTime, beatVoltage, 90);
    APD90s = [APD90s, APD90]; % Store the APD90
end
APD90s


%% Plotting membrane potential and calcium transient
figure(1),set(gcf,'color','w'); hold on
%plot(currents.time, currents.V,'LineWidth', 1.2, 'Color', 'b');
plot(currents.time, currents.V,'Color', colors,'Linewidth',3);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
xlim([-50 param.bcl]);
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
set(gca, 'FontSize', 28); % Increase font size for tick labels
set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
hold off;

figure(2),set(gcf,'color','w');
hold on
plot(currents.time, currents.Cai*1e6,'Color', colors,'Linewidth',3);
xlabel('Time (ms)');
ylabel('[Ca^{2+}]_i (nM)');
xlim([-50 param.bcl]);
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
set(gca, 'FontSize', 28); % Increase font size for tick labels
set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
hold off;

%% Plot INaL,IKr or sub for other current
% figure(3),set(gcf,'color','w');
% hold on
% plot(currents.time-CL, currents.INaL,'Linewidth',3);
% xlabel('Time (ms)');
% ylabel('I_{NaL} (nM)');
% xlim([-50 1000]);
% ylim([-1 0]);
% set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
% set(gca, 'FontSize', 28); % Increase font size for tick labels
% set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
% hold off;
% 
% figure(4),set(gcf,'color','w');
% hold on
% plot(currents.time/1000, currents.IKr,'Color', colors,'Linewidth',3);
% xlabel('Time (ms)');
% ylabel('IKr (uA/uF)');
% xlim([-0.05 1.5]);
% ylim([0 1.2]);
% set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
% set(gca, 'FontSize', 22); % Increase font size for tick labels
% set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
% hold off;

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
