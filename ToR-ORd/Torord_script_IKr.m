%% Setting parameters
clear 
% param is the default model parametrization here
param.cellType = 0; %0 endo 1 epi
param.sexType = 1; %male 1 female 2
param.bcl = 1000; %Change BCL here
param.model = @model_Torord; %Baseline ToR-ORd model
param.IKr_Multiplier = 1; %Initalize IKr multiplier

%% Change IKr Multipliers Here 
% For a list of multipliers
%ikrMultipliers = 0.8:0.1:1.2;

% For a Single Mulitplier
ikrMultipliers=0.1;

params(1:length(ikrMultipliers)) = param; 


for iParam = 1:length(ikrMultipliers)
    params(iParam).IKr_Multiplier = ikrMultipliers(iParam); 
end


options = [];
beats = 200; %Set number of beats
ignoreFirst = beats - 4;% Save last 4 beats

%% Simulation and output extraction

% Sequential loop to run the models
for i = 1:length(params) 
    X0 = getStartingState('Torord_endo');
    [time{i}, X{i}] = modelRunner(X0, options, params(i), beats, ignoreFirst);
    currents{i} = getCurrentsStructure(time{i}, X{i}, params(i), 0);
    currents= currents{i};

    beatDuration = param.bcl;
    numBeats = floor(currents.time(end) / beatDuration); 
    APD90s = []; 
    
        for j = 1:numBeats
        beatStart = (j-1) * beatDuration; 
        beatEnd = j * beatDuration; 
        indices = find(currents.time >= beatStart & currents.time < beatEnd);
    
        beatTime = currents.time(indices) - beatStart; 
        beatVoltage = currents.V(indices);

        APD90 = calculateAPD(beatTime, beatVoltage, 90);
        APD90s = [APD90s, APD90]; 
        end
        APD90s
end

female_color = [23/255, 190/255, 187/255]; 
male_color = [239/255, 62/255, 54/255]; 

%% Plotting membrane potential and calcium transient
figure(1); set(gcf, 'color', 'w'); hold on;
if param.sexType == 1
    colors=male_color;
elseif param.sexType == 2
    colors=female_color;
end

for i = 1:length(params)
    figure(1), set(gcf, 'color', 'w');
    hold on
    plot(currents.time-param.bcl, currents.V, 'Color', colors(i, :), 'LineWidth', 3);%, 'LineStyle', '--'); 
    ylabel('Voltage (mV)'); 
    xlabel('Time (ms)')
    xlim([-100 param.bcl]);
    set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
    set(gca, 'FontSize', 28); % Increase font size for tick labels
    set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
    hold off

    figure(2), set(gcf, 'color', 'w');
    hold on
    plot(currents.time-param.bcl, currents.Cai * 1e6, 'Color', colors(i, :), 'LineWidth', 3);%, 'LineStyle', '--'); 
    ylabel('[Ca]_i (nM)'); 
    xlabel('Time (ms)')
    %ylim([0 250])
    xlim([-100 param.bcl]);
    set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
    set(gca, 'FontSize', 28); % Increase font size for tick labels
    set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
    hold off
end

% title('Exploration of I_{Kr} multiplier');
% legend('0.8', '0.9', '1.0', '1.1', '1.2');
% xlabel('Time (ms)');
% ylabel('Membrane potential (mV)');
% xlim([0 1000]);
% 
% figure(1); set(gcf, 'color', 'w'); hold on;
% set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 12)
% colors = lines(length(params)); % Get a set of colors for plotting
% for i = 1:length(params)
%     hold on
%     subplot(2,1,1)
%     plot(currents{i}.time, currents{i}.V, 'Color', colors(i, :));
%     ylim([-90 40]);
%     xlim([-10 2000]); % Set xlim for the first subplot
%     title('Exploration of I_{Kr} Multiplier'); % Title for the first subplot
%     ylabel('Membrane potential (mV)'); % Y label for the first subplot
%     subplot(2,1,2)
%     plot(currents{i}.time, currents{i}.IKr, 'Color', colors(i, :));
%     hold off
% end


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
