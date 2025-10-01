clear 
param.cellType = 0; %0 endo 1 epi
param.sexType = 1; %male 1 female 2
param.bcl = 1000; %Change BCL here
param.model = @model_Torord; %Change base ToR-ORd model here
param.ICaL_Multiplier = 1; %initialize Multiplier

%%Implement ICaL Changes Here
vCaLChanges = 6;
icalMultipliers= 1.45;

%Loop here also allows you to test a wide range of changes in same
%simulation
params = repmat(param, [length(icalMultipliers), 1]);
for iParam = 1:length(icalMultipliers)
    params(iParam).ICaL_Multiplier = icalMultipliers(iParam);
    params(iParam).VCaL_Change = vCaLChanges;
end
options = [];
beats = 200;%Change number of beats here
ignoreFirst = beats - 4; %Save last 4 (or change) beats here
APD90s = []; % Initialize array for APD90 values


for i = 1:length(params) 
    X0 = getStartingState('Torord_endo');
    [time{i}, X{i}] = modelRunner(X0, options, params(i), beats, ignoreFirst);
    currents{i} = getCurrentsStructure(time{i}, X{i}, params(i), 0);
    currents= currents{i};

    beatDuration = param.bcl; % Duration of a single beat
    numBeats = floor(currents.time(end) / beatDuration); % Number of full beats in the data
    
        for j = 1:numBeats
        beatStart = (j-1) * beatDuration; % Start time of the current beat
        beatEnd = j * beatDuration; % End time of the current beat
    
        % Find indices corresponding to this beat
        indices = find(currents.time >= beatStart & currents.time < beatEnd);
    
        % Extract the time and voltage for this beat
        beatTime = currents.time(indices) - beatStart; % Shift time to start at 0
        beatVoltage = currents.V(indices);

        % Calculate APD90 for this beat
        APD90 = calculateAPD(beatTime, beatVoltage, 90);
        APD90s = [APD90s, APD90]; % Store the APD90
  
    end
end

female_color = [23/255, 190/255, 187/255]; 
male_color = [239/255, 62/255, 54/255]; 


%% Define Colors
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
    xlim([-50 param.bcl]);
    set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
    set(gca, 'FontSize', 28); % Increase font size for tick labels
    set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
    hold off

    figure(2), set(gcf, 'color', 'w');
    hold on
    plot(currents.time-param.bcl, currents.Cai * 1e6, 'Color', colors(i, :), 'LineWidth', 3);%, 'LineStyle', '--'); 
    ylabel('[Ca]_i (nM)'); 
    xlabel('Time (ms)')
    %ylim([0 1600])
    xlim([-50 param.bcl]);
    set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
    set(gca, 'FontSize', 28); % Increase font size for tick labels
    set(gca, 'LineWidth', 1.5);  % Set the thickness of the axis lines
    hold off
end



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
