
gendertype = 2;%1 and 2 in TORORD; 0 and 1 in Grandi
celltype = 0;
CL= 3000;%2000 ICaL; 3000 IKr
dshift= 6;
gCaL_inc= 1.45;
gKr_dec = 0.10;
gCaL_grandi = 0.05;
dshift_grandi = 3;
gKr_ord=0.80;

%subfolder = 'LQT8_CL2000_d6_gcal45';
subfolder = 'LQT2_CL2000';
%subfolder = 'Grandi_INaL_LQT8';
%subfolder = 'ORd_LQT2';
APD90s = [];
APD90_EAD = [];
APD90_no_EAD = [];
female_color = [23/255, 190/255, 187/255]; 
male_color = [239/255, 62/255, 54/255];
dark_female_color = [15/255, 120/255, 118/255]; % Dark blue for EAD in females
dark_male_color = [139/255, 0, 0];  % Dark red for EAD in males

%% Plotting colors
if gendertype == 1
    colors = male_color;
    dark_color = dark_male_color;
elseif gendertype == 2
    colors = female_color;
    dark_color = dark_female_color;
end



for ii = 1:500 %Can change to limit how much of population is plotted
    %% ICaL TORORD
    % Vname = sprintf('%s/Vm_gen%d_cell%d_CL%d_V%d_G%.2f_%d.mat',subfolder,gendertype, celltype, CL, dshift, gCaL_inc, ii); load(Vname);
    % tname = sprintf('%s/time_gen%d_cell%d_CL%d_V%d_G%.2f_%d.mat',subfolder,gendertype,celltype, CL, dshift, gCaL_inc, ii); load(tname);
    % currents_name = sprintf('%s/currents_matrix_gen%d_cell%d_CL%d_V%d_G%.2f_%d.mat',subfolder,gendertype,celltype, CL, dshift, gCaL_inc, ii); load(currents_name);
    % Cai = currents_matrix.Cai; 

    %% IKr TORORD
    Vname = sprintf('%s/Vm_gen%d_cell%d_CL%d_G%.2f_%d.mat',subfolder, gendertype, celltype, CL, gKr_dec, ii); load(Vname);
    tname = sprintf('%s/time_gen%d_cell%d_CL%d_G%.2f_%d.mat',subfolder, gendertype,celltype, CL, gKr_dec, ii); load(tname);
    currents_name = sprintf('%s/currents_matrix_gen%d_cell%d_CL%d_G%.2f_%d.mat',subfolder,gendertype,celltype, CL, gKr_dec, ii); load(currents_name);
    Cai = currents_matrix.Cai; 

    %% ICaL Grandi
    % Vname = sprintf('%s/Vm_gen%d_cell%d_CL%d_V%d_G%.2f_%d.mat',subfolder,gendertype, celltype, CL, dshift_grandi, gCaL_grandi, ii); load(Vname);
    % tname = sprintf('%s/time_gen%d_cell%d_CL%d_V%d_G%.2f_%d.mat',subfolder,gendertype,celltype, CL, dshift_grandi, gCaL_grandi, ii); load(tname);
    % currents_name = sprintf('%s/Cai_gen%d_cell%d_CL%d_V%d_G%.2f_%d.mat',subfolder,gendertype,celltype, CL, dshift_grandi, gCaL_grandi, ii); load(currents_name);

    %% IKr ORd
    % Vname = sprintf('%s/Vm_gen%d_cell%d_CL%d_G%.2f_%d.mat',subfolder,gendertype, celltype, CL, gKr_ord, ii); load(Vname);
    % tname = sprintf('%s/time_gen%d_cell%d_CL%d_G%.2f_%d.mat',subfolder,gendertype,celltype, CL, gKr_ord, ii); load(tname);
    % currents_name = sprintf('%s/Cai_gen%d_cell%d_CL%d_G%.2f_%d.mat',subfolder,gendertype,celltype, CL, gKr_ord, ii); load(currents_name);

    second_last_beat_start = time(end) - 2 * CL;
    second_last_beat_idx = find(time >= second_last_beat_start, 1, 'first');
    time_last_two_beats = time(second_last_beat_idx:end);
    Vm_last_two_beats = Vm(second_last_beat_idx:end);
    Cai_last_two_beats = Cai(second_last_beat_idx:end);

    plot_time = time_last_two_beats - (time(end) - CL); 

    % Extract indices for the plot range (-100 ms to CL)
    plot_idx = find(plot_time >= -100 & plot_time <= CL);
    time_plot = plot_time(plot_idx);
    Vm_plot = Vm_last_two_beats(plot_idx);
    Cai_plot = Cai_last_two_beats(plot_idx);

    % Detect EAD presence
    EAD_present = function_EAD_occurrence(time_plot, Vm_plot, time_plot(1), time_plot(end) - time_plot(1));

    %% Plot Vm
    %figure(1), set(gcf, 'color', 'w'); hold on;
    if EAD_present
        % Plot with darker color for EAD
        plot(time_last_two_beats - time(end) + CL, Vm_last_two_beats, 'LineWidth', 2.5, 'Color', [dark_color 0.6]); 
        APD90_value = calculateAPD(time_plot, Vm_plot, 90);
        APD90_EAD = [APD90_EAD, APD90_value];
        
    else
        % Plot with normal color for no EAD
        plot(time_last_two_beats - time(end) + CL, Vm_last_two_beats, 'LineWidth', 2.5, 'Color', [colors 0.3]);
        APD90_value = calculateAPD(time_plot, Vm_plot, 90);
        APD90_no_EAD = [APD90_no_EAD, APD90_value];
    end
    

end

% Calculates statistics for APD90 without EADs
mean_no_EAD = mean(APD90_no_EAD);
std_no_EAD = std(APD90_no_EAD);
median_no_EAD = median(APD90_no_EAD);
iqr_no_EAD = iqr(APD90_no_EAD); % Interquartile range
q1_no_EAD = quantile(APD90_no_EAD, 0.25); % First quartile
q3_no_EAD = quantile(APD90_no_EAD, 0.75); % Third quartile
range_no_EAD = range(APD90_no_EAD);

% Calculates statistics for APD90 with EADs
mean_EAD = mean(APD90_EAD);
std_EAD = std(APD90_EAD);
median_EAD = median(APD90_EAD);
iqr_EAD = iqr(APD90_EAD); % Interquartile range
q1_EAD = quantile(APD90_EAD, 0.25); % First quartile
q3_EAD = quantile(APD90_EAD, 0.75); % Third quartile
range_EAD = range(APD90_EAD);

% Display results
disp('Statistics for APD90 without EADs:');
disp(['Mean: ', num2str(mean_no_EAD)]);
disp(['Standard Deviation: ', num2str(std_no_EAD)]);
disp(['Median: ', num2str(median_no_EAD)]);
disp(['Q1: ', num2str(q1_no_EAD)]);
disp(['Q3: ', num2str(q3_no_EAD)]);
disp(['Range: ', num2str(range_no_EAD)]);
disp(['IQR: ', num2str(iqr_no_EAD)]);

disp('Statistics for APD90 with EADs:');
disp(['Mean: ', num2str(mean_EAD)]);
disp(['Standard Deviation: ', num2str(std_EAD)]);
disp(['Median: ', num2str(median_EAD)]);
disp(['Q1: ', num2str(q1_EAD)]);
disp(['Q3: ', num2str(q3_EAD)]);
disp(['Range: ', num2str(range_EAD)]);
disp(['IQR: ', num2str(iqr_EAD)]);
% [h_var, p_var] = vartest2(APD90_no_EAD, APD90_EAD);
% disp(['Variance Test P-Value: ', num2str(p_var)]);
% [h_ttest, p_ttest] = ttest2(APD90_no_EAD, APD90_EAD);
% disp(['Two-Sample T-Test P-Value: ', num2str(p_ttest)]);
% [p_mwu, h_mwu] = ranksum(APD90_no_EAD, APD90_EAD);
% disp(['Mann-Whitney U Test P-Value: ', num2str(p_mwu)])

%% Plotting APD90 Histogram
figure(3), set(gcf, 'color', 'w'); hold on;
binEdges = 0:75:2000; %Change binEdges based on population bounds
%binEdges = 0:25:1300;
histogram(APD90_EAD, 'BinEdges', binEdges, 'FaceColor', dark_color, 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.7);
histogram(APD90_no_EAD, 'BinEdges', binEdges, 'FaceColor', colors, 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.7);
%histogram(APD90_no_EAD, 'FaceColor', colors, 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.7);
%histogram(APD90_EAD, 'FaceColor', dark_color, 'EdgeColor', 'k', 'LineWidth', 2, 'FaceAlpha', 0.7);
xlabel('APD_{90} (ms)');
%xlim([250 900]);
%ylim([0 150]);
%xlim([500 1300]);
xlim([150 2000]);
ylim([0 60]);
%xlim([600 1800])
ylabel('Frequency');
set(gca, 'TickDir', 'out', 'TickLength', [0.02, 0.02]); % Makes ticks point outwards
set(gca, 'FontSize', 28); % Increase font size for tick labels
set(gca, 'LineWidth', 2); % Set the thickness of the axis lines
hold off;
