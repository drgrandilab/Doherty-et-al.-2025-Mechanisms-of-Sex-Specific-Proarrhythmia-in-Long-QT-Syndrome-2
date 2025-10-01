clear all;
close all;
clc;

%% Load EAD Output Matrices for Both Genders
% Female data
%load EAD_outputs_matrix_LQT8F;
load EAD_outputs_matrix_LQT2F;
%load EAD_outputs_matrix_LQT8F_grandi
%load EAD_outputs_matrix_LQT2F_ord 
all_outputs_ead_female = all_outputs;
% Male data
%load EAD_outputs_matrix_LQT8M;
load EAD_outputs_matrix_LQT2M ;
%load EAD_outputs_matrix_LQT8M_grandi
%load EAD_outputs_matrix_LQT2M_ord 
all_outputs_ead_male = all_outputs;

%% Analysis EAD/no EAD
aos_female = size(all_outputs_ead_female);
N_beat_female = aos_female(2);
N_cell_female = aos_female(1);
all_outputs_ead_sum_female = sum(all_outputs_ead_female, 2);
ead_presence_female = (all_outputs_ead_sum_female > 1/2);
fraction_ead_female = 100 * sum(ead_presence_female) / N_cell_female;
num_ead_female=sum(ead_presence_female);

aos_male = size(all_outputs_ead_male);
N_beat_male = aos_male(2);
N_cell_male = aos_male(1);
all_outputs_ead_sum_male = sum(all_outputs_ead_male, 2);
ead_presence_male = (all_outputs_ead_sum_male > 1/2);
fraction_ead_male = 100 * sum(ead_presence_male) / N_cell_male;
num_ead_male=sum(ead_presence_male);
%% Load parameters
%load parameter_matrix_500_grandi.mat
%load parameter_matrix_500_ord.mat
load parameter_matrix_500.mat

[N_trials, N_pars] = size(all_parameters);

%% Logistic Regression

%TORORD
parameter_names_b0 = {
     'G_{Na}' 'G_{CaL}' 'G_{Ito}' 'G_{INaL}' 'G_{Kr}'...
     'G_{Ks}' 'G_{K1}' 'G_{Kb}' 'G_{NaCa}' 'G_{NaK}'...
     'G_{NaB}' 'G_{CaB}' 'G_{pCa}' 'G_{CaCl}' 'G_{ClB}'...
     'J_{rel}' 'J_{up}'};
%Grandi
% parameter_names_b0 = {
%     'G_{Na}' 'G_{CaL}' 'G_{Itof}' 'G_{Itos}' 'G_{Kr}'...
%     'G_{Ks}' 'G_{Kp}' 'G_{K1}' 'G_{ClCa}' 'G_{Clb}'...
%     'G_{NaB}' 'G_{CaB}' 'G_{NaK}' 'G_{NaCa}' 'G_{pmCa}'...
%     'J_{up}' 'J_{rel}' 'J_{leak}' 'G_{NaL}'};
%ORD
% parameter_names_b0 = {
%     'G_{Na}' 'G_{NaL}' 'G_{to}' 'P_{Ca}' 'G_{Kr}'...
%     'G_{Ks}' 'G_{K1}' 'G_{ncx}' 'P_{naK}' 'G_{Kb}'...
%     'P_{Nab}' 'P_{CaB}' 'G_{pCa}' 'f_{SERCA}' 'f_{RyR}'...
%     'f_{Ileak}'};
X_LOG = log(all_parameters);
for ii = 1:N_pars % z-score
    X_LOG(:,ii) = (X_LOG(:,ii) - mean(X_LOG(:,ii))) / std(X_LOG(:,ii));
end

%% Logistic Regression for Female
Y_LOGISTIC_female = 1 - (ead_presence_female - 1);
[B_LOGISTIC_female, dev_female, stats_female] = mnrfit(X_LOG, Y_LOGISTIC_female);

%% Logistic Regression for Male
Y_LOGISTIC_male = 1 - (ead_presence_male - 1);
[B_LOGISTIC_male, dev_male, stats_male] = mnrfit(X_LOG, Y_LOGISTIC_male);

%% Calculate EAD Probabilities and PEAD
% Female
B0_female = B_LOGISTIC_female(1);
P_ead_female = 1 ./ (1 + exp(-(B0_female + X_LOG * B_LOGISTIC_female(2:end))));
P_ead_mean_female = mean(P_ead_female);
P_ead_std_female = std(P_ead_female);

% Male
B0_male = B_LOGISTIC_male(1);
P_ead_male = 1 ./ (1 + exp(-(B0_male + X_LOG * B_LOGISTIC_male(2:end))));
P_ead_mean_male = mean(P_ead_male);
P_ead_std_male = std(P_ead_male);

%% Tjur R2 Calculation
array_1_female = P_ead_female(ead_presence_female > 0.5);
array_1_mean_female = mean(array_1_female);
array_0_female = P_ead_female(ead_presence_female < 0.5);
array_0_mean_female = mean(array_0_female);
R2_logistic_female = array_1_mean_female - array_0_mean_female;

array_1_male = P_ead_male(ead_presence_male > 0.5);
array_1_mean_male = mean(array_1_male);
array_0_male = P_ead_male(ead_presence_male < 0.5);
array_0_mean_male = mean(array_0_male);
R2_logistic_male = array_1_mean_male - array_0_mean_male;

%% Remove b0 value
B_LOGISTIC_female = B_LOGISTIC_female(2:end);
B_LOGISTIC_male = B_LOGISTIC_male(2:end);

%% Split coefficients into positive and negative
positive_coefficients_female = B_LOGISTIC_female(B_LOGISTIC_female >= 0);
negative_coefficients_female = B_LOGISTIC_female(B_LOGISTIC_female < 0);

positive_names_female = parameter_names_b0(B_LOGISTIC_female >= 0);
negative_names_female = parameter_names_b0(B_LOGISTIC_female < 0);

% Sort negative coefficients in descending order
[negative_coefficients_female, neg_idx_female] = sort(negative_coefficients_female, 'descend');
negative_names_female = negative_names_female(neg_idx_female);

% Sort positive coefficients in ascending order
[positive_coefficients_female, pos_idx_female] = sort(positive_coefficients_female, 'ascend');
positive_names_female = positive_names_female(pos_idx_female);

%% Plot Positive and Negative Coefficients for Female
figure(1);
set(gcf, 'color', 'w');
hold on;

% Plot negative coefficients
for i = 1:length(negative_coefficients_female)
    barh(i, negative_coefficients_female(i), 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');
end

% Plot positive coefficients
%offset = length(negative_coefficients_female);
for i = 1:length(positive_coefficients_female)
    barh(i, positive_coefficients_female(i), 'FaceColor', 'k', 'EdgeColor', 'none');
end

% Update axes settings
set(gca, 'YGrid', 'off', 'XGrid', 'off');
set(gca, 'YColor', 'none'); % Remove Y-axis line
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 10);
title('Female Logistic Regression Coefficients');
set(gca, 'YTick', []); % Remove Y-axis labels
set(gca, 'FontSize', 28); % Increase font size for tick labels
set(gca, 'LineWidth', 2);  % Set the thickness of the axis lines
xlabel('Regression Coefficients');
hold off;

%% Repeat for Male
positive_coefficients_male = B_LOGISTIC_male(B_LOGISTIC_male >= 0);
negative_coefficients_male = B_LOGISTIC_male(B_LOGISTIC_male < 0);

positive_names_male = parameter_names_b0(B_LOGISTIC_male >= 0);
negative_names_male = parameter_names_b0(B_LOGISTIC_male < 0);

% Sort negative coefficients in descending order
[negative_coefficients_male, neg_idx_male] = sort(negative_coefficients_male, 'descend');
negative_names_male = negative_names_male(neg_idx_male);

% Sort positive coefficients in ascending order
[positive_coefficients_male, pos_idx_male] = sort(positive_coefficients_male, 'ascend');
positive_names_male = positive_names_male(pos_idx_male);

%% Plot Positive and Negative Coefficients for Male
figure(2);
set(gcf, 'color', 'w');
hold on;

% Plot negative coefficients
for i = 1:length(negative_coefficients_male)
    barh(i, negative_coefficients_male(i), 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');
end

% Plot positive coefficients
%offset = length(negative_coefficients_male);
for i = 1:length(positive_coefficients_male)
    barh(i, positive_coefficients_male(i), 'FaceColor', 'k', 'EdgeColor', 'none');
end

% Update axes settings
set(gca, 'YGrid', 'off', 'XGrid', 'off');
set(gca, 'YColor', 'none'); % Remove Y-axis line
set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 10);
title('Male Logistic Regression Coefficients');
set(gca, 'YTick', []); % Remove Y-axis labels
set(gca, 'FontSize', 28); % Increase font size for tick labels
set(gca, 'LineWidth', 2);  % Set the thickness of the axis lines
xlabel('Regression Coefficients');
hold off;

%% Print labels and coefficients
fprintf('Female Negative Coefficients:\n');
for i = 1:length(negative_coefficients_female)
    fprintf('%s: %.4f\n', negative_names_female{i}, negative_coefficients_female(i));
end

fprintf('Female Positive Coefficients:\n');
for i = 1:length(positive_coefficients_female)
    fprintf('%s: %.4f\n', positive_names_female{i}, positive_coefficients_female(i));
end

fprintf('Male Negative Coefficients:\n');
for i = 1:length(negative_coefficients_male)
    fprintf('%s: %.4f\n', negative_names_male{i}, negative_coefficients_male(i));
end

fprintf('Male Positive Coefficients:\n');
for i = 1:length(positive_coefficients_male)
    fprintf('%s: %.4f\n', positive_names_male{i}, positive_coefficients_male(i));
end

%% Display Results
disp('PEAD in the baseline model (Female):');
P_ead_mean_B0_female = 1 / (1 + exp(-B0_female))
disp('PEAD in the baseline model (Male):');
P_ead_mean_B0_male = 1 / (1 + exp(-B0_male))

disp('Fraction of simulations with at least 1 EAD (Female):');
fraction_ead_female
disp('Fraction of simulations with at least 1 EAD (Male):');
fraction_ead_male

disp('R2 (Female):');
R2_logistic_female
disp('R2 (Male):');
R2_logistic_male

%% Plot EAD Probabilities
% figure; set(gcf, 'color', 'w');
% hold on;
% plot(1:length(P_ead_female), P_ead_female, '*', 'Color', female_color, 'DisplayName', 'Female');
% plot(1:length(P_ead_male), P_ead_male, '*', 'Color', male_color, 'DisplayName', 'Male');
% set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 10);
% title('Probability of EAD Development');
% xlabel('Trial');
% ylabel('Probability');
% legend('Female', 'Male');
% hold off;
% 
% figure; set(gcf, 'color', 'w');
% hold on;
% plot(1:length(P_ead_female), P_ead_female, '*', 'Color', female_color, 'DisplayName', 'Female');
% plot(1:length(P_ead_male), P_ead_male, '*', 'Color', male_color, 'DisplayName', 'Male');
% set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 10);
% title('Probability of EAD Development');
% xlabel('Trial');
% ylabel('Probability');
% legend('Female', 'Male');
% hold off;

