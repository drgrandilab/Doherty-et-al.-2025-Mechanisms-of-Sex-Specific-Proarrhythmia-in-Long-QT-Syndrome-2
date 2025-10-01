% This file performs the regression analysis and plots the results.
% It loads the matrix with the perturbations of model parameters

%Change the EAD output matrix for LQT type, sex, and model
%Then change which parameter names get used
clear all
close all
clc

%% Analysis EAD/no EAD
% all_outputs_ead: 500 x 5 (1-5 EAD index)

load EAD_outputs_matrix_LQT8F
%load EAD_outputs_matrix_LQT8M
%load EAD_outputs_matrix_LQT2F 
%load EAD_outputs_matrix_LQT2M 
%load EAD_outputs_matrix_LQT8F_grandi
%load EAD_outputs_matrix_LQT8M_grandi
%load EAD_outputs_matrix_LQT2M_ord 
%load EAD_outputs_matrix_LQT2F_ord 

all_outputs_ead = all_outputs;
aos = size(all_outputs_ead); % 500 x 5
N_beat = aos(2);
N_cell = aos(1);

all_outputs_ead_sum = sum(all_outputs_ead'); % 1 x 1000
ead_presence = (all_outputs_ead_sum>1/2); % 1 for EAD occurrence, 0 for no EAD
disp('Fraction of simulations with at least 1 EAD (%):');
fraction_ead = 100*sum(ead_presence)/N_cell

%% Load parameters
load parameter_matrix_500.mat
%load parameter_matrix_500_grandi.mat
%load parameter_matrix_500_ord.mat
[N_trials N_pars] = size(all_parameters);


%% Logistic regression - EAD/no EAD
parameter_names_b0 = {'b0',...
    'GNa' 'GCaL' 'Gto' 'GNaL' 'GKr'...
    'GKs' 'GK1' 'GKb' 'GNaCa' 'GNaK'...
    'GNaB' 'GCaB' 'GpCa' 'GCaCl' 'GClb'...
    'Jrel' 'Jup'};
%Grandi
% parameter_names_b0 = {'b0',...
%     'GNa' 'GCaL' 'Gtof' 'Gtos' 'GKr'...
%     'GKs' 'GKp' 'GK1' 'GClCa' 'GClb'...
%     'GNaB' 'GCaB' 'GNaK' 'GNaCa' 'GpmCa'...
%     'Jup' 'Jrel' 'Jleak' 'GNaL'};
%ORD
% parameter_names_b0 = {'b0',...
%     'GNa' 'GNaL' 'Gto' 'GCaL' 'GKr'...
%     'GKs' 'GK1' 'Gncx' 'PnaK' 'GKb'...
%     'PNab' 'PCaB' 'GpCa' 'JUp' 'JRel'...
%     'Jleak'};

allpars_LOGISTIC = all_parameters;
X_LOG = log(allpars_LOGISTIC) ;
for ii=1:N_pars % z-score
    X_LOGISTIC(:,ii)=(X_LOG(:,ii)-mean(X_LOG(:,ii)))/std(X_LOG(:,ii));
end

Y_LOGISTIC = 1-(ead_presence-1); % positive integer!
% ead_presence: 0 with no EAD, 1 with EADs 
% Y_LOGISTIC: 1 with EADs, 2 with no EADs % positive integer!
Y_LOGISTIC = Y_LOGISTIC';
[B_LOGISTIC,dev,stats] = mnrfit(X_LOGISTIC,Y_LOGISTIC);

%% Plot regression coeffcients
color = [0 0 0];

figure; set(gcf,'color','w') % with b0
bar(B_LOGISTIC,'FaceColor',color)
set(gca,'box','off','tickdir','out','fontsize',10)
title('Probability EAD development')
set(gca,'XTick',1:N_pars+1)
set(gca,'XTickLabel',parameter_names_b0)
set(gca,'XLim',[0 N_pars+1+1])
rotateXLabels( gca(), 90)

%% EAD Probability
% PEAD in the baseline model
B0 = B_LOGISTIC(1);
pval_LOGISTIC = stats.p;
disp('PEAD in the baseline model:');
P_ead_mean_B0 = 1/(1+exp(-(B0)))

% PEAD in each model of the population
P_ead_array = zeros(1,N_trials);
for iii=1:N_trials
    P_ead_array(iii) = 1/(1+exp(-(B0+sum(B_LOGISTIC(2:end).*X_LOGISTIC(iii,:)'))));
end

figure,set(gcf,'color','w')
plot((1:N_trials),P_ead_array,'*','Color',color)
set(gca,'box','off','tickdir','out','fontsize',10)
title('Probability EAD development')
xlabel('Trial')
ylabel('Probability (-)')
P_ead_mean = mean(P_ead_array);
P_ead_std = std(P_ead_array);

array_1 = P_ead_array(ead_presence>0.5); % EAD
array_1_mean = mean(array_1);
array_1_std = std(array_1);
array_0 = P_ead_array(ead_presence<0.5); % no EAD
array_0_mean = mean(array_0);
array_0_std = std(array_0);

% Tjur (2009)
R2logistic = array_1_mean - array_0_mean

%% Plot PEAD as function of modulation in model parameters (Fig. 4B)
plot_effects = 1;

if plot_effects == 1
    
    rangeM = (0.5:0.001:1.55); % modulation range

    % GKb - index in B 9
    indexGKb = 9;
    muGKb = mean(X_LOG(:, indexGKb - 1));
    sigmaGKb = std(X_LOG(:, indexGKb - 1));
    rangeGKb = (log(rangeM) - muGKb) / sigmaGKb;
    PGKb = 1 ./ (1 + exp(-(B0 + B_LOGISTIC(indexGKb) * rangeGKb)));
    
    % Jup - index in B 18
    indexJup = 18;
    muJup = mean(X_LOG(:, indexJup - 1));
    sigmaJup = std(X_LOG(:, indexJup - 1));
    rangeJup = (log(rangeM) - muJup) / sigmaJup;
    PJup = 1 ./ (1 + exp(-(B0 + B_LOGISTIC(indexJup) * rangeJup)));
    
    % GKr - index in B 6
    indexGKr = 6;
    muGKr = mean(X_LOG(:, indexGKr - 1));
    sigmaGKr = std(X_LOG(:, indexGKr - 1));
    rangeGKr = (log(rangeM) - muGKr) / sigmaGKr;
    PGKr = 1 ./ (1 + exp(-(B0 + B_LOGISTIC(indexGKr) * rangeGKr)));
    
    % GCaL - index in B 3
    indexGCaL = 3;
    muGCaL = mean(X_LOG(:, indexGCaL - 1));
    sigmaGCaL = std(X_LOG(:, indexGCaL - 1));
    rangeGCaL = (log(rangeM) - muGCaL) / sigmaGCaL;
    PGCaL = 1 ./ (1 + exp(-(B0 + B_LOGISTIC(indexGCaL) * rangeGCaL)));
    
    % GNCX - index in B 10
    indexGNCX = 10;
    muGNCX = mean(X_LOG(:, indexGNCX - 1));
    sigmaGNCX = std(X_LOG(:, indexGNCX - 1));
    rangeGNCX = (log(rangeM) - muGNCX) / sigmaGNCX;
    PGNCX = 1 ./ (1 + exp(-(B0 + B_LOGISTIC(indexGNCX) * rangeGNCX)));
    
    % GNaL - index in B 5
    indexGNaL = 5;
    muGNaL = mean(X_LOG(:, indexGNaL - 1));
    sigmaGNaL = std(X_LOG(:, indexGNaL - 1));
    rangeGNaL = (log(rangeM) - muGNaL) / sigmaGNaL;
    PGNaL = 1 ./ (1 + exp(-(B0 + B_LOGISTIC(indexGNaL) * rangeGNaL)));
    
    % Figure
    figure, set(gcf, 'color', 'w'), hold on,
    plot(rangeM, PGKb, rangeM, PJup, rangeM, PGKr, rangeM, PGCaL, rangeM, PGNCX, 'LineWidth', 4);
    % Customize plot
    ylabel('PEAD (-)'), xlabel('Scale Factor (-)');
    title('Probability EAD Development');
    legend('GKb','JUp', 'GKr','GCaL','GNCX','GNaL');
    set(gca, 'box', 'off', 'tickdir', 'out', 'fontsize', 12);
    xlim([rangeM(1), rangeM(end)]);
    set(gca, 'FontSize', 28); 
    set(gca, 'LineWidth', 2);  
end

%% Plot X, Y, B matrices
plot_matrices = 0;

if plot_matrices == 1
    % X - z-score
    figure; set(gcf,'color','w')
    imagesc(X_LOGISTIC); %colormap jet;
    set(gca,'box','off','tickdir','out','fontsize',10)
    set(gca,'YDir','normal')
    title('Parameters (X) - z-score');
    xlabel('Parameters');
    ylabel('Trials');
    %set(gca,'YTick',(1:N_trials))
    set(gca,'XTick',(1:N_pars))
    set(gca,'XTickLabel',parameter_names)
    rotateXLabels( gca(), 90)
    colorbar

    % Y (EAD)
    figure; set(gcf,'color','w')
    imagesc(ead_presence'); %colormap jet;
    set(gca,'box','off','tickdir','out','fontsize',10)
    set(gca,'YDir','normal')
    title('Outputs (Y)');
    xlabel('Output (yes/no)');
    ylabel('Trials');
    %set(gca,'YTick',(1:N_trials))
    set(gca,'XTick',1)
    set(gca,'XTickLabel','EAD')
    %rotateXLabels( gca(), 90)
    colorbar

    % B (EAD)
    figure; set(gcf,'color','w')
    imagesc(B_LOGISTIC); %colormap jet;
    set(gca,'box','off','tickdir','out','fontsize',10)
    set(gca,'YDir','normal')
    title('Regression coefficients (B)');
    xlabel('Output (yes/no)');
    ylabel('Coefficients');
    set(gca,'YTick',(1:N_pars+1))
    set(gca,'YTickLabel',parameter_names_b0)
    set(gca,'XTick',1)
    set(gca,'XTickLabel','EAD')
    %rotateXLabels( gca(), 90)
    colorbar
end
% Heatmap 1: Parameter Matrix (X_LOGISTIC)
figure; set(gcf, 'color', 'w')
imagesc(X_LOGISTIC);
colormap('jet'); colorbar;
set(gca, 'YDir', 'normal')
title('Parameters (X\_LOGISTIC, z-score)');
xlabel('Parameters'); ylabel('Trials');
set(gca, 'XTick', 1:N_pars, 'XTickLabel', parameter_names_b0(2:end));
rotateXLabels(gca(), 90);

% Heatmap 2: EAD Output (Variants on Y-axis, Binary 0/1)
figure; set(gcf, 'color', 'w');
imagesc(ead_presence');
colormap([0.4 0.76 0.9; 0.95 0.4 0.4]); 
colorbar('Ticks',[0.25,0.75],'TickLabels',{'No EAD','EAD'});
set(gca, 'YDir', 'normal');
title('EAD Occurrence (ead\_presence)');
ylabel('Variants'); xlabel('EAD Occurrence');
set(gca, 'XTick', 1, 'XTickLabel', {'EAD Presence'});

% Heatmap 3: Regression Coefficients (Parameters on Y-axis)
figure; set(gcf, 'color', 'w')
imagesc(B_LOGISTIC);
colormap('parula'); colorbar;
set(gca, 'YDir', 'normal')
title('Regression Coefficients (B\_LOGISTIC)');
ylabel('Parameters'); xlabel('Coefficient Value');
set(gca, 'YTick', 1:(N_pars+1), 'YTickLabel', parameter_names_b0);
