
% script used to plot figure 2 in Fennell & Primeau, 2024
% user needs to have run 'driver.m' and have all 26 output est.mat files
% accessible in 'output_mat_files'

% load and parse data
load data.mat
[in] = data;
nD = length(in);

[TCobs]     = in(5,:)';     % TC (umol/kg)
[TAobs]     = in(6,:)';     % TA (umol/kg)
[pco2obs]   = in(10,:)';    % pco2 (uatm)
[co3obs]    = in(11,:)';    % co3 (umol/kg)
[phobs]     = in(9,:)';     % ph
[siobs]     = in(8,:)';     % TSi (umol/kg)
[tpobs]     = in(7,:)';     % TP (umol/kg)

% load 10 pairs to workspace
load output_mat_files/est01.mat;
load output_mat_files/est02.mat;
load output_mat_files/est03.mat;
load output_mat_files/est04.mat;
load output_mat_files/est05.mat;
load output_mat_files/est06.mat;
load output_mat_files/est07.mat;
load output_mat_files/est08.mat;
load output_mat_files/est09.mat;
load output_mat_files/est10.mat;

nD = length(est01);

for i = 1:nD
    input(1).est(i) = est01(i);
    input(2).est(i) = est02(i);
    input(3).est(i) = est03(i);
    input(4).est(i) = est04(i);
    input(5).est(i) = est05(i);
    input(6).est(i) = est06(i);
    input(7).est(i) = est07(i);
    input(8).est(i) = est08(i);
    input(9).est(i) = est09(i);
    input(10).est(i) = est10(i);
end
 
for k = 1:10
    est = input(k).est;
    for i = 1:nD
        pco2(i,k)   = est(i).tp(2).pco2; % posterior pco2
        sig(i,k)    = est(i).tp(2).upco2; % posterior uncertainty
        dpco2(i,k)  = pco2obs(i) - pco2(i,k); % meas - calc
        sig_in(i,k) = 0.01*pco2obs(i);
        perc_sig_out(i,k) = 100*sig(i,k)/pco2(i,k);
        % Z-score = ( meas - calc ) / sigma_meas
        zpco2(i,k)  = dpco2(i,k)/(sig_in(i,k)); % sigma_meas, 1% uncert
        zsig2(i,k)  = dpco2(i,k)/(sqrt(sig(i,k)^2 + sig_in(i,k)^2));
    end
end


% AXES PROPERTIES
set(groot,'defaultAxesFontName','Perpetua',...
    'defaultAxesFontSize',10,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','off',...
    'defaultAxesYMinorTick','off');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Perpetua',...
    'defaultTextInterpreter','latex');

lbl = {'C$_T$, A$_T$',...       % 1
    'C$_T$, pH',...             % 2
    'C$_T$, $p$CO$_{2}$',...    % 3
    'C$_T$, CO$_{3}$', ...      % 4
    'A$_T$, pH',...             % 5
    'A$_T$, $p$CO$_{2}$',...    % 6
    'A$_T$, CO$_{3}$',...       % 7
    'pH, CO$_{3}$',...          % 8
    'pH, $p$CO$_{2}$',...       % 9
    '$p$CO$_{2}$, CO$_{3}$'};   % 10

ddgreen = [0, 102/255, 0];
clr = ddgreen;
%% 

tiledlayout(1,4)

nexttile % delta pco2
b1 = boxchart(dpco2);
b1.JitterOutliers = 'on';
b1.MarkerStyle = '.';
b1.MarkerSize = 7;
b1.BoxWidth = 0.7;
b1.BoxFaceColor = clr;
b1.BoxEdgeColor = clr;
b1.MarkerColor = clr;
b1.WhiskerLineColor = clr;
b1.LineWidth = 1.4;

ylim([-100 100]) % dpco2

ax = gca;
ax.FontSize = 9;
xlabel('QUODcarb Input') 
ylabel('Delta $p$CO$_2$ = Meas - Calc ($\mu$atm)') % pco2
ax.XTickMode = 'auto';
xticklabels(lbl); 
grid on

nexttile % sigma_posterior uatm
% b2 = boxchart(sig); % quodcarb
b2 = boxchart(perc_sig_out); % quodcarb
b2.JitterOutliers = 'on';
b2.MarkerStyle = '.';
b2.MarkerSize = 7;
b2.BoxWidth = 0.7;
b2.BoxFaceColor = clr;
b2.BoxEdgeColor = clr;
b2.MarkerColor = clr;
b2.WhiskerLineColor = clr;
b2.LineWidth = 1.4;

% ylim([0 65])

ax = gca;
ax.FontSize = 9;
xlabel('QUODcarb Input') % quodcarb
ylabel('$p$CO$_2$ $u_{posterior}$ ($\mu$atm)') % upco2
ax.XTickMode = 'auto';
xticklabels(lbl);
grid on

nexttile 
b3 = boxchart(zpco2);
b3.JitterOutliers = 'on';
b3.MarkerStyle = '.';
b3.MarkerSize = 7;
b3.BoxWidth = 0.7;
b3.BoxFaceColor = clr;
b3.BoxEdgeColor = clr;
b3.MarkerColor = clr;
b3.WhiskerLineColor = clr;
b3.LineWidth = 1.4;

ylim([-20 20]) 

ax = gca;
ax.FontSize = 9;
xlabel('QUODcarb Input') % quodcarb
ylabel('Z-scores: ( Meas - Calc ) / $u_{meas}$')
ax.XTickMode = 'auto';
xticklabels(lbl);
grid on

nexttile 
b4 = boxchart(zsig2);
b4.JitterOutliers = 'on';
b4.MarkerStyle = '.';
b4.MarkerSize = 7;
b4.BoxWidth = 0.7;
b4.BoxFaceColor = clr;
b4.BoxEdgeColor = clr;
b4.MarkerColor = clr;
b4.WhiskerLineColor = clr;
b4.LineWidth = 1.4;

% ylim([-20 20]) 

ax = gca;
ax.FontSize = 9;
xlabel('QUODcarb Input') % quodcarb
ylabel('Z-scores: ( Meas - Calc ) / ($u_{meas}^2 + u_{calc}^2$)$^{1/2}$')
ax.XTickMode = 'auto';
xticklabels(lbl);
grid on

h = gcf; 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);

print(h,'f3_pco2_pairs.pdf','-dpdf','-r0');
















