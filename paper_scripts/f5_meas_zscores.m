
% script used to plot figure 5 in Fennell & Primeau, 2024
% user needs to have run 'driver.m' and have output est26.mat file
% accessible in 'output_mat_files'

% load and parse data
load data.mat
[in] = data;
nD = length(in);

[TCobs]     = in(5,:)';     % TC
[TAobs]     = in(6,:)';     % TA
[pco2obs]   = in(10,:)';    % pco2
[co3obs]    = in(11,:)';    % co3
[phobs]     = in(9,:)';     % ph

load output_mat_files/est26.mat;
% calculate Z-scores = (meas - calc)/sigma_meas
for i = 1:nD
    TC(i)       = est26(i).TC;  
    zscore(i,1) = (TCobs(i) - TC(i))/2.00; 

    TA(i)       = est26(i).TA; 
    zscore(i,2) = (TAobs(i) - TA(i))/2.00;

    ph(i)       = est26(i).tp(1).ph; 
    zscore(i,3) = (phobs(i) - ph(i))/0.010;

    pco2(i)     = est26(i).tp(2).pco2; 
    zscore(i,4) = (pco2obs(i) - pco2(i))/(0.01*pco2obs(i)); % 1%

    co3(i)      = est26(i).tp(3).co3;    
    zscore(i,5) = (co3obs(i) - co3(i))/(0.02*co3obs(i)); % 2%
end

% AXES PROPERTIES
set(groot,'defaultAxesFontName','Perpetua',...
    'defaultAxesFontSize',18,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','off',...
    'defaultAxesYMinorTick','off');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Perpetua',...
    'defaultTextInterpreter','latex');


lbl = {'$C_T$', ... 
        '$A_T$',...
        'pH',...
        '$p$CO$_2$',...
        '[CO$_{3}^{2-}$]$_T$'};

ddgreen = [0, 102/255, 0];
clr = ddgreen;

b = boxchart(zscore(:,1:5));
b.JitterOutliers = 'on';
b.MarkerStyle = '.';
b.MarkerSize = 8;
b.LineWidth = 2.0;
b.BoxWidth = 0.85; 
b.BoxFaceColor = clr; % green
b.BoxEdgeColor = clr;
b.MarkerColor = clr;
b.WhiskerLineColor = clr;

ax = gca;
ax.FontSize = 18;
ylabel('Z-scores: ( Meas - Calc ) / ${u_{meas}}$')
ax.XTickMode = 'auto';
ax.TickDir = 'out';
xticklabels(lbl); 
ax.YTick = [-5, -3, -2, 0, 2, 3,  5];
ylim([-5 5])
grid on

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);

print(h,'f5_meas_zscores.pdf','-dpdf','-r0');





