
% scripts used to plot figure 6 of Fennell & Primeau, 2024
% user needs to have run 'driver.m' and have output est26.mat file
% accessible in 'output_mat_files'

% load and parse data
load data.mat
[in] = data;
nD = length(in);

load output_mat_files/est26.mat;

% choose options for opt structure
opt.K1K2 = 10; % option for K1K2 formulation
opt.KSO4 = 1;  % option for KSO4 formulation
opt.KF   = 2;  % option for KF formulation
opt.TB   = 2;  % option for TB formulation
opt.phscale  = 1; % 1 = tot, 2 = sws, 3 = free, 4 = NBS
opt.printcsv = 0; % print est to CSV? 1 = on , 0 = off
opt.co2press = 1; % 1 = on, 0 = off
opt.Revelle  = 0; % 1 = on, 0 = off 
opt.printmes = 0; % 1 = on, 0 = off

opt = check_opt(opt);

% functions
p = @(x) -log10(x); % converts to log space
q = @(x) 10^(-x); % converts from log space back to regular
% ebar: converts log space error to absolute error in regular
ebar = @(px,pe) (0.5 * ( q( px - pe ) - q( px + pe ) ) ); 

% calculate pK's using formulations and posterior T, S, P
% then, calculate Z-score: ( (meas-calc)/sigma_expm )
for i = 1:nD
    sal = est26(i).sal; % posterior salinity

    % tp(1) = pH @ 25C, 0dbar
    % calc_pK to use formulation for pK0, pK1, pK2, pKb
    temp_ph = est26(i).tp(1).T;
    press_ph = est26(i).tp(1).P;
    [pKi,~,upKi] = calc_pK(opt,temp_ph,sal,press_ph); % temp = 25 is tp(2)
    % pK0
    pK0_ph(i)   = pKi(1);     
    upK0_ph(i)  = upKi(1);
    zpK0_ph(i)  = (pK0_ph(i) - est26(i).tp(1).pK0)/upK0_ph(i);
    % pK1
    pK1_ph(i)   = pKi(2);     
    upK1_ph(i)  = upKi(2);
    zpK1_ph(i)  = (pK1_ph(i) - est26(i).tp(1).pK1)/upK1_ph(i);
    % pK2
    pK2_ph(i)   = pKi(3);     
    upK2_ph(i)  = upKi(3);
    zpK2_ph(i)  = (pK2_ph(i) - est26(i).tp(1).pK2)/upK2_ph(i);
    % pKb
    pKb_ph(i)   = pKi(4);     
    upKb_ph(i)  = upKi(4);
    zpKb_ph(i)  = (pKb_ph(i) - est26(i).tp(1).pKb)/upKb_ph(i);

    % tp(2) = pCO2 @ 20C, 0dbar
    % calc_pK to use formulation for pK0, pK1, pK2, pKb
    temp_pc = est26(i).tp(2).T;
    press_pc = est26(i).tp(2).P;
    [pKi,~,upKi] = calc_pK(opt,temp_pc,sal,press_pc); % temp = 20 is tp(3)
    % pK0
    pK0_pc(i)   = pKi(1); 
    upK0_pc(i)  = upKi(1);
    zpK0_pc(i)  = (pK0_pc(i) - est26(i).tp(2).pK0)/upK0_pc(i);
    % pK1
    pK1_pc(i)   = pKi(2); 
    upK1_pc(i)  = upKi(2);
    zpK1_pc(i)  = (pK1_pc(i) - est26(i).tp(2).pK1)/upK1_pc(i);
    % pK2
    pK2_pc(i)   = pKi(3); 
    upK2_pc(i)  = upKi(3);
    zpK2_pc(i)  = (pK2_pc(i) - est26(i).tp(2).pK2)/upK2_pc(i);
    % pKb
    pKb_pc(i)   = pKi(4); 
    upKb_pc(i)  = upKi(4);
    zpKb_pc(i)  = (pKb_pc(i) - est26(i).tp(2).pKb)/upKb_pc(i);

    % tp(3) = CO3 @ 25C, 0dbar
    % calc_pK to use formulation for pK0, pK1, pK2, pKb
    temp_co = est26(i).tp(3).T;
    press_co = est26(i).tp(3).P;
    [pKi,~,upKi] = calc_pK(opt,temp_co,sal,press_co); % temp = 20 is tp(3)
    % pK0
    pK0_co(i)   = pKi(1); 
    upK0_co(i)  = upKi(1);
    zpK0_co(i)  = (pK0_co(i) - est26(i).tp(3).pK0)/upK0_co(i);
    % pK1
    pK1_co(i)   = pKi(2); 
    upK1_co(i)  = upKi(2);
    zpK1_co(i)  = (pK1_co(i) - est26(i).tp(3).pK1)/upK1_co(i);
    % pK2
    pK2_co(i)   = pKi(3); 
    upK2_co(i)  = upKi(3);
    zpK2_co(i)  = (pK2_co(i) - est26(i).tp(3).pK2)/upK2_co(i);
    % pKb
    pKb_co(i)   = pKi(4); 
    upKb_co(i)  = upKi(4);
    zpKb_co(i)  = (pKb_co(i) - est26(i).tp(3).pKb)/upKb_co(i);
end

% combine tp(1) and tp(3), both at 25C & 0dbar
zpK0_25 = [zpK0_ph, zpK0_co];
zpK1_25 = [zpK1_ph, zpK1_co];
zpK2_25 = [zpK2_ph, zpK2_co];
zpKb_25 = [zpKb_ph, zpKb_co];

zscore_pK25(:,1) = zpK0_25;
zscore_pK25(:,2) = zpK1_25;
zscore_pK25(:,3) = zpK2_25;
zscore_pK25(:,4) = zpKb_25;

zscore_pKpc(:,1) = zpK0_pc; % pco2
zscore_pKpc(:,2) = zpK1_pc;
zscore_pKpc(:,3) = zpK2_pc;
zscore_pKpc(:,4) = zpKb_pc;

% next, plot the pK Z-scores

% AXES PROPERTIES
set(groot,'defaultAxesFontName','Perpetua',...
    'defaultAxesFontSize',12,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','off',...
    'defaultAxesYMinorTick','off');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Perpetua',...
    'defaultTextInterpreter','latex');

lbl = {'$pK_0$','$pK_1$','$pK_2$','$pK_B$'};

ddgreen = [0, 102/255, 0];
clr = ddgreen;

tiledlayout(1,2)

nexttile
 
b2 = boxchart(zscore_pK25);
b2.JitterOutliers = 'on';
b2.MarkerStyle = '.';
b2.MarkerSize = 8;
b2.LineWidth = 2.0;
b2.BoxWidth = 0.8;
b2.BoxFaceColor = clr; 
b2.BoxEdgeColor = clr;
b2.MarkerColor = clr;
b2.WhiskerLineColor = clr;

grid on

t = title('$25^{\circ}$C, $101,325 Pa$');
t.FontSize = 15;
xticklabels(lbl);
ylabel('Z-scores: ( Prior - Posterior ) / ${u_{expm}}$');

nexttile
 
b3 = boxchart(zscore_pKpc);
b3.JitterOutliers = 'on';
b3.MarkerStyle = '.';
b3.MarkerSize = 8;
b3.LineWidth = 2.0;
b3.BoxWidth = 0.8;
b3.BoxFaceColor = clr;
b3.BoxEdgeColor = clr;
b3.MarkerColor = clr;
b3.WhiskerLineColor = clr;

grid on

t = title('$20^{\circ}$C, $101,325 Pa$');
t.FontSize = 15;
xticklabels(lbl);
ylabel('Z-scores: ( Prior - Posterior ) / ${u_{expm}}$');


h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);


print(h,'f6_pK_zscores.pdf','-dpdf','-r0');









% subfunctions

% ----------------------------------------------------------------------
% calc_pK
% ----------------------------------------------------------------------

function [pK,gpK,upK] = calc_pK(opt,T,S,P)
% base equations COPIED FROM co2sys.m Orr et al. (2018)  Github
% Originally from  van Heuven et al. (2011)
% Original co2sys is from Lewis and Wallace (1998)

% INPUT:
%   T  = Temp (deg C)
%   S  = Salinity
%   P  = pressure (dbar)

% OUTPUT:
%    pK  = [pK0;pK1;pK2;pKb;pKw;pKs;pKf;pKp1;pKp2;pKp3;pKsi;pKnh4;pKh2s;pp2f;pKar;pKca];
%           -log10 values of equilibrium constants
%   gpK  = [pK_T, pK_S, pK_P]; 
%           first derivatives (gradient of pK wrt T, S, P) 
%           (size: length(pK) x 3 )
%   epK  = [upK0;upK1;upK2;upKb;upKw;upKs;upKf;upKp1;upKp2;upKp3;upKsi;upKnh4;upKh2s;upp2f;upKar;upKca];
%           uncertainties of pK (1 standard deviation)

    TK   = T + 273.15; % convert to Kelvin
    Rgas = 83.14462618; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
    RT   = Rgas * TK;
    RT_T = Rgas;
    
    Pbar    = P / 10; % convert from dbar to bar
    Pbar_P  = 1 / 10;
    A       = 19.924; B = 1000; C = 1.005;
    ions    = @(S) A * S / ( B - C * S); % from DOE handbook
    ions_S  = @(S) ( A * B) / ( B - C * S)^2;
    LOG10   = log(10);
    p       = @(x) -log10(x);
    q       = @(x) 10.^(-x);  % inverse p, i.e., a backward p
    dpdx    = @(x) -1 / (x * LOG10);        % p'
    dqdx    = @(x) -LOG10 * 10.^( -x );     % q'
    
    % corrections for pressure---------------------------------------------
    % Millero 1995, 1992, 1982, 1979; Takahashi et al. 1982;
    %   Culberson & Pytkowicz 1968; Edmond & Gieskes 1970.
    dV      = @(T,a) a(1) + a(2) * T +     a(3) * T^2; 
    dV_T    = @(T,a)        a(2)     + 2 * a(3) * T;
    Ka      = @(T,b) ( b(1) + b(2) * T ) / 1000;
    Ka_T    = @(T,b) b(2) / 1000;
    ppfac   = @(T,Pbar,a,b)... 
              -(( -dV(T,a)    * Pbar  + 0.5 * Ka(T,b)   * Pbar^2 ) / RT )   / ( LOG10 );
    ppfac_T = @(T,Pbar,a,b)...
              -(( -dV_T(T,a) * Pbar   + 0.5 * Ka_T(T,b) * Pbar^2 ) / RT )   / ( LOG10 ) ...
              -(-( -dV(T,a)   * Pbar  + 0.5 * Ka(T,b)   * Pbar^2 ) * RT_T / RT^2 ) / ( LOG10 );        
    ppfac_P = @(T,Pbar,Pbar_P,a,b)... 
              -(( -dV(T,a)         +    Ka(T,b)  * Pbar ) * Pbar_P / RT) / ( LOG10 ) ;
    
    % compute the pK's and their derivatives w.r.t. T,P,and S -------------
    [pp2f    , gpp2f , upp2f ] = calc_p2f(opt,T,RT,RT_T,Pbar,Pbar_P); 
    [pKs     , gpKs  , upKs  ] = calc_pKs(opt,T,S,Pbar,Pbar_P); 
    [pKf     , gpKf  , upKf  ] = calc_pKf(opt,T,S,Pbar,Pbar_P); 
    [pSWS2tot, gpSWS2tot     ] = calc_pSWS2tot(opt,S,pKs,gpKs,pKf,gpKf);
    [pfH     , gpfH  , upfH  ] = calc_pfH(opt,T,S);
    [pK0     , gpK0  , upK0  ] = calc_pK0(opt,T,RT,RT_T,S,Pbar,Pbar_P); 
    [pKb     , gpKb  , upKb  ] = calc_pKb(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot,pfH,gpfH); 
    [pKw     , gpKw  , upKw  ] = calc_pKw(opt,T,S,Pbar,Pbar_P); 
    [pKp1    , gpKp1 , upKp1 ] = calc_pKp1(opt,T,S,Pbar,Pbar_P,pfH,gpfH); 
    [pKp2    , gpKp2 , upKp2 ] = calc_pKp2(opt,T,S,Pbar,Pbar_P,pfH,gpfH); 
    [pKp3    , gpKp3 , upKp3 ] = calc_pKp3(opt,T,S,Pbar,Pbar_P,pfH,gpfH); 
    [pKsi    , gpKsi , upKsi ] = calc_pKsi(opt,T,S,Pbar,Pbar_P,pfH,gpfH); 
    [pK1     , gpK1  , upK1  ] = calc_pK1(opt,T,S,Pbar,Pbar_P,pfH,gpfH,pSWS2tot,gpSWS2tot); 
    [pK2     , gpK2  , upK2  ] = calc_pK2(opt,T,S,Pbar,Pbar_P,pfH,gpfH,pSWS2tot,gpSWS2tot);
    [pKnh4   , gpKnh4, upKnh4] = calc_pKnh4(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot);
    [pKh2s   , gpKh2s, upKh2s] = calc_pKh2s(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot);
    [pKar    , gpKar , upKar ] = calc_pKar(opt,T,S,Pbar,Pbar_P,pfH,gpfH);
    [pKca    , gpKca , upKca ] = calc_pKca(opt,T,S,Pbar,Pbar_P,pfH,gpfH);

    % pressure correction for Ks (Millero, 1995) --------------------------
    a = [ -18.03; 0.0466; 0.000316 ];
    b = [-4.53; 0.09 ];
    pKs     = pKs + ppfac(T,Pbar,a,b);
    pKs_T   = gpKs(1); % isn't this needed?
    pKs_P   = gpKs(3);
    pKs_T   = pKs_T + ppfac_T(T,Pbar,a,b);
    pKs_P   = pKs_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKs(1) = pKs_T;
    gpKs(3) = pKs_P;

    % pressure correction for Kf (Millero, 1995) --------------------------
    a = [ -9.78; -0.009; -0.000942 ];
    b = [ -3.91; 0.054 ];       
    pKf     = pKf + ppfac(T,Pbar,a,b);
    pKf_T   = gpKf(1);
    pKf_P   = gpKf(3);
    pKf_T   = pKf_T + ppfac_T(T,Pbar,a,b);
    pKf_P   = pKf_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKf(1) = pKf_T;
    gpKf(3) = pKf_P;

    % pressure correction for Kb (Millero, 1979) --------------------------
    a = [ -29.48; 0.1622; -0.002608 ];
    b = [  -2.84;   0.0 ];     
    pKb     = pKb + ppfac(T,Pbar,a,b);
    pKb_T   = gpKb(1);
    pKb_P   = gpKb(3);
    pKb_T   = pKb_T + ppfac_T(T,Pbar,a,b);
    pKb_P   = pKb_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKb(1) = pKb_T;
    gpKb(3) = pKb_P;

    % pressure correction for Kw (Millero, 1983) --------------------------
    a = [ -20.02; 0.1119; -0.001409];
    b = [ -5.13; 0.0794 ];
    pKw     = pKw + ppfac(T,Pbar,a,b);
    pKw_T   = gpKw(1);
    pKw_P   = gpKw(3);
    pKw_T   = pKw_T + ppfac_T(T,Pbar,a,b);
    pKw_P   = pKw_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKw(1) = pKw_T;
    gpKw(3) = pKw_P;

    % pressure correction for Kp1 (Millero, 1995; same as Millero, 1983) --
    a = [ -14.51; 0.1211; -0.000321 ];
    b = [  -2.67; 0.0427 ];
    pKp1        = pKp1 + ppfac(T,Pbar,a,b);
    pKp1_T      = gpKp1(1);
    pKp1_P      = gpKp1(3);
    pKp1_T      = pKp1_T + ppfac_T(T,Pbar,a,b);
    pKp1_P      = pKp1_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKp1(1)    = pKp1_T;
    gpKp1(3)    = pKp1_P;
        
    % pressure correction for Kp2 (Millero, 1995; same as Millero, 1983) --
    a = [ -23.12; 0.1758; -0.002647 ];
    b = [ -5.15; 0.09 ]; 
    pKp2        = pKp2 + ppfac(T,Pbar,a,b);
    pKp2_T      = gpKp2(1);
    pKp2_P      = gpKp2(3);
    pKp2_T      = pKp2_T + ppfac_T(T,Pbar,a,b);
    pKp2_P      = pKp2_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKp2(1)    = pKp2_T;
    gpKp2(3)    = pKp2_P;
        
    % pressure correction for Kp3 (Millero, 1995; same as Millero, 1983) --
    a = [ -26.57; 0.202; -0.003042 ];
    b = [ -4.08; 0.0714 ];
    pKp3        = pKp3 + ppfac(T,Pbar,a,b);
    pKp3_T      = gpKp3(1);
    pKp3_P      = gpKp3(3);
    pKp3_T      = pKp3_T + ppfac_T(T,Pbar,a,b);
    pKp3_P      = pKp3_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKp3(1)    = pKp3_T;
    gpKp3(3)    = pKp3_P;

    % pressure correction for Ksi 
    % (Millero, 1995; used the values from boric acid)
    a = [ -29.48; 0.1622; -0.002608 ];
    b =[ -2.84; 0];     
    pKsi        = pKsi + ppfac(T,Pbar,a,b);
    pKsi_T      = gpKsi(1);
    pKsi_P      = gpKsi(3);
    pKsi_T      = pKsi_T + ppfac_T(T,Pbar,a,b);
    pKsi_P      = pKsi_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKsi(1)    = pKsi_T;
    gpKsi(3)    = pKsi_P;

    % pressure correction for K1 (Millero, 1995) --------------------------
    % only for opt.cK1K2 ~=6 & ~=7 & ~=8
    a = [-25.5; 0.1271; 0];
    b = [ -3.08; 0.0877 ];
    pK1     = pK1 + ppfac(T,Pbar,a,b);
    pK1_T   = gpK1(1);
    pK1_P   = gpK1(3);
    pK1_T   = pK1_T + ppfac_T(T,Pbar,a,b);
    pK1_P   = pK1_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpK1(1) = pK1_T;
    gpK1(3) = pK1_P;

    % pressure correction for K2 (Millero, 1995) --------------------------
    % only for opt.cK1K2 ~=6 & ~=7 & ~=8
    a = [ -15.82; -0.0219; 0 ];
    b = [ 1.13; -0.1475 ];
    pK2     = pK2 + ppfac(T,Pbar,a,b);
    pK2_T   = gpK2(1);
    pK2_P   = gpK2(3);
    pK2_T   = pK2_T + ppfac_T(T,Pbar,a,b);
    pK2_P   = pK2_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpK2(1) = pK2_T;
    gpK2(3) = pK2_P;

    % pressure correction for Knh4 (added to CO2SYSv3 by J. Sharp) --------
    a = [ -26.43; 0.0889; -0.000905 ];
    b = [ -5.03; 0.0814 ];
    pKnh4       = pKnh4 + ppfac(T,Pbar,a,b);
    pKnh4_T     = gpKnh4(1);
    pKnh4_P     = gpKnh4(3);
    pKnh4_T     = pKnh4_T + ppfac_T(T,Pbar,a,b);
    pKnh4_P     = pKnh4_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKnh4(1)   = pKnh4_T;
    gpKnh4(3)   = pKnh4_P;
        
    % pressure correction for Kh2s (added to CO2SYSv3 by J. Sharp) --------
    a = [ -11.07; -0.009; -0.000942 ];
    b = [ -2.89;  0.054 ]; 
    pKh2s       = pKh2s + ppfac(T,Pbar,a,b);
    pKh2s_T     = gpKh2s(1);
    pKh2s_P     = gpKh2s(3);
    pKh2s_T     = pKh2s_T + ppfac_T(T,Pbar,a,b);
    pKh2s_P     = pKh2s_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKh2s(1)   = pKh2s_T;
    gpKh2s(3)   = pKh2s_P;
        
    % correct pH scale conversion factors for pressure --------------------
    [pSWS2tot, gpSWS2tot, pFREE2tot, gpFREE2tot] = calc_pSWS2tot(opt,S,pKs,gpKs,pKf,gpKf);

    % find pH scale conversion factor -------------------------------------
    % (pressure corrected)
    if opt.phscale == 1 % pH_total
        phfac = pSWS2tot;
        gphfac = gpSWS2tot;
    elseif opt.phscale == 2 % pH_SWS, they are all on this now
        phfac = 0;
        gphfac = 0;
    elseif opt.phscale == 3 % pH_free
        phfac = pSWS2tot - pFREE2tot;
        gphfac = gpSWS2tot - gpFREE2tot;
    elseif opt.phscale == 4 % pH_NBS
        phfac = pfH;
        gphfac = gpfH;
    else 
        error('Need to input valid pH scale 1-4');
    end

    % convert from SWS to chosen pH scale ---------------------------------    
    pK1   = pK1   + phfac;  gpK1   = gpK1   + gphfac;
    pK2   = pK2   + phfac;  gpK2   = gpK2   + gphfac;
    pKw   = pKw   + phfac;  gpKw   = gpKw   + gphfac;
    pKb   = pKb   + phfac;  gpKb   = gpKb   + gphfac; 
    pKp1  = pKp1  + phfac;  gpKp1  = gpKp1  + gphfac;
    pKp2  = pKp2  + phfac;  gpKp2  = gpKp2  + gphfac;
    pKp3  = pKp3  + phfac;  gpKp3  = gpKp3  + gphfac;
    pKsi  = pKsi  + phfac;  gpKsi  = gpKsi  + gphfac;
    pKnh4 = pKnh4 + phfac;  gpKnh4 = gpKnh4 + gphfac;
    pKh2s = pKh2s + phfac;  gpKh2s = gpKh2s + gphfac;
    % pKar, pKca, pKs, and pKf do not need the conversion
   
    % ---------------------------------------------------------------------
    % output
    % ---------------------------------------------------------------------
    pK   = [  pK0;   pK1;   pK2;    pKb;    pKw;   pKs;   pKf;  pKp1; ...
             pKp2;  pKp3;  pKsi;  pKnh4;  pKh2s;  pp2f;  pKar;  pKca; pfH];

    pK = pK .* opt.mpk(T,S,P);
    
    gpK  = [ gpK0;  gpK1;  gpK2;   gpKb;   gpKw;  gpKs;  gpKf; gpKp1; ...
            gpKp2; gpKp3; gpKsi; gpKnh4; gpKh2s; gpp2f; gpKar; gpKca; gpfH];

    gpK = pK .* opt.gmpk(T,S,P) + gpK.*opt.mpk(T,S,P);
    
    upK  = [ upK0;  upK1;  upK2;   upKb;   upKw;  upKs;  upKf; upKp1; ...
            upKp2; upKp3; upKsi; upKnh4; upKh2s; upp2f; upKar; upKca; upfH];
    upK = upK.*opt.umpk;
    
    % pK0   = 1,  pK1  = 2,  pK2  = 3,  pKb  = 4,  pKw  = 5,  pKs   = 6, 
    % pKf   = 7,  pKp1 = 8,  pKp2 = 9,  pKp3 = 10, pKsi = 11, pKnh4 = 12, 
    % pKh2s = 13, pp2f = 14, pKar = 15, pKca = 16, pfH  = 17 

    % ---------------------------------------------------------------------
    % subfunctions
    % ---------------------------------------------------------------------

    function [pp2f,gpp2f,upp2f] = calc_p2f(opt,T,RT,RT_T,Pbar,Pbar_P)
    % pCO2 to fCO2 conversion, Weiss 1974 Marine Chemistry, 2: 203-215
    % pg. 207 -> valid to within 0.1% (assuming 1 sigma -MF)
        TK = T + 273.15; % convert to Kelvin
        Pstd = 1.01325;
        delC = (57.7 - 0.118.*TK);
        delC_T = -0.118;

        a = [ -1636.75; 12.0408; -0.0327957; 3.16528e-5 ]; 
        
        f  = @(T,a) a(1) + a(2) * T +     a(3) * T^2 +     a(4) * T^3;
        df = @(T,a)        a(2)     + 2 * a(3) * T   + 3 * a(4) * T^2;
  
        if opt.co2press == 0
            % FugFac at 1 atm
            pp2f   = -( ( f(TK,a) + 2 * delC ) * Pstd / RT ) / LOG10;
            pp2f_T = -( ( df(TK,a) + 2 * delC_T ) * Pstd / RT - ...
                    ( f(TK,a) + 2 * delC ) * Pstd * RT_T / RT^2 ) / LOG10;
            pp2f_S = 0;
            pp2f_P = 0;
        elseif opt.co2press == 1
           % FugFac at in situ pressure
           pp2f   = -( ( f(TK,a) + 2 * delC ) * (Pstd+Pbar) / RT ) / LOG10;
           pp2f_T = -( ( df(TK,a) + 2 * delC_T ) * (Pstd+Pbar) / RT - ...
                   ( f(TK,a) + 2 * delC ) * (Pstd+Pbar) * RT_T / RT^2 ) / LOG10;
           pp2f_S = 0;
           pp2f_P = -( ( f(TK,a) + 2 * delC )*Pbar_P / RT ) / LOG10;
        end
        gpp2f = [pp2f_T, pp2f_S, pp2f_P]; % gradient of p(p2f);

        p2f = q(pp2f);
        up2f = 0.001 * p2f; % 0.1% relative uncertainty
        my_abs = @(x) sqrt(x*x);
        upp2f = my_abs( p(p2f + up2f) - pp2f );
    end
    
    function [pK0,gpK0,upK0] = calc_pK0(opt,T,RT,RT_T,S,Pbar,Pbar_P)
    % calculate K0, Weiss 1974, Marine Chemistry, 2: 203-215
    % "data show a root-mean-square deviation from the final fitted
    % equation of1.4*10^-4 mols/l*atm in K0, or about 0.3%"
    % (assuming 0.3% for 1 sigma -MF)
        TK = T + 273.15; % convert to Kelvin
        TK100 = TK./100;
        TK100_T = 1/100;

        a = [  -60.2409;   93.4517;   23.3585  ];
        b = [  0.023517; -0.023656;  0.0047036 ];

        f  = @(T,a) a(1)  + a(2) / T   + a(3) * log(T);
        df = @(T,a)       - a(2) / T^2 + a(3) / T;
        
        g  = @(T,b) b(1) + b(2) * T +     b(3) * T^2;
        dg = @(T,b)        b(2)     + 2 * b(3) * T; 
        
        pK0   = -(  f(TK100,a) +  g(TK100,b) * S ) / LOG10;
        pK0_T = -( df(TK100,a) + dg(TK100,b) * S ) * TK100_T / LOG10;
        pK0_S = -(                g(TK100,b)     ) / LOG10;
        
        pK0_P = 0; % no pressure correction at 1atm

        if opt.co2press == 1   % 1.01325 Bar = 1 atm
            vCO2 = 32.3; % partial molal volume of CO2 (cm3 / mol)
                         % from Weiss (1974, Appendix, paragraph 3)
            % pressure correction at in situ pressure
            ppfacK0   = - ( ((-Pbar)*vCO2) / RT ) / LOG10 ;
            ppfacK0_T = - ( ((-Pbar)*vCO2) * (-RT_T) / RT^2 ) / LOG10 ;
            ppfacK0_P = - ( (-vCO2 *Pbar_P / RT ) ) / LOG10;

            pK0   = pK0 + ppfacK0 ;
            pK0_T = pK0_T + ppfacK0_T;
            pK0_P = pK0_P + ppfacK0_P ;
        end
        % derivatives aka gradient
        gpK0 = [pK0_T, pK0_S, pK0_P];

        % Weiss (1974) reports 0.2 - 0.3% uncertainty on K0
        K0 = q(pK0);
        uK0 = K0 * 0.003; % 0.3% relative on K0
        my_abs = @(x) sqrt(x*x);
        upK0 = my_abs( p(K0 + uK0) - pK0 );
    end
    
    function [pKs,gpKs,upKs] = calc_pKs(opt,T,S,Pbar,Pbar_P)
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS);

        % all calculated on free pH scale
        % stay on free pH scale (no conversion to SWS or total)
        if opt.KSO4 == 1
            % calculate Ks (Dickson 1990a)------------------------------
            % "goodness of fit: 0.021" (assuming 1 sigma -MF)
            a1 = [  -4276.1;  141.328; -23.093 ]; 
            a2 = [ -13856.0;   324.57; -47.986 ];
            a3 = [  35474.0;  -771.54; 114.723 ];
            a4 = [  -2698.0;      0.0;    0.0  ];
            a5 = [   1776.0;      0.0;    0.0  ];

            f  = @(T,a)  a(1) / T    + a(2) + a(3) * log(T);
            df = @(T,a) -a(1) / T^2         + a(3) / T;
        
            pKs = -( f(TK,a1) + ...
                 f(TK,a2) * sqrt(IonS) + ...
                 f(TK,a3) * IonS + ...
                 f(TK,a4) * IonS^(1.5) + ...
                 f(TK,a5) * IonS^2 ) / LOG10 + p(1 - 0.001005 * S );
            pKs_T = -( df(TK,a1) + ...
                   df(TK,a2) * sqrt(IonS) + ...
                   df(TK,a3) * IonS + ...
                   df(TK,a4) * IonS^(1.5) + ...
                   df(TK,a5) * IonS^2 ) / LOG10;
            pKs_S = -( f(TK,a2) * IonS_S * 0.5 / sqrt(IonS) + ...
                   f(TK,a3) * IonS_S + ...
                   f(TK,a4) * IonS_S * 1.5 * IonS^(0.5) + ...
                   f(TK,a5) * IonS_S * 2.0 * IonS ) / ...
                    LOG10 - 0.001005 * dpdx(1 - 0.001005 * S );
        
            ulnKs = 0.021; % from Dickson 1990a pg. 123
            my_abs = @(x) sqrt(x*x);
            upKs = my_abs( -ulnKs/LOG10 ); % 0.021 on lnKs, -lnK/LOG10 converts to epK
        elseif opt.KSO4 == 2
            % calculate Ks (Khoo et al 1977)----------------------------
            % pg. 33 "the standard deviation from regression is 0.0021 in
            % log Beta_HSO4"
            a1 = [   647.59;  -6.3451;  0.0 ]; 
            a2 = [ 0.019085;  -0.5208;  0.0 ];

            f  = @(T,a)  a(1) / T    + a(2) + a(3) * log(T);
            df = @(T,a) -a(1) / T^2         + a(3) / T;
        
            pKs = ( f(TK,a1) + ...
                    a2(1) * TK + ...
                    ( a2(2) * sqrt(IonS)) ) + ...
                    p(1 - 0.001005 * S );
            pKs_T = -( df(TK,a1) + ...
                    a2(1) );
            pKs_S = ( a2(2) * IonS_S * 0.5 / sqrt(IonS) ) - ...
                   0.001005 * dpdx(1 - 0.001005 * S );

            upKs = 0.0021; % given in CO2SYS from Khoo et al 1977
        elseif opt.KSO4 == 3
            % calculate Ks (Waters and Millero, 2013) ---------------------
            % log(K^(.)_HSO4) = logKS0
            % log((K^(*)_HSO4)/(K^(.)_HSO4)) = logKSK0
            % KS = (10^logKSK0))*(10^logKS0)
            % (^ from Sharp's code)
            % eKs = 0.007, assuming 95% confidence (2 sigma) -MF
            a1 = [  0.0         ;    562.69486;  -102.5154 ]; 
            a2 = [ -0.0001117033;    0.2477538;  -13273.76 ];
            a3 = [  4.24666     ;    -0.152671;  0.0267059 ];
            a4 = [ -0.000042128 ;          0.0;        0.0 ];
            a5 = [  0.2542181   ;  -0.00509534; 0.00071589 ];
            a6 = [ -0.00291179  ; 0.0000209968;        0.0 ];
            a7 =   -0.0000403724;

            f1  = @(T,a)  a(1) / T    + a(2) + a(3) * log(T);
            df1 = @(T,a) -a(1) / T^2         + a(3) / T;
            f2  = @(T,a)  a(1) * (T^2) + a(2) * T + a(3) / T ;
            df2 = @(T,a) 2*a(1)*T + a(2) - a(3) / T^2 ;
            f3  = @(T,a)  a(1) + a(2) * T + a(3) * T * log(T) ;
            df3 = @(T,a)  a(2) + a(3)*log(T) + a(3) ;

            logKS0  =   f1(TK,a1) + f2(TK,a2) ;
            logKSK0 = ( f3(TK,a3) + f2(TK,a4) ) * sqrt(S) + ...
                      ( f3(TK,a5) ) * S + ...
                      ( f3(TK,a6) ) * S^(1.5) + ...
                      ( a7 * S^2 ) ;
            
            pKs = (-logKS0 - logKSK0) + p(1 - 0.001005 * S ); % our 'p' form is -log10
         
            logKS0_T  = df1(TK,a1) + df2(TK,a2) ;
            logKSK0_T = ( df3(TK,a3) + df2(TK,a4) ) + ...
                          df3(TK,a5) + ...
                          df3(TK,a6) ;
            logKSK0_S = ( f3(TK,a3) + f2(TK,a4) ) * (-0.5/sqrt(S)) + ...
                        ( f3(TK,a6) * 0.5 * sqrt(S)) + ...
                        ( a7 * 2 * S );

            pKs_T = (-logKS0_T - logKSK0_T) ;
            pKs_S = (-logKSK0_S) - 0.001005 * dpdx(1-0.001005*S);

            Ks = q(pKs);
            uKs = 0.007/2; % QUODcarb uses 1sigma
            my_abs = @(x) sqrt(x*x);
            upKs = my_abs( p(Ks + uKs) - pKs );
        end
        pKs_P = 0.0;
        gpKs = [pKs_T, pKs_S, pKs_P];
    end
    
    function [pKf,gpKf,upKf] = calc_pKf(opt,T,S,Pbar,Pbar_P)
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS);

        % all calculated on free pH scale
        % stay on free pH scale (no conversion to SWS or total)
        if opt.KF == 1
        % calculate Kf (Dickson and Riley 1979)----------------------------   
            a = 1590.2; b = -12.641; c = 1.525;

            pKf = -( a/TK + b + c * sqrt(IonS) ) ./LOG10 ...
                + p(1 - 0.001005 * S) ; % free pH scale
            pKf_T = -(-a/(TK.^2) ) / LOG10; 
            pKf_S =  -c * sqrtIonS_S / LOG10 + dpdx(1 - 0.001005 * S)...
                * (-0.001005) ;  
        elseif opt.KF == 2
        % calculate Perez and Fraga (1987)---------------------------------
        % (to be used for S: 10-40, T: 9-33) (from Sharp's code)
        % "observed experimental error was not greater than pm 0.05 units
        % of lnBeta_HF" pg. 163 (assume 1sigma -MF)
            a = 874;  b = -9.68;  c = 0.111;

            pKf = -( a/TK + b + c * sqrt(S) ) ./LOG10 ; % free pH scale
            pKf_T = -(-a/(TK.^2) ) / LOG10; 
            pKf_S = ( -c * (0.5/sqrt(S)) )/ LOG10 ;  
        end
        pKf_P = 0.0;
        gpKf = [pKf_T, pKf_S, pKf_P];

        ulnKf = 0.05; % 0.05 on lnKf
        my_abs = @(x) sqrt(x*x);
        upKf = my_abs( -ulnKf/ LOG10 ) ; %  -/LOG10 converts to epK 
        % none found in Dickson's, so take Perez and Fraga's
    end
           
    function [pKb,gpKb,upKb] = calc_pKb(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot,pfH,gpfH)   
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 6 || opt.K1K2 == 7
        % GEOSECS Peng et al.
        % Lyman, John, UCLA Thesis, 1957
        % fit by Li et al, JGR 74:5507-5525, 1969
        % (copied from Orr's code)
            pKb = -( -9.26 + 0.00886 * S + 0.01*T) ; % our 'p' is -log10
            pKb_T = -0.01;
            pKb_S = -0.00886;
            pKb_P = 0;

        % convert pKb from NBS scale to SWS scale
            pKb = pKb - pfH ;
            pKb_T = pKb_T - gpfH(1);
            pKb_S = pKb_S - gpfH(2);
        else
        % calculate Kb (Dickson 1990b)--------------------------------------
            a1 = [  -8966.9; -2890.53; -77.942; 1.728; -0.0996 ];
            a2 = [ 148.0248; 137.1942; 1.62142;  0.0 ;    0.0  ];
            a3 = [ -24.4344;  -25.085; -0.2474;  0.0 ;    0.0  ];
            a4 = [   0.0   ; 0.053105;   0.0  ;  0.0 ;    0.0  ];

            f  = @(S,a) a(1) + (a(2) .* sqrt(S)) + (a(3) .* S) + ...
                (a(4) .* S.^(1.5)) + (a(5) .* (S.^2));
            df = @(S,a)  0.5 .* a(2) .* S^(-0.5) + a(3) + ...
                1.5 .* a(4) .* sqrt(S) + 2 .* a(5) .* S;

            pKb   = -(  (f(S,a1) ./ TK) + ...
                f(S,a2) + (f(S,a3) .* log(TK)) + ...
                (f(S,a4) .* TK) ) ./ LOG10;
            pKb   = pKb - pSWS2tot; % convert from total to SWS pH scale
            pKb_T = -( (-f(S,a1) ./ (TK.^2) ) + ...
                ( f(S,a3) ./ TK ) +  f(S,a4) ) ./ LOG10;
            pKb_T = pKb_T - gpSWS2tot(1);
            pKb_S = -( ( df(S,a1) ./ TK ) + df(S,a2) + ...
                (df(S,a3) .* log(TK)) + (df(S,a4) .* TK) ) ./ LOG10;
            pKb_S = pKb_S - gpSWS2tot(2);
            pKb_P = 0;
        end
        gpKb = [pKb_T, pKb_S, pKb_P];

        ulnKb = 0.004; % pg 764 Dickson (1990b) (assume 1 sigma -MF)
        my_abs = @(x) sqrt(x*x);
        upKb = my_abs( -ulnKb/LOG10 ) ; % convert from lnKb to pKb with -/LOG10
        % none found in Li et al's paper
    end
    
    function [pKw,gpKw,upKw] = calc_pKw(opt,T,S,Pbar,Pbar_P)
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 7
        % calculate Kw (Millero, 1979)-------------------------------------
        % temperature dependent parameters only
            a1 = [ 148.9802; -13847.26; -23.6521 ];
            a2 = [ -79.2447;   3298.72;  12.0408 ];
            a3 = [ -0.019813;      0.0;    0.0   ];

            f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
            df = @(T,a)        - a(2) / T^2  + a(3) / T;

            pKw   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pKw_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pKw_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;
            pKw_P = 0;

            ulnKw = 0.0214; % combine Harned and Owen (1958) 0.0014 (table 1)
                            % plus Culberson and Pytkowicz (1973) 0.020 (table 3)
            my_abs = @(x) sqrt(x*x);
            upKw = my_abs( -ulnKw/LOG10 ) ; % convert to epK with -/LOG10
        elseif opt.K1K2 == 8
        % calculate Kw (Millero 1979)--------------------------------------
        % refit data of Harned and Owen, 1958
            a1 = [ 148.9802; -13847.26; -23.6521 ];

            f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
            df = @(T,a)        - a(2) / T^2  + a(3) / T;

            pKw   = -(  f(TK,a1) ) / LOG10;
            pKw_T = -( df(TK,a1) ) / LOG10;
            pKw_S = -(   0   ) / LOG10;
            pKw_P = 0;
            ulnKw = 0.0014; % table 1
            my_abs = @(x) sqrt(x*x);
            upKw = my_abs( -ulnKw/LOG10 ) ; % convert to epK with -/LOG10            
        else
        % calculate Kw (Millero, 1995)-------------------------------------
            a1 = [ 148.9802; -13847.26; -23.6521 ];
            a2 = [   -5.977;    118.67;   1.0495 ];
            a3 = [ -0.01615;     0.0  ;     0.0  ];

            f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
            df = @(T,a)        - a(2) / T^2  + a(3) / T;

            pKw   = -( f(TK,a1) + ...
                       f(TK,a2) * sqrt(S) + ...
                       f(TK,a3) * S ) / LOG10;
            % pKw   =  f(TK,a2) * (S)  ;
            pKw_T = -( df(TK,a1) + ...
                       df(TK,a2) * sqrt(S) + ...
                       df(TK,a3) * S ) / LOG10;
            pKw_S = -( ( f(TK,a2) .* ( 0.5 ./ sqrt(S) ) ) + ...
                       f(TK,a3) ) ./ LOG10;
            % pKw_S = (  f(TK,a2)  );
            pKw_P = 0;

            ulnKw = 0.01; % pg 670 Millero, 1995
            my_abs = @(x) sqrt(x*x);
            upKw = my_abs( -ulnKw/LOG10 ) ; % convert to epK with -/LOG10
        end 

        gpKw = [pKw_T,pKw_S,pKw_P];
    end
    
    function [pKp1,gpKp1,upKp1] = calc_pKp1(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 7
        % calculate Kp1----------------------------------------------------
        % Peng et al don't include the contribution from this term,
        % but it is so small it doesn't contribute. It needs to be kept
        % so that the routines work ok. (From Orr's code)
            pKp1 = p(0.02) - pfH; % Kp1 = 0.02, convert NBS to SWS pH scale
            pKp1_T = 0 - gpfH(1); 
            pKp1_S = 0 - gpfH(2);        
            pKp1_P = 0;
        else
        % calculate Kp1 (Yao and Millero 1995)-----------------------------       
            a1 = [ -4576.752;    115.54;  -18.453 ];
            a2 = [  -106.736;   0.69171;     0.0  ];
            a3 = [  -0.65643;  -0.01844;     0.0  ];

            f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
            df = @(T,a) - a(1) / T^2        + a(3) / T;

            pKp1   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pKp1_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pKp1_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;
            pKp1_P = 0;
        end
        gpKp1 = [pKp1_T,pKp1_S,pKp1_P];

        upKp1 = 0.09; % pg 84 Yao and Millero, 1995 (assume 1 sigma -MF)
    end
    
    function [pKp2,gpKp2,upKp2] = calc_pKp2(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
        TK = T + 273.15; % convert to Kelvin
        if opt.K1K2 == 7

        % calculate Kp2 (Kester and Pytkowicz, 1967)-----------------------
            pKp2 = -(-9.039 - 1450/TK) / LOG10  ; 
            pKp2 = pKp2 - pfH ; % convert from NBS to SWS pH scale

            pKp2_T = -( 1450 / TK^2 ) / LOG10 - gpfH(1) ;
            pKp2_S = -gpfH(2) ;
            pKp2_P = -gpfH(3) ;
        else
        % calculate Kp2 (Yao and Millero 1995)-----------------------------
            a1 = [ -8814.715;  172.1033; -27.927 ]; 
            a2 = [   -160.34;    1.3566;    0.0  ];
            a3 = [   0.37335;  -0.05778;    0.0  ];
        
            f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
            df = @(T,a) - a(1) / T^2        + a(3) / T;

            pKp2   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pKp2_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pKp2_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;
            pKp2_P = 0 ;
        end
        gpKp2 = [pKp2_T, pKp2_S,pKp2_P];

        upKp2 = 0.03; % pg 84 Yao and Millero, 1995 (assume 1 sigma -MF)
    end
    
    function [pKp3,gpKp3,upKp3] = calc_pKp3(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 7
        % calculate Kp3 (Kester and Pytkowicz, 1967)-----------------------
            pKp3 = -(4.466 - 7276/TK) / LOG10;
            pKp3 = pKp3 - pfH ; % convert from NBS to SWS pH scale

            pKp3_T = -( 7276 / TK^2 ) / LOG10 - gpfH(1) ;
            pKp3_S = -gpfH(2) ;
            pKp3_P = -gpfH(3) ;    
        else
        % calculate Kp3 (Yao and Millero 1995)-----------------------------
            a1 = [    -3070.75; -18.126  ];
            a2 = [    17.27039;  2.81197 ];
            a3 = [   -44.99486; -0.09984 ];
                
            f  = @(T,a)   a(1) / T   + a(2);
            df = @(T,a) - a(1) / T^2;
        
            pKp3   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S) +  f(TK,a3) * S ) / LOG10;
            pKp3_T = -( df(TK,a1) + df(TK,a2) * sqrt(S) + df(TK,a3) * S ) / LOG10;
            pKp3_S = -(              f(TK,a2) * 0.5 / sqrt(S) + f(TK,a3) ) / LOG10;
            pKp3_P = 0 ;   
        end
        gpKp3 = [pKp3_T, pKp3_S, pKp3_P]; 

        upKp3 = 0.2; % pg 84 Yao and Millero, 1995 (assume 1 sigma -MF)
    end

    
    function [pKsi,gpKsi,upKsi] = calc_pKsi(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5 * IonS_S / sqrt(IonS);

        if opt.K1K2 == 7
        % calculate Ksi (Sillen, Martell, and Bjerrum, 1964)---------------
            % Ksi = 4e-10; from CO2SYS
            pKsi = p(4e-10) - pfH ; % also convert NBS scale to SWS pH scale
            pKsi_T = -gpfH(1);
            pKsi_S = -gpfH(2);
            pKsi_P = 0;
        else
        % calculate Ksi (Yao and Millero 1995)-----------------------------         
            a1 = [  -8904.2;    117.4; -19.334 ];
            a2 = [  -458.79;   3.5913;   0.0   ]; 
            a3 = [   188.74;  -1.5998;   0.0   ];
            a4 = [ -12.1652;  0.07871;   0.0   ];                       
            f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
            df = @(T,a) - a(1) / T^2        + a(3) / T;
        
            pKsi = -( f(TK,a1) + ...
                  f(TK,a2) * sqrt(IonS) + ...
                  f(TK,a3) * IonS + ...
                  f(TK,a4) * IonS^2  ) / LOG10 + p(1 - 0.001005 * S);
            pKsi_T = -( df(TK,a1) + ...
                    df(TK,a2) * sqrt(IonS) + ...
                    df(TK,a3) * IonS + ...
                    df(TK,a4) * IonS^2 ) / LOG10;
            pKsi_S = -( f(TK,a2) * sqrtIonS_S + ...
                    f(TK,a3) * IonS_S + ...
                    f(TK,a4) * 2 * IonS_S * IonS) / LOG10 ...
                    -0.001005 * dpdx(1 - 0.001005 * S);
            pKsi_P = 0;
        end
        gpKsi = [pKsi_T, pKsi_S, pKsi_P];

        upKsi = 0.02; % pg 84 Yao and Millero, 1995 (assume 1 sigma -MF)
    end
    
    function [pK1,gpK1,upK1] = calc_pK1(opt,T,S,Pbar,Pbar_P,pfH,gpfH,pSWS2tot,gpSWS2tot)
    % calculate pK1 based on user choice input with opt.cK1K2
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 1
        % Roy et al, Marine Chemistry, 44:249-267, 1993 -------------------
            % (see also: Erratum, Marine Chemistry 45:337, 1994
            % and Erratum, Marine Chemistry 52:183, 1996)
            % Typo: in the abstract on p. 249: in the eq. for lnK1* the
            % last term should have S raised to the power 1.5.
            % They claim standard deviations (p. 254) of the fits as
            % .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
            % They also claim (p. 258) 2s precisions of .004 in pK1 and
            % .006 in pK2. These are consistent, but Andrew Dickson
            % (personal communication) obtained an rms deviation of about
            % .004 in pK1 and .003 in pK2. This would be a 2s precision
            % of about 2% in K1 and 1.5% in K2.
            % T:  0-45  S:  5-45. Total Scale. Artificial seawater.
            % This is eq. 29 on p. 254 and what they use in their abstract:    
            a1 = [ 2.83655     ; -2307.1266; -1.5529413];
            a2 = [ -0.20760841 ; -4.0484   ;  0.0      ];
            a3 = [  0.08468345 ;    0.0    ;  0.0      ];
            a4 = [ -0.00654208 ;    0.0    ;  0.0      ];

            f  = @(T,a) a(1) +  a(2) / T   + a(3) * log(T);
            df = @(T,a)       - a(2) / T^2 + a(3) / T;
            
            pK1 = - (f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                f(TK,a3) .* S +  f(TK,a4) .* (S.^(3/2)) )  ... % total scale
                / LOG10 + p(1 - 0.001005 * S) ; % convert to mol/kg-SW
            pK1 = pK1 - pSWS2tot; % convert from pH total to SWS scale
            pK1_T = -( df(TK,a1) + df(TK,a2) .* sqrt(S) + ...
                df(TK,a3) .* S + df(TK,a4) .* (S.^(3/2)) ) / LOG10;
            pK1_T = pK1_T - gpSWS2tot(1); % convert tot to SWS scale
            pK1_S = -(  f(TK,a2) .* 0.5 ./ sqrt(S) +  f(TK,a3) + ...
                (3/2) .* f(TK,a4) .* sqrt(S) ) / LOG10 ...
                -0.001005 * dpdx(1 - 0.001005 * S) ;
            pK1_S = pK1_S - gpSWS2tot(2); % convert tot to SWS scale
            pK1_P = 0;
            % pass all pK1 out of this function as SWS scale

            upK1 = 0.004/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 2 
        % Goyet and Poisson, Deep-Sea Research, 36(11):1635-1654, 1989 ----
            % The 2s precision in pK1 is .011, or 2.5% in K1.
            % The 2s precision in pK2 is .02, or 4.5% in K2.
            % This is in Table 5 on p. 1652 and what they use in the abstract:
            a = 812.27; b = 3.356; c = -0.00171; d = 0.000091;
            pK1 = a ./ TK + b + c .* S .* log(TK) + d .* (S.^2) ; % SWS scale
            pK1_T = - a ./ (TK.^2) + c .* S ./ TK ;
            pK1_S = c .* log(TK) + 2 .* d .* S ;
            pK1_P = 0;
            upK1 = 0.011/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 3
        % Hansson refit by Dickson and Millero, 1987 ----------------------
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
            % and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
            % on the SWS pH scale in mol/kg-SW.
            % Hansson gave his results on the Total scale (he called it
            % the seawater scale) and in mol/kg-SW.
            % Typo in DM on p. 1739 in Table 4: the equation for pK2*
            % for Hansson should have a .000132 *S^2
            % instead of a .000116 *S^2.
            % The 2s precision in pK1 is .013, or 3% in K1.
            % The 2s precision in pK2 is .017, or 4.1% in K2.
            % This is from Table 4 on p. 1739.
            a = 851.4; b = 3.237; c = -0.0106; d = 0.000105;
            pK1 = a ./ TK + b + c .* S + d .* (S.^2) ; % SWS scale
            pK1_T = - a ./ (TK.^2) ;
            pK1_S = c + 2 .* d .* S ;
            pK1_P = 0;
            upK1 = 0.013/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 4
        % Mehrbach refit by Dickson and Millero, 1987 ---------------------
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Mehrbach et al gave results on the NBS scale.
            % The 2s precision in pK1 is .011, or 2.6% in K1.
            % The 2s precision in pK2 is .020, or 4.6% in K2.
	        % Valid for salinity 20-40.
            % This is in Table 4 on p. 1739.
            a = 3670.7; b = -62.008; c = 9.7944; d = -0.0118; g = 0.000116;
            pK1   =  a ./ TK   + b + c .* log(TK) + d .* S + g .* S^2;
            pK1_T = -a ./ TK^2     + c ./ TK;
            pK1_S = d + 2.*g.*S;
            pK1_P = 0;
            upK1 = 0.011/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 5
        % Hansson and Mehrbach refit by Dickson and Millero, 1987 ---------
            % Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
            % Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
            % and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Typo in DM on p. 1740 in Table 5: the second equation
            % should be pK2* =, not pK1* =.
            % The 2s precision in pK1 is .017, or 4% in K1.
            % The 2s precision in pK2 is .026, or 6% in K2.
	        % Valid for salinity 20-40.
            % This is in Table 5 on p. 1740.    
            a = 845; b = 3.248; c = -0.0098; d = 0.000087;
            pK1   = a ./ TK + b + c .* S + d .* (S.^2) ;
            pK1_T = -a ./ (TK.^2) ;
            pK1_S = c + 2 .* d .* S ;
            pK1_P = 0;
            upK1 = 0.017/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 6 || opt.K1K2 == 7
        % GEOSECS and Peng et al use K1, K2 from Mehrbach et al, ----------
            % Limnology and Oceanography, 18(6):897-907, 1973.
	        % I.e., these are the original Mehrbach dissociation constants.
            % The 2s precision in pK1 is .005, or 1.2% in K1.
            % The 2s precision in pK2 is .008, or 2% in K2.

            a = -13.7201; b = 0.031334; c = 3235.76; 
            d = 1.3e-5; g = -0.1032;

            pK1 = a + b .* TK + c ./ TK + d .* S .* TK + g .* sqrt(S) ; % NBS scale
            pK1 = pK1 - pfH ; % convert from NBS to SWS pH scale
            pK1_T = ( b - c ./ (TK.^2) + d .* S ) - gpfH(1) ; % SWS scale
            pK1_S = ( d .* TK + 0.5 .* g ./ sqrt(S) ) - gpfH(2) ; % SWS scale
            pK1_P = 0;
            upK1 = 0.005/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 8
        % PURE WATER CASE -------------------------------------------------
            % Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
            % K1 from refit data from Harned and Davis,
            % J American Chemical Society, 65:2030-2037, 1943.
            % K2 from refit data from Harned and Scholes,
            % J American Chemical Society, 43:1706-1709, 1941.
	            % This is only to be used for Sal=0 water 
                % (note the absence of S in the below formulations)
            % (0-50 C)
            a = 290.9097; b = -14554.21; c = -45.0575;
            pK1 = -(a + b ./ TK + c .* log(TK) ) / LOG10;
            pK1_T = -( -b ./ (TK.^2) + c ./ TK ) / LOG10;
            pK1_S = 0;
            pK1_P = 0;
            ulnK1 = 0.0024;
            my_abs = @(x) sqrt(x*x);
            upK1 = my_abs( -ulnK1/LOG10 ) ; % convert lnK to epK with -/LOG10

        elseif opt.K1K2 == 9
        % Cai and Wang 1998, for estuarine use ----------------------------
            % Data used in this work is from:
	        % K1: Mehrbach (1973) for S>15, for S<15: Mook and Keone (1975)
	        % K2: Mehrbach (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	        % Sigma of residuals between fits and above data: ±0.015, 
            % +0.040 for K1 and K2, respectively.
	        % Sal 0-40, Temp 0.2-30
            % Limnol. Oceanogr. 43(4) (1998) 657-668    
            a = 3404.71; b = 0.032786; c = -14.8435; d = -0.071692; 
            g = 0.0021487; h = 200.1; l = 0.3220;
            f1 = h ./ TK + l;

            pK1 = a ./ TK + b .* TK + c + d .* f1 .* sqrt(S) + g .* S ;
            pK1 = pK1 - pfH; % convert from NBS scale to SWS scale
            pK1_T = ( -a ./ (TK.^2) + b + d .* sqrt(S) .* (-h ./ (TK.^2)) )...
                - gpfH(1) ; % Temp derivative convert to sws scale
            pK1_S = ( 0.5 .* d .* f1 ./ sqrt(S) + g ) - gpfH(2) ;
            pK1_P = 0;
            upK1 = 0.015;

        elseif opt.K1K2 == 10
        % Leuker, Dickson, Keeling, 2000 ----------------------------------
            % This is Mehrbach's data refit after conversion to the 
            % total scale, for comparison with their equilibrator work. 
            % Mar. Chem. 70 (2000) 105-119
            % rms deviation is 0.0055 in pK1 and 0.0100 in pK2
            a = 3633.86; b = -61.2172; c = 9.6777;
            d = -0.011555; g = 0.0001152;

            pK1 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2);
            pK1 = pK1 - pSWS2tot; % convert from total scale to SWS scale
            pK1_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK1_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);
            pK1_P = 0;
            upK1 = 0.0055;

        elseif opt.K1K2 == 11
        % Mojica Prieto and Millero, 2002 ---------------------------------
            % Geochim. et Cosmochim. Acta. 66(14) 2529-2540
            % sigma for pK1 is reported to be 0.0056
	        % sigma for pK2 is reported to be 0.010
	        % This is from the abstract and pages 2536-2537
            a = -43.6977; b = -0.0129037; c = 1.364e-4;
            d = 2885.378; g = 7.045159;

            pK1 = a + b .* S + c .* (S.^2) + d ./ TK + g .* log(TK) ;
            pK1_T = -d ./ (TK.^2) + g ./ TK ;
            pK1_S = b + 2 .* c .* S ;
            pK1_P = 0;
            upK1 = 0.0056;

        elseif opt.K1K2 == 12
        % Millero et al, 2002 ---------------------------------------------
            % Deep-Sea Res. I (49) 1705-1723.
	        % Calculated from overdetermined WOCE-era field measurements 
	        % sigma for pK1 is reported to be 0.005
	        % sigma for pK2 is reported to be 0.008
	        % This is from page 1716
            a = 6.359; b = -0.00664; c = -0.01322; d = 4.989e-5;
            pK1 = a + b .* S + c .* T + d .* (T.^2) ; % tempC
            pK1_T = c + 2 .* d .* T ;
            pK1_S = b ;
            pK1_P = 0;
            upK1 = 0.005;

        elseif opt.K1K2 == 13
        % Millero et al (2006) --------------------------------------------
            % Millero, Graham, Huang, Bustos-Serrano, Pierrot, 2006
            % Mar.Chem. 100 (2006) 80-94.
            % S=1 to 50, T=0 to 50. On seawater scale (SWS). 
            % From titrations in Gulf Stream seawater.
            % sigma pK1 = 0.0054, sigma pK2 = 0.011, from abstract (-MF)
            a1 = [-126.34048;  13.4191; 0.0331; -5.33e-5] ;
            a2 = [  6320.813; -530.123; -6.103;    0.0  ] ;
            a3 = [ 19.568224; -2.06950;   0.0 ;    0.0  ] ;
            
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK1_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK1_P = 0;
            upK1 = 0.0054; 

        elseif opt.K1K2 == 14
        % Millero 2010, for estuarine use ---------------------------------
            % Marine and Freshwater Research, v. 61, p. 139-142.
	        % Fits through compilation of real seawater titration results:
	        % Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), 
            % Millero et al. (2006)
	        % Constants for K's on the SWS; This is from page 141
            % sigma pK1 = 0.005 & sigma pK2 = 0.010 (-MF)
            a1 = [-126.34048;  13.4038; 0.03206; -5.242e-5] ; 
            a2 = [  6320.813; -530.659; -5.8210;    0.0   ] ;
            a3 = [ 19.568224;  -2.0664;   0.0  ;    0.0   ] ;
                        
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK1_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK1_P = 0;
            upK1 = 0.005; % pg 141

        elseif opt.K1K2 == 15
        % Waters, Millero, and Woosley, 2014 ------------------------------
            % Mar. Chem., 165, 66-67, 2014
            % Corrigendum to "The free proton concentration scale for seawater pH".
	        % Effectively, this is an update of Millero (2010) formulation (WhichKs==14)
            % sigma pK1 = 0.0055 & sigma pK2 = 0.0110 (-MF)
            a1 = [-126.34048;  13.409160; 0.031646; -5.1895e-5] ;
            a2 = [  6320.813;  -531.3642;   -5.713;    0.0    ] ;
            a3 = [ 19.568224; -2.0669166;   0.0   ;    0.0    ] ;

            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK1_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK1_P = 0;
            upK1 = 0.0055;

        elseif opt.K1K2 == 16
        % Sulpis et al, 2020 ----------------------------------------------
            % Ocean Science Discussions, 16, 847-862
            % This study uses overdeterminations of the carbonate system to
            % iteratively fit K1 and K2
            % Relative overall uncertainties ~2.5% (sigK/K) for both
            a = 8510.63; b = -172.4493; c = 26.32996; 
            d = -0.011555; g = 0.0001152; 

            pK1 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2) ;
            pK1 = pK1 - pSWS2tot; % convert from tot to SWS scale
            pK1_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK1_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);
            pK1_P = 0;
            K1 = q(pK1);
            eK1 = 0.025 * K1; % ~2.5 % uncertainty on K, pg 854
            my_abs = @(x) sqrt(x*x);
            upK1 = my_abs( p(K1 + eK1) - pK1 ); 

        elseif opt.K1K2 == 17
        % Schockman and Byrne, 2021 ---------------------------------------
            % Geochimica et Cosmochimica Acta, in press
            % This study uses spectrophotometric pH measurements to determine
            % K1*K2 with unprecedented precision, and presents a new
            % parameterization for K2 based on these determinations
            % K1 is taken from Waters, Millero, and Woosley, 2014, on the total pH scale:
            a1 = [-126.34048;  13.568513; 0.031645; -5.3834e-5] ;
            a2 = [  6320.813;  -539.2304;   -5.635;    0.0    ] ;
            a3 = [ 19.568224; -2.0901396;    0.0  ;    0.0    ] ;
                        
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;
            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1 = pK1 - pSWS2tot ; % convert from pHtot to SWS
            pK1_T = ( -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ) ...
                - gpSWS2tot(1) ;
            pK1_S = ( df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK) ) ...
                - gpSWS2tot(2);
            pK1_P = 0;
            upK1 = 0.0055; % same as Waters and Millero formulation

        end

        gpK1 = [pK1_T, pK1_S, pK1_P];     
    end
    
    function [pK2,gpK2,upK2] = calc_pK2(opt,T,S,Pbar,Pbar_P,pfH,gpfH,pSWS2tot,gpSWS2tot)
    % calculate pK2 based on user choice input with opt.cK1K2
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 1
        % Roy et al, Marine Chemistry, 44:249-267, 1993 ------------------
            % (see also: Erratum, Marine Chemistry 45:337, 1994
            % and Erratum, Marine Chemistry 52:183, 1996)
            % Typo: in the abstract on p. 249: in the eq. for lnK1* the
            % last term should have S raised to the power 1.5.
            % They claim standard deviations (p. 254) of the fits as
            % .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
            % They also claim (p. 258) 2s precisions of .004 in pK1 and
            % .006 in pK2. These are consistent, but Andrew Dickson
            % (personal communication) obtained an rms deviation of about
            % .004 in pK1 and .003 in pK2. This would be a 2s precision
            % of about 2% in K1 and 1.5% in K2.
            % T:  0-45  S:  5-45. Total Scale. Artificial seawater.
            % This is eq. 29 on p. 254 and what they use in their abstract:
            a1 = [  -9.226508 ; -3351.6106 ; -0.2005743 ];
            a2 = [-0.106901773;  -23.9722  ;   0.0      ];
            a3 = [  0.1130822 ;    0.0     ;   0.0      ];
            a4 = [ -0.00846934;    0.0     ;   0.0      ];

            f  = @(T,a) a(1) +  a(2) / T   + a(3) * log(T);
            df = @(T,a)       - a(2) / T^2 + a(3) / T;
            
            pK2 = - (f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                f(TK,a3) .* S +  f(TK,a4) .* (S.^(3/2)) )  ... % total scale
                / LOG10 + p(1 - 0.001005 * S) ; % convert to mol/kg-SW
            pK2 = pK2 - pSWS2tot; % convert from pH total to SWS scale

            pK2_T = -( df(TK,a1) + df(TK,a2) .* sqrt(S) + ...
                df(TK,a3) .* S + df(TK,a4) .* (S.^(3/2)) ) / LOG10;
            pK2_T = pK2_T - gpSWS2tot(1); % convert tot to SWS scale
            pK2_S = -(  f(TK,a2) .* 0.5 ./ sqrt(S) +  f(TK,a3) + ...
                (3/2) .* f(TK,a4) .* sqrt(S) ) / LOG10 ...
                -0.001005 * dpdx(1 - 0.001005 * S) ;
            pK2_S = pK2_S - gpSWS2tot(2); % convert tot to SWS scale
            pK2_P = 0;
            % pass all pK2 out of this function as SWS scale

            upK2 = 0.003/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 2 
        % Goyet and Poisson, Deep-Sea Research, 36(11):1635-1654, 1989 ----
            % The 2s precision in pK1 is .011, or 2.5% in K1.
            % The 2s precision in pK2 is .02, or 4.5% in K2.
            % This is in Table 5 on p. 1652 and what they use in the abstract:
            a = 1450.87; b = 4.604; c = -0.00385; d = 0.000182;
            pK2 = a ./ TK + b + c .* S .* log(TK) + d .* (S.^2) ; % SWS scale
            pK2_T = - a ./ (TK.^2) + c .* S ./ TK ;
            pK2_S = c .* log(TK) + 2 .* d .* S ;
            pK2_P = 0;
            upK2 = 0.02/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 3
        % Hansson refit by Dickson and Millero, 1987 ----------------------
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
            % and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
            % on the SWS pH scale in mol/kg-SW.
            % Hansson gave his results on the Total scale (he called it
            % the seawater scale) and in mol/kg-SW.
            % Typo in DM on p. 1739 in Table 4: the equation for pK2*
            % for Hansson should have a .000132 *S^2
            % instead of a .000116 *S^2.
            % The 2s precision in pK1 is .013, or 3% in K1.
            % The 2s precision in pK2 is .017, or 4.1% in K2.
            % This is from Table 4 on p. 1739.
            a = -3885.4; b = 125.844; c = -18.141; 
            d = -0.0192; g = 0.000132;

            pK2 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2) ; % SWS scale

            pK2_T = - a ./ (TK.^2) + c ./ TK ;
            pK2_S = d + 2 .* g .* S ;
            pK2_P = 0;
            upK2 = 0.017/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 4
        % Mehrbach refit by Dickson and Millero, 1987 ---------------------
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Mehrbach et al gave results on the NBS scale.
            % The 2s precision in pK1 is .011, or 2.6% in K1.
            % The 2s precision in pK2 is .020, or 4.6% in K2.
	        % Valid for salinity 20-40.
            % This is in Table 4 on p. 1739.
            a = 1394.7; b = 4.777; c = -0.0184; d = 0.000118;
            pK2 = a ./ TK + b + c .* S + d .* (S.^2);
            pK2_T = -a ./ (TK.^2);
            pK2_S = c + 2 .* d .* S;
            pK2_P = 0;
            upK2 = 0.020/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 5
        % Hansson and Mehrbach refit by Dickson and Millero, 1987 ---------
            % Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
            % Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
            % and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Typo in DM on p. 1740 in Table 5: the second equation
            % should be pK2* =, not pK1* =.
            % The 2s precision in pK1 is .017, or 4% in K1.
            % The 2s precision in pK2 is .026, or 6% in K2.
	        % Valid for salinity 20-40.
            % This is in Table 5 on p. 1740.    
            a = 1377.3; b = 4.824; c = -0.0185; d = 0.000122;
            pK2   = a ./ TK + b + c .* S + d .* (S.^2) ;
            pK2_T = -a ./ (TK.^2) ;
            pK2_S = c + 2 .* d .* S ;
            pK2_P = 0;
            upK2 = 0.026/2; % QUODcarb uses 1sigma
        
        elseif opt.K1K2 == 6 || opt.K1K2 == 7
        % GEOSECS and Peng et al use K1, K2 from Mehrbach et al, ----------
            % Limnology and Oceanography, 18(6):897-907, 1973.
	        % I.e., these are the original Mehrbach dissociation constants.
            % The 2s precision in pK1 is .005, or 1.2% in K1.
            % The 2s precision in pK2 is .008, or 2% in K2.
            a1 = [ 5371.9654; -128375.28;  1.671221 ] ;
            a2 = [  0.22913 ;    2.136  ; -8.0944e-4] ;
            a3 = [  18.3802 ;  -5617.11 ;     0.0   ] ;
            b = -2194.3055;

            f  = @(T,a) a(1) +  a(2) ./ T   + a(3) .* T ;
            df = @(T,a)    - a(2) ./ (T.^2) + a(3) ./ T ;
            
            pK2 = ( f(TK,a1) + f(TK,a2) .* S + f(TK,a3) .* log10(S) + ...
                b .* log10(TK) ) - pfH ; % convert from NBS to SWS scale

            pK2_T = ( df(TK,a1) + df(TK,a2) .* S + df(TK,a3) .* log10(S) + ...
                b ./ (TK .* LOG10) ) - gpfH(1) ;
            pK2_S = ( f(TK,a2) + f(TK,a3) ./ (S .* LOG10) ) - gpfH(2) ;
            pK2_P = 0;
            upK2 = 0.008/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 8
        % PURE WATER CASE -------------------------------------------------
            % Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
            % K1 from refit data from Harned and Davis,
            % J American Chemical Society, 65:2030-2037, 1943.
            % K2 from refit data from Harned and Scholes,
            % J American Chemical Society, 43:1706-1709, 1941.
	            % This is only to be used for Sal=0 water 
                    % (note the absence of S in the below formulations)
            a = 207.6548; b = -11843.79; c = -33.6485;
            pK2 = -(a + b ./ TK + c .* log(TK) ) / LOG10;
            pK2_T = -( -b ./ (TK.^2) + c ./ TK ) / LOG10;
            pK2_S = 0;
            pK2_P = 0;
            ulnK2 = 0.0033;
            my_abs = @(x) sqrt(x*x);
            upK2 = my_abs( -ulnK2/LOG10 ) ; % convert lnK to epK with -/LOG10

        elseif opt.K1K2 == 9
        % Cai and Wang 1998, for estuarine use ----------------------------
            % Data used in this work is from:
	        % K1: Mehrbach (1973) for S>15, for S<15: Mook and Keone (1975)
	        % K2: Mehrbach (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	        % Sigma of residuals between fits and above data: ±0.015, 
            % +0.040 for K1 and K2, respectively.
	        % Sal 0-40, Temp 0.2-30
            % Limnol. Oceanogr. 43(4) (1998) 657-668    
            a = 2902.39; b = 0.02379; c = -6.4980; d = -0.3191; 
            g = 0.0198; h = -129.24; l = 1.4381;
            f1 = h ./ TK + l;

            pK2 = a ./ TK + b .* TK + c + d .* f1 .* sqrt(S) + g .* S ;
            pK2 = pK2 - pfH; % convert from NBS scale to SWS scale

            pK2_T = ( -a ./ (TK.^2) + b + d .* sqrt(S) .* (-h ./ (TK.^2)) )...
                - gpfH(1) ; % Temp derivative convert to sws scale
            pK2_S = ( 0.5 .* d .* f1 ./ sqrt(S) + g ) - gpfH(2) ;
            pK2_P = 0;
            upK2 = 0.040;

        elseif opt.K1K2 == 10
        % Leuker, Dickson, Keeling, 2000 ----------------------------------
            % This is Mehrbach's data refit after conversion to the 
            % total scale, for comparison with their equilibrator work. 
            % Mar. Chem. 70 (2000) 105-119
            % rms deviation is 0.0055 in pK1 and 0.0100 in pK2 (-MF)
            a = 471.78; b = 25.929; c = -3.16967;
            d = -0.01781; g = 0.0001122;

            pK2 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2);
            pK2 = pK2 - pSWS2tot; % convert from total scale to SWS scale

            pK2_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK2_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);
            pK2_P = 0;
            upK2 = 0.0100;

        elseif opt.K1K2 == 11
        % Mojica Prieto and Millero, 2002 ---------------------------------
            % Geochim. et Cosmochim. Acta. 66(14) 2529-2540
            % sigma for pK1 is reported to be 0.0056
	        % sigma for pK2 is reported to be 0.010
	        % This is from the abstract and pages 2536-2537
            a1 = [ -452.0940;  21263.61; 68.483143 ] ;
            a2 = [ 13.142162; -581.4428; -1.967035 ] ;
            a3 = [ -8.101e-4;  0.259601;    0.0    ] ;

            f = @(T,a) a(1) + a(2) ./ TK + a(3) .* log(TK) ;
            df = @(T,a) -a(2) ./ (TK.^2) + a(3) ./ TK ;

            pK2 = f(TK,a1) + f(TK,a2) .* S + f(TK,a3) .* (S.^2) ;

            pK2_T = df(TK,a1) + df(TK,a2) .* S + df(TK,a3) .* (S.^2) ;
            pK2_S = f(TK,a2) + 2 .* f(TK,a3) .* S ;
            pK2_P = 0;
            upK2 = 0.010;

        elseif opt.K1K2 == 12
        % Millero et al, 2002 ---------------------------------------------
            % Deep-Sea Res. I (49) 1705-1723.
	        % Calculated from overdetermined WOCE-era field measurements 
	        % sigma for pK1 is reported to be 0.005
	        % sigma for pK2 is reported to be 0.008
	        % This is from page 1715-1716
            a = 9.867; b = -0.01314; c = -0.01904; d = 2.448e-5;
            pK2 = a + b .* S + c .* T + d .* (T.^2) ; % tempC
            pK2_T = c + 2 .* d .* T ;
            pK2_S = b ;
            pK2_P = 0;
            upK2 = 0.008;

        elseif opt.K1K2 == 13
        % Millero et al (2006) --------------------------------------------
            % Millero, Graham, Huang, Bustos-Serrano, Pierrot, 2006
            % Mar.Chem. 100 (2006) 80-94.
            % S=1 to 50, T=0 to 50. On seawater scale (SWS). 
            % From titrations in Gulf Stream seawater.
            % sigma pK1 = 0.0054, sigma pK2 = 0.011, from abstract (-MF)
            a1 = [ -90.18333;  21.0894;  0.1248; -3.687e-4 ] ;
            a2 = [  5143.692; -772.483; -20.051;    0.0    ] ;
            a3 = [ 14.613358;  -3.3336;    0.0 ;    0.0    ] ;
            
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK2 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;

            pK2_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK2_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK2_P = 0;
            upK2 = 0.011;

        elseif opt.K1K2 == 14
        % Millero 2010, for estuarine use ---------------------------------
            % Marine and Freshwater Research, v. 61, p. 139-142.
	        % Fits through compilation of real seawater titration results:
	        % Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), 
            % Millero et al. (2006)
	        % Constants for K's on the SWS; This is from page 141
            % sigma pK1 = 0.005 & sigma pK2 = 0.010 (-MF)
            a1 = [ -90.18333;  21.3728;  0.1218; -3.688e-4] ; 
            a2 = [  5143.692; -788.289; -19.189;    0.0   ] ;
            a3 = [ 14.613358;   -3.374;   0.0  ;    0.0   ] ;
                        
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK2 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK2_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK2_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK2_P = 0;
            upK2 = 0.010; % pg 141

        elseif opt.K1K2 == 15
        % Waters, Millero, and Woosley, 2014 ------------------------------
            % Mar. Chem., 165, 66-67, 2014
            % Corrigendum to "The free proton concentration scale for seawater pH".
	        % Effectively, this is an update of Millero (2010) formulation (WhichKs==14)
            % sigma pK1 = 0.0055 & sigma pK2 = 0.0110 (-MF)
            a1 = [  -90.18333;  21.225890; 0.12450870; -3.7243e-4] ;
            a2 = [  5143.692;  -779.3444;  -19.91739;    0.0    ] ;
            a3 = [ 14.613358; -3.3534679;     0.0   ;    0.0    ] ;

            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK2 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK2_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK2_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK2_P = 0;
            upK2 = 0.0110;

        elseif opt.K1K2 == 16
        % Sulpis et al, 2020 ----------------------------------------------
            % Ocean Science Discussions, 16, 847-862
            % This study uses overdeterminations of the carbonate system to
            % iteratively fit K1 and K2
            % Relative overall uncertainties ~2.5% (sigK/K) for both (-MF)
            a = 4226.23; b = -59.4636; c = 9.60817; 
            d = -0.01781; g = 0.0001122; 

            pK2 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2) ;
            pK2 = pK2 - pSWS2tot; % convert from tot to SWS scale

            pK2_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK2_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);
            pK2_P = 0;
            K2 = q(pK2);
            uK2 = 0.025 * K2; % ~2.5 % uncertainty in K, pg 854
            my_abs = @(x) sqrt(x*x);
            upK2 = my_abs( p(K2 + uK2) - pK2 ); 
            
        elseif opt.K1K2 == 17
        % Schockman and Byrne, 2021 ---------------------------------------
            % Geochimica et Cosmochimica Acta, in press
            % This study uses spectrophotometric pH measurements to determine
            % K1*K2 with unprecedented precision, and presents a new
            % parameterization for K2 based on these determinations
            % K1 is taken from Waters, Millero, and Woosley, 2014, on the total pH scale:
            a = 116.8067; b = -3655.02; c = -16.45817; 
            d = 0.04523; g = -0.615; h = -0.0002799; l = 4.969;

            pK2 = a + b ./ TK + c .* log(TK) + d .* S + g .* sqrt(S) + ...
                h .* (S.^2) + l .* (S ./ TK) ;
            pK2 = pK2 - pSWS2tot ; % convert from pHtot to SWS

            pK2_T = ( -b ./ (TK.^2) + c ./ TK - l .* S ./ (TK.^2) ) ...
                - gpSWS2tot(1) ;
            pK2_S = ( d + 0.5 .* g ./ sqrt(S) + 2 .* h .* S + l ./ TK ) ...
                - gpSWS2tot(2);
            pK2_P = 0;
            upK2 = 0.010; % in abstract

        end
        
        gpK2 = [pK2_T, pK2_S, pK2_P]; 
    end
    
    function [pKnh4, gpKnh4,upKnh4] = calc_pKnh4(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot)
    % calcaulate pKnh4
    TK = T + 273.15; % convert to Kelvin
        if (opt.K1K2 == 6 || opt.K1K2 == 7 || opt.K1K2 == 8)
            % Knh4 = 0;
            pKnh4 = Inf;
            pKnh4_T = 0;
            pKnh4_S = 0;
            pKnh4_P = 0;
        else
        % Ammonia, added by Sharp et al 2021, from Clegg and Whitfield (1995)
        a = 9.244605; b = -2729.33; c = 1/298.15; d = 0.04203362; g = -11.24742; 

        f  = @(T,a) a(1)  + a(2) * sqrt(T) + a(3) * T + a(4) / T;
        df = @(T,a)    0.5 *a(2) / sqrt(T) + a(3)     - a(4) / T^2;

        a1 = [ -13.6416     ; 1.176949     ; -0.02860785  ; 545.4834    ];
        a2 = [ -0.1462507   ; 0.0090226468 ; -0.0001471361; 10.5425     ];
        a3 = [ 0.004669309  ; -0.0001691742;   0.0        ; -0.5677934  ];
        a4 = [ -2.354039e-05;   0.0        ;   0.0        ; 0.009698623 ];

        pKnh4 = a + b * (c - (1 / TK)) + (d + (g / TK)) * S^(0.25) + ...
                f(TK,a1) * S^0.5 + ...
                f(TK,a2) * S^1.5 + ...
                f(TK,a3) * S^2   + ...
                f(TK,a4) * S^2.5 + ...
                + p(1 - 0.001005 * S) ...
                - pSWS2tot; % convert from total scale to SWS pH scale 
        
        pKnh4_T = b / (TK^2) - g / (TK^2) * S^(0.25) + ...
                  df(TK,a1) * S^0.5 + ...
                  df(TK,a2) * S^1.5 + ...
                  df(TK,a3) * S^2   + ...
                  df(TK,a4) * S^2.5 ...
                  - gpSWS2tot(1); % convert from total scale to SWS pH scale

        pKnh4_S = 0.25 * (d + g / TK) * S^(-0.75) + ...
                  f(TK,a1) * 0.5 * S^-0.5 + ...
                  f(TK,a2) * 1.5 * S^0.5 + ...
                  f(TK,a3) * 2.0 * S   + ...
                  f(TK,a4) * 2.5 * S^1.5 - ...
                  0.001005 * dpdx(1 - 0.001005 * S) ...
                  - gpSWS2tot(2); % convert from total scale to SWS pH scale 
        pKnh4_P = 0;
        end
        gpKnh4 = [pKnh4_T, pKnh4_S, pKnh4_P];

        upKnh4 = 0.00017; % pg 2416 of Clegg and Whitefield (1995)
    end
    
    function [pKh2s,gpKh2s,upKh2s] = calc_pKh2s(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot)
    % calculate pKh2s
    TK = T + 273.15; % convert to Kelvin
        if (opt.K1K2 == 6 || opt.K1K2 == 7 || opt.K1K2 == 8)
            % Kh2s = 0;
            pKh2s = Inf;
            pKh2s_T = 0;
            pKh2s_S = 0;
            pKh2s_P = 0;
        else
        % Millero et al (1988)
            a = 225.838; b = -13275.3; c = -34.6435;
            d = 0.3449; h = -0.0274;
        
            pKh2s = (- (a + b./TK + c.*log(TK) + d.*sqrt(S) + h.*S) ./ LOG10 ) ...
                   - pSWS2tot; % convert from total scale to SWS pH scale
            pKh2s_T = (- (-b./(TK.^2) + c./TK ) ./ LOG10 ) ...
                   - gpSWS2tot(1); % convert from total scale to SWS pH scale
            pKh2s_S = (- (0.5.*d./sqrt(S) + h) ./ LOG10 ) ...
                      - gpSWS2tot(2); % convert from total scale to SWS pH scale
            pKh2s_P = 0;
        end
        gpKh2s = [pKh2s_T, pKh2s_S, pKh2s_P];

        upKh2s = 0.033; % from Millero et al (1988), in abstract
    end
    
    function [pKar, gpKar,upKar] = calc_pKar(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
    % Aragonite solubility
        TK = T + 273.15;
        Rgas = 83.14462618; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
        RT   = Rgas * TK;
        RT_T = Rgas;

        if opt.K1K2 == 6 || opt.K1K2 == 7
        % calculate Kar-Aragonite for GEOSECS -----------------------------
            % Berner, R. A., American Journal of Science 276:713-730, 1976:
            % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
            a = 0.0000001; b = -34.452;      c = -39.866; 
            d = 110.21;    g = 0.0000075752;
            
            % Berner (p.722) states that he uses 1.48.
            % It appears that 1.45 was used in the GEOSECS calculations.
            Kar = 1.45 .* a .* ( b + c .* S^(1/3) + d .* log10(S) ...
                    + g .* (TK.^2) ) ;
            Kar_T = 1.45 .* a .* 2 .* g .* TK ; 
            Kar_S = 1.45 .* a .* ( (1/3) .* c .* S ^ (-2/3) + ...
                    d ./ (S .* LOG10) ) ;

            pKar = p(Kar); % - pfH ; % 'p' it and convert to SWS pH scale            
            pKar_T = dpdx(Kar) .* Kar_T ; %- gpfH(1) ; 
            pKar_S = dpdx(Kar) .* Kar_S ; %- gpfH(2) ;

            % pressure correction
            pKar = pKar - ((33.3 - 0.22 .* T) .* Pbar ./ RT ) ./ LOG10 ; % T = tempC
            pKar_T = pKar_T - (Pbar .* Rgas .* ...
                    (59.8730 .* T - 9155.988) ./ (RT)^2 ) ./ LOG10;
            pKar_P = - ((33.3 - 0.22 .* T) * Pbar_P ./ RT ) ./ LOG10;

        else
        % calculate Kar-Aragonite (Mucci, 1983) ---------------------------
            a1 = [ -171.945; -0.077993; 2903.293; 71.595] ;
            a2 = [-0.068393; 0.0017276;   88.135;   0.0 ] ;
            b = -0.10018; c = 0.0059415;

            f  = @(T,a) a(1) + a(2) .* T + a(3) ./ T + a(4) .* log10(T) ;
            df = @(T,a)      a(2) - a(3)./ (T.^2) + a(4) ./ (T .* LOG10); 

            log10Kar = f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                        b .* S + c .* S ^(3/2) ;
            pKar = -log10Kar ; % pK = -log10(K);

            pKar_T = - ( df(TK,a1) + df(TK,a2) .* sqrt(S) ) ;
            pKar_S = - ( 0.5 .* f(TK,a2) ./ sqrt(S) + ...
                        b + (3/2) .* c .* sqrt(S) ) ;

            % pressure correction
            d1 = [(-48.76 + 2.8); 0.5304; 0.0 ] ;
            d2 = [-11.76; 0.3692] ;

            pKar = pKar + ppfac(T,Pbar,d1,d2);
            pKar_T = pKar_T + ppfac_T(T,Pbar,d1,d2);
            pKar_P = ppfac_P(T,Pbar,Pbar_P,d1,d2);

            % pKa std = 0.03 (Mucci 1983)
            % std for Sal part = 0.009 (as given in CO2SYSv3)
        end
        gpKar = [pKar_T, pKar_S, pKar_P];

        upKar = 0.009; % from Mucci table 7
    end

    function [pKca, gpKca,upKca] = calc_pKca(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
    % Calcite solubility
        TK = T + 273.15;
        Rgas = 83.14462618; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
        RT   = Rgas * TK;
        RT_T = Rgas;

        if opt.K1K2 == 6 || opt.K1K2 == 7
        % calculate Kca-Calcite (Berner, 1976) ----------------------------
            a = 0.0000001; b = -34.452;      c = -39.866; 
            d = 110.21;    g = 0.0000075752;

            Kca = a .* ( b + c .* S^(1/3) + d .* log10(S) ...
                    + g .* (TK.^2) ) ;
            Kca_T = a .* 2 .* g .* TK ; 
            Kca_S = a .* ( (1/3) .* c .* S ^ (-2/3) + ...
                    d ./ (S .* LOG10) ) ;

            pKca = p(Kca) ;            
            pKca_T = dpdx(Kca) .* Kca_T ; 
            pKca_S = dpdx(Kca) .* Kca_S ; 

            % pressure correction
            pKca = pKca - ((36 - 0.2 .* T) .* Pbar ./ RT ) ./ LOG10 ; % T = tempC
            pKca_T = pKca_T - (Pbar .* Rgas .* ...
                    (54.43 .* T - 9888.03) ./ (RT)^2 ) ./ LOG10;
            pKca_P = - ((36 - 0.2 .* T) *Pbar_P ./ RT ) ./ LOG10;

        else
        % calculate Kca-Calcite (Mucci, 1983) -----------------------------
            a1 = [-171.9065; -0.077993; 2839.319; 71.595] ;
            a2 = [ -0.77712; 0.0028426;   178.34;   0.0 ] ;
            b = -0.07711; c = 0.0041249;

            f  = @(T,a) a(1) + a(2) .* T + a(3) ./ T + a(4) .* log10(T) ;
            df = @(T,a)      a(2) - a(3)./ (T.^2) + a(4) ./ (T .* LOG10); 

            log10Kca = f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                        b .* S + c .* S ^(3/2) ;
            pKca = -log10Kca ; % pK = -log10(K);

            pKca_T = - ( df(TK,a1) + df(TK,a2) .* sqrt(S) );
            pKca_S = - ( 0.5 .* f(TK,a2) ./ sqrt(S) + ...
                        b + (3/2) .* c .* sqrt(S) ) ;

            % pressure correction
            d1 = [-48.76; 0.5304; 0.0 ] ;
            d2 = [-11.76; 0.3692] ;

            pKca = pKca + ppfac(T,Pbar,d1,d2);
            pKca_T = pKca_T + ppfac_T(T,Pbar,d1,d2);
            pKca_P = ppfac_P(T,Pbar,Pbar_P,d1,d2);
            
            % pKca std = 0.03 (Mucci 1983)
            % std for Sal part = 0.01 (as given in CO2SYSv3)
        end
        gpKca = [pKca_T, pKca_S, pKca_P];

        upKca = 0.010; % from Mucci 1983, table 7
    end
    
    function [pSWS2tot,gpSWS2tot,pFREE2tot,gpFREE2tot] = calc_pSWS2tot(opt,S,pKs,gpKs,pKf,gpKf)
    % pH scale conversion factors (not pressure corrected)----------- 

    % calculate TF (Riley 1965)--------------------------------------------
        TF = (0.000067./18.998).*(S./1.80655); % mol/kg-SW
        TF_S = (0.000067/18.998)*(1/1.80655); % derivative wrt S
        pTF = p(TF);
        pTF_S = dpdx(TF).*TF_S;
        
    % calculate TS (Morris & Riley 1966)-----------------------------------
        TS = (0.14./96.062).*(S./1.80655); % mol/kg-SW
        TS_S = (0.14/96.062)*(1/1.80655); % derivative wrt S
        pTS = p(TS);
        pTS_S = dpdx(TS).*TS_S;
        
        top = 1 + q(pTS - pKs) ;
        top_T = dqdx(pTS - pKs) * (-gpKs(1)) ;
        top_S = dqdx(pTS - pKs) * (pTS_S - gpKs(2)) ;
        top_P = dqdx(pTS - pKs) * (-gpKs(3)) ;
        bot = top + q(pTF - pKf) ;
        bot_T = top_T + dqdx(pTF - pKf) * (-gpKf(1)) ;
        bot_S = top_S + dqdx(pTF - pKf) * (pTF_S - gpKf(2));
        bot_P = top_P + dqdx(pTF - pKf) * ( -gpKf(3) );
        pSWS2tot = p(top./bot);
        pSWS2tot_T = (-top_T / (top*LOG10) ) + (bot_T / (bot*LOG10) );
        pSWS2tot_S = (-top_S / (top*LOG10) ) + (bot_S / (bot*LOG10) );
        pSWS2tot_P = (-top_P / (top*LOG10) ) + (bot_P / (bot*LOG10) );
        gpSWS2tot = [pSWS2tot_T, pSWS2tot_S, pSWS2tot_P];

     % FREE2tot = 1 + TS./Ks; ---------------------------------------------
        pFREE2tot = p(top);
        pFREE2tot_T = (-top_T / (top*LOG10) ) ;
        pFREE2tot_S = (-top_S / (top*LOG10) ) ;
        pFREE2tot_P = (-top_P / (top*LOG10) ) ;
        gpFREE2tot = [pFREE2tot_T, pFREE2tot_S, pFREE2tot_P];

        if (S==0)
            pSWS2tot = 0;
            pSWS2tot_T = 0;
            pSWS2tot_S = 0;
            pSWS2tot_P = 0;
            gpSWS2tot = [0,0,0];
                     
            pFREE2tot = 0;
            pFREE2tot_T = 0;
            pFREE2tot_S = 0;
            pFREE2tot_P = 0;
            gpFREE2tot = [0,0,0];
        end
    end

    function [pfH,gpfH,upfH] = calc_pfH(opt,T,S)
        % fH = [H]/(1 + TS/Ks)
        TK = T + 273.15; % convert to Kelvin
        if opt.K1K2 == 8
            fH = 1; 
            pfH = p(fH);
            gpfH_T = 0;
            gpfH_S = 0;
        elseif opt.K1K2 == 7
        % fH def'n: Peng et al, Tellus 39B: 439-458, 1987: ----------------
            % They reference the GEOSECS report, but round the value
            % given there off so that it is about .008 (1%) lower. It
            % doesn't agree with the check value they give on p. 456.
            fH   = ( 1.29 - 0.00204.*(TK) + ...
                (0.00046 - (0.00000148 * (TK))) * (S).^2 ) ;
            pfH = p(fH);
            gpfH_T = dpdx(pfH) - 0.00204 - 0.00000148 * (S.^2)  ;
            gpfH_S = dpdx(pfH) * 2 * (0.00046 - (0.00000148 * (TK))) * S ;
        else
        % fH def'n: Takahashi et al, Ch 3 in GEOSECS v.3, 1982 ------------
            fH = 1.2948 - 0.002036*TK + ...
                ( 0.0004607 - 0.000001475 * (TK) ) .* (S).^2 ;
            fH_T = -0.002036 + ( -0.000001475 ) .* (S).^2  ;
            fH_S = 2 * ( 0.0004607 - ( 0.000001475 .* (TK) ) ) .* (S) ;
            pfH = p(fH);
            gpfH_T = dpdx(fH) * fH_T;
            gpfH_S = dpdx(fH) * fH_S;
        end
        gpfH = [gpfH_T, gpfH_S, 0];
        % assumed independent of pressure
        
        ufH = 0.005; % ± 0.005 on fH from Culberson, Pytkowicz, 
                     % and Hawley 1970 Journal of Marine Research
                     % assume 1 sigma -MF
        my_abs = @(x) sqrt(x*x);
        upfH = my_abs( p(fH + ufH) - pfH ) ; 
    end 
    
end


% ----------------------------------------------------------------------
% check_opt
% ----------------------------------------------------------------------

function [opt] = check_opt(opt)
    % check opt input
    isbad = @(thing) (isempty(thing) & sum(isnan(thing)));

    % opt.printmes
    if ~isfield(opt,'printmes') || isbad(opt.printmes)
        opt.printmes = 1; % default on
    end
    % opt.K1K2
    if ~isfield(opt,'K1K2') || isbad(opt.K1K2)
        if opt.printmes ~= 0
            fprintf('No K1K2 formulation chosen. Assuming opt.K1K2 = 10\n');
        end
        opt.K1K2 = 10; % default K1K2 setting
    elseif opt.K1K2 > 18 || opt.K1K2 < 1 || ...
                opt.K1K2 == 6 || opt.K1K2 == 7 
        if opt.printmes ~= 0
            error('Invalid K1K2 formulation chosen. ');
        end
    end
    % opt.TB
    if ~isfield(opt,'TB') || isbad(opt.TB)
        if opt.printmes ~= 0
            fprintf('No TB formulation chosen. Assuming opt.TB = 2\n');
        end
        opt.TB = 2;
    elseif opt.TB > 2 || opt.TB < 1
        if opt.printmes ~= 0
            fprintf(['Invalid TB formulation chosen. ' ...
                     'Assuming opt.TB = 2\n']);
        end
        opt.TB = 2;
    end
    % opt.KSO4
    if ~isfield(opt,'KSO4') || isbad(opt.KSO4)
        if opt.printmes ~= 0
            fprintf('No KSO4 formulation chosen. Assuming opt.KSO4 = 1\n');
        end
        opt.KSO4 = 1; % default opt.KSO4 setting
    elseif opt.KSO4 > 3 || opt.KSO4 < 1
        if opt.printmes ~= 0
            fprintf(['Invalid KSO4 formulation chosen. ' ...
                     'Assuming opt.KSO4 = 1\n']);
        end
        opt.KSO4 = 1; % default opt.KSO4 setting
    end
    % opt.KF
    if ~isfield(opt,'KF') || isbad(opt.KF)
        if opt.printmes ~= 0
            fprintf('No KF formulation chosen. Assuming opt.KF = 2 \n');
        end
        opt.KF = 2; % default KF
    elseif opt.KF > 2 || opt.KF < 1 
        if opt.printmes ~= 0
            fprintf(['Invalid KF formulation chosen. ' ...
                     'Assuming opt.KF = 2 \n']);
        end
        opt.KF = 2;
    end
    % opt.phscale
    if ~isfield(opt,'phscale') || isbad(opt.phscale)
        error(['No opt.phscale chosen, must choose 1 = tot, ' ...
               '2 = sws, 3 = free, 4 = NBS \n'])
    elseif opt.phscale > 4 || opt.phscale < 1
        eror(['Invalid opt.phscale chosen, must choose 1 = tot, ' ...
              '2 = sws, 3 = free, 4 = NBS \n'])
    end
    % opt.printcsv and opt.fname
    if ~isfield(opt,'printcsv') || isbad(opt.printcsv)
        opt.printcsv = 0; % default off
    elseif opt.printcsv > 1 || opt.printcsv < 0
        if opt.printmes ~= 0
            fprintf('Invalid CSV opt chosen. Assuming opt.csv = 1\n');
        end
    else
        if ~isfield(opt,'fname') || isbad(opt.fname)
            opt.fname = 'QUODcarb_output.csv';
            if opt.printmes ~= 0
                fprintf(['Invalid CSV filename. Assuming opt.fname' ...
                    ' = ''QUODcarb_output.csv'' \n']);
            end
        end
    end
    % opt.co2press
    if ~isfield(opt,'co2press') || isbad(opt.co2press)
        opt.co2press = 1; % on
        if opt.printmes ~=0
            fprintf('No opt.co2press chosen. Assuming opt.co2press = 1 (on). \n');
        end
    end
    % opt.Revelle
    if ~isfield(opt,'Revelle') || isbad(opt.Revelle)
        opt.Revelle = 0;
        if opt.printmes ~= 0
            fprintf('No opt.Revelle chosen. Assuming opt.Revelle = 0 (off). \n');
        end
    end
    % organic alkalinity
    if ~isfield(opt,'pKalpha') || isbad(opt.pKalpha)
        opt.pKalpha = 0; % off
    end
    if ~isfield(opt,'pKbeta') || isbad(opt.pKbeta)
        opt.pKbeta = 0; % off
    end
    % opt.turnoff 
    if ~isfield(opt,'turnoff')
        opt.turnoff.TB = 0;
        opt.turnoff.pK1 = 0;
        opt.turnoff.pK2 = 0;
    end
    if ~isfield(opt.turnoff,'TB') || isbad(opt.turnoff.TB)
        opt.turnoff.TB = 0; % default = not turned off
    end
    if ~isfield(opt.turnoff,'pK1') || isbad(opt.turnoff.pK1)
        opt.turnoff.pK1 = 0; % default = not turned off
    end
    if ~isfield(opt.turnoff,'pK2') || isbad(opt.turnoff.pK2)
        opt.turnoff.pK2 = 0; % default = not turned off
    end
    if opt.turnoff.pK1 == 1 && opt.turnoff.pK2 == 1
        error('opt.turnoff can only turn off pK1 or pK2, not both')
    end
    % mpk, gmpk, umpk
    if ~isfield(opt,'mpk')
        mpk0 = @(T,S,P) 1;
        mpk1 = @(T,S,P) 1;
        mpk2 = @(T,S,P) 1;
        mpk = @(T,S,P) [mpk0(T,S,P);mpk1(T,S,P);mpk2(T,S,P);1;1;1;1;1;1;1;1;1;1;1;1;1;1];
        opt.mpk = mpk;
    end
    if ~isfield(opt,'gmpk')        
        gmpk0 = @(T,S,P) [0 0 0];
        gmpk1 = @(T,S,P) [0 0 0];
        gmpk2 = @(T,S,P) [0 0 0];
        gmpk = @(T,S,P) [gmpk0(T,S,P);...
                         gmpk1(T,S,P);...
                         gmpk2(T,S,P);...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0];
        opt.gmpk = gmpk;
    end
    if ~isfield(opt,'umpk')
        umpk0 = 1;
        umpk1 = 1;
        umpk2 = 1;
        umpk = [umpk0;umpk1;umpk2;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
        opt.umpk = umpk;
    end
    
end




