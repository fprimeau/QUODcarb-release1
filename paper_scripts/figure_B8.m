
% figure B8, pTOT Z-scores plot

% script used to plot figure B8 of Fennell & Primeau, 2024
% user needs to have run 'driver.m' and have output est26.mat file
% accessible in 'output_mat_files'


% load and parse data
load data.mat
[in] = data;
nD = length(in);

[siobs]     = in(8,:)';     % TSi (umol/kg)
[tpobs]     = in(7,:)';     % TP (umol/kg)

% reset zero to 1e-3 umol/kg, this is necessary for QUODcarb
for i = 1:nD
    if tpobs(i) == 0
        tpobs(i) = 1e-3;
    end
    if siobs(i) == 0
        siobs(i) = 1e-3;
    end
end

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


% calculate pTOT's using formulations and posterior S
% then, calculate Z-score: ( (meas-calc)/sigma_expm )
for i = 1:nD
    sal = est26(i).sal;

    % calc_pTOT to use formulation for prior TB, TS, TF, TCa
    [pT,~,~,upT] = calc_pTOT(opt,sal);
    % TB
    TB(i)   = q(pT(1))*1e6; % convert from log space, then to umol/kg
    uTB(i)  = ebar(pT(1),upT(1))*1e6; % convert from log space, then to umol/kg
    zTB(i)  = (TB(i) - est26(i).TB)/uTB(i); 
    % TS
    TS(i)   = q(pT(2))*1e6; 
    uTS(i)  = ebar(pT(2),upT(2))*1e6;
    zTS(i)  = (TS(i) - est26(i).TS)/uTS(i);
    % TF
    TF(i)   = q(pT(3))*1e6; 
    uTF(i)  = ebar(pT(3),upT(3))*1e6;
    zTF(i)  = (TF(i) - est26(i).TF)/uTF(i);
    % TCa
    TCa(i)  = q(pT(4))*1e6; 
    uTCa(i) = ebar(pT(4),upT(4))*1e6;
    zTCa(i) = (TCa(i) - est26(i).TCa)/uTCa(i);

    % measured TP and TSi
    zTP(i)= (tpobs(i) - est26(i).TP)/(0.0040); %  median 2% 
    zTSi(i) = (siobs(i) - est26(i).TSi)/(0.0620); % median 2%
end

% measured totals: TP, TSi
zscore_meas(:,1) = zTP; % input sig
zscore_meas(:,2) = zTSi;

% formulated totals: TB, TS, TF, TCa
zscore_form(:,1) = zTB;
zscore_form(:,2) = zTS;
zscore_form(:,3) = zTF;
zscore_form(:,4) = zTCa;

% AXES PROPERTIES
set(groot,'defaultAxesFontName','Perpetua',...
    'defaultAxesFontSize',18,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultAxesXMinorTick','off',...
    'defaultAxesYMinorTick','off');
% TEXT PROPERTIES
set(groot,'defaultTextFontName','Perpetua',...
    'defaultTextInterpreter','latex');

lbl1 = {'P$_T$', 'Si$_T$'};
lbl2 = {'B$_T$', 'S$_T$', 'F$_T$', 'Ca$_T$'};

ddgreen = [0, 102/255, 0];
clr = ddgreen;

tiledlayout(1,2)

nexttile
 
b2 = boxchart(zscore_meas);
b2.JitterOutliers = 'on';
b2.MarkerStyle = '.';
b2.MarkerSize = 8;
b2.LineWidth = 2.0;
b2.BoxWidth = 0.8;
b2.BoxFaceColor = clr; % was dblue
b2.BoxEdgeColor = clr;
b2.MarkerColor = clr;
b2.WhiskerLineColor = clr;

grid on

t = title('Measured Totals');
t.FontSize = 18;
xticklabels(lbl1);
ylabel('Z-scores: ( Meas - Calc ) / ${u_{meas}}$');

nexttile
 
b3 = boxchart(zscore_form);
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

t = title('Formulated Totals');
t.FontSize = 18;
xticklabels(lbl2);
ylabel('Z-scores: ( Prior - Posterior ) / ${u_{expm}}$');

h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);


print(h,'figure_B8.pdf','-dpdf','-r0');









% subfunctions

% ----------------------------------------------------------------------
% calc_pTOT
% ----------------------------------------------------------------------

function [pT,gpT,ggpT,upT] = calc_pTOT(opt,S)
% base equations COPIED from co2sys.m Orr et al. (2018) Github
% Originally from van Heuven et al. (2011)
% Original co2sys is from Lewis and Wallace (1998)

% INPUT:
%   S   = Salinity

% OUTPUT:
%   pT  = [   pTB;   pTS;   pTF;   pTCa ]
%           -log10 values of totals (pT)
%  gpT  = [  gpTB;  gpTS;  gpTF;  gpTCa ]
%           first derivatives (gradient of pT wrt S) 
% ggpT  = [ ggpTB; ggpTS; ggpTF; ggpTCa ]
%           second derivatives (Hessian of pT)
%  upT  = [  upTB;  upTS;  upTF;  upTCa ]
%           precisions of pT (w of errors)

    % utility functions
    LOG10   = log(10);
    p       = @(x) -log10(x);
    q       = @(x) 10.^(-x);  % inverse p, i.e., a backward p
    dpdx    = @(x) -1 / (x * LOG10);        % p'
    d2pdx2  = @(x) 1 / (x^2 * LOG10);       % p''
    dqdx    = @(x) -LOG10 * 10.^( -x );     % q'
    d2qdx2  = @(x) LOG10^2 * 10.^( -x );    % q''
    my_abs  = @(x) sqrt(x*x);

    % compute the totals and their derivatives
    [pTB  , gpTB  , ggpTB  , upTB  ] = calc_pTB(opt,S);
    [pTS  , gpTS  , ggpTS  , upTS  ] = calc_pTS(opt,S);
    [pTF  , gpTF  , ggpTF  , upTF  ] = calc_pTF(opt,S);
    [pTCa , gpTCa , ggpTCa , upTCa ] = calc_pTCa(opt,S);

    % ---------------------------------------------------------------------
    % output
    % ---------------------------------------------------------------------
    pT   = [   pTB;   pTS;   pTF;   pTCa ];
    gpT  = [  gpTB;  gpTS;  gpTF;  gpTCa ];
    ggpT = [ ggpTB; ggpTS; ggpTF; ggpTCa ];
    upT  = [  upTB;  upTS;  upTF;  upTCa ];

    % ---------------------------------------------------------------------
    % subfunctions
    % ---------------------------------------------------------------------

    function [pTB,gpTB,ggpTB,upTB] = calc_pTB(opt,S)
        if (opt.TB == 1)
            % Uppstrom, L., Deep-Sea Research 21:161-162, 1974
            % ( copied from Orr's code )
            % TB = ( 0.000232/ 10.811) * (sal/1.80655)
            TB     = 0.0004157 * S / 35;
            pTB    =  p( TB );
            gpTB   = dpdx( TB )   *  0.0004157 / 35 ;
            ggpTB  = d2pdx2( TB ) * (0.0004157 / 35 ) ^ 2;

            % std 5e-6 on avg 2.32e-4 for (B mg kg^-1)/(Cl o/oo)
            TBu     = ( ( (2.32e-4 + 5e-6)/10.811) * S/1.80655 ); 
            TBl     = ( ( (2.32e-4 - 5e-6)/10.811) * S/1.80655 );
            uTB     = (TBu - TBl) /2 ;
            upTB    = my_abs( p(TB + uTB) - pTB ); % mol/kg
            
        elseif (opt.TB == 2)
            % Lee, Kim, Myrne, Millero, Feely, Yong-Ming Liu. 2010.
            % Geochemica Et Cosmochimica Acta 74 (6): 1801-1811.
            % ( copied from Sharp's code )
            % TB = (0.0002414/ 10.811) * (sal/1.80655)
            TB      = 0.0004326 * S / 35;
            pTB     = p(TB);
            gpTB    = dpdx( TB )   *   0.0004326 / 35;
            ggpTB   = d2pdx2( TB ) * ( 0.0004326 / 35 ) ^ 2;
            
            % std 9e-7 on avg 2.414e-4
            TBu     = ( ( (2.414e-4 + 9e-7)/10.811) * S/1.80655);
            TBl     = ( ( (2.414e-4 - 9e-7)/10.811) * S/1.80655);
            uTB     = (TBu - TBl) /2;
            upTB    = my_abs( p(TB + uTB) - pTB ); % mol/kg
            
        elseif (opt.K1K2 == 6) || (opt.K1K2 == 7)
            % this is about 1% lower than Uppstrom's value
            % Culkin, F., in Chemical Oceanography,
            % ed. Riley and Skirrow, 1965: GEOSECS references this
            % (copied from Orr's code)
            TB      = 0.0004106 * S / 35;
            pTB     = p( TB );
            gpTB    = dpdx( TB )      *    0.0004106 / 35 ;
            ggpTB   = sys.d2pdx2( TB ) * ( 0.0004106 / 35 ) ^ 2;

            % can't find paper, assume same as Uppstrom
            % std 5e-6 on avg 2.32e-4
            TBu     = ( ( (2.32e-4 + 5e-6)/10.811) * S/1.80655 );
            TBl     = ( ( (2.32e-4 - 5e-6)/10.811) * S/1.80655 );
            uTB     = (TBu - TBl) /2 ;
            upTB    = my_abs( p(TB + uTB) - pTB ); % mol/kg
        end 
    end

    function [pTS,gpTS,ggpTS,upTS] = calc_pTS(opt,S)
        % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
        % copied from Orr's code
        TS      = ( 0.14 / 96.062 ) * ( S / 1.80655 );
        pTS     = p( TS );
        gpTS    = dpdx( TS )   *   (0.14 / 96.062 ) / 1.80655 ;
        ggpTS   = d2pdx2( TS ) * ( (0.14 / 96.062 ) / 1.80655 ) ^ 2 ;
        
        % 0.14000 ± 0.00023
        TSu     = ( ( (0.14+0.00023)/96.062 ) * S/ 1.80655 );
        TSl     = ( ( (0.14-0.00023)/96.062 ) * S/ 1.80655 );
        uTS     = (TSu - TSl) / 2;
        my_abs  = @(x) sqrt(x*x);
        upTS    = my_abs( p(TS + uTS) - pTS );
    end

    function [pTF,gpTF,ggpTF,upTF] = calc_pTF(opt,S)
        % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
        % this is .000068.*Sali./35. = .00000195.*Sali   
        TF      = ( 0.000067 / 18.998 ) * ( S / 1.80655 ); 
        pTF     = p( TF );
        gpTF    = dpdx( TF )   *   (0.000067 / 18.998 ) / 1.80655 ;
        ggpTF   = d2pdx2( TF ) * ( (0.000067 / 18.998 ) / 1.80655 ) ^ 2 ;

        % 6.7 ± 0.1 e-5
        TFu     = ( ( (6.7e-5 + 0.1e-5)/18.998) * S/1.80655 );
        TFl     = ( ( (6.7e-5 - 0.1e-5)/18.998) * S/1.80655 );
        uTF     = (TFu - TFl) / 2;
        my_abs  = @(x) sqrt(x*x);
        upTF    = my_abs( p(TF + uTF) - pTF );
    end

    function [pTCa,gpTCa,ggpTCa,upTCa] = calc_pTCa(opt,S)
        if opt.K1K2 == 6 || opt.K1K2 == 7
            % Calculate Ca for GEOSECS, Riley and Skirrow 1965
            TCa     =  ( 0.01026 * S / 35) ;
            pTCa    = p( TCa );
            gpTCa   = dpdx( TCa )   *   0.01026 / 35 ;
            ggpTCa  = d2pdx2( TCa ) * ( 0.01026 / 35 ) ^ 2 ;
        else
            % Calculate Ca, Riley and Tongudai 1967
            % this is 0.010285.*obs.sal./35;
            TCa     = ( 0.02128 / 40.087 * ( S / 1.80655 ) );
            pTCa    = p( TCa ); 
            gpTCa   = dpdx( TCa )   *   ( 0.02128 / 40.087 ) / 1.80655 ; 
            ggpTCa  = d2pdx2( TCa ) * ( ( 0.02128 / 40.087 ) / 1.80655 ) ^ 2 ; 
        end
        % mean 0.02128 ± 0.00006 Ca/Cl ratio (g/kg)/(o/oo)
        TCau    = ( (0.02128 + 6e-5)/ 40.087 * ( S / 1.80655 ) );
        TCal    = ( (0.02128 - 6e-5)/ 40.087 * ( S / 1.80655 ) );
        uTCa    = (TCau - TCal) / 2;
        my_abs  = @(x) sqrt(x*x);
        upTCa   = my_abs( p(TCa + uTCa) - pTCa );
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







