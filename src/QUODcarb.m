function [est,obs,sys,iflag,opt] = QUODcarb(obs,opt)
%
% OUTPUT:
%   est := posterior estimates of co2-system variables, equilibrium constants, and precisions
%   obs := same as input except that the variable names are all there
%   sys := indexing scheme as created in 'mksys', useful for debugging
% iflag := 0 if solver converged to specified accuracy 
%          1 after reaching maximum number of iterations without converging
%          2 if it was one order of magnitude away from converging at
%               maximum iteration
%   opt := same as input except that all fields are now filled in
%   CSV := optional CSV output, if turned on in opt.printcsv
%
% INPUT:
%   obs  := co2-system measured quantities with overall uncertainties
%   opt  := solver options input by user 
%
% -------------------------------------------------------------------------
%
% SYNTAX example:
%   obs.sal         = salinity;     (PSU)           
%   obs.usal        = sal_uncert;   (±sigma) (default 0.002)     
%   obs.TC          = total_c;      (umol/kg-SW)    
%   obs.uTC         = TC_uncert;    (±sigma) (default 2 umol/kg)
%   obs.TA          = total_alk;    (umol/kg-SW)    
%   obs.uTA         = alk_uncert;   (±sigma) (default 2 umol/kg)
%   obs.tp(1).T     = temp;         (deg C)         
%   obs.tp(1).uT    = temp_uncert;  (±sigma) (default 0.1 degC)
%   obs.tp(1).P     = pressure;     (dbar, 0 = surface)          
%   obs.tp(1).uP    = pres_uncert;  (±sigma) (default 0.1 dbar)
%   obs.tp(1).ph    = ph_meas;      
%   obs.tp(1).uph   = ph_uncert;    (±sigma) (default 0.010)
%
%   opt.K1K2        = 10;           % (Lueker et al 2000)
%   opt.KSO4        = 1;            % (Dickson et al 1990a) 
%   opt.KF          = 2;            % (Perez and Fraga 1987
%   opt.TB          = 2;            % (Lee et al. 2010)
%   opt.phscale     = 1;            % (1=tot, 2=free, 3=sws, 4=nbs)
%   opt.printcsv    = 1;            % (1=on, 0=off)
%   opt.fname       = 'output.csv'; % (CSV filename)
%   opt.printmes    = 1;            % (1=on, 0=off)
%   opt.co2press    = 1;            % (1=on, 0=off)
%   opt.Revelle     = 1;            % (1=on, 0=off)
%
%--------------------------------------------------------------------------
% 
% INPUT OPTIONS:
%
%   opt.K1K2  -> choice of K1 and K2 formulation    (T-range)   (S-range)
%           1 = Roy et al, 1993                     (T:0-45)    (S:5-45)
%           2 = Goyet & Poisson, 1989               (T:-1-40)   (S:10-50)
%           3 = Hansson, 1973          REFIT by Dickson & Millero, 1987
%                                                   (T:2-35)    (S:20-40)
%           4 = Mehrbach et al., 1973  REFIT by Dickson & Millero, 1987
%                                                   (T:2-35)    (S:20-40)
%           5 = Hansson, 1973 and Mehrbach, 1973 
%                                      REFIT by Dickson & Millero, 1987
%                                                   (T:2-35)    (S:20-40)
%           x = x(GEOSECS)            ~NOT AVAILABLE IN QUODcarb~
%           x = x(Peng)               ~NOT AVAILABLE IN QUODcarb~
%           8 = Millero, 1979 PURE WATER ONLY       (T:0-50)    (S: 0.)
%           9 = Cai and Wang, 1998                  (T:2-35)    (S:0-49)
%          10 = Lueker et al., 2000    (DEFAULT)    (T:2-35)    (S:19-43)
%          11 = Mojica Prieto and Millero, 2002     (T:0-45)    (S:5-42)
%          12 = Millero et al., 2002                (T:-1.6-35) (S:34-37)
%          13 = Millero et al., 2006                (T:0-50)    (S:1-50)
%          14 = Millero et al., 2010                (T:0-50)    (S:1-50)
%          15 = Waters, Millero, and Woosley, 2014  (T:0-50)    (S:1-50)
%          16 = Sulpis et al., 2020                 (T:-1.7-32) (S:31-38)
%          17 = Schockman and Byrne, 2021           (T:15-35)   (S:19-41)
%
%   opt.KSO4  -> choice of KSO4 formulation
%           1 = Dickson (1990a)         (DEFAULT)
%           2 = Khoo et al., 1977
%           3 = Waters and Millero, 2013
%
%   opt.KF    -> choice of KF formulation
%           1 = Dickson and Riley, 1979
%           2 = Perez and Fraga, 1987   (DEFAULT)
%
%   opt.TB    -> choice of total borate formulation
%           1 = Uppstrom, 1979
%           2 = Lee et al., 2010        (DEFAULT)
%
%   opt.phscale -> pH scale of choice for calculations, the equilibrium
%                   constants are returned on this scale, REQUIRED
%           1 = total scale
%           2 = free scale
%           3 = seawater scale (SWS)
%           4 = National Bureau of Standards scale (NBS)
%
%   opt.printcsv -> option to print the output est structure to a csv file
%           1 = on
%           0 = off     (DEFAULT)
%
%   opt.fname -> chosen filename for CSV, if opt.printcsv = 1
%                   must be a .csv file with no spaces
%           opt.fname = 'QUODcarb_output.csv' (DEFAULT)
%
%   opt.printmes -> should QUODcarb print messages to your command line?
%                   recommended 'on' for new users
%           1 = on      (DEFAULT)
%           0 = off
%
%   opt.co2press -> turn on or off the pressure dependencies for K0 and
%           pCO2 to fCO2 fugacity factor (p2f)
%           1 = on      (DEFAULT)
%           0 = off
%
%   opt.Revelle -> calculate the Revelle factor? reminder: QUODcarb cannot
%           calculate uRevelle
%           1 = on
%           0 = off     (DEFAULT)
%
%   opt.tol -> change default tolerance to some other value, recommended to
%               stay within 1e-5 and 1e-8
%           opt.tol = 1e-6; (DEFAULT)
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%   est ->  'est' structure with best estimate contains:
%               1. p(value) and p(error) where p(x) = -log10(x)
%               2. value and average error about the value in 'q'
%                   where q(x) = x^(-10)
%               3. upper and lower bounds in 'q' space, not symmetric
%                   about the value in 'q' space
%                   - as given in 'parse_output'
%   obs -> 'obs' structure with updated values as given in 'parse_input'
%               so all variable names have a value, nan or otherwise
%                   - restates what user input and adds omitted fields or
%                   default values
%   sys -> 'sys' structure with the system's variable names and associated
%               indexing scheme, as created in 'mksys' 
%                   -very useful when debugging and/or checking est.C 
%                       output values
%   iflag -> 0 if solver converged to specified tolerance (opt.tol) 
%            1 after reaching maximum number of iterations without converging
%                   - MAXIT in newtn.m
%            2 if it was one order of magnitude away from converging at
%               maximum iteration, consider keeping or rerunning with
%               larger tolerance
%   opt -> 'opt' structure with all variable names as given in 'check_opt'
%                   -restates what user input and adds omitted fields
%
%   csv ->  if opt.printcsv = 1, csv file with most of est populated in a 
%               spreadsheet, contains column headers with labels and units                 
%                    -does not include upper and lower errors
%
%--------------------------------------------------------------------------
%
% Changes? -> the only things you may want to change are: 
%               1. Max Iteration number -> MAXIT in newtn.m
%               2. Print iteration number and F0 value to screen 
%                       -> iprint in newtn.m
%
%--------------------------------------------------------------------------

    opt = check_opt(opt);                   % check opt structure

    sys = mksys(obs(1),opt.phscale,opt.pKalpha,opt.pKbeta);
    nD  = length(obs);  
    nv  = size(sys.K,2);
    nlam = size(sys.K,1)+size(sys.M,1);

    % populate obs, yobs, wobs at each datapoint
    [obs,yobs,wobs,sys] = parse_input(obs,sys,opt,nD);

    for i = 1:nD % loop over the full data set

        z0 = init(opt,yobs(i,:),sys);       % initialize
        tol = opt.tol;                      % tolerance, 1e-6 default
        % negative of the log of the posterior 
        % aka the log improbability (limp for short)
        fun = @(z) limp(z,yobs(i,:),wobs(i,:),obs(i),sys,opt);

        % find the maximum of the posterior probability 
        if (sys.freshwater)
            nlam = sys.nMr+sys.nKr;
            z0 = [z0(sys.ic); zeros(nlam,1)];...
        end
        [zhat,J,iflag(i)] = newtn(z0,fun,tol);
        zhat_saved = zhat;
        % residual f value, to tack onto est
        [~,~,f] = limp(zhat,yobs(i,:),wobs(i,:),obs(i),sys,opt);
        if (sys.freshwater)
            y = zeros(nv,1);
            lamM = zeros(sys.nMr_,1);
            lamK = zeros(sys.nKr_,1);
            y(sys.ic) = zhat(1:sys.nKc);
            lamM(sys.iMr) = zhat(sys.nKc+(1:sys.nMr));
            lamK(sys.iKr) = zhat((sys.nKc+sys.nMr)+(1:sys.nKr));
            zhat = [y;lamM;lamK];
        end
        if (iflag(i) ~=0) && (opt.printmes ~= 0)
            fprintf('Newton''s method iflag = %i at i = %i \n',iflag(i),i);
        end

        % calculate the marginalized posterior uncertainty using Laplace's approximation
        C = inv(J);
        
        if (sys.freshwater)
            C = C(1:sys.nKc,1:sys.nKc);
            bigC = zeros(nv,nv);
            bigC(sys.ic,sys.ic) = C;
            C = bigC;
        else
            C = C(1:nv,1:nv);
        end

        sigx = sqrt(full(diag(C)));
        if (opt.printmes ~= 0)
            if (sum(isnan(sigx)) > 0) || (sum(isinf(sigx)) > 0) 
            fprintf('NaN found in output means faulty run. i = %i\n',i)
            end
        end

        % populate est
        [est(i)] = parse_output(zhat,sigx,sys,f,C);    
        
        % calculate the Revelle factor if opt.Revelle = 1 ('on')
        if opt.Revelle == 1
            for j = 1:length(sys.tp(:))
                % Revelle
                ifree   = sys.tp(j).ifree;

                ei      = zeros(length(ifree),1);
                ei(1)   = 1;
                jac     = sys.tp(j).dcdx_pTAfixed(zhat(ifree));
                z       = ei - ( jac.' ) * ( ( jac * jac.' ) \ ( jac*ei ) );
                est(i).tp(j).Revelle = z(2)/z(1);

                % dpfCO2dpTA (similar to Revelle but TC held fixed)
                jfree   = sys.tp(j).jfree;
                ej      = zeros(length(jfree),1);

                ej(1)   = 1;
                jac     = sys.tp(j).dcdx_pTCfixed(zhat(jfree)) ;
                z       = ej - ( jac.' ) * ( ( jac * jac.' ) \ ( jac*ej ) );
                est(i).tp(j).dpfco2dpTA = z(2)/z(1);
            end
        end
    end
    % PrintCSV if opt.printcsv = 1 using filename opt.fname
    PrintCSV(est,obs,iflag,opt);
end

% -------------------------------------------------------------------------
% Subfunctions
% -------------------------------------------------------------------------

function [g,H,f] = limp(z,y,w,obs,sys,opt)
% [g,H,f] = limp(z,y,w,obs,sys,opt)
%
% Negative log probability for the co2 system  a.k.a. log improbability i.e., limp!
%
% INPUT:
%
%   y  := measured components, non-measured components set to NaN
%   w  := measurement precisions, same size and order as y, 
%           non-measured components set to anything, they are ignored
% gpK  := first order derivatives of pK wrt T,S,P
% ggpK := second order derivatives of pK wrt T,S,P

% OUTPUT:
%
%   f := limp
%   g := grad f w.r.t. x
%   h := hessian of f w.r.t. x

    p   = sys.p;    
    q   = sys.q;
    M   = sys.M;    
    K   = sys.K;
    if (sys.freshwater)
        % remove the relationships for the unused mass totals
        M = M(sys.iMr,:);
        M = M(:,sys.ic);
        % remove the relationships for the unused mass action laws 
        K = K(sys.iKr,:);
        K = K(:,sys.ic);        
    end
        
    nrk     = size(K,1);
    nTP     = length(sys.tp); 
    nv      = size(M,2);
    x       = z(1:nv);      % thermodynamic system state, PP*x is modeled measured vars

    % fill ypK, gypK, ggypK, and with associated calculated pK, gpK, and ggpK values
    % update the pK values based on the new estimate of (T,P)
    if (sys.freshwater)
        big_x = zeros(sys.nKc_,1);
        big_x(sys.ic) = x;
    else
        big_x = x;
    end
    [y, gy, ggy ] = update_y(y,big_x,obs,sys,opt);
    if (sys.freshwater)
        % remove the unused system state variables
        y = y(sys.ic);
        gy = gy(sys.ic,:);
        gy = gy(:,sys.ic);
        ggy = ggy(sys.ic,:,:);
        ggy = ggy(:,sys.ic,:);
        ggy = ggy(:,:,sys.ic);
    end

    % Make a vector of measured quantities
    id  = find(~isnan(y));
    y   = y(id).';
    gy  = gy(id,:);
    ggy = ggy(id,:,:);
    
    % Make a precision matrix
    if (sys.freshwater)
        % remove the precisions for the unused state variables
        w = w(sys.ic);
    end       
    W   = diag(w(id));

    % Build a matrix that Picks out the measured components of x
    if (sys.freshwater)
        I   = eye(sys.nKc);  % for chain rule
    else
        I = eye(nv);
    end
    PP  = I(id,:);  % picking/pick out the measured ones
    e   = PP*x - y; % calculated - measured (minus)
    ge  = PP - gy;
    
    nlam    = size(M,1) + size(K,1);
    lam     = z(nv+1:end);  % Lagrange multipliers 

    % constraint equations
    c   = [  M * q( x ); ...
             K * x ] ;
    f   = 0.5 * e.' * W * e  + lam.' * c ;  % limp, method of lagrange multipliers    
    % -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    dcdx = [ M * diag( sys.dqdx( x ) ); ...
             K  ];

    %iTSP = [[sys.tp(:).iT], sys.isal, [sys.tp(:).iP]];
    g    = [ e.' * W * ge + lam.' * dcdx,  c.' ];

    ddq     =  diag( sys.d2qdx2( x ) ); % q"
        
    [nr,nc] = size(M);
    gg      = zeros(nc);
    for row = (1:nr)
        gg  = gg + lam(row)*diag(M(row,:))*ddq;
    end
    tmp = zeros(length(x));
    eW  = e.'*W;
    for jj = 1:size(ggy,1)
        tmp = tmp + eW(jj)*(squeeze(ggy(jj,:,:)));
    end
    H   = [  ge.'*W*ge-tmp+gg,  dcdx.'    ; ... % derivatives wrt lambdas
             dcdx            ,  zeros(nlam)  ]; % derivatives wrt var's
    g = g(:); % make sure g is returned as a column vector
end

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
% mksys
% ----------------------------------------------------------------------

function sys = mksys(obs,phscale,optpKalpha,optpKbeta) % ORG ALK
%
% Private function for QUODcarb.m
% it creates the K and M matrices
%
% utility functions and constants
    LOG10   = log(10);
    p       = @(x) -log10( x );  % inverse of q    
    q       = @(x)  10.^( -x );  % inverse of p
    sys.p   = p;
    sys.q   = q;
    dqdx    = @(x) - LOG10 * 10.^( -x );  % q'
    d2qdx2  = @(x) LOG10^2 * 10.^(-x ); % q"
    dpdx    = @(x) -1 / (x .* LOG10);     % p'
    d2pdx2  = @(x) 1 / (x.^2 * LOG10);  % p"
    sys.dqdx    = dqdx;
    sys.dpdx    = dpdx;
    sys.d2pdx2  = d2pdx2;
    sys.d2qdx2  = d2qdx2;
    isgood = @(thing) (~sum(isnan(thing)) & ~isempty(thing));
    if (~isfield(obs,'tp'))
        error('Need to provide temperature and pressure measurement.')
    end
    field_names = union( fieldnames(obs), fieldnames(obs.tp) );
    if ~ismember('sal',field_names)
        error('Error no salinity specified (obs.sal is missing). Use obs.sal = 0 for freshwater case.')
    end
    if ~ismember('T',field_names)
        error('Error no temperature specified (obs.tp(:).T is missing).')
    end
    if ~ismember('P',field_names)
        error('Error no pressure specified (obs.tp(:).P is missing).')
    end
    if (obs.sal == 0) %freshwater case
        freshwater = true;
        iremove_vars = [];
        iremove_Krows = [];
        iremove_Mrows = [];
    else
        freshwater = false;
    end
        
    % carbonate:
    isal = 1;  sys.isal  = isal;   % salinity
    ipTC = 2;  sys.ipTC  = ipTC;   % p(total carbonate)
    ipTA = 3;  sys.ipTA  = ipTA;   % p(total alkalinity)

    i = 3;

    % borate: Kb = [h][boh4]/[boh3]
    i = i + 1;
    ipTB = i;  sys.ipTB = ipTB; % p(total borate)
        
    % sulfate: Ks = [hf][so4]/[hso4]
    i = i + 1;
    ipTS = i;  sys.ipTS = ipTS; % p(total sulfate)
    
    % fluoride: Kf = [h][F]/[HF]
    i = i + 1;
    ipTF = i;  sys.ipTF = ipTF; % p(total fluoride)
    
    % phosphate: Kp1 = [h][h2po4]/[h3po4], Kp2 = [h][hpo4]/[h2po4], Kp3 = [h][po4]/[hpo4]
    i = i + 1;
    ipTP = i; sys.ipTP = ipTP; % p(total phosphate)

    % KSi = [h][siooh3]/[sioh4]
    i = i + 1;
    ipTSi = i; sys.ipTSi = ipTSi; % p(total silicate)
    
    % Knh4 = [h][nh3]/[nh4]
    i = i + 1;
    ipTNH4 = i; sys.ipTNH4 = ipTNH4; % p(total amonia)
    
    % Kh2s = [h][hs]/[h2s]
    i = i + 1;
    ipTH2S = i; sys.ipTH2S = ipTH2S; % p(total sulfide)
    
    % Kar = [co3][ca]/OmegaAr
    % Kca = [co3][ca]/OmegaCa
    i = i + 1;
    ipTCa = i; sys.ipTCa = ipTCa ; % p(total calcium)

    if optpKalpha == 1
        % Kalpha = [alpha][H]/[Halpha] unknown alkalinity
        i = i + 1;
        ipTAlpha = i; sys.ipTAlpha = ipTAlpha;
    end

    if optpKbeta == 1
        % Kbeta = [beta][H]/[Hbeta] unknown alkalinity
        i = i + 1;
        ipTBeta = i; sys.ipTBeta = ipTBeta;
    end
    
    if (freshwater)
        iremove_vars = union(iremove_vars,[isal,ipTB,ipTH2S,ipTNH4,ipTSi,ipTP,ipTS,ipTF,ipTCa]);
    end

    nTP = length(obs.tp); % number of different (T,P) sub systems
    for j = 1:nTP % loop over (T,P) sub systems
        i = i + 1;    iT = i;     tp(j).iT = iT;  % temperature
        i = i + 1;    iP = i;     tp(j).iP = iP;  % pressure
            
        % K0 = [co2st]/fco2  
        % K1 = [h][hco3]/[co2st]
        % K2 = [h][co3]/[hco3]
        nrk = 3;              i = i + 1; 
        tp(j).ipK0      = i;  i = i + 1; 
        tp(j).ipK1      = i;  i = i + 1;
        tp(j).ipK2      = i;  i = i + 1;
        tp(j).ipfco2    = i;  i = i + 1;
        tp(j).ipco2st   = i;  i = i + 1;
        tp(j).iphco3    = i;  i = i + 1;
        tp(j).ipco3     = i;  i = i + 1;
        tp(j).iph       = i; 
        
        % Kb = [h][boh4]/[boh3]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKb   = i;     i = i + 1;
        tp(j).ipboh4 = i;     i = i + 1;
        tp(j).ipboh3 = i;

        % Kw = [h][oh] 
        nrk = nrk + 1;       i = i + 1;
        tp(j).ipKw = i;      i = i + 1;
        tp(j).ipoh = i;      
     
        % Ks  = [hf][so4]/[hso4]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKs     = i;   i = i + 1;
        tp(j).ipso4    = i;   i = i + 1;
        tp(j).iphso4   = i;
        
        % Kf = [h][F]/[HF]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKf = i;       i = i + 1;
        tp(j).ipF  = i;       i = i + 1;
        tp(j).ipHF = i;
        % Kp1 = [h][h2po4]/[h3po4]
        % Kp2 = [h][hpo4]/[h2po4]
        % Kp3 = [h][po4]/[hpo4]
        nrk = nrk + 3;        i = i + 1;
        tp(j).ipKp1   = i;    i = i + 1;
        tp(j).ipKp2   = i;    i = i + 1;
        tp(j).ipKp3   = i;    i = i + 1;
        tp(j).iph3po4 = i;    i = i + 1;
        tp(j).iph2po4 = i;    i = i + 1;
        tp(j).iphpo4  = i;    i = i + 1;
        tp(j).ippo4   = i;
        
        % KSi = [h][siooh3]/[sioh4]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKsi    = i;   i = i + 1;
        tp(j).ipsiooh3 = i;   i = i + 1;
        tp(j).ipsioh4  = i;
        
        % Knh4 = [h][nh3]/[nh4]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKnh4 = i;     i = i + 1;
        tp(j).ipnh3  = i;     i = i + 1;
        tp(j).ipnh4  = i;
        
        % Kh2s = [h][hs]/[h2s]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKh2s = i;     i = i + 1;
        tp(j).ipHS   = i;     i = i + 1;
        tp(j).ipH2S  = i;    

        % fco2 = pco2 * p2f;
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipp2f     = i;  i = i + 1;
        tp(j).ippco2    = i;  
        
        % Kar = [co3][ca]/OmegaAr
        % Kca = [co3][ca]/OmegaCa
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKar     = i;  i = i + 1;
        tp(j).ipca      = i;  i = i + 1;
        tp(j).ipOmegaAr = i;
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKca     = i;  i = i + 1;
        tp(j).ipOmegaCa = i;  i = i+1;
        tp(j).ipfH      = i; i = i + 1;
        
        %  ph scales 
        tp(j).iph_tot  = i;   i = i + 1;
        tp(j).iph_free = i;   i = i + 1;
        tp(j).iph_sws  = i;   i = i + 1;
        tp(j).iph_nbs  = i;
        % changing the order in this for loop messes up the code -MF 5/2/24

        if optpKalpha == 1
            % Kalpha = [alpha][H]/[Halpha] unknown alkalinity
            nrk = nrk + 1;          i = i + 1;
            tp(j).ipKalpha = i;     i = i + 1; % K(alpha)
            tp(j).ipalpha  = i;     i = i + 1; % p(alpha)
            tp(j).iphalpha = i;                % H-alpha
        end

        if optpKbeta == 1
            % Kbeta = [beta][H]/[Hbeta] unknown alkalinity
            nrk = nrk + 1;          i = i + 1;
            tp(j).ipKbeta  = i;     i = i + 1; % K(beta)
            tp(j).ipbeta   = i;     i = i + 1; % p(beta)
            tp(j).iphbeta  = i;                % H-beta
        end
        
        if (freshwater)
            iremove_vars = ...
                union(iremove_vars,...
                      [tp(j).ipKb,    tp(j).ipboh4,   tp(j).ipboh3,...
                       tp(j).ipKp1,   tp(j).ipKp2,    tp(j).ipKp3,...
                       tp(j).iph3po4, tp(j).iph2po4,  tp(j).iphpo4, tp(j).ippo4,...
                       tp(j).ipKsi,   tp(j).ipsiooh3, tp(j).ipsioh4,...
                       tp(j).ipKnh4,  tp(j).ipnh3,    tp(j).ipnh4,...
                       tp(j).ipKh2s,  tp(j).ipHS,     tp(j).ipH2S,...
                       tp(j).ipKf,    tp(j).ipF,      tp(j).ipHF,...
                       tp(j).ipKs,    tp(j).ipso4,    tp(j).iphso4,...
                       tp(j).ipKca,   tp(j).ipKar,    tp(j).ipca,...
                       tp(j).ipOmegaAr, tp(j).ipOmegaCa...
                      ]);
        end
    end

    nv = i;
    K = sparse(nTP*(nrk+2),nv);
    row = 0;
    for j = 1:nTP
        kr = [];
        kc = [];
        % K0 = [CO2*]/fCO2 ==> -pK0 + pco2st - pfco2 = 0         
        row = row + 1;
        K(row,[ tp(j).ipK0, tp(j).ipco2st, tp(j).ipfco2]) = [ -1, 1, -1 ];
        kc = union(kc,[tp(j).ipK0,tp(j).ipco2st, tp(j).ipfco2]);
        kr = [kr, row];
        
        % K1 = [HCO3][H]/[CO2*] ==> -pK1 + phco3 + ph - pco2st = 0
        row = row + 1;
        K(row,[ tp(j).ipK1, tp(j).iph, tp(j).iphco3, tp(j).ipco2st]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[tp(j).ipK1, tp(j).iph, tp(j).iphco3, tp(j).ipco2st]);
        kr = [kr, row];
        
        % K2 = [CO3][H]/[HCO3]
        row = row + 1;
        K(row,[ tp(j).ipK2, tp(j).iph, tp(j).ipco3, tp(j).iphco3 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[tp(j).ipK2, tp(j).iph, tp(j).ipco3, tp(j).iphco3]);
        kr = [kr, row];
        
        % Kw = [OH][H]
        row = row + 1;
        K(row,[ tp(j).ipKw, tp(j).iph, tp(j).ipoh ]) = [ -1, 1, 1 ];
        kc = union(kc,[tp(j).ipKw, tp(j).iph, tp(j).ipoh]);
        kr = [kr, row];
        
        % Kb = [H][BOH4]/[BOH3]
        row = row+1;
        K(row,[ tp(j).ipKb, tp(j).iph, tp(j).ipboh4, tp(j).ipboh3 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[tp(j).ipKb, tp(j).iph, tp(j).ipboh4, tp(j).ipboh3]);
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end

        % Ks  = [H]free[SO4]/[HSO4] 
        row = row+1;
        K(row,[ tp(j).ipKs, tp(j).iph_free, tp(j).ipso4, tp(j).iphso4 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc, [tp(j).ipKs, tp(j).iph_free, tp(j).ipso4, tp(j).iphso4]);
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end
        
        % Kf = [H]free[F]/[HF]     
        row = row+1;
        K(row,[ tp(j).ipKf, tp(j).iph_free, tp(j).ipF, tp(j).ipHF ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKf, tp(j).iph_free, tp(j).ipF, tp(j).ipHF]);
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end
        
        % Kp1 = [H][H2PO4]/[H3PO4]
        row = row+1;
        K(row,[ tp(j).ipKp1, tp(j).iph, tp(j).iph2po4, tp(j).iph3po4]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKp1, tp(j).iph, tp(j).iph2po4, tp(j).iph3po4]);
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end

        % Kp2 = [H][HPO4]/[H2PO4]
        row = row + 1;
        K(row,[ tp(j).ipKp2, tp(j).iph, tp(j).iphpo4, tp(j).iph2po4 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKp2, tp(j).iph, tp(j).iphpo4, tp(j).iph2po4] );
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end
        
        % Kp3 = [H][PO4]/[HPO4]        
        row = row + 1;
        K(row,[ tp(j).ipKp3, tp(j).iph, tp(j).ippo4,tp(j).iphpo4 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKp3, tp(j).iph, tp(j).ippo4, tp(j).iphpo4]);
        kr = [kr, row];       
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end
            
        % KSi = [H][SiO(OH)3]/[Si(OH)4]
        row = row + 1;
        K(row,[ tp(j).ipKsi, tp(j).iph, tp(j).ipsiooh3, tp(j).ipsioh4 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKsi, tp(j).iph, tp(j).ipsiooh3, tp(j).ipsioh4]);
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end
        
        % Knh4 = [H][NH3]/[NH4+]
        row = row + 1;
        K(row,[ tp(j).ipKnh4, tp(j).iph, tp(j).ipnh3, tp(j).ipnh4]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKnh4, tp(j).iph, tp(j).ipnh3, tp(j).ipnh4]);
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end
        
        % Kh2s = [H][HS]/[H2S]
        row = row + 1;
        K(row,[ tp(j).ipKh2s, tp(j).iph, tp(j).ipHS, tp(j).ipH2S]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKh2s, tp(j).iph, tp(j).ipHS, tp(j).ipH2S]);
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end
            
        % fco2 = pco2 * p2f;
        row = row + 1;
        K(row,[ tp(j).ipfco2, tp(j).ippco2, tp(j).ipp2f ]) = [ -1, 1, 1 ];
        kc = union(kc,[ tp(j).ipfco2, tp(j).ippco2, tp(j).ipp2f]);
        kr = [kr, row];
          
        % Kar = [co3][ca]/OmegaAr ==> -pKar + pco3 + pca - pOmegaAr = 0
        row = row + 1;
        K(row,[ tp(j).ipKar, tp(j).ipco3, tp(j).ipca, tp(j).ipOmegaAr]) = [ -1, 1, 1, -1 ];
        kc = union(kc, [ tp(j).ipKar, tp(j).ipco3, tp(j).ipca, tp(j).ipOmegaAr]);
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end
        
        % Kca = [co3][ca]/OmegaCa ==> -pKca + pco3 + pca - pOmegaCa = 0
        row = row + 1;
        K(row, [ tp(j).ipKca, tp(j).ipco3, tp(j).ipca, tp(j).ipOmegaCa]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKca, tp(j).ipco3, tp(j).ipca, tp(j).ipOmegaCa]);
        kr = [kr, row];
        if (freshwater)
            iremove_Krows = union(iremove_Krows,row);
        end
        
        row = row + 1;
        switch phscale % working ph to phscale
          case 1
            K(row,[ tp(j).iph  tp(j).iph_tot]) = [ -1, 1 ];
            kc = union(kc,[ tp(j).iph, tp(j).iph_tot]);
          case 2
            K(row,[ tp(j).iph tp(j).iph_sws] ) = [ -1, 1 ];
            kc = union(kc,[ tp(j).iph, tp(j).iph_sws]);
          case 3
            K(row,[ tp(j).iph tp(j).iph_free] ) = [ -1, 1 ];
            kc = union(kc,[ tp(j).iph, tp(j).iph_free]);
          case 4
            K(row, [tp(j).iph tp(j).iph_nbs] ) = [ -1, 1 ];
            kc = union(kc,[ tp(j).iph, tp(j).iph_nbs]);
        end        
        tp(j).kphscale = row;
        kr = [kr, row];

        % def'n for nbs
        row = row + 1;
        K(row,[ tp(j).iph_nbs, tp(j).iph_sws, tp(j).ipfH ]) = [ 1, -1, -1 ];
        kc = union(kc, [ tp(j).iph_nbs, tp(j).iph_sws, tp(j).ipfH]);            
        kr = [kr, row];

        nr = 12; % TA, TC, TB, TS, TF, TP, TSi, TNH4, TH2S, TCa, 
                 % ph_tot to ph_free, ph_sws to ph_free

        if optpKalpha == 1
            % Kalpha = [alpha][H]/[Halpha] -> -pKalpha + palpha + pH - pHalpha = 0
            row = row + 1;
            K(row, [tp(j).ipKalpha, tp(j).ipalpha, tp(j).iph, tp(j).iphalpha]) = [-1, 1, 1, -1];
            kc = union(kc,[ tp(j).ipKalpha, tp(j).ipalpha, tp(j).iph, tp(j).iphalpha]);
            kr = [kr, row];
            nr = nr + 1;
        end

        if optpKbeta == 1
            % Kbeta = [beta][H]/[Hbeta] -> -pKbeta + pbeta + pH - pHbeta = 0
            row = row + 1;
            K(row, [tp(j).ipKbeta, tp(j).ipbeta, tp(j).iph, tp(j).iphbeta]) = [-1, 1, 1, -1];
            kc = union(kc,[ tp(j).ipKbeta, tp(j).ipbeta, tp(j).iph, tp(j).iphbeta]);
            kr = [kr, row];
            nr = nr + 1;
        end

        tp(j).kr = kr;
        tp(j).kc = kc;
    end

    % "mass conservation" equations
    M = sparse(nTP*nr,nv);
    row = 0;
    for j = 1:nTP
        mr = [];
        mc = [];

        % Total alkalinity: TA - [HCO3] - 2[CO3] ( - [OH] )...
        row = row + 1;
        row_alk = row;
        tp(j).row_alk = row_alk;
        mr = [mr, row];

        % carbonate alkalinity
        M(row_alk,[ ipTA, tp(j).iphco3, tp(j).ipco3]) = [ 1, -1, -2 ];
        mc = union(mc,[ipTA, tp(j).iphco3, tp(j).ipco3]);
        
        % Total carbonate: TC - [CO2*] - [HCO3] - [CO3] = 0
        row = row + 1;
        M(row, [ ipTC, tp(j).ipco2st, tp(j).iphco3, tp(j).ipco3 ]) = [ 1, -1, -1, -1 ];
        mc = union(mc,[ipTC, tp(j).ipco2st, tp(j).iphco3, tp(j).ipco3]);
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e2;
        
        % Total borate
        row = row + 1;
        M(row,[ ipTB, tp(j).ipboh3, tp(j).ipboh4 ])   =  [ 1, -1, -1 ];
        mc = union(mc,[ipTB, tp(j).ipboh3, tp(j).ipboh4]);
        M(row_alk, tp(j).ipboh4) = -1;
        M(row_alk, tp(j).ipoh) = -1;        
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e3;
        if (freshwater)
            iremove_Mrows = union(iremove_Mrows,row);
        end
        
        % Total sulfate
        row = row + 1;
        M(row, [ ipTS, tp(j).iphso4, tp(j).ipso4 ])   =  [ 1, -1, -1 ];
        M(row_alk,[ tp(j).iph_free, tp(j).iphso4 ])   =  [1, 1]; 
        mc = union(mc,[ipTS, tp(j).iphso4, tp(j).ipso4]);
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e1;
        if (freshwater)
            iremove_Mrows = union(iremove_Mrows,row);
        end
        
        % Total fluoride
        row = row + 1;
        M(row, [ ipTF, tp(j).ipHF, tp(j).ipF ]) =  [ 1, -1, -1 ];
        mc = union(mc,[ipTF, tp(j).ipHF, tp(j).ipF]);
        M(row_alk, tp(j).ipHF) =  1;
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e3;
        if (freshwater)
            iremove_Mrows = union(iremove_Mrows,row);
        end
        
        % Total phosphate
        row = row + 1;
        M(row, [ ipTP, tp(j).iph3po4, tp(j).iph2po4, tp(j).iphpo4, tp(j).ippo4 ]) = [ 1, -1, -1, -1, -1 ];
        mc = union(mc,[ipTP, tp(j).iph3po4, tp(j).iph2po4, tp(j).iphpo4, tp(j).ippo4]);
        M(row_alk, [ tp(j).iphpo4, tp(j).ippo4, tp(j).iph3po4 ]) = [ -1, -2, 1 ]; 
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e7;
        if (freshwater)
            iremove_Mrows = union(iremove_Mrows,row);
        end
        
        % Total silicate
        row = row + 1;
        M(row, [ ipTSi, tp(j).ipsioh4, tp(j).ipsiooh3]) = [ 1, -1, -1 ];
        mc = union(mc,[ipTSi, tp(j).ipsioh4, tp(j).ipsiooh3]);
        M(row_alk, tp(j).ipsiooh3) = -1; 
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e6;
        if (freshwater)
            iremove_Mrows = union(iremove_Mrows,row);
        end
        
        % Total amonia
        row = row + 1;
        M(row, [ ipTNH4, tp(j).ipnh4, tp(j).ipnh3 ]) = [ 1, -1, -1 ];
        mc = union(mc,[ipTNH4, tp(j).ipnh4, tp(j).ipnh3]);
        M(row_alk,tp(j).ipnh3) = -1; 
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e8;
        if (freshwater)
            iremove_Mrows = union(iremove_Mrows,row);
        end
        
        % Total sulfide
        row = row + 1;
        M(row,[ ipTH2S, tp(j).ipH2S, tp(j).ipHS ]) = [ 1, -1, -1 ];
        mc = union(mc,[ipTH2S, tp(j).ipH2S, tp(j).ipHS]);
        M(row_alk,tp(j).ipHS) = -1; 
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e8;
        M(row_alk,:) = M(row_alk,:)*1e2;
        if (freshwater)
            iremove_Mrows = union(iremove_Mrows,row);
        end
        
        % total Ca
        row = row + 1;
        M(row, [ ipTCa, tp(j).ipca]) = [ 1, -1 ];
        mc = union(mc,[ipTCa, tp(j).ipca]);
        mr = [mr,row];
        M(row,:) = M(row,:)*1e2;
        if (freshwater)
            iremove_Mrows = union(iremove_Mrows,row);
        end
        
        % ph_tot and ph_free relationship
        row = row + 1;
        M(row, [ tp(j).iph_tot, tp(j).iph_free, tp(j).iphso4 ] ) = [ 1, -1, -1 ];
        mc = union(mc,[tp(j).iph_tot, tp(j).iph_free, tp(j).iphso4]);
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e8;

        % ph_sws and ph_free relationship
        row = row + 1;
        M(row,[ tp(j).iph_sws, tp(j).iph_free, tp(j).iphso4, tp(j).ipHF ]) = [ 1, -1, -1, -1 ];
        mc = union(mc,[tp(j).iph_sws, tp(j).iph_free, tp(j).iphso4, tp(j).ipHF]);
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e8;

        if optpKalpha == 1
            % Kalpha = [alpha][H]/[Halpha] unknown alkalinity TAlpha
            row = row + 1;
            M(row,[ ipTAlpha, tp(j).ipalpha, tp(j).iphalpha]) = [1, -1, -1];
            mc = union(mc, [ipTAlpha, tp(j).ipalpha, tp(j).iphalpha]);
            mr = [mr,row];
            M(row_alk, [tp(j).ipalpha,tp(j).iphalpha] ) = [-1, +1]; % or -, +?
            % rescale
            M(row,:) = M(row,:)*1e5; % 1e5 for 1umol order
            M(row_alk,:) = M(row_alk,:)*1e1; % 1e2 for TS 0.1 order
        end

        if optpKbeta == 1
            % Kbeta = [beta][H]/[Hbeta] unknown alkalinity TBeta
            row = row + 1;
            M(row,[ ipTBeta, tp(j).ipbeta, tp(j).iphbeta]) = [1, -1, -1];
            mc = union(mc, [ipTBeta, tp(j).ipbeta, tp(j).iphbeta]);
            mr = [mr,row];
            M(row_alk, [tp(j).ipbeta,tp(j).iphbeta] ) = [-1, +1]; 
            % rescale
            M(row,:) = M(row,:)*1e5; % 1e5 for 1umol order
            M(row_alk,:) = M(row_alk,:)*1e1; % 1e2 for TS 0.1 order
        end

        tp(j).mr = mr;
        tp(j).mc = mc.';
    end
    for j = 1:nTP
        % STUFF needed to compute the Revelle buffer factor
        ifixed = [ ipTA, ipTB, ipTS, ipTF, ipTP, ipTSi, ipTNH4, ipTH2S, ipTCa, isal ];
        ifixed = [ifixed, tp(j).ipK0,  tp(j).ipK1,  tp(j).ipK2,  tp(j).ipKb,  tp(j).ipKw, tp(j).ipKs,  tp(j).ipKf, ...
                 tp(j).ipKp1, tp(j).ipKp2, tp(j).ipKp3, tp(j).ipKsi, tp(j).ipKnh4, tp(j).ipKh2s,...
                 tp(j).ipp2f, tp(j).ipKar, tp(j).ipKca, tp(j).ipfH,  tp(j).iP, tp(j).iT];
        if optpKalpha == 1
            ifixed(end) = ipTAlpha;
            ifixed(end) = tp(j).ipKalpha;
        end
        if optpKbeta == 1
            ifixed(end) = ipTBeta;
            ifixed(end) = tp(j).ipKbeta;
        end
        tp(j).ifixed = ifixed;
        tp(j).ifree = setdiff(union(tp(j).mc,tp(j).kc),ifixed);

        ML = M(tp(j).mr,:); KL = K(tp(j).kr,:);
        ML = ML(:,tp(j).ifree); KL = KL(:,tp(j).ifree);
        
        %tp(j).c = @(xfree,zhat) [ML*q(xfree); KL*xfree] + [MR*q(zhat(ifixed)); KR*zhat(ifixed)];
        tp(j).dcdx_pTAfixed = @(xfree) [ML*diag(dqdx(xfree));KL];
                
        % STUFF needed to compute a dpfco2/dpTA holding pTC fixed 
        jfixed = [ ipTC, ipTB, ipTS, ipTF, ipTP, ipTSi, ipTNH4, ipTH2S, ipTCa, isal ];
        jfixed = [jfixed, tp(j).ipK0,  tp(j).ipK1,  tp(j).ipK2,  tp(j).ipKb,  tp(j).ipKw, tp(j).ipKs,  tp(j).ipKf, ...
                 tp(j).ipKp1, tp(j).ipKp2, tp(j).ipKp3, tp(j).ipKsi, tp(j).ipKnh4, tp(j).ipKh2s,...
                 tp(j).ipp2f, tp(j).ipKar, tp(j).ipKca, tp(j).ipfH,  tp(j).iP, tp(j).iT];
        if optpKalpha == 1
            jfixed(end) = ipTAlpha;
            jfixed(end) = tp(j).ipKalpha;
        end
        if optpKbeta == 1
            jfixed(end) = ipTBeta;
            jfixed(end) = tp(j).ipKbeta;
        end
        tp(j).jfixed = jfixed;
        tp(j).jfree = setdiff(union(tp(j).mc,tp(j).kc),jfixed);

        ML = M(tp(j).mr,:); KL = K(tp(j).kr,:);
        ML = ML(:,tp(j).jfree); KL = KL(:,tp(j).jfree);
        tp(j).dcdx_pTCfixed = @(xfree) [ML*diag(dqdx(xfree));KL];
    end
    sys.M  = M;
    sys.K  = K;
    sys.tp = tp;
    sys.freshwater = freshwater;
    if (freshwater)
        [nKr_,nKc_] = size(K);
        [nMr_,nMc_] = size(M);

        ic = setdiff(1:nMc_,iremove_vars);
        iMr = setdiff(1:nMr_,iremove_Mrows);
        M = M(iMr,:);
        M = M(:,ic);
        
        iKr = setdiff(1:nKr_,iremove_Krows);
        K = K(iKr,:);
        K = K(:,ic);        
    
        [nKr,nKc] = size(K);
        [nMr,nMc] = size(M);

        sys.iremove_vars = iremove_vars;
        sys.iremove_Krows = iremove_Krows;
        sys.iremove_Mrows = iremove_Mrows;
        sys.ic = ic;
        sys.iMr = iMr;
        sys.iKr = iKr;
        sys.nKr_ = nKr_; 
        sys.nKc_ = nKc_;
        sys.nMr_ = nMr_;
        sys.nMc_ = nMc_;
        sys.nKr = nKr;   
        sys.nKc = nKc;
        sys.nMr = nMr;
        sys.nMc = nMc;        
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
    % opt.tol
    if ~isfield(opt,'tol') || isbad(opt.tol)
        opt.tol = 1e-6; % default tolerance
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

% ----------------------------------------------------------------------
% parse_input
% ----------------------------------------------------------------------

function [obs,yobs,wobs,sys] = parse_input(obs,sys,opt,nD)
    % ORG ALK
    isgood  = @(thing) (~isempty(thing) & ~sum(isnan(thing)));
    iszero  = @(thing) (isequal(thing,0));
    p       = sys.p;
    q       = sys.q;
    % convert x+/-e into precision for p(x) (precision = 1/variance)
    w       = @(x,e) abs( p(1 + e./x) ).^(-2);

    nv      = size(sys.K,2);
    yobs    = nan(nD,nv);
    wobs    = nan(nD,nv);

    if (~isfield(obs,'tp'))
        error('Need to provide temperature and pressure measurement.')
    end
    
    nTP = length(obs(1).tp); 
    for i = 1:nD % loop over all the stations
        % make sure all the required fields in the obs struct exist
        if (~isfield(obs(i), 'sal'))
            error('Need to provide salinity measurement.');
        else
            yobs(i,sys.isal) = obs(i).sal;
        end
        if (~isfield(obs(i), 'usal')) || (~isgood(obs(i).usal)) || (iszero(obs(i).usal))
            obs(i).usal         = 0.002; % 1std = 0.002 PSU
            wobs(i,sys.isal)    = (obs(i).usal)^(-2);
            if opt.printmes ~= 0
                fprintf('Warning: Assuming salinity uncertainty is 0.002 PSU \n');
            end
        else
            wobs(i,sys.isal)    = (obs(i).usal)^(-2); % std u -> w
        end
        if (~isfield(obs(i),'TC')) || (~isgood(obs(i).TC))
            obs(i).TC        = nan; 
            yobs(i,sys.ipTC) = nan;
        else
            yobs(i,sys.ipTC) = p((obs(i).TC)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'uTC')) || (~isgood(obs(i).uTC)) || (iszero(obs(i).uTC))
            obs(i).uTC       = nan;
            wobs(i,sys.ipTC) = nan;
            if (isgood(obs(i).TC))
                if opt.printmes ~= 0
                fprintf('Warning: Assuming TC uncertainty is 2.00 umol/kg \n');
                end
                wobs(i,sys.ipTC) = w(obs(i).TC,2.00);
            end
        else
            wobs(i,sys.ipTC) = w(obs(i).TC,obs(i).uTC); % std u -> w
        end
        if(~isfield(obs(i),'TA'))  || (~isgood(obs(i).TA))
            obs(i).TA        = nan; %[]
            yobs(i,sys.ipTA) = nan;
        else
            yobs(i,sys.ipTA) = p((obs(i).TA)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'uTA'))  || (~isgood(obs(i).uTA)) || (iszero(obs(i).uTA))
            obs(i).uTA       = nan;
            wobs(i,sys.ipTA) = nan;
            if (isgood(obs(i).TA))
                if opt.printmes ~= 0
                fprintf('Warning: Assuming TA uncertainty is 2.00 umol/kg \n');
                end
                wobs(i,sys.ipTA) = w(obs(i).TA,2.00);
            end
        else
            wobs(i,sys.ipTA) = w(obs(i).TA,obs(i).uTA); % std u -> w
        end
        % calculate totals that are a function of salinity
        [pT,~,~,upT]  =  calc_pTOT(opt,obs(i).sal);
        % (see Ref's within calc_pTOT)
        pTB  = pT(1); TB  = q(pTB)*1e6;  upTB  = upT(1);
        pTS  = pT(2); TS  = q(pTS)*1e6;  upTS  = upT(2);
        pTF  = pT(3); TF  = q(pTF)*1e6;  upTF  = upT(3);
        pTCa = pT(4); TCa = q(pTCa)*1e6; upTCa = upT(4); 
        % total borate
        if (~isfield(obs(i),'TB')) || (~isgood(obs(i).TB))
            obs(i).TB = nan;
            yobs(i,sys.ipTB) = pTB;
        else
            if ((obs(i).TB) == 0)
                obs(i).TB   = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTB) = p(obs(i).TB*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'uTB')) || (~isgood(obs(i).uTB))
            obs(i).uTB       = nan;
            wobs(i,sys.ipTB) = (upTB)^(-2); % convert to precision
        else
            wobs(i,sys.ipTB) = w(obs(i).TB,obs(i).uTB); % µmol/kg
        end
        
        % total sulfate
        if (~isfield(obs(i), 'TS')) || (~isgood(obs(i).TS))
            obs(i).TS        = nan;
            yobs(i,sys.ipTS) = pTS;
        else
            if ((obs(i).TS) == 0)
                obs(i).TS = 1e-3; % µmol/kg reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTS)     = p(obs(i).TS*1e-6); % mol/kg
        end
        if (~isfield(obs(i), 'uTS')) || (~isgood(obs(i).uTS))
            obs(i).uTS       = nan ;
            wobs(i,sys.ipTS) = (upTS)^(-2); % convert to precision
        else
            wobs(i,sys.ipTS)    = w(obs(i).TS,obs(i).uTS);
        end

        % total fluoride
        if (~isfield(obs(i), 'TF')) || (~isgood(obs(i).TF))
            obs(i).TF        = nan; 
            yobs(i,sys.ipTF) = pTF;
        else
            if ((obs(i).TF) == 0)
                obs(i).TF       = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTF)    = p(obs(i).TF*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'uTF'))  || (~isgood(obs(i).uTF))
            obs(i).uTF       = nan;
            wobs(i,sys.ipTF) = (upTF)^(-2);
        else
            wobs(i,sys.ipTF)    = w(obs(i).TF,obs(i).uTF);
        end        
        
        % total phosphate
        if (~isfield(obs(i), 'TP'))  || (~isgood(obs(i).TP))
            obs(i).TP           = nan;
            yobs(i,sys.ipTP)    = p(1e-9); % resest minimum (mol/kg)
        else
            if ((obs(i).TP) == 0) % zero po4 is very unlikely and breaks the code
                obs(i).TP       = 1e-3; % umol/kg-SW, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTP)     = p(obs(i).TP*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'uTP'))  || (~isgood(obs(i).uTP))
            obs(i).uTP          = nan;
            wobs(i,sys.ipTP)    = w(1e-3,1e-3); 
        else
            if ((obs(i).uTP) == 0)
                obs(i).uTP      = 1e-3; % umol/kg, reset minimum if zero
            end
            wobs(i,sys.ipTP)    = w(obs(i).TP,obs(i).uTP); % mol/kg
        end
        
        % total silicate
        if (~isfield(obs(i), 'TSi'))  || (~isgood(obs(i).TSi))
            obs(i).TSi          = nan;
            yobs(i,sys.ipTSi)   = p(1e-9); % reset minimum (mol/kg)
        else
            if ((obs(i).TSi) == 0) % zero silicate very unlikely and breaks code
                obs(i).TSi      = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTSi)   = p(obs(i).TSi*1e-6);
        end
        if (~isfield(obs(i), 'uTSi'))  || (~isgood(obs(i).uTSi))
            obs(i).uTSi         = nan;            
            wobs(i,sys.ipTSi)   = w(1e-3,1e-3); % mol/kg
            if (isgood(obs(i).TSi)) && opt.printmes ~= 0
                fprintf('Warning, no obs.uTSi input with obs.uTSi. Assuming 1 nanomolar. \n' )
            end
        else
            if ((obs(i).uTSi) == 0)
                obs(i).uTSi     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            wobs(i,sys.ipTSi)   = w(obs(i).TSi,obs(i).uTSi); %  mol/kg
        end

        % total amonia
        if (~isfield(obs(i), 'TNH4'))  || (~isgood(obs(i).TNH4))
            obs(i).TNH4         = nan;
            yobs(i,sys.ipTNH4)  = p(1e-9); % reset minimum (mol/kg)
        else
            if ((obs(i).TNH4) == 0)
                obs(i).TNH4     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTNH4)  = p(obs(i).TNH4*1e-6);
        end
        if (~isfield(obs(i), 'uTNH4'))  || (~isgood(obs(i).uTNH4))
            uTNH4               = 5e-4; % µmol/kg
            wobs(i,sys.ipTNH4)  = w(1e-3,uTNH4); % mol/kg
            obs(i).uTNH4        = nan;
            if (isgood(obs(i).TNH4)) && opt.printmes ~= 0
                fprintf('Warning, no obs.uTNH4 input with obs.uTNH4. Assuming 5e-4 umol/kg. \n' )
            end
        else
            wobs(i,sys.ipTNH4)  = w(obs(i).TNH4,obs(i).uTNH4);
        end

        % total sulfide
        if (~isfield(obs(i), 'TH2S'))  || (~isgood(obs(i).TH2S))
            obs(i).TH2S         = nan; 
            yobs(i,sys.ipTH2S)  = p(1e-9); % reset minimum (mol/kg)
        else
            if ((obs(i).TH2S) == 0)
                obs(i).TH2S     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTH2S)  = p(obs(i).TH2S*1e-6);
        end
        if (~isfield(obs(i), 'uTH2S'))  || (~isgood(obs(i).uTH2S))
            uTH2S               = 5e-4; % µmol/kg
            wobs(i,sys.ipTH2S)  = w(1e-3,uTH2S); % mol/kg
            obs(i).uTH2S        = nan;
            if (isgood(obs(i).TH2S)) && opt.printmes ~= 0
                fprintf(' Warning, no obs.uTH2S input with obs.uTNH4. Assuming 5e-4 umol/kg. \n' )
            end
        else
            wobs(i,sys.ipTH2S)  = w(obs(i).TH2S,obs(i).uTH2S);
        end
            
        % total calcium
        if (~isfield(obs(i), 'TCa'))  || (~isgood(obs(i).TCa))
            obs(i).TCa        = nan;
            yobs(i,sys.ipTCa) = pTCa;
        else
            if ((obs(i).TCa) == 0)
                obs(i).TCa      = 1e-9; % mol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTCa)   = p(obs(i).TCa); % assume user input of mol/kg
        end
        if (~isfield(obs(i), 'uTCa'))  || (~isgood(obs(i).uTCa))
            obs(i).uTCa       = nan;
            wobs(i,sys.ipTCa) = (upTCa)^(-2);
            if (isgood(obs(i).TCa)) && opt.printmes ~= 0
                fprintf(' Warning, no obs.uTCa input with obs.uTCa. Assuming 6e-5 mol/kg. \n' )
            end
        else
            wobs(i,sys.ipTCa)   = w(obs(i).TCa,obs(i).uTCa);
        end

        if opt.pKalpha == 1
            % Organic Alkalinity Alpha
            if nTP > 1
                error(['opt.pKalpha only works with one tp, not 2+.'...
                        'Must input at tp(1) only.']);
            end
            if (~isfield(obs(i), 'TAlpha')) || (~isgood(obs(i).TAlpha))
                obs(i).TAlpha           = nan; 
                yobs(i,sys.ipTAlpha)    = nan;
            else
                if ((obs(i).TAlpha) == 0)
                    obs(i).TAlpha       = 1e-9;
                end
                yobs(i,sys.ipTAlpha)    = p(obs(i).TAlpha*1e-6);
            end
            if (~isfield(obs(i), 'uTAlpha')) || (~isgood(obs(i).uTAlpha))
                obs(i).uTAlpha          = nan; 
                wobs(i,sys.ipTAlpha)    = nan; 
            else
                wobs(i,sys.ipTAlpha)    = w(obs(i).TAlpha,obs(i).uTAlpha);
            end
        end

        if opt.pKbeta == 1
            % Organic Alkalinity Beta
            if nTP > 1
                error(['opt.pKbeta only works with one tp, not 2+.'...
                        'Must input at tp(1) only.']);
            end
            if (~isfield(obs(i),'TBeta')) || (~isgood(obs(i).TBeta))
                obs(i).TBeta        = nan;
                yobs(i,sys.ipTBeta) = nan;
            else
                if ((obs(i).TBeta) == 0)
                    obs(i).TBeta    = 1e-9;
                end
                yobs(i,sys.ipTBeta) = p(obs(i).TBeta*1e-6);
            end
            if (~isfield(obs(i),'uTBeta')) || (~isgood(obs(i).uTBeta))
                obs(i).uTBeta       = nan;
                wobs(i,sys.ipTBeta) = nan;
            else
                wobs(i,sys.ipTBeta) = w(obs(i).TBeta,obs(i).uTBeta);
            end
        end

        for j = 1:nTP % loop over (T,P) systems
            if (~isfield(obs(i).tp(j),'T')) || (~isgood(obs(i).tp(j).T))
                error('Must supply temperature at obs.tp(%d).T',j)
            end
            yobs(i,sys.tp(j).iT)   = obs(i).tp(j).T;
            if (~isfield(obs(i).tp(j),'uT')) || (~isgood(obs(i).tp(j).uT)) ...
                    || (iszero(obs(i).tp(j).uT))
                if opt.printmes ~= 0
                fprintf('Warning: Assuming obs.tp(%d).uT 0.1 deg Celsius \n',j);
                end
                obs(i).tp(j).uT = 0.1;
            end
            wobs(i,sys.tp(j).iT)   = (obs(i).tp(j).uT)^(-2);

            if (~isfield(obs(i).tp(j),'P')) || (~isgood(obs(i).tp(j).P))
                error('Must supply pressure at obs.tp(%d).P',j)
            end
            yobs(i,sys.tp(j).iP)   = obs(i).tp(j).P;
            if (~isfield(obs(i).tp(j),'uP')) || (~isgood(obs(i).tp(j).uP)) ...
                    || (iszero(obs(i).tp(j).uP))
                if opt.printmes ~= 0
                fprintf('Warning: Assuming obs.tp(%d).uP is 0.1 dbar \n',j);
                end
                obs(i).tp(j).uP = 0.1;
            end
            wobs(i,sys.tp(j).iP)   = (obs(i).tp(j).uP)^(-2);

            % calculate the equilibrium constants
            [pK,~,upK] = calc_pK(opt,obs(i).tp(j).T,obs(i).sal,obs(i).tp(j).P); % T, S, P

            pK0   = pK(1);      pK1  = pK(2);     pK2   = pK(3);  
            pKb   = pK(4);      pKw  = pK(5);     pKs   = pK(6);  
            pKf   = pK(7);      pKp1 = pK(8);     pKp2  = pK(9);  
            pKp3  = pK(10);     pKsi = pK(11);    pKnh4 = pK(12);
            pKh2s = pK(13);     pp2f = pK(14);    pKar  = pK(15); 
            pKca  = pK(16);     pfH  = pK(17);
            
            upK0   = upK(1);    upK1  = upK(2);   upK2   = upK(3);  
            upKb   = upK(4);    upKw  = upK(5);   upKs   = upK(6);  
            upKf   = upK(7);    upKp1 = upK(8);   upKp2  = upK(9);  
            upKp3  = upK(10);   upKsi = upK(11);  upKnh4 = upK(12);
            upKh2s = upK(13);   upp2f = upK(14);  upKar  = upK(15); 
            upKca  = upK(16);   upfH  = upK(17);
            
            % add "observations" for the equilibrium constants
            % and transfer from obs struct to yobs and wobs
            
            % co2 solubility and fugacity
            if (~isfield(obs(i).tp(j),'pK0')) || (~isgood(obs(i).tp(j).pK0))
                obs(i).tp(j).pK0        = nan;
                yobs(i,sys.tp(j).ipK0)  = pK0;
            else
                yobs(i,sys.tp(j).ipK0)  = obs(i).tp(j).pK0;
            end
            if (~isfield(obs(i).tp(j),'upK0')) || (~isgood(obs(i).tp(j).upK0))
                obs(i).tp(j).upK0       = nan;
                wobs(i,sys.tp(j).ipK0)  = (upK0)^(-2);
            else
                wobs(i,sys.tp(j).ipK0)  = (obs(i).tp(j).pK0)^(-2);
            end
            if (~isfield(obs(i).tp(j),'co2st')) || (~isgood(obs(i).tp(j).co2st))
                obs(i).tp(j).co2st          = nan;
                yobs(i,sys.tp(j).ipco2st)   = nan;
            else
                yobs(i,sys.tp(j).ipco2st)   = p(obs(i).tp(j).co2st*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uco2st')) || (~isgood(obs(i).tp(j).uco2st))
                obs(i).tp(j).uco2st         = nan;
                wobs(i,sys.tp(j).ipco2st)   = nan;
            else
                wobs(i,sys.tp(j).ipco2st) = w(obs(i).tp(j).co2st,obs(i).tp(j).uco2st);
            end
            if (~isfield(obs(i).tp(j),'fco2')) || (~isgood(obs(i).tp(j).fco2))
                obs(i).tp(j).fco2           = nan;
                yobs(i,sys.tp(j).ipfco2)    = nan;
            else
                yobs(i,sys.tp(j).ipfco2)    = p(obs(i).tp(j).fco2*1e-6); % convt µatm to atm
            end
            if (~isfield(obs(i).tp(j),'ufco2')) || (~isgood(obs(i).tp(j).ufco2)) ...
                    || (iszero(obs(i).tp(j).ufco2))
                obs(i).tp(j).ufco2          = nan;
                wobs(i,sys.tp(j).ipfco2)    = nan;
                if (isgood(obs(i).tp(j).fco2))
                    if opt.printmes ~= 0
                        fprintf('Warning: Assuming obs.tp(%d).ufco2 is 1%% \n',j);
                    end
                    wobs(i,sys.tp(j).ipfco2) = w(obs(i).tp(j).fco2,0.01*obs(i).tp(j).fco2);
                end
            else
                wobs(i,sys.tp(j).ifco2) = w(obs(i).tp(j).fco2, obs(i).tp(j).ufco2);
            end
            if (~isfield(obs(i).tp(j),'pp2f')) || (~isgood(obs(i).tp(j).pp2f))
                obs(i).tp(j).pp2f           = nan;
                yobs(i,sys.tp(j).ipp2f)     = pp2f;
            else
                yobs(i,sys.tp(j).ipp2f)     = obs(i).tp(j).pp2f;
            end
            if (~isfield(obs(i).tp(j),'upp2f')) || (~isgood(obs(i).tp(j).upp2f))
                obs(i).tp(j).upp2f         = nan;
                wobs(i,sys.tp(j).ipp2f)    = (upp2f)^(-2);
            else
                wobs(i,sys.tp(j).ipp2f) = (obs(i).tp(j).upp2f)^(-2);
            end
            if (~isfield(obs(i).tp(j),'pco2')) || (~isgood(obs(i).tp(j).pco2))
                obs(i).tp(j).pco2          = nan;
                yobs(i,sys.tp(j).ippco2)   = nan;
            else
                yobs(i,sys.tp(j).ippco2)   = p(obs(i).tp(j).pco2*1e-6); % convt µatm to atm
            end
            if (~isfield(obs(i).tp(j),'upco2')) || (~isgood(obs(i).tp(j).upco2)) ...
                    || (iszero(obs(i).tp(j).upco2))
                obs(i).tp(j).upco2         = nan;
                wobs(i,sys.tp(j).ippco2)   = nan;
                if (isgood(obs(i).tp(j).pco2))
                    if opt.printmes ~= 0
                        fprintf('Warning: Assuming obs.tp(%d).upco2 is 1%% \n',j);
                    end
                    wobs(i,sys.tp(j).ippco2) = w(obs(i).tp(j).pco2,0.01*obs(i).tp(j).pco2);
                end
            else
                wobs(i,sys.tp(j).ippco2) = w(obs(i).tp(j).pco2,obs(i).tp(j).upco2);
            end

            % carbonate system
            if (~isfield(obs(i).tp(j),'pK1')) || (~isgood(obs(i).tp(j).pK1))
                obs(i).tp(j).pK1       = nan;
                if opt.turnoff.pK1 == 1 % do not input pK1
                    if nTP > 1
                        error(['opt.turnoff.pK1 only works with one tp, not 2+. ' ...
                                'Must input at tp(1) only.']);
                    else
                        yobs(i,sys.tp(j).ipK1) = nan;
                    end
                else % do input pK1
                    yobs(i,sys.tp(j).ipK1) = pK1;
                end
            else
                yobs(i,sys.tp(j).ipK1) = obs(i).tp(j).pK1;
            end
            if (~isfield(obs(i).tp(j),'upK1')) || (~isgood(obs(i).tp(j).upK1))
                obs(i).tp(j).upK1      = nan;
                if opt.turnoff.pK1 == 1
                    wobs(i,sys.tp(j).ipK1) = nan;
                else
                    wobs(i,sys.tp(j).ipK1) = (upK1)^(-2);
                end
            else
                wobs(i,sys.tp(j).ipK1) = (obs(i).tp(j).upK1)^(-2);
            end
            if (~isfield(obs(i).tp(j),'pK2')) || (~isgood(obs(i).tp(j).pK2))
                obs(i).tp(j).pK2       = nan;
                if opt.turnoff.pK2 == 1
                    if nTP > 1
                        error(['opt.turnoff.pK2 only works with one tp, not 2+. ' ...
                                'Must input at tp(1) only.']);
                    else
                        yobs(i,sys.tp(j).ipK2) = nan;
                    end
                else
                    yobs(i,sys.tp(j).ipK2) = pK2;
                end
            else
                yobs(i,sys.tp(j).ipK2) = obs(i).tp(j).pK2;
            end
            if (~isfield(obs(i).tp(j),'upK2')) || (~isgood(obs(i).tp(j).upK2))
                obs(i).tp(j).upK2      = nan;
                if opt.turnoff.pK2 == 1
                    wobs(i,sys.tp(j).ipK2) = nan;
                else
                    wobs(i,sys.tp(j).ipK2) = (upK2)^(-2);
                end
            else
                wobs(i,sys.tp(j).ipK2) = (obs(i).tp(j).upK2)^(-2);
            end
            if (~isfield(obs(i).tp(j),'hco3')) || (~isgood(obs(i).tp(j).hco3))
                obs(i).tp(j).hco3          = nan;
                yobs(i,sys.tp(j).iphco3)   = nan;   
            else
                yobs(i,sys.tp(j).iphco3) = p(obs(i).tp(j).hco3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uhco3')) || (~isgood(obs(i).tp(j).uhco3))
                obs(i).tp(j).uhco3         = nan;
                wobs(i,sys.tp(j).iphco3)   = nan;
            else
                wobs(i,sys.tp(j).iphco3) = w(obs(i).tp(j).hco3,obs(i).tp(j).uhco3);
            end
            if (~isfield(obs(i).tp(j),'co3')) || (~isgood(obs(i).tp(j).co3))
                obs(i).tp(j).co3           = nan;
                yobs(i,sys.tp(j).ipco3)    = nan;
            else
                yobs(i,sys.tp(j).ipco3) = p(obs(i).tp(j).co3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uco3')) || (~isgood(obs(i).tp(j).uco3)) ...
                    || (iszero(obs(i).tp(j).uco3))
                obs(i).tp(j).uco3          = nan;
                wobs(i,sys.tp(j).ipco3)    = nan;
                if (isgood(obs(i).tp(j).co3))
                    if opt.printmes ~= 0
                        fprintf('Warning: Assuming obs.tp(%d).uco3 is 2%% \n',j);
                    end
                    wobs(i,sys.tp(j).ipco3) = w(obs(i).tp(j).co3,0.02*obs(i).tp(j).co3);
                end
            else
                wobs(i,sys.tp(j).ipco3) = w(obs(i).tp(j).co3,obs(i).tp(j).uco3);
            end
            if (~isfield(obs(i).tp(j),'ph')) || (~isgood(obs(i).tp(j).ph))
                obs(i).tp(j).ph        = nan;
                yobs(i,sys.tp(j).iph)  = nan;
            else
                yobs(i,sys.tp(j).iph)  = obs(i).tp(j).ph;
            end
            if (~isfield(obs(i).tp(j),'uph')) || (~isgood(obs(i).tp(j).uph)) ...
                    || (iszero(obs(i).tp(j).uph))
                obs(i).tp(j).uph       = nan;
                wobs(i,sys.tp(j).iph)  = nan;
                if (isgood(obs(i).tp(j).ph))
                    if opt.printmes ~= 0
                        fprintf('Warning: Assuming obs.tp(%d).uph is 0.010\n',j);
                    end
                    wobs(i,sys.tp(j).iph) = w(obs(i).tp(j).ph,0.010);
                end
            else
                wobs(i,sys.tp(j).iph)  = (obs(i).tp(j).uph).^(-2);
            end
            if (~isfield(obs(i).tp(j),'ph_tot')) || (~isgood(obs(i).tp(j).ph_tot))
                obs(i).tp(j).ph_tot       = nan;
                yobs(i,sys.tp(j).iph_tot) = nan;
            else
                yobs(i,sys.tp(j).iph_tot) = obs(i).tp(j).ph_tot;
            end
            if (~isfield(obs(i).tp(j),'uph_tot')) || (~isgood(obs(i).tp(j).uph_tot))
                obs(i).tp(j).uph_tot      = nan;
                wobs(i,sys.tp(j).iph_tot) = nan;
            else
                wobs(i,sys.tp(j).iph_tot) = obs(i).tp(j).uph_tot.^(-2);
            end
            if (~isfield(obs(i).tp(j),'ph_free')) || (~isgood(obs(i).tp(j).ph_free))
                obs(i).tp(j).ph_free       = nan;
                yobs(i,sys.tp(j).iph_free) = nan;
            else
                yobs(i,sys.tp(j).iph_free) = obs(i).tp(j).ph_free ;
            end
            if (~isfield(obs(i).tp(j),'uph_free')) || (~isgood(obs(i).tp(j).uph_free))
                obs(i).tp(j).uph_free      = nan;
                wobs(i,sys.tp(j).iph_free) = nan;
            else
                wobs(i,sys.tp(j).iph_free) = obs(i).tp(j).uph_free.^(-2);
            end
            if (~isfield(obs(i).tp(j),'ph_sws')) || (~isgood(obs(i).tp(j).ph_sws))
                obs(i).tp(j).ph_sws       = nan;
                yobs(i,sys.tp(j).iph_sws) = nan;
            else
                yobs(i,sys.tp(j).iph_sws) = obs(i).tp(j).ph_sws;
            end
            if (~isfield(obs(i).tp(j),'uph_sws')) || (~isgood(obs(i).tp(j).uph_sws))
                obs(i).tp(j).uph_sws      = nan;
                wobs(i,sys.tp(j).iph_sws) = nan;
            else
                wobs(i,sys.tp(j).iph_sws) = obs(i).tp(j).uph_sws.^(-2);
            end
            if (~isfield(obs(i).tp(j),'ph_nbs')) || (~isgood(obs(i).tp(j).ph_nbs))
                obs(i).tp(j).ph_nbs       = nan;
                yobs(i,sys.tp(j).iph_nbs) = nan;
            else
                yobs(i,sys.tp(j).iph_nbs) = obs(i).tp(j).ph_nbs;
            end
            if (~isfield(obs(i).tp(j),'uph_nbs')) || (~isgood(obs(i).tp(j).uph_nbs))
                obs(i).tp(j).uph_nbs      = nan;
                wobs(i,sys.tp(j).iph_nbs) = nan;
            else
                wobs(i,sys.tp(j).iph_nbs) = obs(i).tp(j).uph_nbs.^(-2);
            end
            if (~isfield(obs(i).tp(j),'pfH')) || (~isgood(obs(i).tp(j).pfH))
                obs(i).tp(j).pfH           = nan;
                yobs(i,sys.tp(j).ipfH)     = pfH;
            else
                yobs(i,sys.tp(j).ipfH)     = obs(i).tp(j).pfH;
            end
            if (~isfield(obs(i).tp(j),'upfH')) || (~isgood(obs(i).tp(j).upfH))
                obs(i).tp(j).upfH          = nan;
                wobs(i,sys.tp(j).ipfH)     = (upfH).^(-2);
            else
                wobs(i,sys.tp(j).ipfH)     = (obs(i).tp(j).upfH).^(-2) ;
            end
            
            % water dissociation
            if (~isfield(obs(i).tp(j),'pKw')) || (~isgood(obs(i).tp(j).pKw))
                obs(i).tp(j).pKw       = nan;
                yobs(i,sys.tp(j).ipKw) = pKw;
            else
                yobs(i,sys.tp(j).ipKw) = obs(i).tp(j).pKw;
            end
            if (~isfield(obs(i).tp(j),'upKw')) || (~isgood(obs(i).tp(j).upKw))
                obs(i).tp(j).upKw      = nan;
                wobs(i,sys.tp(j).ipKw) = (upKw).^(-2);
            else
                wobs(i,sys.tp(j).ipKw) = (obs(i).tp(j).upKw).^(-2);
            end
            if (~isfield(obs(i).tp(j),'oh')) || (~isgood(obs(i).tp(j).oh))
                obs(i).tp(j).oh        = nan;
                yobs(i,sys.tp(j).ipoh) = nan;
            else
                yobs(i,sys.tp(j).ipoh) = p(obs(i).tp(j).oh*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uoh')) || (~isgood(obs(i).tp(j).uoh))
                obs(i).tp(j).uoh       = nan;
                wobs(i,sys.tp(j).ipoh) = nan;
            else
                wobs(i,sys.tp(j).ipoh) = w(obs(i).tp(j).oh,obs(i).tp(j).uoh);
            end

            % borate system 
            if (~isfield(obs(i).tp(j),'pKb')) || (~isgood(obs(i).tp(j).pKb))
                obs(i).tp(j).pKb       = nan;
                yobs(i,sys.tp(j).ipKb) = pKb; 
            else
                yobs(i,sys.tp(j).ipKb) = obs(i).tp(j).pKb;
            end
            if (~isfield(obs(i).tp(j),'upKb')) || (~isgood(obs(i).tp(j).upKb))
                obs(i).tp(j).upKb      = nan;
                wobs(i,sys.tp(j).ipKb) = (upKb).^(-2);
            else
                wobs(i,sys.tp(j).ipKb) = (obs(i).tp(j).upKb).^(-2);
            end
            if (~isfield(obs(i).tp(j),'boh3')) || (~isgood(obs(i).tp(j).boh3))
                obs(i).tp(j).boh3          = nan;
                yobs(i,sys.tp(j).ipboh3)   = nan;
            else
                yobs(i,sys.tp(j).ipboh3) = p(obs(i).tp(j).boh3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uboh3')) || (~isgood(obs(i).tp(j).uboh3))
                obs(i).tp(j).uboh3         = nan;
                wobs(i,sys.tp(j).ipboh3)   = nan;
            else
                wobs(i,sys.tp(j).ipboh3) = w(obs(i).tp(j).boh3,obs(i).tp(j).uboh3);
            end
            if (~isfield(obs(i).tp(j),'boh4')) || (~isgood(obs(i).tp(j).boh4))
                obs(i).tp(j).boh4          = nan;
                yobs(i,sys.tp(j).ipboh4)   = nan;
            else
                yobs(i,sys.tp(j).ipboh4) = p(obs(i).tp(j).boh4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uboh4')) || (~isgood(obs(i).tp(j).uboh4))
                obs(i).tp(j).uboh4         = nan;
                wobs(i,sys.tp(j).ipboh4)   = nan;
            else
                wobs(i,sys.tp(j).ipboh4) = w(obs(i).tp(j).boh4,obs(i).tp(j).uboh4);
            end
            
            % sulfate system
            if (~isfield(obs(i).tp(j),'pKs')) || (~isgood(obs(i).tp(j).pKs))
                obs(i).tp(j).pKs       = nan;
                yobs(i,sys.tp(j).ipKs) = pKs;
            else
                yobs(i,sys.tp(j).ipKs) = obs(i).tp(j).pKs;
            end
            if (~isfield(obs(i).tp(j),'upKs')) || (~isgood(obs(i).tp(j).upKs))
                obs(i).tp(j).upKs      = nan;
                wobs(i,sys.tp(j).ipKs) = (upKs).^(-2);
            else
                wobs(i,sys.tp(j).ipKs) = (obs(i).tp(j).upKs).^(-2);
            end
            if (~isfield(obs(i).tp(j),'hso4')) || (~isgood(obs(i).tp(j).hso4))
                obs(i).tp(j).hso4          = nan;
                yobs(i,sys.tp(j).iphso4)   = nan;
            else
                yobs(i,sys.tp(j).iphso4) = p(obs(i).tp(j).hso4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uhso4')) || (~isgood(obs(i).tp(j).uhso4))
                obs(i).tp(j).uhso4         = nan;
                wobs(i,sys.tp(j).iphso4)   = nan;
            else
                wobs(i,sys.tp(j).iphso4) = w(obs(i).tp(j).hso4,obs(i).tp(j).uhso4);
            end
            if (~isfield(obs(i).tp(j),'so4')) || (~isgood(obs(i).tp(j).so4))
                obs(i).tp(j).so4           = nan;
                yobs(i,sys.tp(j).ipso4)    = nan;
            else
                yobs(i,sys.tp(j).ipso4) = p(obs(i).tp(j).so4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uso4')) || (~isgood(obs(i).tp(j).uso4))
                obs(i).tp(j).uso4          = nan;
                wobs(i,sys.tp(j).ipso4)    = nan;
            else
                wobs(i,sys.tp(j).ipso4) = w(obs(i).tp(j).so4,obs(i).tp(j).uso4);
            end

            % fluoride system
            if (~isfield(obs(i).tp(j),'pKf')) || (~isgood(obs(i).tp(j).pKf))
                obs(i).tp(j).pKf       = nan;
                yobs(i,sys.tp(j).ipKf) = pKf;
            else
                yobs(i,sys.tp(j).ipKf) = obs(i).tp(j).pKf;
            end
            if (~isfield(obs(i).tp(j),'upKf')) || (~isgood(obs(i).tp(j).upKf))
                obs(i).tp(j).upKf      = nan;
                wobs(i,sys.tp(j).ipKf) = (upKf).^(-2);
            else
                wobs(i,sys.tp(j).ipKf) = (obs(i).tp(j).upKf).^(-2);
            end
            if (~isfield(obs(i).tp(j),'HF')) || (~isgood(obs(i).tp(j).HF))
                obs(i).tp(j).HF        = nan;
                yobs(i,sys.tp(j).ipHF) = nan;
            else
                yobs(i,sys.tp(j).ipHF) = p(obs(i).tp(j).HF*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uHF')) || (~isgood(obs(i).tp(j).uHF))
                obs(i).tp(j).uHF       = nan;
                wobs(i,sys.tp(j).ipHF) = nan;
            else
                wobs(i,sys.tp(j).ipHF) = w(obs(i).tp(j).HF,obs(i).tp(j).uHF);
            end
            if (~isfield(obs(i).tp(j),'F')) || (~isgood(obs(i).tp(j).F))
                obs(i).tp(j).F         = nan;
                yobs(i,sys.tp(j).ipF)  = nan;
            else
                yobs(i,sys.tp(j).ipF) = p(obs(i).tp(j).F*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uF')) || (~isgood(obs(i).tp(j).uF))
                obs(i).tp(j).uF        = nan;
                wobs(i,sys.tp(j).ipF)  = nan;
            else
                wobs(i,sys.tp(j).ipF) = w(obs(i).tp(j).F,obs(i).tp(j).uF);
            end
            
            % phosphate system
            if (~isfield(obs(i).tp(j),'pKp1')) || (~isgood(obs(i).tp(j).pKp1))
                obs(i).tp(j).pKp1          = nan;
                yobs(i,sys.tp(j).ipKp1)    = pKp1;
            else
               yobs(i,sys.tp(j).ipKp1)    = obs(i).tp(j).pKp1; 
            end
            if (~isfield(obs(i).tp(j),'upKp1')) || (~isgood(obs(i).tp(j).upKp1))
                obs(i).tp(j).upKp1         = nan;
                wobs(i,sys.tp(j).ipKp1)    = (upKp1).^(-2);
            else
                wobs(i,sys.tp(j).ipKp1) = (obs(i).tp(j).upKp1).^(-2);
            end
            if (~isfield(obs(i).tp(j),'pKp2')) || (~isgood(obs(i).tp(j).pKp2))
               obs(i).tp(j).pKp2          = nan;
                yobs(i,sys.tp(j).ipKp2)    = pKp2; 
            else
                yobs(i,sys.tp(j).ipKp2)    = obs(i).tp(j).pKp2;
            end
            if (~isfield(obs(i).tp(j),'upKp2')) || (~isgood(obs(i).tp(j).upKp2))
                obs(i).tp(j).upKp2         = nan;
                wobs(i,sys.tp(j).ipKp2)    = (upKp2).^(-2);
            else
                wobs(i,sys.tp(j).ipKp2) = (obs(i).tp(j).upKp2).^(-2);
            end
            if (~isfield(obs(i).tp(j),'pKp3')) || (~isgood(obs(i).tp(j).pKp3))
                obs(i).tp(j).pKp3          = nan;
                yobs(i,sys.tp(j).ipKp3)    = pKp3;
            else
                yobs(i,sys.tp(j).ipKp3)    = obs(i).tp(j).pKp3;
            end
            if (~isfield(obs(i).tp(j),'upKp3')) || (~isgood(obs(i).tp(j).upKp3))
                obs(i).tp(j).upKp3         = nan;
                wobs(i,sys.tp(j).ipKp3)    = (upKp3).^(-2);
            else
                wobs(i,sys.tp(j).ipKp3) = (obs(i).tp(j).upKp3).^(-2);
            end        
            if (~isfield(obs(i).tp(j),'h3po4')) || (~isgood(obs(i).tp(j).h3po4))
                obs(i).tp(j).h3po4         = nan;
                yobs(i,sys.tp(j).iph3po4)  = nan;
            else
                yobs(i,sys.tp(j).iph3po4) = p(obs(i).tp(j).h3po4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uh3po4')) || (~isgood(obs(i).tp(j).uh3po4))
                obs(i).tp(j).uh3po4        = nan;
                wobs(i,sys.tp(j).iph3po4)  = nan;
            else
                wobs(i,sys.tp(j).iph3po4) = w(obs(i).tp(j).h3po4,obs(i).tp(j).uh3po4);
            end
            if (~isfield(obs(i).tp(j),'h2po4')) || (~isgood(obs(i).tp(j).h2po4))
                obs(i).tp(j).h2po4         = nan;
                yobs(i,sys.tp(j).iph2po4)  = nan;
            else
                yobs(i,sys.tp(j).iph2po4) = p(obs(i).tp(j).h2po4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uh2po4')) || (~isgood(obs(i).tp(j).uh2po4))
                obs(i).tp(j).uh2po4        = nan;
                wobs(i,sys.tp(j).iph2po4)  = nan;
            else
                wobs(i,sys.tp(j).iph2po4) = w(obs(i).tp(j).h2po4,obs(i).tp(j).uh2po4);
            end
            if (~isfield(obs(i).tp(j),'hpo4')) || (~isgood(obs(i).tp(j).hpo4))
                obs(i).tp(j).hpo4          = nan;
                yobs(i,sys.tp(j).iphpo4)   = nan;
            else
               yobs(i,sys.tp(j).iphpo4) = p(obs(i).tp(j).hpo4*1e-6); % convt µmol/kg to mol/kg 
            end
            if (~isfield(obs(i).tp(j),'uhpo4')) || (~isgood(obs(i).tp(j).uhpo4))
                obs(i).tp(j).uhpo4         = nan;
                wobs(i,sys.tp(j).iphpo4)   = nan;
            else
                wobs(i,sys.tp(j).iphpo4) = w(obs(i).tp(j).hpo4,obs(i).tp(j).uhpo4);
            end
            if (~isfield(obs(i).tp(j),'po4')) || (~isgood(obs(i).tp(j).po4))
                obs(i).tp(j).po4           = nan;
                yobs(i,sys.tp(j).ippo4)    = nan;
            else
                yobs(i,sys.tp(j).ippo4) = p(obs(i).tp(j).po4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'upo4')) || (~isgood(obs(i).tp(j).upo4))
                obs(i).tp(j).upo4          = nan;
                wobs(i,sys.tp(j).ippo4)    = nan;
            else
                wobs(i,sys.tp(j).ippo4) = w(obs(i).tp(j).po4,obs(i).tp(j).upo4);
            end

            % silicate system
            if (~isfield(obs(i).tp(j),'pKsi')) || (~isgood(obs(i).tp(j).pKsi))
                obs(i).tp(j).pKsi          = nan;
                yobs(i,sys.tp(j).ipKsi)    = pKsi;
            else
                yobs(i,sys.tp(j).ipKsi)    = obs(i).tp(j).pKsi;
            end
            if (~isfield(obs(i).tp(j),'upKsi')) || (~isgood(obs(i).tp(j).upKsi))
                obs(i).tp(j).upKsi         = nan;
                wobs(i,sys.tp(j).ipKsi)    = (upKsi).^(-2); 
            else
                wobs(i,sys.tp(j).ipKsi) = (obs(i).tp(j).upKsi).^(-2);
            end
            if (~isfield(obs(i).tp(j),'siooh3')) || (~isgood(obs(i).tp(j).siooh3))
                obs(i).tp(j).siooh3        = nan;
                yobs(i,sys.tp(j).ipsiooh3) = nan;
            else
                yobs(i,sys.tp(j).ipsiooh3) = p(obs(i).tp(j).siooh3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'usiooh3')) || (~isgood(obs(i).tp(j).usiooh3))
                obs(i).tp(j).usiooh3       = nan;
                wobs(i,sys.tp(j).ipsiooh3) = nan;
            else
                wobs(i,sys.tp(j).ipsiooh3) = w(obs(i).tp(j).siooh3,obs(i).tp(j).usiooh3);
            end
            if (~isfield(obs(i).tp(j),'sioh4')) || (~isgood(obs(i).tp(j).sioh4))
                obs(i).tp(j).sioh4         = nan;
                yobs(i,sys.tp(j).ipsioh4)  = nan;
            else
                yobs(i,sys.tp(j).ipsioh4) = p(obs(i).tp(j).sioh4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'usioh4')) || (~isgood(obs(i).tp(j).usioh4))
                obs(i).tp(j).usioh4        = nan;
                wobs(i,sys.tp(j).ipsioh4)  = nan;
            else
                wobs(i,sys.tp(j).ipsioh4) = w(obs(i).tp(j).sioh4,obs(i).tp(j).usioh4);
            end

            % ammonia system
            if (~isfield(obs(i).tp(j),'pKnh4')) || (~isgood(obs(i).tp(j).pKnh4))
                obs(i).tp(j).pKnh4         = nan;
                yobs(i,sys.tp(j).ipKnh4)   = pKnh4;
            else
                yobs(i,sys.tp(j).ipKnh4)   = obs(i).tp(j).pKnh4;
            end
            if (~isfield(obs(i).tp(j),'upKnh4')) || (~isgood(obs(i).tp(j).upKnh4))
                obs(i).tp(j).upKnh4        = nan;
                wobs(i,sys.tp(j).ipKnh4)   = (upKnh4).^(-2);
            else
                wobs(i,sys.tp(j).ipKnh4)   = (obs(i).tp(j).upKnh4).^(-2);
            end
            if (~isfield(obs(i).tp(j),'nh4')) || (~isgood(obs(i).tp(j).nh4))
                obs(i).tp(j).nh4           = nan;
                yobs(i,sys.tp(j).ipnh4)    = nan;
            else
                yobs(i,sys.tp(j).ipnh4)    = p(obs(i).tp(j).nh4*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'unh4')) || (~isgood(obs(i).tp(j).unh4))
                obs(i).tp(j).unh4          = nan;
                wobs(i,sys.tp(j).ipnh4)    = nan;
            else
                wobs(i,sys.tp(j).ipnh4) = w(obs(i).tp(j).nh4,obs(i).tp(j).unh4);
            end
            if (~isfield(obs(i).tp(j),'nh3')) || (~isgood(obs(i).tp(j).nh3))
                obs(i).tp(j).nh3           = nan;
                yobs(i,sys.tp(j).ipnh3)    = nan;
            else
                yobs(i,sys.tp(j).ipnh3)    = p(obs(i).tp(j).nh3*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'unh3')) || (~isgood(obs(i).tp(j).unh3))
                obs(i).tp(j).unh3          = nan;
                wobs(i,sys.tp(j).ipnh3)    = nan;
            else
                wobs(i,sys.tp(j).ipnh3) = w(obs(i).tp(j).nh3,obs(i).tp(j).unh3);
            end

            % sulfide system
            if (~isfield(obs(i).tp(j),'pKh2s')) || (~isgood(obs(i).tp(j).pKh2s))
                obs(i).tp(j).pKh2s         = nan;
                yobs(i,sys.tp(j).ipKh2s)   = pKh2s;
            else
                yobs(i,sys.tp(j).ipKh2s)   = obs(i).tp(j).pKh2s;
            end
            if (~isfield(obs(i).tp(j),'upKh2s')) || (~isgood(obs(i).tp(j).upKh2s))
                obs(i).tp(j).upKh2s        = nan;
                wobs(i,sys.tp(j).ipKh2s)   = (upKh2s).^(-2);
            else
                wobs(i,sys.tp(j).ipKh2s)   = (obs(i).tp(j).upKh2s).^(-2);
            end
            if (~isfield(obs(i).tp(j),'H2S')) || (~isgood(obs(i).tp(j).H2S))
                obs(i).tp(j).H2S           = nan;
                yobs(i,sys.tp(j).ipH2S)    = nan;
            else
                yobs(i,sys.tp(j).ipH2S) = p(obs(i).tp(j).H2S*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uH2S')) || (~isgood(obs(i).tp(j).uH2S))
                obs(i).tp(j).uH2S          = nan;
                wobs(i,sys.tp(j).ipH2S)    = nan;
            else
                wobs(i,sys.tp(j).ipH2S) = w(obs(i).tp(j).H2S,obs(i).tp(j).uH2S);
            end            
            if (~isfield(obs(i).tp(j),'HS')) || (~isgood(obs(i).tp(j).HS))
                obs(i).tp(j).HS            = nan;
                yobs(i,sys.tp(j).ipHS)     = nan;
            else
                yobs(i,sys.tp(j).ipHS)     = p(obs(i).tp(j).HS*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uHS')) || (~isgood(obs(i).tp(j).uHS))
                obs(i).tp(j).uHS           = nan;
                wobs(i,sys.tp(j).ipHS)     = nan;
            else
                wobs(i,sys.tp(j).ipHS) = w(obs(i).tp(j).HS,obs(i).tp(j).uHS);
            end
            
            % calcium carbonate solubility system
            if (~isfield(obs(i).tp(j),'pKar')) || (~isgood(obs(i).tp(j).pKar))
                obs(i).tp(j).pKar          = nan;
                yobs(i,sys.tp(j).ipKar)    = pKar;
            else
                yobs(i,sys.tp(j).ipKar)    = obs(i).tp(j).pKar;
            end
            if (~isfield(obs(i).tp(j),'upKar')) || (~isgood(obs(i).tp(j).upKar))
                obs(i).tp(j).upKar         = nan;
                wobs(i,sys.tp(j).ipKar)    = (upKar).^(-2);
            else
                wobs(i,sys.tp(j).ipKar)    = (obs(i).tp(j).upKar).^(-2);
            end
            if (~isfield(obs(i).tp(j),'pKca')) || (~isgood(obs(i).tp(j).pKca))
                obs(i).tp(j).pKca          = nan;
                yobs(i,sys.tp(j).ipKca)    = pKca;
            else
                yobs(i,sys.tp(j).ipKca)    = obs(i).tp(j).pKca;
            end
            if (~isfield(obs(i).tp(j),'upKca')) || (~isgood(obs(i).tp(j).upKca))
                obs(i).tp(j).upKca         = nan;
                wobs(i,sys.tp(j).ipKca)    = (upKca).^(-2);
            else
                wobs(i,sys.tp(j).ipKca)    = (obs(i).tp(j).upKca).^(-2);
            end
            if (~isfield(obs(i).tp(j),'OmegaAr')) || (~isgood(obs(i).tp(j).OmegaAr))
                obs(i).tp(j).OmegaAr           = nan;
                yobs(i,sys.tp(j).ipOmegaAr)    = nan;
            else
                yobs(i,sys.tp(j).ipOmegaAr) = p(obs(i).tp(j).OmegaAr); % Omega is dimensionless
            end
            if (~isfield(obs(i).tp(j),'uOmegaAr')) || (~isgood(obs(i).tp(j).uOmegaAr))
                obs(i).tp(j).uOmegaAr          = nan;
                wobs(i,sys.tp(j).ipOmegaAr)    = nan;
            else
                wobs(i,sys.tp(j).ipOmegaAr) = w(obs(i).tp(j).OmegaAr,obs(i).tp(j).uOmegaAr);
            end
            if (~isfield(obs(i).tp(j),'OmegaCa')) || (~isgood(obs(i).tp(j).OmegaCa))
                obs(i).tp(j).OmegaCa           = nan;
                yobs(i,sys.tp(j).ipOmegaCa)    = nan;
            else
                yobs(i,sys.tp(j).ipOmegaCa) = p(obs(i).tp(j).OmegaCa); % Omega is dimensionless
            end
            if (~isfield(obs(i).tp(j),'uOmegaCa')) || (~isgood(obs(i).tp(j).uOmegaCa))
                obs(i).tp(j).uOmegaCa          = nan;
                wobs(i,sys.tp(j).ipOmegaCa)    = nan;
            else
                wobs(i,sys.tp(j).ipOmegaCa) = w(obs(i).tp(j).OmegaCa,obs(i).tp(j).uOmegaCa);
            end
            if (~isfield(obs(i).tp(j),'ca')) || (~isgood(obs(i).tp(j).ca))
               obs(i).tp(j).ca        = nan;
                yobs(i,sys.tp(j).ipca) = nan; 
            else
                yobs(i,sys.tp(j).ipca) = p(obs(i).tp(j).ca*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i).tp(j),'uca')) || (~isgood(obs(i).tp(j).uca))
                obs(i).tp(j).uca       = nan;
                wobs(i,sys.tp(j).ipca) = nan;
            else
                wobs(i,sys.tp(j).ipca)  = w(obs(i).tp(j).ca,obs(i).tp(j).uca);
            end

            if opt.pKalpha == 1
                % pKalpha system
                if (~isfield(obs(i).tp(j),'pKalpha')) || (~isgood(obs(i).tp(j).pKalpha))
                    obs(i).tp(j).pKalpha        = nan;
                    yobs(i,sys.tp(j).ipKalpha)  = nan; 
                else
                    pKalpha = obs(i).tp(j).pKalpha;
                    yobs(i,sys.tp(j).ipKalpha)  = obs(i).tp(j).pKalpha;
                end
                if (~isfield(obs(i).tp(j),'upKalpha')) || (~isgood(obs(i).tp(j).upKalpha))
                    obs(i).tp(j).upKalpha       = nan;
                    wobs(i,sys.tp(j).ipKalpha)  = nan;
                else
                    wobs(i,sys.tp(j).ipKalpha)  = (obs(i).tp(j).upKalpha)^(-2);
                end
                if (pKalpha > 4.5) % > from Kerr, < from Humphreys
                    sys.M(sys.tp(j).row_alk, sys.tp(j).ipalpha) = 0;
                else
                    sys.M(sys.tp(j).row_alk, sys.tp(j).iphalpha) = 0;
                end
                % other unknowns in system
                obs(i).tp(j).palpha         = nan; % p(alpha)
                yobs(i,sys.tp(j).ipalpha)   = nan;
                obs(i).tp(j).upalpha        = nan;
                wobs(i,sys.tp(j).ipalpha)   = nan;
                obs(i).tp(j).phalpha        = nan; % p(H-alpha)
                yobs(i,sys.tp(j).iphalpha)  = nan;
                obs(i).tp(j).uphalpha       = nan;
                wobs(i,sys.tp(j).iphalpha)  = nan;
            end

            if opt.pKbeta == 1
                % pKbeta system
                if (~isfield(obs(i).tp(j),'pKbeta')) || (~isgood(obs(i).tp(j).pKbeta))
                    obs(i).tp(j).pKbeta         = nan;
                    yobs(i,sys.tp(j).ipKbeta)   = nan; 
                else
                    pKbeta                      = obs(i).tp(j).pKbeta;
                    yobs(i,sys.tp(j).ipKbeta)   = obs(i).tp(j).pKbeta;
                end
                if (~isfield(obs(i).tp(j),'upKbeta')) || (~isgood(obs(i).tp(j).upKbeta))
                    obs(i).tp(j).upKbeta        = nan;
                    wobs(i,sys.tp(j).ipKbeta)   = nan;
                else
                    wobs(i,sys.tp(j).ipKbeta)   = (obs(i).tp(j).upKbeta)^(-2);
                end
                if (pKbeta > 4.5)
                    sys.M(sys.tp(j).row_alk, sys.tp(j).ipbeta) = 0;
                else
                    sys.M(sys.tp(j).row_alk, sys.tp(j).iphbeta) = 0;
                end
                % other unknowns in system
                obs(i).tp(j).pbeta          = nan; % p(beta)
                yobs(i,sys.tp(j).ipbeta)    = nan;
                obs(i).tp(j).upbeta         = nan;
                wobs(i,sys.tp(j).ipbeta)    = nan;
                obs(i).tp(j).phbeta         = nan; % p(H-beta)
                yobs(i,sys.tp(j).iphbeta)   = nan;
                obs(i).tp(j).uphbeta        = nan;
                wobs(i,sys.tp(j).iphbeta)   = nan;
            end

        end % for j = 1:nTP

    end
end

% ----------------------------------------------------------------------
% update_y
% ----------------------------------------------------------------------

function [y,gy,ggy] = update_y(y,x,obs,sys,opt)
    nTP         = length(sys.tp);
    M           = sys.M;
    K           = sys.K;
    nv          = size(M,2);
    gy          = zeros(nv,length(y));
    ggy         = zeros(nv,length(y),length(y));
    sal         = x(sys.isal);
    e3          = eps^3;
    ie3         = sqrt(-1)*e3;
    % calculate totals (see Ref's in calc_pTOT)
    [pT,gpT,ggpT,~] = calc_pTOT(opt,sal);
    pTB  = pT(1);  gpTB = gpT(1);  ggpTB = ggpT(1);
    pTS  = pT(2);  gpTS = gpT(2);  ggpTS = ggpT(2);
    pTF  = pT(3);  gpTF = gpT(3);  ggpTF = ggpT(3);
    pTCa = pT(4); gpTCa = gpT(4); ggpTCa = ggpT(4);
    % update y
    if (isnan(obs.TB))
        if opt.turnoff.TB ~= 1
            y(sys.ipTB)                     = pTB;
            gy(sys.ipTB,sys.isal)           = gpTB;
            ggy(sys.ipTB,sys.isal,sys.isal) = ggpTB;
        end
    end
    if (isnan(obs.TS))
        y(sys.ipTS)                     = pTS;
        gy(sys.ipTS,sys.isal)           = gpTS;
        ggy(sys.ipTS,sys.isal,sys.isal) = ggpTS;
    end
    if (isnan(obs.TF))
        y(sys.ipTF)                     = pTF;
        gy(sys.ipTF,sys.isal)           = gpTF;
        ggy(sys.ipTF,sys.isal,sys.isal) = ggpTF;
    end
    if (isnan(obs.TCa))
        y(sys.ipTCa)                        = pTCa;
        gy(sys.ipTCa,sys.isal)              = gpTCa;
        ggy(sys.ipTCa,sys.isal,sys.isal)    = ggpTCa;
    end

    for i = 1:nTP
        % use complex step method to get ∂T, ∂S, ∂P
        [pK,gpK] = calc_pK(opt,x(sys.tp(i).iT), x(sys.isal), x(sys.tp(i).iP));
        
        [pK_T,gpK_T] = calc_pK(opt, x(sys.tp(i).iT) + ie3  , x(sys.isal)       , x(sys.tp(i).iP));
        [pK_S,gpK_S] = calc_pK(opt, x(sys.tp(i).iT)        , x(sys.isal) + ie3 , x(sys.tp(i).iP));
        [pK_P,gpK_P] = calc_pK(opt, x(sys.tp(i).iT)        , x(sys.isal)       , x(sys.tp(i).iP) + ie3);
        %
        ggpK = zeros(length(pK),3,3);
        ggpK(:,1,:) = imag(gpK_T)/e3;
        ggpK(:,2,:) = imag(gpK_S)/e3;
        ggpK(:,3,:) = imag(gpK_P)/e3;
        
        iTSP        = [ sys.tp(i).iT, sys.isal, sys.tp(i).iP];
        if (isnan(obs.tp(i).pK0))
            y(sys.tp(i).ipK0)             = pK(1); 
            gy(sys.tp(i).ipK0,iTSP)       = gpK(1,:);
            ggy(sys.tp(i).ipK0,iTSP,iTSP) = ggpK(1,:,:);
        end
        if (isnan(obs.tp(i).pK1))
            if opt.turnoff.pK1 ~= 1
                y(sys.tp(i).ipK1)             = pK(2);
                gy(sys.tp(i).ipK1,iTSP)       = gpK(2,:);
                ggy(sys.tp(i).ipK1,iTSP,iTSP) = ggpK(2,:,:);
            end
        end
        if (isnan(obs.tp(i).pK2))
            if opt.turnoff.pK2 ~=1
                y(sys.tp(i).ipK2)             = pK(3);
                gy(sys.tp(i).ipK2,iTSP)       = gpK(3,:);
                ggy(sys.tp(i).ipK2,iTSP,iTSP) = ggpK(3,:,:);
            end
        end   
        if (isnan(obs.tp(i).pKb))
            y(sys.tp(i).ipKb)             = pK(4);
            gy(sys.tp(i).ipKb,iTSP)       = gpK(4,:);
            ggy(sys.tp(i).ipKb,iTSP,iTSP) = ggpK(4,:,:);
        end
        if (isnan(obs.tp(i).pKw))
            y(sys.tp(i).ipKw)             = pK(5);
            gy(sys.tp(i).ipKw,iTSP)       = gpK(5,:);
            ggy(sys.tp(i).ipKw,iTSP,iTSP) = ggpK(5,:,:);
        end
        if (isnan(obs.tp(i).pKs))
            y(sys.tp(i).ipKs)             = pK(6);
            gy(sys.tp(i).ipKs,iTSP)       = gpK(6,:);
            ggy(sys.tp(i).ipKs,iTSP,iTSP) = ggpK(6,:,:);
        end
        if (isnan(obs.tp(i).pKf))
            y(sys.tp(i).ipKf)             = pK(7);
            gy(sys.tp(i).ipKf,iTSP)       = gpK(7,:);
            ggy(sys.tp(i).ipKf,iTSP,iTSP) = ggpK(7,:,:);
        end
        if (isnan(obs.tp(i).pKp1))
            y(sys.tp(i).ipKp1)             = pK(8); 
            gy(sys.tp(i).ipKp1,iTSP)       = gpK(8,:);
            ggy(sys.tp(i).ipKp1,iTSP,iTSP) = ggpK(8,:,:);
        end
        if (isnan(obs.tp(i).pKp2))
            y(sys.tp(i).ipKp2)             = pK(9);
            gy(sys.tp(i).ipKp2,iTSP)       = gpK(9,:);
            ggy(sys.tp(i).ipKp2,iTSP,iTSP) = ggpK(9,:,:);
        end
        if (isnan(obs.tp(i).pKp3))
            y(sys.tp(i).ipKp3)             = pK(10); 
            gy(sys.tp(i).ipKp3,iTSP)       = gpK(10,:);
            ggy(sys.tp(i).ipKp3,iTSP,iTSP) = ggpK(10,:,:);
        end  
        if (isnan(obs.tp(i).pKsi))
            y(sys.tp(i).ipKsi)             = pK(11);  
            gy(sys.tp(i).ipKsi,iTSP)       = gpK(11,:);
            ggy(sys.tp(i).ipKsi,iTSP,iTSP) = ggpK(11,:,:);
        end        
        if (isnan(obs.tp(i).pKnh4))
            y(sys.tp(i).ipKnh4)             = pK(12); 
            gy(sys.tp(i).ipKnh4,iTSP)       = gpK(12,:);
            ggy(sys.tp(i).ipKnh4,iTSP,iTSP) = ggpK(12,:,:);
        end   
        if (isnan(obs.tp(i).pKh2s))
            y(sys.tp(i).ipKh2s)             = pK(13); 
            gy(sys.tp(i).ipKh2s,iTSP)       = gpK(13,:); 
            ggy(sys.tp(i).ipKh2s,iTSP,iTSP) = ggpK(13,:,:); 
        end
        if (isnan(obs.tp(i).pp2f))
            y(sys.tp(i).ipp2f)             = pK(14);
            gy(sys.tp(i).ipp2f,iTSP)       = gpK(14,:);
            ggy(sys.tp(i).ipp2f,iTSP,iTSP) = ggpK(14,:,:);
        end
        if (isnan(obs.tp(i).pKar))
            y(sys.tp(i).ipKar)             = pK(15);
            gy(sys.tp(i).ipKar,iTSP)       = gpK(15,:);
            ggy(sys.tp(i).ipKar,iTSP,iTSP) = ggpK(15,:,:);
        end
        if (isnan(obs.tp(i).pKca))
            y(sys.tp(i).ipKca)             = pK(16);
            gy(sys.tp(i).ipKca,iTSP)       = gpK(16,:);
            ggy(sys.tp(i).ipKca,iTSP,iTSP) = ggpK(16,:,:);
        end

    end
end

% ----------------------------------------------------------------------
% init
% ----------------------------------------------------------------------

function z0 = init(opt,yobs,sys)
    q = sys.q;
    p = sys.p;
    
    y0  = yobs; 
    dic = q(yobs(sys.ipTC));
    alk = q(yobs(sys.ipTA));
    if (isnan(dic))
        dic         = 2200e-6;
        y0(sys.ipTC) = p(dic);
    end
    if (isnan(alk))
        alk         = 2200e-6;
        y0(sys.ipTA) = p(alk);
    end
    % calculate totals
    S = yobs(sys.isal);
    [pT,~,~,~] = calc_pTOT(opt,S);
    pTB  = pT(1);

    % calculate tp-independent vars
    gam     = dic/alk;

    if opt.turnoff.TB == 1
        TB = q(pTB);
        y0(sys.ipTB) = pTB;
    else
        TB = q(yobs(sys.ipTB));
        y0(sys.ipTB) = p(TB);
    end

    TS      = q(yobs(sys.ipTS));
    y0(sys.ipTS) = p(TS);

    TF      = q(yobs(sys.ipTF));
    y0(sys.ipTF) = p(TF);

    TP      = q(yobs(sys.ipTP));
    y0(sys.ipTP) = p(TP);

    TSi     = q(yobs(sys.ipTSi));
    y0(sys.ipTSi) = p(TSi);

    TNH4    = q(yobs(sys.ipTNH4));
    y0(sys.ipTNH4) = p(TNH4);

    TH2S    = q(yobs(sys.ipTH2S));
    y0(sys.ipTH2S) = p(TH2S);

    TCa     = q(yobs(sys.ipTCa));
    y0(sys.ipTCa) = p(TCa);

    if opt.pKalpha == 1 
        TAlpha = q(yobs(sys.ipTAlpha));
        y0(sys.ipTAlpha) = p(TAlpha);
    end
    if opt.pKbeta == 1
        TBeta = q(yobs(sys.ipTBeta));
        y0(sys.ipTBeta) = p(TBeta);
    end

    nTP = length(sys.tp);
    for i = 1:nTP
        % calculate pK1
        T = y0(sys.tp(i).iT); P = y0(sys.tp(i).iP);
        [pK,~,~] = calc_pK(opt,T,S,P);
        if opt.turnoff.pK1 == 1
            pK1 = pK(2);
            K1  = q(pK1);
            y0(sys.tp(i).ipK1) = pK1;
        else
            K1  = q(y0(sys.tp(i).ipK1));
        end
        if opt.turnoff.pK2 == 1
            pK2 = pK(3);
            K2  = q(pK2);
            y0(sys.tp(i).ipK2) = pK2;
        else
            K2  = q(y0(sys.tp(i).ipK2));
        end

        % solve for the [H+] using only the carbonate alkalinity
        K0      = q(y0(sys.tp(i).ipK0));
        % K1      = q(y0(sys.tp(i).ipK1));
        % K2      = q(y0(sys.tp(i).ipK2));
        h       = 0.5*( ( gam - 1 ) * K1 + ( ( 1 - gam )^2 * ...
                                             K1^2 - 4 * K1 * K2 * ( 1 - 2 * gam ) ).^0.5 ) ;
        hco3    = h * alk / (h + 2 * K2 );
        co2st   = h * hco3 / K1 ;
        co3     = dic*K1*K2/(K1*h + h*h + K1*K2) ;
        fco2    = co2st/K0;
        
        y0(sys.tp(i).iph)     = p(h);
        y0(sys.tp(i).iphco3)  = p(hco3);
        y0(sys.tp(i).ipco2st) = p(co2st);
        y0(sys.tp(i).ipco3)   = p(co3);
        y0(sys.tp(i).ipfco2)  = p(fco2);
        
        Kb      = q(y0(sys.tp(i).ipKb));
        boh4    = TB * Kb / (Kb + h) ;
        boh3    = TB - boh4;
        y0(sys.tp(i).ipboh3) = p(boh3);
        y0(sys.tp(i).ipboh4) = p(boh4);
        
        Kw      = q(y0(sys.tp(i).ipKw));
        oh      = Kw / h;
        y0(sys.tp(i).ipoh)   = p(oh);
        
        Ks      = q(y0(sys.tp(i).ipKs));
        
        h_tot   = h * ( 1 + TS / Ks );
        h_free  = h;
        hso4    = h_tot - h_free;
        so4     = Ks * hso4 / h;
        y0(sys.tp(i).iphso4)    = p(hso4);
        y0(sys.tp(i).ipso4)     = p(so4);
        
        Kf      = q(y0(sys.tp(i).ipKf));
        HF      = TF / ( 1 + Kf / h_free );
        F       = Kf * HF / h_free;
        h_sws   = h_tot + HF;
        y0(sys.tp(i).ipF)    = p(F);
        y0(sys.tp(i).ipHF)   = p(HF);
        
        Kp1     = q(y0(sys.tp(i).ipKp1));
        Kp2     = q(y0(sys.tp(i).ipKp2));
        Kp3     = q(y0(sys.tp(i).ipKp3));
        d       = ( h^3 + Kp1 * h^2 + Kp1 * Kp2 * h + Kp1 * Kp2 * Kp3);
        h3po4   = TP * h^3 / d;
        h2po4   = TP * Kp1 * h^2 / d;
        hpo4    = TP * Kp1 * Kp2 * h / d;
        po4     = TP * Kp1 * Kp2 * Kp3 / d;
        y0(sys.tp(i).iph3po4)    = p(h3po4);
        y0(sys.tp(i).iph2po4)    = p(h2po4);
        y0(sys.tp(i).iphpo4)     = p(hpo4);
        y0(sys.tp(i).ippo4)      = p(po4);
        
        Ksi     = q(y0(sys.tp(i).ipKsi));
        siooh3  = TSi / ( 1 + h / Ksi );
        sioh4   = TSi - siooh3;
        y0(sys.tp(i).ipsiooh3)   = p(siooh3);
        y0(sys.tp(i).ipsioh4)    = p(sioh4);
        
        Knh4    = q(y0(sys.tp(i).ipKnh4));
        nh3     = TNH4 / ( 1 + h / Knh4 );
        nh4     = TNH4 - nh3 ;
        y0(sys.tp(i).ipnh3)  = p(nh3);
        y0(sys.tp(i).ipnh4)  = p(nh4);
        
        Kh2s    = q(y0(sys.tp(i).ipKh2s));
        hs      = TH2S / ( 1 + h / Kh2s );
        h2s     = TH2S - hs ;
        y0(sys.tp(i).ipHS)   = p(hs);
        y0(sys.tp(i).ipH2S)  = p(h2s);
        
        p2f     = q(y0(sys.tp(i).ipp2f));
        pco2    = fco2/p2f;
        y0(sys.tp(i).ippco2)  = p(pco2);
        
        Kar     = q(y0(sys.tp(i).ipKar));
        OmegaAr = co3 * TCa / Kar;
        Kca     = q(y0(sys.tp(i).ipKca));
        OmegaCa = co3 * TCa / Kca ;
        y0(sys.tp(i).ipca)       = p(TCa);
        y0(sys.tp(i).ipOmegaAr)  = p(OmegaAr);
        y0(sys.tp(i).ipOmegaCa)  = p(OmegaCa);
                
        y0(sys.tp(i).iph_tot)   = p(h_tot);
        y0(sys.tp(i).iph_sws)   = p(h_sws);
        y0(sys.tp(i).iph_free)  = p(h_free);
        fH      = q(y0(sys.tp(i).ipfH));
        h_nbs   = h_free*fH;
        y0(sys.tp(i).iph_nbs)   = p(h_nbs);

        if opt.pKalpha == 1
            Kalpha = q(y0(sys.tp(i).ipKalpha));   
            alpha = (Kalpha*TAlpha)/(h+Kalpha);
            halpha = TAlpha - alpha;
            % alpha = (1/2)*Talpha; % make a function of pH
            % halpha = (1/2)*Talpha;
            y0(sys.tp(i).ipalpha) = p(alpha);
            y0(sys.tp(i).iphalpha) = p(halpha);
        end
        if opt.pKbeta == 1
            Kbeta = q(y0(sys.tp(i).ipKbeta));
            beta = (Kbeta*TBeta)/(h+Kbeta);
            hbeta = TBeta - beta;
            % beta = (1/2)*Tbeta;
            % hbeta = (1/2)*Tbeta;
            y0(sys.tp(i).ipbeta) = p(beta);
            y0(sys.tp(i).iphbeta) = p(hbeta);
        end
    end
    nlam    = size(sys.M,1) + size(sys.K,1);
    lam     = zeros(nlam,1);
    z0      = [y0(:);lam(:)];
end

% ----------------------------------------------------------------------
% parse_output
% ----------------------------------------------------------------------

function [est] = parse_output(z,sigx,sys,f,C)
    % populate est, output structure with best estimates
    %
    % INPUT:
    %   z  := system state including equilibrium constants and lagrange multipliers
    
    p       = sys.p;
    q       = sys.q;

    ebar    = @(j) (0.5 * ( q( z(j) - sigx(j) ) - q( z(j) + sigx(j) ) ) );
    ebar_l  = @(j) ( q( z(j)) - ( q( z(j) + sigx(j) ) )  ); % lower sigma
    ebar_u  = @(j) ( q( z(j) - sigx(j) ) - q(z(j)) ); % upper sigma
        
    % populate 'est' structure with best estimate:
    %   1. p(value) and p(error) where p(x) = -log10(x)
    %   2. value and average uncertainty about the value in 'q' 
    %           where q(x) = x^(-10)
    %   3. upper and lower bounds in 'q' space, not symmetric
    %           about the value in 'q' space

    % PLEASE ADD VARIABLES AT END
    % CHANGING THIS ORDER BREAKS PrintCSV
    est.sal     = z(sys.isal);
    est.usal    = sigx(sys.isal);
    
    % TC (DIC)
    est.pTC     = z(sys.ipTC);               
    est.upTC    = sigx(sys.ipTC);    
    est.TC      = q(z(sys.ipTC))*1e6; % 1e6  converts mol/kg to µmol/kg
    est.uTC     = ebar(sys.ipTC)*1e6;      
    est.uTC_l   = ebar_l(sys.ipTC)*1e6;    
    est.uTC_u   = ebar_u(sys.ipTC)*1e6;

    % TA Alkalinity
    est.pTA     = z(sys.ipTA);               
    est.upTA    = sigx(sys.ipTA);
    est.TA      = q(z(sys.ipTA))*1e6; % convt 
    est.uTA     = ebar(sys.ipTA)*1e6;      
    est.uTA_l   = ebar_l(sys.ipTA)*1e6;    
    est.uTA_u   = ebar_u(sys.ipTA)*1e6;

    % TB borate
    est.pTB     = z(sys.ipTB);
    est.upTB    = sigx(sys.ipTB);
    if (sys.freshwater)
        est.TB      = 0;
        est.uTB     = 0;
        est.uTB_l   = 0;
        est.uTB_u   = 0;
    else
        est.TB      = q(z(sys.ipTB))*1e6; % convt mol/kg to µmol/kg
        est.uTB     = ebar(sys.ipTB)*1e6;
        est.uTB_l   = ebar_l(sys.ipTB)*1e6;
        est.uTB_u   = ebar_u(sys.ipTB)*1e6;
    end
    
    % TS sulfate
    est.pTS     = z(sys.ipTS);
    est.upTS    = sigx(sys.ipTS);
    if (sys.freshwater)
        est.TS      = 0;
        est.uTS     = 0;
        est.uTS_l   = 0;
        est.uTS_u   = 0;
    else
        est.TS      = q(z(sys.ipTS))*1e6;  % convt mol/kg to µmol/kg
        est.uTS     = ebar(sys.ipTS)*1e6;
        est.uTS_l   = ebar_l(sys.ipTS)*1e6;
        est.uTS_u   = ebar_u(sys.ipTS)*1e6;
    end
    
    % TF fluoride
    est.pTF     = z(sys.ipTF);
    est.upTF        = sigx(sys.ipTF);
    if (sys.freshwater)
        est.TF      = 0;
        est.uTF     = 0;
        est.uTF_l   = 0;
        est.uTF_u   = 0;
    else
        est.TF      = q(z(sys.ipTF))*1e6; % convt mol/kg to µmol/kg
        est.uTF     = ebar(sys.ipTF)*1e6;
        est.uTF_l   = ebar_l(sys.ipTF)*1e6;
        est.uTF_u   = ebar_u(sys.ipTF)*1e6;
    end

    % TP Phosphate
    est.pTP     = z(sys.ipTP);           
    est.upTP    = sigx(sys.ipTP);
    if (sys.freshwater)
        est.TP      = 0;
        est.uTP     = 0;
        est.uTP_l   = 0;
        est.uTP_u   = 0;
    else
        est.TP      = q(z(sys.ipTP))*1e6;   % convt mol/kg to µmol/kg
        est.uTP     = ebar(sys.ipTP)*1e6;
        est.uTP_l   = ebar_l(sys.ipTP)*1e6;
        est.uTP_u   = ebar_u(sys.ipTP)*1e6;
    end

    % TSi silicate
    est.pTSi    = z(sys.ipTSi);         
    est.upTSi   = sigx(sys.ipTSi);
    if (sys.freshwater)
        est.TSi      = 0;
        est.uTSi     = 0;
        est.uTSi_l   = 0;
        est.uTSi_u   = 0;
    else
        est.TSi      = q(z(sys.ipTSi))*1e6;  % convt mol/kg to µmol/kg
        est.uTSi     = ebar(sys.ipTSi)*1e6;
        est.uTSi_l   = ebar_l(sys.ipTSi)*1e6;
        est.uTSi_u   = ebar_u(sys.ipTSi)*1e6;
    end
    
    % TNH4 nitrate
    est.pTNH4   = z(sys.ipTNH4);       
    est.upTNH4  = sigx(sys.ipTNH4);
    if (sys.freshwater)
        est.TNH4     = 0;
        est.uTNH4    = 0;
        est.uTNH4_l  = 0;
        est.uTNH4_u  = 0;
    else
        est.TNH4     = q(z(sys.ipTNH4))*1e6;   % convt mol/kg to µmol/kg
        est.uTNH4    = ebar(sys.ipTNH4)*1e6;
        est.uTNH4_l  = ebar_l(sys.ipTNH4)*1e6;
        est.uTNH4_u  = ebar_u(sys.ipTNH4)*1e6;
    end
    
    % TH2S sulfide
    est.pTH2S   = z(sys.ipTH2S);       
    est.upTH2S  = sigx(sys.ipTH2S);
    if (sys.freshwater)
        est.TH2S      = 0;
        est.uTH2S     = 0;
        est.uTH2S_l   = 0;
        est.uTH2S_u   = 0;
    else
        est.TH2S      = q(z(sys.ipTH2S))*1e6; % convt mol/kg to µmol/kg
        est.uTH2S     = ebar(sys.ipTH2S)*1e6;
        est.uTH2S_l   = ebar_l(sys.ipTH2S)*1e6;
        est.uTH2S_u   = ebar_u(sys.ipTH2S)*1e6;
    end

    % TCa calcium
    est.pTCa   = z(sys.ipTCa);       
    est.upTCa  = sigx(sys.ipTCa);
    if (sys.freshwater)
        est.TCa      = 0;
        est.uTCa     = 0;
        est.uTCa_l   = 0;
        est.uTCa_u   = 0;
    else
        est.TCa      = q(z(sys.ipTCa))*1e6; % convt mol/kg to µmol/kg
        est.uTCa     = ebar(sys.ipTCa)*1e6;
        est.uTCa_l   = ebar_l(sys.ipTCa)*1e6;
        est.uTCa_u   = ebar_u(sys.ipTCa)*1e6;
    end
    
    if (isfield(sys,'ipTAlpha'))
        est.pTAlpha     = z(sys.ipTAlpha);
        est.upTAlpha    = sigx(sys.ipTAlpha);
        if (sys.freshwater)
            est.TAlpha      = 0;
            est.uTAlpha     = 0;
            est.uTAlpha_l   = 0;
            est.uTAlpha_u   = 0;
        else
            est.TAlpha      = q(z(sys.ipTAlpha))*1e6;
            est.uTAlpha     = ebar(sys.ipTAlpha)*1e6;
            est.uTAlpha_l   = ebar_l(sys.ipTAlpha)*1e6;
            est.uTAlpha_u   = ebar_u(sys.ipTAlpha)*1e6;
        end
    end
    if (isfield(sys,'ipTBeta'))
        est.pTBeta      = z(sys.ipTBeta);
        est.upTBeta     = sigx(sys.ipTBeta);
        if (sys.freshwater)
            est.TBeta      = 0;
            est.uTBeta     = 0;
            est.uTBeta_l   = 0;
            est.uTBeta_u   = 0;
        else
            est.TBeta      = q(z(sys.ipTBeta))*1e6;
            est.uTBeta     = ebar(sys.ipTBeta)*1e6;
            est.uTBeta_l   = ebar_l(sys.ipTBeta)*1e6;
            est.uTBeta_u   = ebar_u(sys.ipTBeta)*1e6;
        end
    end

    nTP = length(sys.tp);
    for i = 1:nTP
        % temp (deg C)
        est.tp(i).T     = z(sys.tp(i).iT);
        est.tp(i).uT    = sigx(sys.tp(i).iT);
        est.tp(i).uT_l  = z(sys.tp(i).iT)-sigx(sys.tp(i).iT);
        est.tp(i).uT_u  = z(sys.tp(i).iT)+sigx(sys.tp(i).iT);
        
        % pressure (dbar)
        est.tp(i).P     = z(sys.tp(i).iP);        
        est.tp(i).uP    = sigx(sys.tp(i).iP);
        est.tp(i).uP_l  = z(sys.tp(i).iP)-sigx(sys.tp(i).iP);
        est.tp(i).uP_u  = z(sys.tp(i).iP)+sigx(sys.tp(i).iP);
        
        % fCO2
        est.tp(i).fco2      = q(z(sys.tp(i).ipfco2)) * 1e6; % convt atm to µatm
        est.tp(i).ufco2     = ebar(sys.tp(i).ipfco2) * 1e6;
        est.tp(i).ufco2_l   = ebar_l(sys.tp(i).ipfco2) * 1e6;
        est.tp(i).ufco2_u   = ebar_u(sys.tp(i).ipfco2) * 1e6;
        est.tp(i).pfco2     = z(sys.tp(i).ipfco2);
        est.tp(i).upfco2    = sigx(sys.tp(i).ipfco2);

        % pCO2
        est.tp(i).pco2      = q(z(sys.tp(i).ippco2)) * 1e6; % convt atm to µatm
        est.tp(i).upco2     = ebar(sys.tp(i).ippco2) * 1e6;
        est.tp(i).upco2_l   = ebar_l(sys.tp(i).ippco2) * 1e6;
        est.tp(i).upco2_u   = ebar_u(sys.tp(i).ippco2) * 1e6;
        est.tp(i).ppco2     = z(sys.tp(i).ippco2);
        est.tp(i).uppco2    = sigx(sys.tp(i).ippco2);
        
        % HCO3
        est.tp(i).hco3      = q(z(sys.tp(i).iphco3))*1e6; % convt mol/kg to µmol/kg
        est.tp(i).uhco3     = ebar(sys.tp(i).iphco3)*1e6;
        est.tp(i).uhco3_l   = ebar_l(sys.tp(i).iphco3)*1e6;
        est.tp(i).uhco3_u   = ebar_u(sys.tp(i).iphco3)*1e6;
        est.tp(i).phco3     = z(sys.tp(i).iphco3);
        est.tp(i).uphco3    = sigx(sys.tp(i).iphco3);

        % CO3
        est.tp(i).co3       = q(z(sys.tp(i).ipco3))*1e6; % convt mol/kg to µmol/kg
        est.tp(i).uco3      = ebar(sys.tp(i).ipco3)*1e6;
        est.tp(i).uco3_l    = ebar_l(sys.tp(i).ipco3)*1e6;
        est.tp(i).uco3_u    = ebar_u(sys.tp(i).ipco3)*1e6;
        est.tp(i).pco3      = z(sys.tp(i).ipco3);
        est.tp(i).upco3     = sigx(sys.tp(i).ipco3);

        % CO2*
        est.tp(i).co2st     = q(z(sys.tp(i).ipco2st)); % convt mol/kg to µmol/kg
        est.tp(i).uco2st    = ebar(sys.tp(i).ipco2st);
        est.tp(i).uco2st_l  = ebar_l(sys.tp(i).ipco2st);
        est.tp(i).uco2st_u  = ebar_u(sys.tp(i).ipco2st);
        est.tp(i).pco2st    = z(sys.tp(i).ipco2st);
        est.tp(i).upco2st   = sigx(sys.tp(i).ipco2st);

        % pH on the scale opt.phscale used to compute the pK values
        est.tp(i).ph        = z(sys.tp(i).iph); % (9)
        est.tp(i).uph       = sigx(sys.tp(i).iph);
        est.tp(i).h         = q(z(sys.tp(i).iph)) * 1e6;
        est.tp(i).uh        = ebar(sys.tp(i).iph) * 1e6;
        est.tp(i).uh_l      = ebar_l(sys.tp(i).iph) * 1e6;
        est.tp(i).uh_u      = ebar_u(sys.tp(i).iph) * 1e6;

        % pH_free
        est.tp(i).ph_free    = z(sys.tp(i).iph_free); % (15)
        est.tp(i).uph_free   = sigx(sys.tp(i).iph_free);
        est.tp(i).h_free     = q(z(sys.tp(i).iph_free)) * 1e6;
        est.tp(i).uh_free    = ebar(sys.tp(i).iph_free) * 1e6;
        est.tp(i).uh_free_l  = ebar_l(sys.tp(i).iph_free) * 1e6;
        est.tp(i).uh_free_u  = ebar_u(sys.tp(i).iph_free) * 1e6;

        % pH_tot
        est.tp(i).ph_tot    = z(sys.tp(i).iph_tot); % (21)
        est.tp(i).uph_tot   = sigx(sys.tp(i).iph_tot);
        est.tp(i).h_tot     = q(z(sys.tp(i).iph_tot)) * 1e6;
        est.tp(i).uh_tot    = ebar(sys.tp(i).iph_tot) * 1e6;
        est.tp(i).uh_tot_l  = ebar_l(sys.tp(i).iph_tot) * 1e6;
        est.tp(i).uh_tot_u  = ebar_u(sys.tp(i).iph_tot) * 1e6;

        % pH_sws
        est.tp(i).ph_sws    = z(sys.tp(i).iph_sws);
        est.tp(i).uph_sws   = sigx(sys.tp(i).iph_sws);
        est.tp(i).h_sws     = q(z(sys.tp(i).iph_sws)) * 1e6;
        est.tp(i).uh_sws    = ebar(sys.tp(i).iph_sws) * 1e6;
        est.tp(i).uh_sws_l  = ebar_l(sys.tp(i).iph_sws) * 1e6;
        est.tp(i).uh_sws_u  = ebar_u(sys.tp(i).iph_sws) * 1e6;

        % pH_nbs
        est.tp(i).ph_nbs    = z(sys.tp(i).iph_nbs);
        est.tp(i).uph_nbs   = sigx(sys.tp(i).iph_nbs);
        est.tp(i).h_nbs     = q(z(sys.tp(i).iph_nbs)) * 1e6;
        est.tp(i).uh_nbs    = ebar(sys.tp(i).iph_nbs) * 1e6;
        est.tp(i).uh_nbs_l  = ebar_l(sys.tp(i).iph_nbs) * 1e6;
        est.tp(i).uh_nbs_u  = ebar_u(sys.tp(i).iph_nbs) * 1e6;

        % fH = activity coefficient
        if (sys.freshwater)
            est.tp(i).fH      = 0;
            est.tp(i).ufH     = 0;
            est.tp(i).ufH_l   = 0;
            est.tp(i).ufH_u   = 0;
        else
            est.tp(i).fH        = q(z(sys.tp(i).ipfH)) ;
            est.tp(i).ufH       = ebar(sys.tp(i).ipfH) ;
            est.tp(i).ufH_l     = ebar_l(sys.tp(i).ipfH) ;
            est.tp(i).ufH_u     = ebar_u(sys.tp(i).ipfH) ;
        end
        est.tp(i).pfH       = z(sys.tp(i).ipfH);
        est.tp(i).upfH      = sigx(sys.tp(i).ipfH);

        % p2f
        est.tp(i).p2f       = q(z(sys.tp(i).ipp2f));
        est.tp(i).up2f      = ebar(sys.tp(i).ipp2f);
        est.tp(i).pp2f      = z(sys.tp(i).ipp2f); 
        est.tp(i).upp2f     = sigx(sys.tp(i).ipp2f);

        % pK0 
        est.tp(i).pK0       = z(sys.tp(i).ipK0); % (79)
        est.tp(i).upK0      = sigx(sys.tp(i).ipK0);
        est.tp(i).K0        = q(z(sys.tp(i).ipK0));
        est.tp(i).uK0       = ebar(sys.tp(i).ipK0);
        est.tp(i).uK0_l     = ebar_l(sys.tp(i).ipK0);
        est.tp(i).uK0_u     = ebar_u(sys.tp(i).ipK0);

        % pK1
        est.tp(i).pK1       = z(sys.tp(i).ipK1);
        est.tp(i).upK1      = sigx(sys.tp(i).ipK1);
        est.tp(i).K1        = q(z(sys.tp(i).ipK1));
        est.tp(i).uK1       = ebar(sys.tp(i).ipK1);
        est.tp(i).uK1_l     = ebar_l(sys.tp(i).ipK1); 
        est.tp(i).uK1_u     = ebar_u(sys.tp(i).ipK1);

        % pK2
        est.tp(i).pK2       = z(sys.tp(i).ipK2);
        est.tp(i).upK2      = sigx(sys.tp(i).ipK2);
        est.tp(i).K2        = q(z(sys.tp(i).ipK2));
        est.tp(i).uK2       = ebar(sys.tp(i).ipK2);
        est.tp(i).uK2_l     = ebar_l(sys.tp(i).ipK2);
        est.tp(i).uK2_u     = ebar_u(sys.tp(i).ipK2);
        
        % OH 
        est.tp(i).oh        = q(z(sys.tp(i).ipoh))*1e6; % convt
        est.tp(i).uoh       = ebar(sys.tp(i).ipoh)*1e6;
        est.tp(i).uoh_l     = ebar_l(sys.tp(i).ipoh)*1e6;
        est.tp(i).uoh_u     = ebar_u(sys.tp(i).ipoh)*1e6;
        est.tp(i).poh       = z(sys.tp(i).ipoh);
        est.tp(i).upoh      = sigx(sys.tp(i).ipoh);

        % pKw 
        est.tp(i).pKw       = z(sys.tp(i).ipKw); % (103)
        est.tp(i).upKw      = sigx(sys.tp(i).ipKw);
        est.tp(i).Kw        = q(z(sys.tp(i).ipKw));
        est.tp(i).uKw       = ebar(sys.tp(i).ipKw);
        est.tp(i).uKw_l     = ebar_l(sys.tp(i).ipKw);
        est.tp(i).uKw_u     = ebar_u(sys.tp(i).ipKw);

        % BOH4 borate
        if (sys.freshwater)
            est.tp(i).boh4      = 0;
            est.tp(i).uboh4     = 0;
            est.tp(i).uboh4_l   = 0;
            est.tp(i).uboh4_u   = 0;
        else
            est.tp(i).boh4      = q(z(sys.tp(i).ipboh4))*1e6; % convt mol/kg to µmol/kg
            est.tp(i).uboh4     = ebar(sys.tp(i).ipboh4)*1e6;
            est.tp(i).uboh4_l   = ebar_l(sys.tp(i).ipboh4)*1e6;
            est.tp(i).uboh4_u   = ebar_u(sys.tp(i).ipboh4)*1e6;
        end
        est.tp(i).pboh4     = z(sys.tp(i).ipboh4);
        est.tp(i).upboh4    = sigx(sys.tp(i).ipboh4);

        % BOH3
        if (sys.freshwater)
            est.tp(i).boh3      = 0;
            est.tp(i).uboh3     = 0;
            est.tp(i).uboh3_l   = 0;
            est.tp(i).uboh3_u   = 0;
        else
            est.tp(i).boh3      = q(z(sys.tp(i).ipboh3))*1e6;
            est.tp(i).uboh3     = ebar(sys.tp(i).ipboh3)*1e6;
            est.tp(i).uboh3_l   = ebar_l(sys.tp(i).ipboh3)*1e6;
            est.tp(i).uboh3_u   = ebar_u(sys.tp(i).ipboh3)*1e6;
        end
        est.tp(i).pboh3     = z(sys.tp(i).ipboh3);
        est.tp(i).upboh3    = sigx(sys.tp(i).ipboh3);
            
        % pKb
        est.tp(i).pKb       = z(sys.tp(i).ipKb); % (121)
        est.tp(i).upKb      = sigx(sys.tp(i).ipKb);
        if (sys.freshwater)
            est.tp(i).Kb      = 0;
            est.tp(i).uKb     = 0;
            est.tp(i).uKb_l   = 0;
            est.tp(i).uKb_u   = 0;
        else
            est.tp(i).Kb      = q(z(sys.tp(i).ipKb));
            est.tp(i).uKb     = ebar(sys.tp(i).ipKb);
            est.tp(i).uKb_l   = ebar_l(sys.tp(i).ipKb);
            est.tp(i).uKb_u   = ebar_u(sys.tp(i).ipKb);
        end

        % SO4 sulfate
        if (sys.freshwater)
            est.tp(i).so4      = 0;
            est.tp(i).uso4     = 0;
            est.tp(i).uso4_l   = 0;
            est.tp(i).uso4_u   = 0;
        else
            est.tp(i).so4      = q(z(sys.tp(i).ipso4))*1e6; % umol/kg
            est.tp(i).uso4     = ebar(sys.tp(i).ipso4)*1e6;
            est.tp(i).uso4_l   = ebar_l(sys.tp(i).ipso4)*1e6;
            est.tp(i).uso4_u   = ebar_u(sys.tp(i).ipso4)*1e6;
        end
        est.tp(i).pso4      = z(sys.tp(i).ipso4);
        est.tp(i).upso4     = sigx(sys.tp(i).ipso4);

        % HSO4
        if (sys.freshwater)
            est.tp(i).hso4      = 0;
            est.tp(i).uhso4     = 0;
            est.tp(i).uhso4_l   = 0;
            est.tp(i).uhso4_u   = 0;
        else
            est.tp(i).hso4      = q(z(sys.tp(i).iphso4))*1e6;
            est.tp(i).uhso4     = ebar(sys.tp(i).iphso4)*1e6;
            est.tp(i).uhso4_l   = ebar_l(sys.tp(i).iphso4)*1e6;
            est.tp(i).uhso4_u   = ebar_u(sys.tp(i).iphso4)*1e6;
        end
        est.tp(i).phso4     = z(sys.tp(i).iphso4);
        est.tp(i).uphso4    = sigx(sys.tp(i).iphso4);

        % pKs
        est.tp(i).pKs       = z(sys.tp(i).ipKs); % (145)
        est.tp(i).upKs      = sigx(sys.tp(i).ipKs);
        if (sys.freshwater)
            est.tp(i).Ks      = 0;
            est.tp(i).uKs     = 0;
            est.tp(i).uKs_l   = 0;
            est.tp(i).uKs_u   = 0;
        else
            est.tp(i).Ks      = q(z(sys.tp(i).ipKs));
            est.tp(i).uKs     = ebar(sys.tp(i).ipKs);
            est.tp(i).uKs_l   = ebar_l(sys.tp(i).ipKs);
            est.tp(i).uKs_u   = ebar_u(sys.tp(i).ipKs);
        end

        % F fluoride 
        if (sys.freshwater)
            est.tp(i).F      = 0;
            est.tp(i).uF     = 0;
            est.tp(i).uF_l   = 0;
            est.tp(i).uF_u   = 0;
        else
            est.tp(i).F      = q(z(sys.tp(i).ipF))*1e6; % convt
            est.tp(i).uF     = ebar(sys.tp(i).ipF)*1e6;
            est.tp(i).uF_l   = ebar_l(sys.tp(i).ipF)*1e6;
            est.tp(i).uF_u   = ebar_u(sys.tp(i).ipF)*1e6;
        end
        est.tp(i).pF    = z(sys.tp(i).ipF);
        est.tp(i).upF   = sigx(sys.tp(i).ipF);

        % HF 
        if (sys.freshwater)
            est.tp(i).HF      = 0;
            est.tp(i).uHF     = 0;
            est.tp(i).uHF_l   = 0;
            est.tp(i).uHF_u   = 0;
        else
            est.tp(i).HF      = q(z(sys.tp(i).ipHF))*1e6;
            est.tp(i).uHF     = ebar(sys.tp(i).ipHF)*1e6;
            est.tp(i).uHF_l   = ebar_l(sys.tp(i).ipHF)*1e6;
            est.tp(i).uHF_u   = ebar_u(sys.tp(i).ipHF)*1e6;
        end
        est.tp(i).pHF   = z(sys.tp(i).ipHF);
        est.tp(i).upHF  = sigx(sys.tp(i).ipHF);

        % pKf
        est.tp(i).pKf   = z(sys.tp(i).ipKf); % (163)
        est.tp(i).upKf  = sigx(sys.tp(i).ipKf);
        if (sys.freshwater)
            est.tp(i).Kf      = 0;
            est.tp(i).uKf     = 0;
            est.tp(i).uKf_l   = 0;
            est.tp(i).uKf_u   = 0;
        else
            est.tp(i).Kf      = q(z(sys.tp(i).ipKf));
            est.tp(i).uKf     = ebar(sys.tp(i).ipKf);
            est.tp(i).uKf_l   = ebar_l(sys.tp(i).ipKf);
            est.tp(i).uKf_u   = ebar_u(sys.tp(i).ipKf);
        end

        % PO4
        if (sys.freshwater)
            est.tp(i).po4      = 0;
            est.tp(i).upo4     = 0;
            est.tp(i).upo4_l   = 0;
            est.tp(i).upo4_u   = 0;
        else
            est.tp(i).po4      = q(z(sys.tp(i).ippo4))*1e6; % convt
            est.tp(i).upo4     = ebar(sys.tp(i).ippo4)*1e6;
            est.tp(i).upo4_l   = ebar_l(sys.tp(i).ippo4)*1e6;
            est.tp(i).upo4_u   = ebar_u(sys.tp(i).ippo4)*1e6;
        end
        est.tp(i).ppo4  = z(sys.tp(i).ippo4);
        est.tp(i).uppo4 = sigx(sys.tp(i).ippo4);

        % HPO4
        if (sys.freshwater)
            est.tp(i).hpo4      = 0;
            est.tp(i).uhpo4     = 0;
            est.tp(i).uhpo4_l   = 0;
            est.tp(i).uhpo4_u   = 0;
        else
            est.tp(i).hpo4      = q(z(sys.tp(i).iphpo4))*1e6;
            est.tp(i).uhpo4     = ebar(sys.tp(i).iphpo4)*1e6;
            est.tp(i).uhpo4_l   = ebar_l(sys.tp(i).iphpo4)*1e6;
            est.tp(i).uhpo4_u   = ebar_u(sys.tp(i).iphpo4)*1e6;
        end
        est.tp(i).phpo4     = z(sys.tp(i).iphpo4);
        est.tp(i).uphpo4    = sigx(sys.tp(i).iphpo4);

        % H2PO4
        if (sys.freshwater)
            est.tp(i).h2po4      = 0;
            est.tp(i).uh2po4     = 0;
            est.tp(i).uh2po4_l   = 0;
            est.tp(i).uh2po4_u   = 0;
        else
            est.tp(i).h2po4      = q(z(sys.tp(i).iph2po4))*1e6;
            est.tp(i).uh2po4     = ebar(sys.tp(i).iph2po4)*1e6;
            est.tp(i).uh2po4_l   = ebar_l(sys.tp(i).iph2po4)*1e6;
            est.tp(i).uh2po4_u   = ebar_u(sys.tp(i).iph2po4)*1e6;
        end
        est.tp(i).ph2po4    = z(sys.tp(i).iph2po4);
        est.tp(i).uph2po4   = sigx(sys.tp(i).iph2po4);            

        % H3PO4
        if (sys.freshwater)
            est.tp(i).h3po4      = 0;
            est.tp(i).uh3po4     = 0;
            est.tp(i).uh3po4_l   = 0;
            est.tp(i).uh3po4_u   = 0;
        else
            est.tp(i).h3po4      = q(z(sys.tp(i).iph3po4))*1e6;
            est.tp(i).uh3po4     = ebar(sys.tp(i).iph3po4)*1e6;
            est.tp(i).uh3po4_l   = ebar_l(sys.tp(i).iph3po4)*1e6;
            est.tp(i).uh3po4_u   = ebar_u(sys.tp(i).iph3po4)*1e6;
        end
        est.tp(i).ph3po4    = z(sys.tp(i).iph3po4);
        est.tp(i).uph3po4   = sigx(sys.tp(i).iph3po4);

        % pKp1
        est.tp(i).pKp1      = z(sys.tp(i).ipKp1);
        est.tp(i).upKp1     = sigx(sys.tp(i).ipKp1);
        if (sys.freshwater)
            est.tp(i).Kp1      = 0;
            est.tp(i).uKp1     = 0;
            est.tp(i).uKp1_l   = 0;
            est.tp(i).uKp1_u   = 0;
        else
            est.tp(i).Kp1      = q(z(sys.tp(i).ipKp1));
            est.tp(i).uKp1     = ebar(sys.tp(i).ipKp1);
            est.tp(i).uKp1_l   = ebar_l(sys.tp(i).ipKp1);
            est.tp(i).uKp1_u   = ebar_u(sys.tp(i).ipKp1);
        end

        % pKp2
        est.tp(i).pKp2      = z(sys.tp(i).ipKp2);
        est.tp(i).upKp2     = sigx(sys.tp(i).ipKp2);
        if (sys.freshwater)
            est.tp(i).Kp2      = 0;
            est.tp(i).uKp2     = 0;
            est.tp(i).uKp2_l   = 0;
            est.tp(i).uKp2_u   = 0;
        else
            est.tp(i).Kp2      = q(z(sys.tp(i).ipKp2));
            est.tp(i).uKp2     = ebar(sys.tp(i).ipKp2);
            est.tp(i).uKp2_l   = ebar_l(sys.tp(i).ipKp2);
            est.tp(i).uKp2_u   = ebar_u(sys.tp(i).ipKp2);
        end

        % pKp3
        est.tp(i).pKp3      = z(sys.tp(i).ipKp3);
        est.tp(i).upKp3     = sigx(sys.tp(i).ipKp3);
        if (sys.freshwater)
            est.tp(i).Kp3      = 0;
            est.tp(i).uKp3     = 0;
            est.tp(i).uKp3_l   = 0;
            est.tp(i).uKp3_u   = 0;
        else
            est.tp(i).Kp3      = q(z(sys.tp(i).ipKp3));
            est.tp(i).uKp3     = ebar(sys.tp(i).ipKp3);
            est.tp(i).uKp3_l   = ebar_l(sys.tp(i).ipKp3);
            est.tp(i).uKp3_u   = ebar_u(sys.tp(i).ipKp3);
        end
        
        % SiOH4
        if (sys.freshwater)
            est.tp(i).sioh4     = 0;
            est.tp(i).usioh4    = 0;
            est.tp(i).usioh4_l  = 0;
            est.tp(i).usioh4_u  = 0;
        else
            est.tp(i).sioh4     = q(z(sys.tp(i).ipsioh4))*1e6; % convt
            est.tp(i).usioh4    = ebar(sys.tp(i).ipsioh4)*1e6;
            est.tp(i).usioh4_l  = ebar_l(sys.tp(i).ipsioh4)*1e6;
            est.tp(i).usioh4_u  = ebar_u(sys.tp(i).ipsioh4)*1e6;
        end
        est.tp(i).psioh4    = z(sys.tp(i).ipsioh4);
        est.tp(i).upsioh4   = sigx(sys.tp(i).ipsioh4);

        % SiOH3
        if (sys.freshwater)
            est.tp(i).siooh3      = 0;
            est.tp(i).usiooh3     = 0;
            est.tp(i).usiooh3_l   = 0;
            est.tp(i).usiooh3_u   = 0;
        else
            est.tp(i).siooh3      = q(z(sys.tp(i).ipsiooh3))*1e6;
            est.tp(i).usiooh3     = ebar(sys.tp(i).ipsiooh3)*1e6;
            est.tp(i).usiooh3_l   = ebar_l(sys.tp(i).ipsiooh3)*1e6;
            est.tp(i).usiooh3_u   = ebar_u(sys.tp(i).ipsiooh3)*1e6;
        end
        est.tp(i).psiooh3   = z(sys.tp(i).ipsiooh3);
        est.tp(i).upsiooh3  = sigx(sys.tp(i).ipsiooh3);

        % pKsi
        est.tp(i).pKsi      = z(sys.tp(i).ipKsi);
        est.tp(i).upKsi     = sigx(sys.tp(i).ipKsi);
        if (sys.freshwater)
            est.tp(i).Ksi      = 0;
            est.tp(i).uKsi     = 0;
            est.tp(i).uKsi_l   = 0;
            est.tp(i).uKsi_u   = 0;
        else
            est.tp(i).Ksi      = q(z(sys.tp(i).ipKsi));
            est.tp(i).uKsi     = ebar(sys.tp(i).ipKsi);
            est.tp(i).uKsi_l   = ebar_l(sys.tp(i).ipKsi);
            est.tp(i).uKsi_u   = ebar_u(sys.tp(i).ipKsi);
        end

        % NH3
        if (sys.freshwater)
            est.tp(i).nh3      = 0;
            est.tp(i).unh3     = 0;
            est.tp(i).unh3_l   = 0;
            est.tp(i).unh3_u   = 0;
        else
            est.tp(i).nh3      = q(z(sys.tp(i).ipnh3))*1e6; % convt
            est.tp(i).unh3     = ebar(sys.tp(i).ipnh3)*1e6;
            est.tp(i).unh3_l   = ebar_l(sys.tp(i).ipnh3)*1e6;
            est.tp(i).unh3_u   = ebar_u(sys.tp(i).ipnh3)*1e6;
        end
        est.tp(i).pnh3      = z(sys.tp(i).ipnh3);
        est.tp(i).upnh3     = sigx(sys.tp(i).ipnh3);

        % NH4
        if (sys.freshwater)
            est.tp(i).nh4      = 0;
            est.tp(i).unh4     = 0;
            est.tp(i).unh4_l   = 0;
            est.tp(i).unh4_u   = 0;
        else
            est.tp(i).nh4      = q(z(sys.tp(i).ipnh4))*1e6;
            est.tp(i).unh4     = ebar(sys.tp(i).ipnh4)*1e6;
            est.tp(i).unh4_l   = ebar_l(sys.tp(i).ipnh4)*1e6;
            est.tp(i).unh4_u   = ebar_u(sys.tp(i).ipnh4)*1e6;
        end
        est.tp(i).pnh4      = z(sys.tp(i).ipnh4);
        est.tp(i).upnh4     = sigx(sys.tp(i).ipnh4);

        % pKNH4
        est.tp(i).pKnh4     = z(sys.tp(i).ipKnh4);
        est.tp(i).upKnh4    = sigx(sys.tp(i).ipKnh4);
        if (sys.freshwater)
            est.tp(i).Knh4      = 0;
            est.tp(i).uKnh4     = 0;
            est.tp(i).uKnh4_l   = 0;
            est.tp(i).uKnh4_u   = 0;
        else
            est.tp(i).Knh4      = q(z(sys.tp(i).ipKnh4));
            est.tp(i).uKnh4     = ebar(sys.tp(i).ipKnh4);
            est.tp(i).uKnh4_l   = ebar_l(sys.tp(i).ipKnh4);
            est.tp(i).uKnh4_u   = ebar_u(sys.tp(i).ipKnh4);
        end

        % HS
        if (sys.freshwater)
            est.tp(i).HS      = 0;
            est.tp(i).uHS     = 0;
            est.tp(i).uHS_l   = 0;
            est.tp(i).uHS_u   = 0;
        else
            est.tp(i).HS      = q(z(sys.tp(i).ipHS))*1e6; % convt
            est.tp(i).uHS     = ebar(sys.tp(i).ipHS)*1e6;
            est.tp(i).uHS_l   = ebar_l(sys.tp(i).ipHS)*1e6;
            est.tp(i).uHS_u   = ebar_u(sys.tp(i).ipHS)*1e6;
        end
        est.tp(i).pHS       = z(sys.tp(i).ipHS);
        est.tp(i).upHS      = sigx(sys.tp(i).ipHS);

        % H2S
        if (sys.freshwater)
            est.tp(i).H2S      = 0;
            est.tp(i).uH2S     = 0;
            est.tp(i).uH2S_l   = 0;
            est.tp(i).uH2S_u   = 0;
        else
            est.tp(i).H2S      = q(z(sys.tp(i).ipH2S))*1e6;
            est.tp(i).uH2S     = ebar(sys.tp(i).ipH2S)*1e6;
            est.tp(i).uH2S_l   = ebar_l(sys.tp(i).ipH2S)*1e6;
            est.tp(i).uHS2_u   = ebar_u(sys.tp(i).ipH2S)*1e6;
        end
        est.tp(i).pH2S      = z(sys.tp(i).ipH2S);
        est.tp(i).upH2S     = sigx(sys.tp(i).ipH2S);

        % pKh2s
        est.tp(i).pKh2s     = z(sys.tp(i).ipKh2s);
        est.tp(i).upKh2s    = sigx(sys.tp(i).ipKh2s);
        if (sys.freshwater)
            est.tp(i).Kh2s      = 0;
            est.tp(i).uKh2s     = 0;
            est.tp(i).uKh2s_l   = 0;
            est.tp(i).uKh2s_u   = 0;
        else
            est.tp(i).Kh2s      = q(z(sys.tp(i).ipKh2s));
            est.tp(i).uKh2s     = ebar(sys.tp(i).ipKh2s);
            est.tp(i).uKh2s_l   = ebar_l(sys.tp(i).ipKh2s);
            est.tp(i).uKh2s_u   = ebar_u(sys.tp(i).ipKh2s);
        end

        % Ca
        if (sys.freshwater)
            est.tp(i).ca      = 0;
            est.tp(i).uca     = 0;
            est.tp(i).uca_l   = 0;
            est.tp(i).uca_u   = 0;
        else
            est.tp(i).ca      = q(z(sys.tp(i).ipca))*1e6;
            est.tp(i).uca     = ebar(sys.tp(i).ipca)*1e6;
            est.tp(i).uca_l   = ebar_l(sys.tp(i).ipca)*1e6;
            est.tp(i).uca_u   = ebar_u(sys.tp(i).ipca)*1e6;
        end
        est.tp(i).pca       = z(sys.tp(i).ipca);
        est.tp(i).upca      = sigx(sys.tp(i).ipca);

        % Omega_Ar
        if (sys.freshwater)
            est.tp(i).OmegaAr      = 0;
            est.tp(i).uOmegaAr     = 0;
            est.tp(i).uOmegaAr_l   = 0;
            est.tp(i).uOmegaAr_u   = 0;
        else
            est.tp(i).OmegaAr      = q(z(sys.tp(i).ipOmegaAr)); % unitless
            est.tp(i).uOmegaAr     = ebar(sys.tp(i).ipOmegaAr);
            est.tp(i).uOmegaAr_l   = ebar_l(sys.tp(i).ipOmegaAr);
            est.tp(i).uOmegaAr_u   = ebar_u(sys.tp(i).ipOmegaAr);
        end
        est.tp(i).pOmegaAr   = z(sys.tp(i).ipOmegaAr);
        est.tp(i).upOmegaAr  = sigx(sys.tp(i).ipOmegaAr);

        % pKar
        est.tp(i).pKar      = z(sys.tp(i).ipKar);
        est.tp(i).upKar     = sigx(sys.tp(i).ipKar);
        if (sys.freshwater)
            est.tp(i).Kar      = 0;
            est.tp(i).uKar     = 0;
            est.tp(i).uKar_l   = 0;
            est.tp(i).uKar_u   = 0;
        else
            est.tp(i).Kar      = q(z(sys.tp(i).ipKar));
            est.tp(i).uKar     = ebar(sys.tp(i).ipKar);
            est.tp(i).uKar_l   = ebar_l(sys.tp(i).ipKar);
            est.tp(i).uKar_u   = ebar_u(sys.tp(i).ipKar);
        end

        % Omega_Ca
        if (sys.freshwater)
            est.tp(i).OmegaCa      = 0;
            est.tp(i).uOmegaCa     = 0;
            est.tp(i).uOmegaCa_l   = 0;
            est.tp(i).uOmegaCa_u   = 0;
        else
            est.tp(i).OmegaCa      = q(z(sys.tp(i).ipOmegaCa));
            est.tp(i).uOmegaCa     = ebar(sys.tp(i).ipOmegaCa);
            est.tp(i).uOmegaCa_l   = ebar_l(sys.tp(i).ipOmegaCa);
            est.tp(i).uOmegaCa_u   = ebar_u(sys.tp(i).ipOmegaCa);
        end
        est.tp(i).pOmegaCa    = z(sys.tp(i).ipOmegaCa);
        est.tp(i).upOmegaCa   = sigx(sys.tp(i).ipOmegaCa);

        % pKca
        est.tp(i).pKca    = z(sys.tp(i).ipKca);
        est.tp(i).upKca   = sigx(sys.tp(i).ipKca);
        if (sys.freshwater)
            est.tp(i).Kca      = 0;
            est.tp(i).uKca     = 0;
            est.tp(i).uKca_l   = 0;
            est.tp(i).uKca_u   = 0;
        else
            est.tp(i).Kca      = q(z(sys.tp(i).ipKca));
            est.tp(i).uKca     = ebar(sys.tp(i).ipKca);
            est.tp(i).uKca_l   = ebar_l(sys.tp(i).ipKca);
            est.tp(i).uKca_u   = ebar_u(sys.tp(i).ipKca); % #288 in CSV
        end

        if (isfield(sys,'ipTAlpha'))
            % alpha
            if (sys.freshwater)
                est.tp(i).alpha     = 0;
                est.tp(i).ualpha    = 0;
                est.tp(i).ualpha_l  = 0;
                est.tp(i).ualpha_u  = 0;
            else
                est.tp(i).alpha     = q(z(sys.tp(i).ipalpha))*1e6;
                est.tp(i).ualpha    = ebar(sys.tp(i).ipalpha)*1e6;
                est.tp(i).ualpha_l  = ebar_l(sys.tp(i).ipalpha)*1e6;
                est.tp(i).ualpha_u  = ebar_u(sys.tp(i).ipalpha)*1e6;
            end
            est.tp(i).palpha    = z(sys.tp(i).ipalpha);
            est.tp(i).upalpha   = sigx(sys.tp(i).ipalpha);

            % halpha
            if (sys.freshwater)
                est.tp(i).halpha     = 0;
                est.tp(i).uhalpha    = 0;
                est.tp(i).uhalpha_l  = 0;
                est.tp(i).uhalpha_u  = 0;
            else
                est.tp(i).halpha     = q(z(sys.tp(i).iphalpha))*1e6;
                est.tp(i).uhalpha    = ebar(sys.tp(i).iphalpha)*1e6;
                est.tp(i).uhalpha_l  = ebar_l(sys.tp(i).iphalpha)*1e6;
                est.tp(i).uhalpha_u  = ebar_u(sys.tp(i).iphalpha)*1e6;
            end
            est.tp(i).phalpha   = z(sys.tp(i).iphalpha);
            est.tp(i).uphalpha  = sigx(sys.tp(i).iphalpha);
            
            % pKalpha
            est.tp(i).pKalpha   = z(sys.tp(i).ipKalpha);
            est.tp(i).upKalpha  = sigx(sys.tp(i).ipKalpha);
            if (sys.freshwater)
                est.tp(i).Kalpha     = 0;
                est.tp(i).uKalpha    = 0;
                est.tp(i).uKalpha_l  = 0;
                est.tp(i).uKalpha_u  = 0;
            else
                est.tp(i).Kalpha     = q(z(sys.tp(i).ipKalpha));
                est.tp(i).uKalpha    = ebar(sys.tp(i).ipKalpha);
                est.tp(i).uKalpha_l  = ebar_l(sys.tp(i).ipKalpha);
                est.tp(i).uKalpha_u  = ebar_u(sys.tp(i).ipKalpha);
            end
        end

        if (isfield(sys,'ipTBeta'))
            % beta
            if (sys.freshwater)
                est.tp(i).beta      = 0;
                est.tp(i).ubeta     = 0;
                est.tp(i).ubeta_l   = 0;
                est.tp(i).ubeta_u   = 0;
            else
                est.tp(i).beta      = q(z(sys.tp(i).ipbeta))*1e6;
                est.tp(i).ubeta     = ebar(sys.tp(i).ipbeta)*1e6;
                est.tp(i).ubeta_l   = ebar_l(sys.tp(i).ipbeta)*1e6;
                est.tp(i).ubeta_u   = ebar_u(sys.tp(i).ipbeta)*1e6;
            end
            est.tp(i).pbeta     = z(sys.tp(i).ipbeta);
            est.tp(i).upbeta    = sigx(sys.tp(i).ipbeta);

            % hbeta
            if (sys.freshwater)
                est.tp(i).hbeta     = 0;
                est.tp(i).uhbeta    = 0;
                est.tp(i).uhbeta_l  = 0;
                est.tp(i).uhbeta_u  = 0;
            else
                est.tp(i).hbeta     = q(z(sys.tp(i).iphbeta))*1e6;
                est.tp(i).uhbeta    = ebar(sys.tp(i).iphbeta)*1e6;
                est.tp(i).uhbeta_l  = ebar_l(sys.tp(i).iphbeta)*1e6;
                est.tp(i).uhbeta_u  = ebar_u(sys.tp(i).iphbeta)*1e6;
            end
            est.tp(i).phbeta    = z(sys.tp(i).iphbeta);
            est.tp(i).uphbeta   = sigx(sys.tp(i).iphbeta);

            % pKbeta
            est.tp(i).pKbeta    = z(sys.tp(i).ipKbeta);
            est.tp(i).upKbeta   = sigx(sys.tp(i).ipKbeta);
            if (sys.freshwater)
                est.tp(i).Kbeta     = 0;
                est.tp(i).uKbeta    = 0;
                est.tp(i).uKbeta_l  = 0;
                est.tp(i).uKbeta_u  = 0;
            else
                est.tp(i).Kbeta     = q(z(sys.tp(i).ipKbeta));
                est.tp(i).uKbeta    = ebar(sys.tp(i).ipKbeta);
                est.tp(i).uKbeta_l  = ebar_l(sys.tp(i).ipKbeta);
                est.tp(i).uKbeta_u  = ebar_u(sys.tp(i).ipKbeta);
            end
        end

        % PLEASE ADD NEW VARIABLES HERE AT END (for PrintCSV's sake)
        est.f       = f; % residual f value, from limp
        est.C       = C;
        nr = size(C,1);
        est.yhat    = z(1:nr);
    end
end

% ----------------------------------------------------------------------
% newtn
% ----------------------------------------------------------------------

function [x,J,iflag] = newtn(x0,F,tol)
%[x,J,iflag] = newt(x0,F,tol);
% Simple function that applies Newton's method to find the root of

%  F(x) = 0
% 
% Input: 
% x0: initial iterate
% F: function for which we seek a zero
% tol: convergence criteria ||F(x)||<tol
%
% Output:
% x: steady state solution
% J: Jacobian = [partial F/partial x]
% iflag: 0 ==> Newton's method converged to the desired
%                      tolerance
%        1 ==> Newton's method did not converge 
    iprint = 0; % 0 = off, 1 = on
    MAXIT = 50;
    x = x0;
    if (nargin==4)
        [F0,iJ] = F(x);
    else
        [F0, J] = F(x); 
    end
    iflag = 0; itno = 0;
    while (((norm(F0) > tol) && (itno<MAXIT)) )
        if (nargin==4)            
            dx = -iJ(F0);
            x = x + dx;
            [F0, iJ] = F(x);
        else
            dx = -J\F0;
            x = x + dx;
            [F0, J] = F(x);
        end
        itno = itno+1;
        if (iprint)
            fprintf('%i: norm(F0) = %e norm(F0) = %e \n',itno,norm(F0),norm(F0(1:end-1)));
        end
    end
    if (itno>=MAXIT)
        if norm(F0(1:end-1)) < tol*1e1
            iflag = 2;
            fprintf('Warning Newton''s Method did not converge.\n ')
            fprintf('But value was only 10x more than tolerance...\n')
            fprintf('so it was close, recommend keeping. \n')
        else
            iflag = 1;
            fprintf('Warning Newton''s Method did not converge.\n')
        end
    end
    if (isnan(norm(F0)))
        iflag = 1;
        fprintf('Warning limp is NaN!!!!\n')
    end
    if (nargin==4)
        J = iJ;
    end
end

% ----------------------------------------------------------------------
% PrintCSV, make_headers, parse_CSV
% ----------------------------------------------------------------------

function PrintCSV(est,obs,iflag,opt)
% subfunction of QUODcarb to print results to a CSV
% uses filename opt.fname 
%               opt.fname = 'QUODcarb_output.csv'(default)

    if opt.printcsv == 1
        nD = length(obs);
        for i = 1:nD
            if i == 1
                fname = fopen(opt.fname,'w');
                % make column headers
                make_headers(est(i),opt,fname); 
            end
            % fill one row with data
            parse_CSV(est(i),obs(i),iflag(i),opt,fname);
        end 
    end
end

function make_headers(est,opt,fid)
% Print the column headers  
    nTP = length(est.tp);
    % row 1
    fprintf(fid, '%s, ', '  ');
    fprintf(fid, 'est = output, ');
    fprintf(fid, 'u = 1sigma, ');

    fn = fieldnames(est);
    fnl = length(fn)-3; 
    for i = 5:6:fnl % temperature independent totals
        fprintf(fid, '%s, ', '  '); % est
        fprintf(fid, '%s, ', '  '); % est.u
    end
    for j = 1:nTP
        fnj = fieldnames(est.tp);
        for i = 1:4:8 % T, P
            fprintf(fid, 'tp(%i), ', j); % est
            fprintf(fid, '%s, ', '  ');  % est.u
        end
        for i = 9:6:39 % fco2, pco2, hco3, co3, co2*, ph
            fprintf(fid, 'tp(%i), ', j); % est
            fprintf(fid, '%s, ', '  ');  % est.u
        end
        for i = 45:6:75 % ph_free, ph_tot, ph_sws, ph_nbs, fH, p2f
            fprintf(fid, 'tp(%i), ', j); % est
            fprintf(fid, '%s, ', '  ');  % est.u
        end
        for i = 79:6:length(fnj) % all the rest
            fprintf(fid, 'tp(%i), ', j); % est
            fprintf(fid, '%s, ', '  ');  % est.u
        end
    end

    fprintf(fid,'\n'); % finish first row

    % row 2
    fprintf(fid, '%s, ','iflag');

    fprintf(fid,'est.sal, '); % sal is first
    fprintf(fid,'est.usal, '); % usal

    for i = 5:6:fnl % temperature independent totals
        fprintf(fid,'est.%s, ',fn{i} );
        fprintf(fid,'est.u%s, ',fn{i} ); % u
    end
    for j = 1:nTP
        fnj = fieldnames(est.tp(j));
        for i = 1:4:8 % T, P
            fprintf(fid,'est.%s, ',  fnj{i});
            fprintf(fid,'est.u%s, ', fnj{i}); % u = uncertainty
        end
        for i = 9:6:39 % fco2, pco2, hco3, co3, co2*, ph
            fprintf(fid,'est.%s, ',  fnj{i});
            fprintf(fid,'est.u%s, ', fnj{i}); % u = uncertainty
        end
        for i = 45:6:75 % ph_free, ph_tot, ph_sws, ph_nbs, fH, p2f
            fprintf(fid,'est.%s, ',  fnj{i});
            fprintf(fid,'est.u%s, ', fnj{i}); % u = uncertainty
        end
        for i = 79:6:288 % all the rest of the regulars
            fprintf(fid,'est.%s, ',  fnj{i});
            fprintf(fid,'est.u%s, ', fnj{i}); % u = uncertainty
        end
        n = 0;
        if opt.pKalpha == 1
            for i = 289:6:306
                fprintf(fid,'est.%s, ',  fnj{i});
                fprintf(fid,'est.u%s, ', fnj{i}); % u = uncertainty
            end
            n = 18; 
        end
        if opt.pKbeta == 1
           for i = n+(289:6:306) 
               fprintf(fid,'est.%s, ',  fnj{i});
                fprintf(fid,'est.u%s, ', fnj{i}); % u = uncertainty
           end
        end
        if opt.Revelle == 1
            fprintf(fid,'est.%s, ', fnj{end-1}); % Revelle
            fprintf(fid,'est.%s, ', fnj{end}); % dpfco2dpTA
        end    
    end

    fprintf(fid,'\n'); % finish second row

    % row 3
    fprintf(fid,'%s, ','(0=good)');
       
    fprintf(fid,'%s, ', '(S_P)'); % est.sal
    fprintf(fid,'%s, ', '(S_P)'); % est.usal
    
    for i = 3:6:fnl % temperature independent totals
        fprintf(fid,'%s, ', 'umol/kg');
        fprintf(fid,'%s, ', 'umol/kg');
    end
    for j = 1:nTP
        fprintf(fid, '%s, ', 'deg C'); % est.T
        fprintf(fid, '%s, ', 'deg C');
        fprintf(fid, '%s, ', 'dbar'); % est.P
        fprintf(fid, '%s, ', 'dbar');
        for i = 9:6:15 % fco2, pco2
            fprintf(fid, '%s, ', 'uatm');
            fprintf(fid, '%s, ', 'uatm');
        end
        for i = 21:6:33 % hco3, co3, co2*;
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', 'umol/kg');
        end
        fprintf(fid, '%s, ', '(-log10)'); % ph, log10 unitless
        fprintf(fid, '%s, ', '  ');
        for i = 45:6:63 % ph_free, ph_tot, ph_sws, ph_nbs
            fprintf(fid, '%s, ', '(-log10)');      % log10 unitless
            fprintf(fid, '%s, ', '  ');
        end
        for i = 69:6:75 % fH, p2f aka FugFac
            fprintf(fid, '%s, ', '(unitless)'); 
            fprintf(fid, '%s, ', '  ');
        end
        for i = 79:6:91 % pK0, pK1, pK2
            fprintf(fid, '%s, ', '(-log10)');      % log10 unitless
            fprintf(fid, '%s, ', '  ');
        end
        fprintf(fid, '%s, ', 'umol/kg'); % oh (97)
        fprintf(fid, '%s, ', 'umol/kg');
        fprintf(fid, '%s, ', '(-log10)'); % pKw (103)
        fprintf(fid, '%s, ', '  ');
        for i = 1:3 % (boh3, boh4, pKb)(so4, hso4, pKs)(F, HF, pKf)
            fprintf(fid, '%s, ', 'umol/kg'); % est 1st
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', 'umol/kg'); % est 2nd
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', '(-log10)'); % est.pK
            fprintf(fid, '%s, ', '  ');
        end
        for i = 1:4 % po4, hpo4, h2po4, h3po4
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', 'umol/kg');
        end
        for i = 1:3 % pKp1, pKp2, pKp3
            fprintf(fid, '%s, ', '(-log10)');
            fprintf(fid, '%s, ', '  ');
        end
        for i = 1:3 % (sioh4, siooh3, pKsi)(nh3, nh4, pKnh4)(HS, H2S, pKh2s)
            fprintf(fid, '%s, ', 'umol/kg'); % est 1st
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', 'umol/kg'); % est 2nd
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', '(-log10)'); % est.pK
            fprintf(fid, '%s, ', '  ');
        end
        fprintf(fid, '%s, ', 'umol/kg'); % ca
        fprintf(fid, '%s, ', 'umol/kg');
        for i = 1:2 % (Omega_Ar, pKar)(Omega_Ca, pKca)
            fprintf(fid, '%s, ', 'umol/kg'); % est 1st
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', '(-log10)'); % est.pK
            fprintf(fid, '%s, ', '  ');
        end
        if opt.pKalpha == 1 % (alpha, halpha, pKalpha)
            fprintf(fid, '%s, ', 'umol/kg'); % est 1st
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', 'umol/kg'); % est 2nd
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', '(-log10)'); % est.pK
            fprintf(fid, '%s, ', '  ');
        end
        if opt.pKbeta == 1 % (beta, hbeta, pKbeta)
            fprintf(fid, '%s, ', 'umol/kg'); % est 1st
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', 'umol/kg'); % est 2nd
            fprintf(fid, '%s, ', 'umol/kg');
            fprintf(fid, '%s, ', '(-log10)'); % est.pK
            fprintf(fid, '%s, ', '  ');
        end
        if opt.Revelle == 1
            fprintf(fid, '%s, ', '(unitless)'); % Revelle
            fprintf(fid, '%s, ', '  ');         % dpfco2dpTA
        end
    end

    fprintf(fid,'\n'); % finish third row

end

function parse_CSV(varargin)
% Print the output
    est     = varargin{1}; % posterior marginal precision
    obs     = varargin{2}; % measurement
    iflag   = varargin{3}; % Newton solver convergence flag: 0 converged, 1 not converged
    opt     = varargin{4}; % options
    fid   = varargin{5}; % filename with fopen command

    nTP     = length(est.tp);

    fprintf(fid,'%i, ',iflag);

    % salinity
    fprintf(fid,'%f, ', est.sal);
    fprintf(fid,'%f, ', est.usal);
    % TC (DIC)
    fprintf(fid,'%f, ', est.TC);
    fprintf(fid,'%f, ', est.uTC);
    % TA
    fprintf(fid,'%f, ', est.TA);
    fprintf(fid,'%f, ', est.uTA);
    % TB borate
    fprintf(fid,'%f, ', est.TB);
    fprintf(fid,'%f, ', est.uTB);
    % TS sulfate
    fprintf(fid,'%f, ', est.TS);
    fprintf(fid,'%f, ', est.uTS);
    % TF fluoride
    fprintf(fid,'%f, ', est.TF);
    fprintf(fid,'%f, ', est.uTF);
    % TP phosphate
    fprintf(fid,'%f, ', est.TP);
    fprintf(fid,'%f, ', est.uTP);
    % TSi silicate
    fprintf(fid,'%f, ', est.TSi);
    fprintf(fid,'%f, ', est.uTSi);
    % TNH4 nitrate
    fprintf(fid,'%f, ', est.TNH4);
    fprintf(fid,'%f, ', est.uTNH4);
    % TH2S sulfide
    fprintf(fid,'%f, ', est.TH2S);
    fprintf(fid,'%f, ', est.uTH2S);
    % TCa calcium solubility
    fprintf(fid,'%f, ', est.TCa);
    fprintf(fid,'%f, ', est.uTCa);

    if opt.pKalpha == 1
        % TAlpha
        fprintf(fid,'%f, ', est.TAlpha);
        fprintf(fid,'%f, ', est.uTAlpha);
    end
    if opt.pKbeta == 1
        % TBeta
        fprintf(fid,'%f, ', est.TBeta);
        fprintf(fid,'%f, ', est.uTBeta);
    end

    for j = 1:nTP
        % Temp
        fprintf(fid,'%f, ', est.tp(j).T);
        fprintf(fid,'%f, ', est.tp(j).uT);
        % pres
        fprintf(fid,'%f, ', est.tp(j).P);
        fprintf(fid,'%f, ', est.tp(j).uP);

        % fco2
        fprintf(fid,'%f, ', est.tp(j).fco2);
        fprintf(fid,'%f, ', est.tp(j).ufco2);
        % pco2
        fprintf(fid,'%f, ', est.tp(j).pco2);
        fprintf(fid,'%f, ', est.tp(j).upco2);
        % hco3
        fprintf(fid,'%f, ', est.tp(j).hco3);
        fprintf(fid,'%f, ', est.tp(j).uhco3);
        % co3
        fprintf(fid,'%f, ', est.tp(j).co3);
        fprintf(fid,'%f, ', est.tp(j).uco3);
        % co2* (co2st)
        fprintf(fid,'%f, ', est.tp(j).co2st);
        fprintf(fid,'%f, ', est.tp(j).uco2st);

        % ph
        fprintf(fid,'%f, ', est.tp(j).ph);
        fprintf(fid,'%f, ', est.tp(j).uph);
        % ph_free
        fprintf(fid,'%f, ', est.tp(j).ph_free);
        fprintf(fid,'%f, ', est.tp(j).uph_free);
        % ph_tot
        fprintf(fid,'%f, ', est.tp(j).ph_tot);
        fprintf(fid,'%f, ', est.tp(j).uph);
        % ph_sws
        fprintf(fid,'%f, ', est.tp(j).ph_sws);
        fprintf(fid,'%f, ', est.tp(j).uph);
        % ph_nbs
        fprintf(fid,'%f, ', est.tp(j).ph_nbs);
        fprintf(fid,'%f, ', est.tp(j).uph);

         % fH = activity coefficient
        fprintf(fid,'%f, ', est.tp(j).fH);
        fprintf(fid,'%f, ', est.tp(j).ufH);
        % pp2f
        fprintf(fid,'%f, ', est.tp(j).pp2f);
        fprintf(fid,'%f, ', est.tp(j).upp2f);

        % pK0
        fprintf(fid,'%f, ', est.tp(j).pK0);
        fprintf(fid,'%f, ', est.tp(j).upK0);
        % pK1
        fprintf(fid,'%f, ', est.tp(j).pK1);
        fprintf(fid,'%f, ', est.tp(j).upK1);
        % pK2
        fprintf(fid,'%f, ', est.tp(j).pK2);
        fprintf(fid,'%f, ', est.tp(j).upK2);

        % oh
        fprintf(fid,'%f, ', est.tp(j).oh);
        fprintf(fid,'%f, ', est.tp(j).uoh);
        % pKw = [h][oh]
        fprintf(fid,'%f, ', est.tp(j).pKw);
        fprintf(fid,'%f, ', est.tp(j).upKw);

        % boh4
        fprintf(fid,'%f, ', est.tp(j).boh4);
        fprintf(fid,'%f, ', est.tp(j).uboh4);
        % boh3
        fprintf(fid,'%f, ', est.tp(j).boh3);
        fprintf(fid,'%f, ', est.tp(j).uboh3);
        % pKb = [h][boh4]/[boh3]
        fprintf(fid,'%f, ', est.tp(j).pKb);
        fprintf(fid,'%f, ', est.tp(j).upKb);

        % so4
        fprintf(fid,'%f, ', est.tp(j).so4);
        fprintf(fid,'%f, ', est.tp(j).uso4);
        % hso4
        fprintf(fid,'%f, ', est.tp(j).hso4);
        fprintf(fid,'%f, ', est.tp(j).uhso4);
        % pKs  = [hf][so4]/[hso4]
        fprintf(fid,'%f, ', est.tp(j).pKs);
        fprintf(fid,'%f, ', est.tp(j).upKs);

        % [F]
        fprintf(fid,'%f, ', est.tp(j).F);
        fprintf(fid,'%f, ', est.tp(j).uF);
        % [HF] hydrogen fluoride
        fprintf(fid,'%f, ', est.tp(j).HF);
        fprintf(fid,'%f, ', est.tp(j).uHF);
        % pKf = [h][F]/[HF]
        fprintf(fid,'%f, ', est.tp(j).pKf);
        fprintf(fid,'%f, ', est.tp(j).upKf);

        % po4
        fprintf(fid,'%f, ', est.tp(j).po4);
        fprintf(fid,'%f, ', est.tp(j).upo4);
        % hpo4
        fprintf(fid,'%f, ', est.tp(j).hpo4);
        fprintf(fid,'%f, ', est.tp(j).uhpo4);
        % h2po4
        fprintf(fid,'%f, ', est.tp(j).h2po4);
        fprintf(fid,'%f, ', est.tp(j).uh2po4);
        % h3po4
        fprintf(fid,'%f, ', est.tp(j).h3po4);
        fprintf(fid,'%f, ', est.tp(j).uh3po4);
        % pKp1 = [h][h2po4]/[h3po4]
        fprintf(fid,'%f, ', est.tp(j).pKp1);
        fprintf(fid,'%f, ', est.tp(j).upKp1);
        % pKp2 = [h][hpo4]/[h2po4]
        fprintf(fid,'%f, ', est.tp(j).pKp2);
        fprintf(fid,'%f, ', est.tp(j).upKp2);
        % pKp3 = [h][po4]/[hpo4]
        fprintf(fid,'%f, ', est.tp(j).pKp3);
        fprintf(fid,'%f, ', est.tp(j).upKp3);

        % sioh4
        fprintf(fid,'%f, ', est.tp(j).sioh4);
        fprintf(fid,'%f, ', est.tp(j).usioh4);
        % siooh3
        fprintf(fid,'%f, ', est.tp(j).siooh3);
        fprintf(fid,'%f, ', est.tp(j).usiooh3);
        % pKSi = [h][siooh3]/[sioh4]
        fprintf(fid,'%f, ', est.tp(j).pKsi);
        fprintf(fid,'%f, ', est.tp(j).upKsi);

        % nh3
        fprintf(fid,'%f, ', est.tp(j).nh3);
        fprintf(fid,'%f, ', est.tp(j).unh3);
        % nh4
        fprintf(fid,'%f, ', est.tp(j).nh4);
        fprintf(fid,'%f, ', est.tp(j).unh4);
        % pKnh4 = [h][nh3]/[nh4]
        fprintf(fid,'%f, ', est.tp(j).pKnh4);
        fprintf(fid,'%f, ', est.tp(j).upKnh4);

        % hs
        fprintf(fid,'%f, ', est.tp(j).HS);
        fprintf(fid,'%f, ', est.tp(j).uHS);
        % h2s
        fprintf(fid,'%f, ', est.tp(j).H2S);
        fprintf(fid,'%f, ', est.tp(j).uH2S);
        % pKh2s = [h][hs]/[h2s]
        fprintf(fid,'%f, ', est.tp(j).pKh2s);
        fprintf(fid,'%f, ', est.tp(j).upKh2s);

        % ca
        fprintf(fid,'%f, ', est.tp(j).ca);
        fprintf(fid,'%f, ', est.tp(j).uca);
        % OmegaAr
        fprintf(fid,'%f, ', est.tp(j).OmegaAr);
        fprintf(fid,'%f, ', est.tp(j).uOmegaAr);
        % OmegaCa
        fprintf(fid,'%f, ', est.tp(j).OmegaCa);
        fprintf(fid,'%f, ', est.tp(j).uOmegaCa);
        % pKar = [ca][co3]/[omegaAr]
        fprintf(fid,'%f, ', est.tp(j).pKar);
        fprintf(fid,'%f, ', est.tp(j).upKar);
        % pKca = [ca][co3]/[omegaCa]
        fprintf(fid,'%f, ', est.tp(j).pKca); 
        fprintf(fid,'%f, ', est.tp(j).upKca);

        if opt.pKalpha == 1
            % alpha
            fprintf(fid,'%f, ', est.tp(j).alpha);
            fprintf(fid,'%f, ', est.tp(j).ualpha);
            % halpha
            fprintf(fid,'%f, ', est.tp(j).halpha);
            fprintf(fid,'%f, ', est.tp(j).uhalpha);
            % pKalpha = [h][alpha]/[halpha]
            fprintf(fid,'%f, ', est.tp(j).pKalpha);
            fprintf(fid,'%f, ', est.tp(j).upKalpha);
        end
        if opt.pKbeta == 1
            % beta
            fprintf(fid,'%f, ', est.tp(j).beta);
            fprintf(fid,'%f, ', est.tp(j).ubeta);
            % hbeta
            fprintf(fid,'%f, ', est.tp(j).hbeta);
            fprintf(fid,'%f, ', est.tp(j).uhbeta);
            % pKbeta = [h][beta]/[hbeta]
            fprintf(fid,'%f, ', est.tp(j).pKbeta);
            fprintf(fid,'%f, ', est.tp(j).upKbeta);
        end
        if opt.Revelle == 1
            fprintf(fid,'%f, ', est.tp(j).Revelle);
            fprintf(fid,'%f, ', est.tp(j).dpfco2dpTA);
        end
    end

fprintf(fid,'\n'); % end row of data
out = [];

end % function




