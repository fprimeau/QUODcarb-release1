
% QUODcarb example 3
% example to show more syntax and possibilities

load fakedata.mat; 
% variable name is 'fakedata'
% sal  in fakedata(1,:);
% TC   in fakedata(2,:);
% TA   in fakedata(3,:);
% ph   in fakedata(4,:); at 25 degC
% pco2 in fakedata(5,:); at 20 degC
% co3  in fakedata(6,:); at 25 degC
% TP   in fakedata(7,:);
% TSi  in fakedata(8,:);
% T insitu in fakedata(9,:);
% P insitu in fakedata(10,:);

% populate opt structure
opt.K1K2 = 10; % option for K1K2 formulation
opt.KSO4 = 1;   % option for KSO4 formulation
opt.KF   = 2;   % option for KF formulation
opt.TB   = 2;   % option for TB formulation
opt.phscale     = 1; % 1 = tot, 2 = sws, 3 = free, 4 = nbs
opt.printcsv    = 1; % print output to CSV? (1=on, 0=off)
opt.fname       = 'example3.csv'; % not used if opt.printcsv = 0
opt.printmes    = 0; % print screen messages? (1=on, 0=off)
opt.co2press    = 1; % pressure correction for p2f and K0 (1=on, 0=off)
opt.Revelle     = 1; % calculate Revelle factor? (1=on, 0=off)

q  = @(x) 10^(-x);
nD = length(fakedata); % number of datapoints

for i = 1:nD
    % temperature and pressure INdependent
    obs(i).sal  = fakedata(1,i); % salinity
    obs(i).usal = 0.001; % salinity uncertainty, 1 sigma
    obs(i).TC   = fakedata(2,i); % TC (umol/kg)
    obs(i).uTC  = 2.00; % TC uncertainty, 1 sigma
    obs(i).TA   = fakedata(3,i); % TA (umol/kg)
    obs(i).uTA  = 2.00; % TA uncertainty, 1 sigma
    
    obs(i).TP   = fakedata(7,i); % total phosphate (umol/kg)
    obs(i).uTP  = 0.01*fakedata(7,i); % 1% uncertainty TP, 1 sigma
    obs(i).TSi  = fakedata(8,i); % total silicate (umol/kg)
    obs(i).uTSi = 0.01*fakedata(8,i); % 1% uncertainty TSi, 1 sigma

    % temperature and pressure dependent 1st system, in-situ
    obs(i).tp(1).T  = fakedata(9,i); % in-situ temperature (Celsius)
    obs(i).tp(1).uT = 0.01; % in-situ temperature uncertainty, 1 sigma
    obs(i).tp(1).P  = fakedata(10,i); % in-situ pressure (dbar)
    obs(i).tp(1).uP = 0.63; % in-situ pressure uncertainty (dbar), 1 sigma
    obs(i).tp(1).upK1 = 0.01; % want to reset the upK1 to 0.01;

    % temperature and pressure dependent 2nd system, for ph & co3 25 degC
    obs(i).tp(2).T  = 25.0; % temperature for ph in Celsius
    obs(i).tp(2).uT = 0.01; % ph temperature uncertainty, 1 sigma
    obs(i).tp(2).P  = 0.0; % pressure (dbar)
    obs(i).tp(2).uP = 0.07; % pressure uncertainty, 1 sigma
    obs(i).tp(2).ph     = fakedata(4,i); % ph measured on total scale
    obs(i).tp(2).uph    = 0.01; % ph uncertainty, 1 sigma
    obs(i).tp(2).co3    = fakedata(6,i); % total carbonate ion (umol/kg)
    obs(i).tp(2).uco3   = 0.02*fakedata(6,i); % 2% CO3 uncertainty (umol/kg), 1 sigma
    obs(i).tp(2).upK1   = 0.01; % want to reset the upK1 to 0.01;

    % temperature and pressure dependent 3rd system, for pco2 20 degC
    obs(i).tp(3).T  = 20.0; % temperature for pCO2 in Celsius
    obs(i).tp(3).uT = 0.02; % pCO2 temperature uncertainty, 1 sigma
    obs(i).tp(3).P  = 0.0; % pressure in dbar for pCO2
    obs(i).tp(3).uP = 0.03; % dbar, pressure for pCO2 uncertainty, 1 sigma
    obs(i).tp(3).pco2   = fakedata(5,i); % partial pressure CO2 (pCO2) (uatm)
    obs(i).tp(3).upco2  = 0.01*fakedata(5,i); % 1% pCO2 uncertainty (uatm), 1 sigma
    obs(i).tp(3).upK1   = 0.01; % want to reset the upK1 to 0.01;

    % temperature and pressure dependend 4th system
    obs(i).tp(4).T  = 15.0; % want an output at 15 deg C
    obs(i).tp(4).uT = 0.01; 
    obs(i).tp(4).P  = 10.0; % want output at 10 dbar
    obs(i).tp(4).uP = 0.03;
    obs(i).tp(4).upK1 = 0.01; % want to reset the upK1 to 0.01;

end

[est,obs,sys,iflag,opt] = QUODcarb(obs,opt);

% save example3.mat est; % save output estimate as a mat file

% output est(1) values should be:
% est(1).TC         = 2194.8;       % posterior TC
% est(1).TA         = 2308.3;       % posterior TA
% est(1).tp(1).T    = 10.0932;      % in-situ temperature
% est(1).tp(1).P    = 650.96;       % in-situ pressure
% est(1).tp(2).T    = 25.0000;      % temperature at which ph and co3 were measured
% est(1).tp(2).ph   = 7.6378;       % posterior ph
% est(1).tp(2).co3  = 99.2912;      % posterior co3
% est(1).tp(3).T    = 19.9983;      % temperature at which pco2 was measured
% est(1).tp(3).pco2 = 938.3595;     % posterior pco2
% est(1).tp(4).T    = 15;           % output system @15degC
% est(1).tp(4).P    = 10;
% est(1).tp(2).OmegaAr = 1.5735;    
% est(1).tp(2).OmegaCa = 2.3864;    
% est(1).tp(2).Revelle = 14.4902;   
% est(1).f          = 13.6064;      % f(yhat) in paper, numerical measure 
                                    % of internal consistency 

% calculate SIR: HCO3/H_free at tp(1) and est(1)
nr = size(est(1).C,1);
v = zeros(nr,1);
v( [ sys.tp(1).iphco3, sys.tp(1).iph_free ] )  = [ 1; -1 ];
sig = sqrt( v' * est(1).C * v);
% crank 1:
SIR = q(v' * est(1).yhat);
% crank 2:
SIRu = q( v' * est(1).yhat - sig );
% crank 3:
SIRl = q( v' * est(1).yhat + sig );
uSIR = 0.5*(SIRu-SIRl);





