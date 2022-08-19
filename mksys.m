function sys = mksys(sys)
%
% Private function for QUODcarb.m
% it creates the K amd M matrices and if necessary the free2tot function
%
% utility functions and constants
    LOG10 = log(10);
    p = @(x) -log10( x );  % inverse of q    
    q = @(x)  10.^( -x );  % inverse of p
    sys.p = p;
    sys.q = q;
    dqdx = @(x) - LOG10 * 10.^( -x );  % q'
    d2qdx2 = @(x) LOG10^2 * 10.^(-x ); % q"
    dpdx = @(x) -1 / (x .* LOG10); % p'
    d2pdx2 = @(x) 1 / (x.^2 * LOG10); % p"

    sys.dqdx = dqdx;
    sys.dpdx = dpdx;
    sys.d2pdx2 = d2pdx2;
    sys.d2qdx2 = d2qdx2;
    % list of acid-base reactions
    if ismember('carbonate',sys.abr)
        % K0 = [co2st]/fco2   
        % K1 = [h][hco3]/[co2st]
        % K2 = [h][co3]/[hco3]
        carbonate =  { 'T','S','P','K0','K1','K2','TC','TA','fco2', ...
                        'co2st', 'hco3', 'co3', 'h', 'p2f', 'pco2'};

        sys.variables = carbonate;
        nrk = 3; % number of rows
        if ismember('water',sys.abr)
            % Kw = [h][oh]
            water = {'Kw','oh' };
            nrk = nrk+1; % number of rows
            sys.variables = {sys.variables{:},water{:}};
        end
    end
    
    if ismember('borate',sys.abr)
        % Kb = [h][boh4]/[boh3]
        borate = {'Kb','TB', 'boh4', 'boh3'};
        nrk = nrk+1; % number of rows
        sys.variables = {sys.variables{:},borate{:}};
    end
    
    if ismember('sulfate',sys.abr)
        % Ks  = [hf][so4]/[hso4]
        sulfate = {'Ks','TS','hf', 'so4','hso4'};
        nrk = nrk+1; % number of rows
        sys.variables = {sys.variables{:},sulfate{:}};
    end
    
    if ismember('fluoride',sys.abr)
        % KF = [h][F]/[HF]
        fluoride = {'KF','TF', 'F', 'HF'};
        nrk = nrk+1; % number of rows
        sys.variables = {sys.variables{:},fluoride{:}};
    end
    
    if (ismember('phosphate',sys.abr)|ismember('K1p',sys.abr)|ismember('K2p',sys.abr)|ismember('K3p',sys.abr))
        % K1p = [h][h2po4]/[h3po4]
        % K2p = [h][hpo4]/[h2po4]
        % K3p = [h][po4]/[hpo4]
        phosphate = {'K1p','K2p','K3p','TP', 'h3po4', 'h2po4', 'hpo4', 'po4' };
        nrk = nrk+3; % number of rows
        sys.variables = {sys.variables{:},phosphate{:}};
    end
    
    if ismember('silicate',sys.abr)
        % KSi = [h][siooh3]/[sioh4]
        silicate = {'Ksi','TSi', 'siooh3', 'sioh4'};
        nrk = nrk+1; % number of rows
        sys.variables = {sys.variables{:},silicate{:}};
    end
    
    if ismember('ammonia',sys.abr)
        % Knh4 = [h][nh3]/[nh4]
        ammonia = {'Knh4','TNH3', 'nh3', 'nh4'};
        nrk = nrk+1; % number of rows
        sys.variables = {sys.variables{:},ammonia{:}};
    end
    
    if ismember('sulfide',sys.abr)
        % Kh2s = [h][hs]/[h2s]
        sulfide = {'Kh2s','TH2S','hs', 'h2s' };
        nrk = nrk+1; % number of rows
        sys.variables = {sys.variables{:},sulfide{:}};
    end
    %
    % assign integer indices to each species
    %
    nc = length(sys.variables); % number of columns
    for i = 1:nc
        cmd = sprintf('i%s = %i;',sys.variables{i},i);
        eval(cmd);
        cmd = sprintf('sys.i%s = %i;',sys.variables{i},i);
        eval(cmd);
    end

    
    %
    K = zeros(2*nrk,nc);
    %
    row = 0;
    if (ismember('K0',sys.variables))        % K0 = [CO2*]/fCO2 ==> -pK0 + pco2st - pfco2 = 0         
        row = row+1;
        K(row,[iK0,ico2st,ifco2]) = [-1, 1, -1];
        sys.system{row} = 'K0';
    
        row = row+1;
        jK0 = row;
        sys.jK0 = jK0;
        K(jK0,iK0) = 1; % K0
        sys.system{jK0} = 'pK0(T,S,P)';  
    end
    
    if (ismember('K1',sys.variables))        % K1 = [HCO3][H]/[CO2*] ==> -pK1 + phco3 + ph - pfco2 = 0
        row = row+1;
        K(row,[iK1,ih,ihco3,ico2st]) = [-1, 1, 1, -1];
        sys.system{row} = 'K1';
    
        row = row+1;
        jK1 = row;
        sys.jK1 = jK1;
        K(jK1,iK1) = 1; % K1
        sys.system{jK1} = 'pK1(T,S,P)';

    end
    
    if (ismember('K2',sys.variables))
        % K2 = [CO3][H]/[HCO3]
        
        row = row+1;
        K(row,[iK2,ih,ico3,ihco3]) = [-1, 1, 1, -1];
        sys.system{row} = 'K2';
    
        row = row+1;
        jK2 = row;
        sys.jK2 = jK2;
        K(jK2,iK2) = 1; % K2
        sys.system{jK2} = 'pK2(T,S,P)';

    end

    if (ismember('p2f', sys.variables))
        %fco2 = pco2*p2f;

        row = row+1;
        K(row,[ifco2,ipco2,ip2f]) = [1 -1 -1];
        sys.system{row} = 'p2f';
        
        row = row+1;
        jp2f = row;
        sys.jp2f = jp2f;
        K(jp2f,ip2f) = 1; % p2f
        sys.system{jp2f} = 'p2f(T,S,P)';
    
    end
    
    if (ismember('Kw',sys.variables))
        % Kw = [OH][H]

        row = row+1;
        K(row,[iKw,ih,ioh]) = [-1, 1, 1];
        sys.system{row} = 'Kw';
        
        row = row+1;
        jKw = row;
        sys.jKw = jKw;
        K(jKw,iKw) = 1; % Kw
        sys.system{jKw} = 'pKw(T,S,P)';

    end
    
    if (ismember('Kb',sys.variables))
        % Kb = [H][BOH4]/[BOH3]

        row = row+1;
        K(row,[iKb,ih,iboh4,iboh3]) = [-1, 1, 1, -1];
        sys.system{row} = 'Kb';
        
        row = row+1;
        jKb = row;
        sys.jKb = jKb;
        K(jKb,iKb) = 1; % Kb
        sys.system{jKb} = 'pKb(T,S,P)';

    
    end
    
    if (ismember('Ks',sys.variables))
        % KS  = [H]F[SO4]/[HSO4]
        
        row = row+1;
        K(row,[iKs,ihf,iso4,ihso4]) = [-1, 1, 1, -1];
        sys.system{row} = 'Ks';
    
        row = row+1;
        jKs = row;
        sys.jKs = jKs;
        K(jKs,iKs) = 1; % Ks
        sys.system{jKs} = 'pKs(T,S,P)';

    end
    
    if (ismember('KF',sys.variables))
        % KF = [H][F]/[HF]

        row = row+1;
        K(row,[iKF, ih, iF, iHF]) = [-1, 1, 1, -1];
        sys.system{row} = 'KF';

        row = row+1;
        jKF = row;
        sys.jKF = jKF;
        K(jKF,iKF) = 1; % KF
        sys.system{jKF} = 'pKF(T,S,P)';

    end
    
    if (ismember('K1p',sys.variables))
        % K1p = [H][H2PO4]/[H3PO4]

        row = row+1;
        K(row,[iK1p,ih,ih2po4,ih3po4]) = [-1, 1, 1, -1];
        sys.system{row} = 'K1p';

        row = row+1;
        jK1p = row;
        sys.jK1p = jK1p;
        K(jK1p,iK1p) = 1; % K1p
        sys.system{jK1p} = 'pK1p(T,S,P)';

    end
    if (ismember('K2p',sys.variables))
        % K2p = [H][HPO4]/[H2PO4]
        K(row,[iK2p,ih,ihpo4,ih2po4]) = [-1, 1, 1, -1];
        sys.system{row} = 'K2p';
        row = row+1;
        jK2p = row;
        sys.jK2p = jK2p;
        K(jK2p,iK2p) = 1; % K2p
        sys.system{jK2p} = 'pK2p(T,S,P)';

    end

    if (ismember('K3p',sys.variables))
        % K3p = [H][PO4]/[HPO4]
        
        row = row+1;
        K(row,[iK3p,ih,ipo4,ihpo4]) = [-1, 1, 1, -1];
        sys.system{row} = 'K3p';

        row = row+1;
        jK3p = row;
        sys.jK3p = jK3p;
        K(jK3p,iK3p) = 1; % K3p
        sys.system{jK3p} = 'pK3p(T,S,P)';

    end

    if (ismember('Ksi',sys.variables))
        % KSi = [H][SiO(OH)3]/[Si(OH)4]
        
        row = row+1;
        K(row,[iKsi,ih,isiooh3,isioh4]) = [-1, 1, 1, -1];
        sys.system{row} = 'Ksi';

        row = row+1;
        jKsi = row;
        sys.jKsi = jKsi;
        K(jKsi,iKsi) = 1; % Ksi
        sys.system{jKsi} = 'pKSi(T,S,P)';
    
    end

    if (ismember('Knh4',sys.variables))
        % Knh4 = [H][NH3]/[NH4+]

        row = row+1;
        K(row,[iKnh4,ih,inh3,inh4]) = [-1, 1, 1, -1];
        sys.system{row} = 'Knh4';

        row = row+1;
        jKnh4 = row;
        sys.jKnh4 = jKnh4;
        K(jKnh4,iKnh4) = 1; % Knh4
        sys.system{jKnh4} = 'pKnh4(T,S,P)';

    end

    if (ismember('Kh2s',sys.variables))
        % Kh2s = [H][HS]/[H2S]

        row = row+1;
        K(row,[iKh2s,ih,ihs,ih2s]) = [-1, 1, 1, -1];
        sys.system{row} = 'Kh2s';

        row = row+1;
        jKh2s = row;
        sys.jKh2s = jKh2s;
        K(jKh2s,iKh2s) = 1; % Kh2s
        sys.system{jKh2s} = 'pKh2s(T,S,P)';
    end

    % "mass" conservation" equations
    nr = 2; %TC and TA
    if (ismember('TB',sys.variables));
        nr = nr + 1;
    end
    if (ismember('TS',sys.variables));
        nr = nr + 1; % add an extra equation relating [H]F to [H]
    end
    if (ismember('TF',sys.variables));
        nr = nr + 1;
    end
    if (ismember('TP',sys.variables));
        nr = nr + 1;
    end
    if (ismember('TSi',sys.variables));
        nr = nr + 1;
    end
    if (ismember('NH3T',sys.variables));
        nr = nr + 1;
    end
    if (ismember('H2ST',sys.variables));
        nr = nr + 1;
    end
    
    M = zeros(nr,nc);

    % Total carbonate: TC - [CO2*] - [HCO3] - [CO3] = 0
    M(1,[iTC,ico2st,ihco3,ico3])    =  [1, -1, -1, -1];
    sys.mass{1} = 'carbon';
    
    % Total alkalinity: TA - [HCO3] - 2[CO3]  ...
    M(2,[iTA,ihco3,ico3])    =  [1, -1, -2];
    sys.mass{2} = 'alkalinity';
    if (ismember('Kw',sys.variables));
        M(2,ioh) = -1;
    end
    row = 3;
    
    % Total borate
    if (ismember('TB',sys.variables));
        sys.mass{row} = 'boron';
        M(row,[iTB,iboh3,iboh4])   =  [1, -1, -1];
        M(2,iboh4) = -1; % alk
        row = row+1;
    end
    
    % Total sulfate
    if (ismember('TS',sys.variables))
        sys.mass{row} = 'sulfur';
        M(row,[iTS,ihso4,iso4])   =  [1, -1, -1];
        M(2,[ihf,ihso4])   =  [1, 1]; % alk
        sys.f2t = @(z) z(ihf)  + p( q( z(iKs) ) + q( z(iTS) ) ) - z(iKs) - z(ih);
        f2t_phf = @(z) 1;
        f2t_ph = @(z) -1;
        f2t_pTS = @(z) dpdx( q( z(iKs) ) + q( z(iTS) ) ) * dqdx( z(iTS) );
        f2t_pKs = @(z) dpdx( q( z(iKs) ) + q( z(iTS) ) ) * dqdx( z(iKs) ) - 1;

        f2t_2pKs = @(z) dpdx( q( z(iKs) ) + q( z(iTS) ) ) * d2qdx2( z(iKs) ) + ...
            d2pdx2( q( z(iKs) ) + q( z(iTS) ) ) * dqdx( z(iKs) ).^2;
        
        f2t_2pTS = @(z) dpdx( q( z(iKs) ) + q( z(iTS) ) ) * d2qdx2( z(iTS) ) + ...
            d2pdx2( q( z(iKs) ) + q( z(iTS) ) ) * dqdx( z(iTS) ).^2;
        
        f2t_pTS_pKs = @(z) d2pdx2( q( z(iKs) ) + q( z(iTS) ) ) * dqdx( z(iTS) ) * dqdx( z(iKs) );

        sys.gf2t = @(z) [ f2t_ph(z), f2t_pKs(z), f2t_pTS(z), f2t_phf(z)];
        sys.ggf2t = @(z) [ [f2t_2pKs(z), f2t_pTS_pKs(z)]; ...
                           [f2t_pTS_pKs(z), f2t_2pTS(z)]  ];
        row = row+1;
    end

    % Total fluoride
    if (ismember('TF',sys.variables))
        sys.mass{row} = 'fluorine';
        M(row,[iTF,iHF,iF]) =  [1, -1, -1];
        M(2,iHF) =  1;  % alk
        row = row+1;
    end

    % Total phosphate
    if (ismember('TP',sys.variables))
        sys.mass{row} = 'phosphorous';             
        M(row,[iTP,ih3po4,ih2po4,ihpo4,ipo4])    =  [1, -1, -1, -1, -1];
        M(2,[ihpo4,ipo4,ih3po4])  = [-1, -2, 1]; % alk
        row = row+1;
    end

    % Total silicate
    if (ismember('TSi',sys.variables))
        sys.mass{row} = 'silicon';
        M(row,[iTSi, isioh4,isiooh3]) = [1, -1, -1];
        M(2,isiooh3) = -1; % alk
        row = row+1;
    end

    % Total ammonia
    if (ismember('TNH3',sys.variables))
        sys.mass{row} = 'ammonia';
        M(row,[iTNH3,inh4,inh3]) = [1, -1, -1];
        M(2,inh3) = -1; % alk
        row = row+1;
    end

    % Total sulfide
    if (ismember('TH2S',sys.variables))
        sys.mass{row} = 'sulfide';
        M(row,[iTH2S,ih2s,ihs]) = [1, -1, -1];
        M(2,ihs) = -1; %alk
        row = row+1;
    end
    
    sys.M = M;
    sys.K = K;
   
end
