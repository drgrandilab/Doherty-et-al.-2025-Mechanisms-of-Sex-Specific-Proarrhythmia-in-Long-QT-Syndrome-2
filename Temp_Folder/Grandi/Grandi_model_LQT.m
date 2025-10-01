%% Grandi model - ODE file

% 2019-Dec-02: corrected equation for J_CaB_cytosol
% 2020-Aug-25: corrected NCX

function output = Grandi_model_LQT(t, y, p, runType)

%% State variables
%   1       2       3       4       5       6       7       8       9       10      11      12      13
%   m       h       j       d       f       fcaBj   fcaBsl  xtos   ytos     xtof    ytof    xkr     xks   
%   14      15      16      17      18      19      20      21      22      23      24
%   25      26      27      28      29      30
%   SRB     SLLj   SLLsl    SLHj    SLHsl  Csqnb
%    31      32      33      34      35      36      37     38     39    
%    Ca_sr   Naj     Nasl    Nai     K  i      Caj    Casl    Cai    Vm
%   40   41
% rtos  ?
%CaMKt=y(58)
%mL=y(61)
%hL=y(60)
%hLp=y(59)

ydot = zeros(size(y));

% Mutation ICaL
mutation_flag_LQT8 = p(3);
mutation_change_GCaL_LQT8 = p(4);
mutation_change_vCaL_LQT8 = p(5);
mutation_flag_LQT2 = p(13);
mutation_change_GKr = p(14);
mutation_flag_LQT3 = p(15);
mutation_change_GNaL = p(16);

% Perturbations Female
female_flag = p(6);
female_change_Gto = p(7);
female_change_GKr = p(8);
female_change_GKs = p(9);
female_change_GK1 = p(10);
female_change_vNCX = p(11);
female_change_vPMCA = p(12);

%% Current multipliers
par_SA = p(17:end);

INa_Multiplier = par_SA(1);
ICaL_Multiplier = par_SA(2)*(1+mutation_change_GCaL_LQT8*mutation_flag_LQT8);
Itof_Multiplier = par_SA(3)*(1+female_change_Gto*female_flag);
Itos_Multiplier = par_SA(4)*(1+female_change_Gto*female_flag);
IKr_Multiplier = par_SA(5)*(1+female_change_GKr*female_flag)*(1-mutation_change_GKr*mutation_flag_LQT2);
IKs_Multiplier = par_SA(6)*(1+female_change_GKs*female_flag);
IKp_Multiplier = par_SA(7);
IK1_Multiplier = par_SA(8)*(1+female_change_GK1*female_flag);
IClCa_Multiplier = par_SA(9);
IClB_Multiplier = par_SA(10);
INaB_Multiplier = par_SA(11);
ICaB_Multiplier = par_SA(12);
INaK_Multiplier = par_SA(13);
INaCa_Multiplier = par_SA(14)*(1+female_change_vNCX*female_flag);
IpmCa_Multiplier = par_SA(15)*(1+female_change_vPMCA*female_flag);
Jup_Multiplier = par_SA(16);
Jrel_Multiplier = par_SA(17);
Jleak_Multiplier = par_SA(18);
INaL_Multiplier = par_SA(19)*(1+mutation_change_GNaL*mutation_flag_LQT3);
%% Input Parameters
% EPI with 1, ENDO otherwise
epi = p(1);
celltype=p(1);
% Stimulation period
stimPeriod = p(2);

%% Model Parameters
% Constants
R = 8314;       % [J/kmol*K]  
Frdy = 96485;   % [C/mol]  
Temp = 310;     % [K]
FoRT = Frdy/R/Temp;
Cmem = 1.3810e-10;   % [F] membrane capacitance
Qpow = (Temp-310)/10;

% Cell geometry
cellLength = 100;     % cell length [um]
cellRadius = 10.25;   % cell radius [um]
junctionLength = 160e-3;  % junc length [um]
junctionRadius = 15e-3;   % junc radius [um]
distSLcyto = 0.45;    % dist. SL to cytosol [um]
distJuncSL = 0.5;  % dist. junc to SL [um]
DcaJuncSL = 1.64e-6;  % Dca junc to SL [cm^2/sec]
DcaSLcyto = 1.22e-6; % Dca SL to cyto [cm^2/sec]
DnaJuncSL = 1.09e-5;  % Dna junc to SL [cm^2/sec]
DnaSLcyto = 1.79e-5;  % Dna SL to cyto [cm^2/sec] 
Vcell = pi*cellRadius^2*cellLength*1e-15;    % [L]
Vmyo = 0.65*Vcell; Vsr = 0.035*Vcell; Vsl = 0.02*Vcell; Vjunc = 0.0539*.01*Vcell; 
SAjunc = 20150*pi*2*junctionLength*junctionRadius;  % [um^2]
SAsl = pi*2*cellRadius*cellLength;          % [um^2]
%J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  % [L/msec] = 2.3056e-11
%J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 1.1074e-13
%J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 1.5714e-12
%J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;% [L/msec] = 7.36e-13
%J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  % [L/msec] = 2.3056e-11
% tau's from c-code, not used here
J_ca_juncsl =1/1.2134e12; % [L/msec] = 8.2413e-13
J_ca_slmyo = 1/2.68510e11; % [L/msec] = 3.2743e-12
J_na_juncsl = 1/(1.6382e12/3*100); % [L/msec] = 6.1043e-13
J_na_slmyo = 1/(1.8308e10/3*100);  % [L/msec] = 5.4621e-11

% Fractional currents in compartments
Fjunc = 0.11;   Fsl = 1-Fjunc;
Fjunc_CaL = 0.9; Fsl_CaL = 1-Fjunc_CaL;

% Fixed ion concentrations     
Cli = 15;   % Intracellular Cl  [mM]
Clo = 150;  % Extracellular Cl  [mM]
Ko = 5.4;   % Extracellular K   [mM]
Nao = 140;  % Extracellular Na  [mM]
Cao = 1.8;  % Extracellular Ca  [mM]1.8
Mgi = 1;    % Intracellular Mg  [mM]

% Nernst Potentials
ena_junc = (1/FoRT)*log(Nao/y(32));     % [mV]
ena_sl = (1/FoRT)*log(Nao/y(33));       % [mV]
ENa= (1/FoRT)*log(Nao/y(34));       % [mV]
% ena_junc = (1/FoRT)*log(Nao/7.65);     % [mV]
% ena_sl = (1/FoRT)*log(Nao/7.65);       % [mV]
ek = (1/FoRT)*log(Ko/y(35));	        % [mV]
eca_junc = (1/FoRT/2)*log(Cao/y(36));   % [mV]
eca_sl = (1/FoRT/2)*log(Cao/y(37));     % [mV]
ecl = (1/FoRT)*log(Cli/Clo);            % [mV]

%% Na transport parameters
GNa = INa_Multiplier*23;
GNaB = INaB_Multiplier*0.597e-3;    % [mS/uF] 0.897e-3
IbarNaK = INaK_Multiplier*1.8;  %1.90719;     % [uA/uF]
KmNaip = 11;        % [mM]
KmKo = 1.5;         % [mM]
Q10NaK = 1.63;  
Q10KmNai = 1.39;

%% K current parameters
pNaK = 0.01833;      
gkp = IKp_Multiplier*(2*0.001);

%% Cl current parameters
GClCa = IClCa_Multiplier*(0.5*0.109625);   % [mS/uF]
GClB = IClB_Multiplier*9e-3;        % [mS/uF]
KdClCa = 100e-3;    % [mM]

%% ICaL parameters
pNa = ICaL_Multiplier*0.50*1.5e-8;       % [cm/sec]
pCa = ICaL_Multiplier*0.50*5.4e-4;       % [cm/sec]
pK = ICaL_Multiplier*0.50*2.7e-7;        % [cm/sec]
Q10CaL = 1.8;       

%% Ca transport parameters
IbarNCX = INaCa_Multiplier*4.5;      % [uA/uF]5.5 before - 9 in rabbit
KmCai = 3.59e-3;    % [mM]
KmCao = 1.3;        % [mM]
KmNai = 12.29;      % [mM]
KmNao = 87.5;       % [mM]
ksat = 0.32;        % [none]  
nu = 0.27;          % [none]
Kdact = 0.150e-3;   % [mM] 
Q10NCX = 1.57;      % [none]
IbarSLCaP = IpmCa_Multiplier*0.0673; % IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
KmPCa = 0.5e-3;     % [mM] 
GCaB = ICaB_Multiplier*5.513e-4;    % [uA/uF] 3
Q10SLCaP = 2.35;    % [none]

% SR flux parameters
Q10SRCaP = 2.6;          % [none]
Vmax_SRCaP = Jup_Multiplier*5.3114e-3;  % [mM/msec] (286 umol/L cytosol/sec)
Kmf = 0.246e-3;          % [mM] default
%Kmf = 0.175e-3;          % [mM]
Kmr = 1.7;               % [mM]L cytosol
hillSRCaP = 1.787;       % [mM]
ks = Jrel_Multiplier*25;                 % [1/ms]      
koCa = 10;               % [mM^-2 1/ms]   %default 10   modified 20
kom = 0.06;              % [1/ms]     
kiCa = 0.5;              % [1/mM/ms]
kim = 0.005;             % [1/ms]
ec50SR = 0.45;           % [mM]

% Buffering parameters
% koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
Bmax_Naj = 7.561;       % [mM] % Bmax_Naj = 3.7; (c-code difference?)  % Na buffering
Bmax_Nasl = 1.65;       % [mM]
koff_na = 1e-3;         % [1/ms]
kon_na = 0.1e-3;        % [1/mM/ms]
Bmax_TnClow = 70e-3;    % [mM]                      % TnC low affinity
koff_tncl = 19.6e-3;    % [1/ms] 
kon_tncl = 32.7;        % [1/mM/ms]
Bmax_TnChigh = 140e-3;  % [mM]                      % TnC high affinity 
koff_tnchca = 0.032e-3; % [1/ms] 
kon_tnchca = 2.37;      % [1/mM/ms]
koff_tnchmg = 3.33e-3;  % [1/ms] 
kon_tnchmg = 3e-3;      % [1/mM/ms]
Bmax_CaM = 24e-3;       % [mM] **? about setting to 0 in c-code**   % CaM buffering
koff_cam = 238e-3;      % [1/ms] 
kon_cam = 34;           % [1/mM/ms]
Bmax_myosin = 140e-3;   % [mM]                      % Myosin buffering
koff_myoca = 0.46e-3;   % [1/ms]
kon_myoca = 13.8;       % [1/mM/ms]
koff_myomg = 0.057e-3;  % [1/ms]
kon_myomg = 0.0157;     % [1/mM/ms]
Bmax_SR = 19*.9e-3;     % [mM] (Bers text says 47e-3) 19e-3
koff_sr = 60e-3;        % [1/ms]
kon_sr = 100;           % [1/mM/ms]
Bmax_SLlowsl = 37.4e-3*Vmyo/Vsl;        % [mM]    % SL buffering
Bmax_SLlowj = 4.6e-3*Vmyo/Vjunc*0.1;    % [mM]    %Fei *0.1!!! junction reduction factor
koff_sll = 1300e-3;     % [1/ms]
kon_sll = 100;          % [1/mM/ms]
Bmax_SLhighsl = 13.4e-3*Vmyo/Vsl;       % [mM] 
Bmax_SLhighj = 1.65e-3*Vmyo/Vjunc*0.1;  % [mM] %Fei *0.1!!! junction reduction factor
koff_slh = 30e-3;       % [1/ms]
kon_slh = 100;          % [1/mM/ms]
Bmax_Csqn = 140e-3*Vmyo/Vsr;            % [mM] % Bmax_Csqn = 2.6;      % Csqn buffering
koff_csqn = 65;         % [1/ms] 
kon_csqn = 100;         % [1/mM/ms] 

% Membrane Currents
%% I_Na: Fast Na Current
% am = 0.32*(y(39)+47.13)/(1-exp(-0.1*(y(39)+47.13)));
% bm = 0.08*exp(-y(39)/11);
% if y(39) >= -40
%     ah = 0; aj = 0;
%     bh = 1/(0.13*(1+exp(-(y(39)+10.66)/11.1)));
%     bj = 0.3*exp(-2.535e-7*y(39))/(1+exp(-0.1*(y(39)+32)));
% else
%     ah = 0.135*exp((80+y(39))/-6.8);
%     bh = 3.56*exp(0.079*y(39))+3.1e5*exp(0.35*y(39));
%     aj = (-1.2714e5*exp(0.2444*y(39))-3.474e-5*exp(-0.04391*y(39)))*(y(39)+37.78)/(1+exp(0.311*(y(39)+79.23)));
%     bj = 0.1212*exp(-0.01052*y(39))/(1+exp(-0.1378*(y(39)+40.14)));
% end
% ydot(1) = am*(1-y(1))-bm*y(1);
% ydot(2) = ah*(1-y(2))-bh*y(2);
% ydot(3) = aj*(1-y(3))-bj*y(3);

mss = 1 / ((1 + exp( -(56.86 + y(39)) / 9.03 ))^2);
taum = 0.1292 * exp(-((y(39)+45.79)/15.54)^2) + 0.06487 * exp(-((y(39)-4.823)/51.12)^2);                 
 
ah = (y(39) >= -40) * (0)... 
   + (y(39) < -40) * (0.057 * exp( -(y(39) + 80) / 6.8 )); 
bh = (y(39) >= -40) * (0.77 / (0.13*(1 + exp( -(y(39) + 10.66) / 11.1 )))) ...
   + (y(39) < -40) * ((2.7 * exp( 0.079 * y(39)) + 3.1*10^5 * exp(0.3485 * y(39)))); 
tauh = 1 / (ah + bh); 
hss = 1 / ((1 + exp( (y(39) + 71.55)/7.43 ))^2);
 
aj = (y(39) >= -40) * (0) ...
    +(y(39) < -40) * (((-2.5428 * 10^4*exp(0.2444*y(39)) - 6.948*10^-6 * exp(-0.04391*y(39))) * (y(39) + 37.78)) / ...
                     (1 + exp( 0.311 * (y(39) + 79.23) )));
bj = (y(39) >= -40) * ((0.6 * exp( 0.057 * y(39))) / (1 + exp( -0.1 * (y(39) + 32) ))) ...
   + (y(39) < -40) * ((0.02424 * exp( -0.01052 * y(39) )) / (1 + exp( -0.1378 * (y(39) + 40.14) ))); 
tauj = 1 / (aj + bj);
jss = 1 / ((1 + exp( (y(39) + 71.55)/7.43 ))^2);         
 
ydot(1) = (mss - y(1)) / taum;
ydot(2) = (hss - y(2)) / tauh;
ydot(3) = (jss - y(3)) / tauj;
    
I_Na_junc = Fjunc*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_junc);
I_Na_sl = Fsl*GNa*y(1)^3*y(2)*y(3)*(y(39)-ena_sl);
I_Na = I_Na_junc+I_Na_sl;

%CaMK constants
%update CaMK
%37-cass
%CaMKt=y(58)
KmCaMK=0.15;
aCaMK=0.05;
bCaMK=0.00068;
CaMKo=0.05;
KmCaM=0.0015;
CaMKb=CaMKo*(1.0-y(58))/(1.0+KmCaM/y(37));
CaMKa=CaMKb+y(58);
ydot(58)=aCaMK*CaMKb*(CaMKb+y(58))-bCaMK*y(58);%dCaMKt

%% I_NaL: Late Sodium
%CaMKt=y(58)
%mL=y(61),44
%hL=y(60),43
%hLp=y(59),42
mLss=1.0/(1.0+exp((-(y(39)+42.85))/5.264));
tmL=taum;
ydot(61)=(mLss-y(61))/tmL;%dmL
hLss=1.0/(1.0+exp((y(39)+87.61)/7.488));
thL=200.0;
ydot(60)=(hLss-y(60))/thL;%dhL
hLssp=1.0/(1.0+exp((y(39)+93.81)/7.488));
thLp=3.0*thL;
ydot(59)=(hLssp-y(59))/thLp;%dhLp
GNaL=0.0075;
if epi==1
    GNaL=GNaL*0.6;
end
% fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
% INaL=INaL_scale*GNaL*(v-ENa)*mL*((1.0-fINaLp)*hL+fINaLp*hLp);
fINaLp=(1.0/(1.0+KmCaMK/CaMKa));
I_nal_junc = INaL_Multiplier*Fjunc*GNaL*(y(39)-ena_junc)*y(61)*((1.0-fINaLp)*y(60)+fINaLp*y(59));
I_nal_sl = INaL_Multiplier*Fsl*GNaL*(y(39)-ena_sl)*y(61)*((1.0-fINaLp)*y(60)+fINaLp*y(59));
I_NaL=I_nal_junc+I_nal_sl;

%% I_nabk: Na Background Current
I_nabk_junc = Fjunc*GNaB*(y(39)-ena_junc);
I_nabk_sl = Fsl*GNaB*(y(39)-ena_sl);
I_nabk = I_nabk_junc+I_nabk_sl;

%% I_nak: Na/K Pump Current

if celltype ==1
    scale_INaK = 0.94;   %1
else
    scale_INaK = 1;
end
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0365*sigma*exp(-y(39)*FoRT));
I_nak_junc = scale_INaK*1*Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/y(32))^4) /(Ko+KmKo);
I_nak_sl = scale_INaK*1*Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/y(33))^4) /(Ko+KmKo);
% I_nak_junc = Fjunc*IbarNaK*fnak*Ko /(1+(KmNaip/9.1)^4) /(Ko+KmKo);
% I_nak_sl = Fsl*IbarNaK*fnak*Ko /(1+(KmNaip/9.1)^4) /(Ko+KmKo);
I_nak = I_nak_junc+I_nak_sl;
% fnak=(1/(1+0.1245*exp(-0.1*y(39)*FoRT)+0.0353*exp(-y(39)*FoRT)));
% I_nak_junc = Fjunc*IbarNaK*fnak*Ko /(Ko+KmKo)*(y(32)/(y(32)+KmNaip));
% I_nak_sl = Fsl*IbarNaK*fnak*Ko /(Ko+KmKo)*(y(33)/(y(33)+KmNaip));
% I_nak = I_nak_junc+I_nak_sl;

%% I_kr: Rapidly Activating K Current
% gkr = IKr_Multiplier*0.035*sqrt(Ko/5.4);
xrss = 1/(1+exp(-(y(39)+10)/5));
tauxr = 550/(1+exp((-22-y(39))/9))*6/(1+exp((y(39)-(-11))/9))+230/(1+exp((y(39)-(-40))/20));
ydot(12) = (xrss-y(12))/tauxr;
rkr = 1/(1+exp((y(39)+74)/24));

if female_flag == 0
    % if celltype == 1; gkr=0.0744;
    % elseif celltype == 0; gkr=0.04; end
    if celltype == 1; gkr_scale=0.0744/0.046;
    elseif celltype == 0; gkr_scale=1; end

elseif female_flag == 1
    if celltype == 1; gkr_scale= (0.0744/0.046)*0.8;
    elseif celltype == 0; gkr_scale=0.8; end
end
gkr = 0.035*sqrt(Ko/5.4)*gkr_scale;
% axr1=450/(1+exp((-45-y(39))/10));
% bxr1=6/(1+exp((y(39)-(-30))/11.5));
% xr1inf=1/(1+exp((-26-y(39))/7));
% tauxr1=axr1*bxr1;
% ydot(40)=(xr1inf-y(40))/tauxr1;
% xr1inf=y(40);
% axr2=3/(1+exp((-60-y(39))/20));
% bxr2=1.12/(1+exp((y(39)-60)/20));
% xr2inf=1/(1+exp((y(39)-(-88))/24));
% tauxr2=axr2*bxr2;
% ydot(41)=(xr2inf-y(41))/tauxr2;
% xr2inf=y(41);

I_kr = IKr_Multiplier*gkr*y(12)*rkr*(y(39)-ek);
% I_kr = gkr*y(40)*y(41)*(y(39)-ek);

%% I_ks: Slowly Activating K Current
markov_iks = 0;
% pcaks_junc = -log10(y(36))+3.0; 
% pcaks_sl = -log10(y(37))+3.0;  
% gks_junc = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
% gks_sl = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6)));     

eks = (1/FoRT)*log((Ko+pNaK*Nao)/(y(35)+pNaK*y(34)));
%scale 
if markov_iks == 0
    if female_flag== 0
        if celltype == 1
            iks_scale=1.04;
            gks_junc = IKs_Multiplier*0.0035*iks_scale;
            gks_sl = IKs_Multiplier*0.0035*iks_scale; 
            xsss = 1 / (1+exp((-y(39) + 3.8)/14.25)*iks_scale);
            tauxs = 990.1/(1+exp(-(y(39)+2.436)/14.12)*iks_scale);
            ydot(13) = (xsss-y(13))/tauxs;
            I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
            I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);                                                                                                                                   
            I_ks = I_ks_junc+I_ks_sl;
        elseif celltype == 0
            iks_scale=1.0;
            gks_junc = IKs_Multiplier*0.0035*iks_scale;
            gks_sl = IKs_Multiplier*0.0035*iks_scale; 
            xsss = 1 / (1+exp((-y(39) + 3.8)/14.25)*iks_scale);
            tauxs = 990.1/(1+exp(-(y(39)+2.436)/14.12)*iks_scale);
            ydot(13) = (xsss-y(13))/tauxs;
            I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
            I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);                                                                                                                                   
            I_ks = I_ks_junc+I_ks_sl;
        end
    elseif female_flag == 1
        if celltype == 1
            iks_scale=0.87;
            gks_junc = IKs_Multiplier*0.0035*iks_scale;
            gks_sl = IKs_Multiplier*0.0035*iks_scale; 
            xsss = 1 / (1+exp((-y(39) + 3.8)/14.25)*iks_scale);
            tauxs = 990.1/(1+exp(-(y(39)+2.436)/14.12)*iks_scale);
            ydot(13) = (xsss-y(13))/tauxs;
            I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
            I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);                                                                                                                                   
            I_ks = I_ks_junc+I_ks_sl;
        elseif celltype == 0
            iks_scale=0.83;
            gks_junc = IKs_Multiplier*0.0035*iks_scale;
            gks_sl = IKs_Multiplier*0.0035*iks_scale; 
            xsss = 1 / (1+exp((-y(39) + 3.8)/14.25)*iks_scale);
            tauxs = 990.1/(1+exp(-(y(39)+2.436)/14.12)*iks_scale);
            ydot(13) = (xsss-y(13))/tauxs;
            I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
            I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);                                                                                                                                   
            I_ks = I_ks_junc+I_ks_sl;
        end
    end
    % gks_junc = IKs_Multiplier*0.0035;%scale both junc and sl same sex difference for both
    % gks_sl = IKs_Multiplier*0.0035; %FRA
    % xsss = 1 / (1+exp(-(y(39) + 3.8)/14.25)); % fitting Fra; after 14.25 scale by 1.04 male epi, male endo no scale
    % tauxs = 990.1/(1+exp(-(y(39)+2.436)/14.12));
    % ydot(13) = (xsss-y(13))/tauxs;
    % I_ks_junc = Fjunc*gks_junc*y(13)^2*(y(39)-eks);
    % I_ks_sl = Fsl*gks_sl*y(13)^2*(y(39)-eks);                                                                                                                                   
    % I_ks = I_ks_junc+I_ks_sl;
else
    gks_junc=IKs_Multiplier*0.0065;
    gks_sl=IKs_Multiplier*0.0065; %FRA
    alpha=3.98e-4*exp(3.61e-1*y(39)*FoRT);
    beta=5.74e-5*exp(-9.23e-2*y(39)*FoRT);
    gamma=3.41e-3*exp(8.68e-1*y(39)*FoRT);
    delta=1.2e-3*exp(-3.3e-1*y(39)*FoRT);
    teta=6.47e-3;
    eta=1.25e-2*exp(-4.81e-1*y(39)*FoRT);
    psi=6.33e-3*exp(1.27*y(39)*FoRT);
    omega=4.91e-3*exp(-6.79e-1*y(39)*FoRT);
    
    ydot(42)=-4*alpha*y(42)+beta*y(43);
    ydot(43)=4*alpha*y(42)-(beta+gamma+3*alpha)*y(43)+2*beta*y(44);
    ydot(44)=3*alpha*y(43)-(2*beta+2*gamma+2*alpha)*y(44)+3*beta*y(45);
    ydot(45)=2*alpha*y(44)-(3*beta+3*gamma+alpha)*y(45)+4*beta*y(46);
    ydot(46)=1*alpha*y(44)-(4*beta+4*gamma)*y(46)+delta*y(50);    
    ydot(47)=gamma*y(43)-(delta+3*alpha)*y(47)+beta*y(48);   
    ydot(48)=2*gamma*y(44)+3*alpha*y(47)-(delta+beta+2*alpha+gamma)*y(48)+2*beta*y(49)+2*delta*y(51);
    ydot(49)=3*gamma*y(45)+2*alpha*y(48)-(delta+2*beta+1*alpha+2*gamma)*y(49)+3*beta*y(50)+2*delta*y(52);
    ydot(50)=4*gamma*y(46)+1*alpha*y(49)-(delta+3*beta+0*alpha+3*gamma)*y(50)+2*delta*y(53);
    ydot(51)=1*gamma*y(48)-(2*delta+2*alpha)*y(51)+beta*y(52);  
    ydot(52)=2*gamma*y(49)+2*alpha*y(51)-(2*delta+beta+1*alpha+gamma)*y(52)+2*beta*y(53)+3*delta*y(54);
    ydot(53)=3*gamma*y(50)+1*alpha*y(52)-(2*delta+2*beta+2*gamma)*y(53)+3*delta*y(55);
    ydot(54)=1*gamma*y(52)-(3*delta+1*alpha)*y(54)+beta*y(55);  
    ydot(55)=2*gamma*y(53)+1*alpha*y(54)-(3*delta+1*beta+1*gamma)*y(55)+4*delta*y(56);
    ydot(56)=1*gamma*y(55)-(4*delta+teta)*y(56)+eta*y(57);
    O2=1-(y(42)+y(43)+y(44)+y(45)+y(46)+y(47)+y(49)+y(48)+y(50)+y(51)+y(52)+y(53)+y(54)+y(55)+y(56)+y(57));
    ydot(57)=1*teta*y(56)-(eta+psi)*y(57)+omega*O2;
    I_ks_junc = Fjunc*gks_junc*(y(57)+O2)*(y(39)-eks);
    I_ks_sl = Fsl*gks_sl*(y(57)+O2)*(y(39)-eks);                                                                                                                                   
    I_ks = I_ks_junc+I_ks_sl;
end

%% I_kp: Plateau K current
kp_kp = 1/(1+exp(7.488-y(39)/5.98));
I_kp_junc = Fjunc*gkp*kp_kp*(y(39)-ek);
I_kp_sl = Fsl*gkp*kp_kp*(y(39)-ek);
I_kp = I_kp_junc+I_kp_sl;

%% I_to: Transient Outward K Current (slow and fast components)
% modified for human myocytes
if female_flag == 0
    if celltype == 1; Itos_Multiplier = 0.6;
    elseif celltype == 0; Itos_Multiplier = 1; end
    
elseif female_flag == 1
    if celltype == 1; Itos_Multiplier = 0.26;
    elseif celltype == 0; Itos_Multiplier = 0.64; end
   
else 
    Itos_Multiplier = 1;
   
end

if celltype == 1
    GtoSlow = Itos_Multiplier*1.0*0.13*0.12; %epi
    GtoFast = Itof_Multiplier*1.0*0.13*0.88; %epi0.88
    %tau_multiplier = 
else
    GtoSlow = Itos_Multiplier*0.13*0.3*0.964; %endo
    GtoFast = Itof_Multiplier*0.13*0.3*0.036; %endo
end

xtoss = 1/(1+exp(-(y(39)-19.0)/13));
ytoss = 1/(1+exp((y(39)+19.5)/5));
% rtoss = 1/(1+exp((y(39)+33.5)/10));
tauxtos = 9/(1+exp((y(39)+3.0)/15))+0.5;
tauytos = 800/(1+exp((y(39)+60.0)/10))+30;
% taurtos = 2.8e3/(1+exp((y(39)+60.0)/10))+220; %Fei changed here!! time-dependent gating variable
ydot(8) = (xtoss-y(8))/tauxtos;
ydot(9) = (ytoss-y(9))/tauytos;
% ydot(40)=0;
I_tos = GtoSlow*y(8)*y(9)*(y(39)-ek);    % [uA/uF]

tauxtof = 8.5*exp(-((y(39)+45)/50)^2)+0.5;
%tauxtof = 3.5*exp(-((y(39)+3)/30)^2)+1.5;
tauytof = 85*exp((-(y(39)+40)^2/220))+7;
%tauytof = 20.0/(1+exp((y(39)+33.5)/10))+20.0;
ydot(10) = (xtoss-y(10))/tauxtof;
ydot(11) = (ytoss-y(11))/tauytof;
I_tof = GtoFast*y(10)*y(11)*(y(39)-ek);
I_to = I_tos + I_tof;

%% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(y(39)-ek-59.215)));
bki =(0.49124*exp(0.08032*(y(39)+5.476-ek)) + exp(0.06175*(y(39)-ek-594.31))) /(1 + exp(-0.5143*(y(39)-ek+4.753)));
% Ak1=0.1/(1+exp(0.06*(y(39)-ek-200)));
% Bk1=(3*exp(0.0002*(y(39)-ek+100))+ exp(0.1*(y(39)-ek-10)))/(1+exp(-0.5*(y(39)-ek)));
% kiss=Ak1/(Ak1+Bk1);
% I_ki = GK1*sqrt(Ko/5.4)*kiss*(y(39)-ek);
kiss = aki/(aki+bki);
if female_flag == 0
    if celltype == 1; GK1=0.35*0.98;
    elseif celltype == 0; GK1=0.35; end
    
elseif female_flag == 1
    if celltype == 1; GK1=0.35*0.74;
    elseif celltype == 0; GK1=0.35*0.86; end
end
%I_ki = IK1_Multiplier*0.35*sqrt(Ko/5.4)*kiss*(y(39)-ek);
I_ki = IK1_Multiplier*GK1*sqrt(Ko/5.4)*kiss*(y(39)-ek);%???? replace 0.35 with GK1 like OHd

%% I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
I_ClCa_junc = Fjunc*GClCa/(1+KdClCa/y(36))*(y(39)-ecl);
I_ClCa_sl = Fsl*GClCa/(1+KdClCa/y(37))*(y(39)-ecl);
I_ClCa = I_ClCa_junc+I_ClCa_sl;

I_Clbk = GClB*(y(39)-ecl);

%% I_Ca: L-type Calcium Current
% dss = 1/(1+exp(-(y(39)+4.8)/6.2));
% Ad=4/(1+exp((-45-y(39))/33));
% Bd=5/(1+exp((y(39)+30)/13));
% Cd=10/(1+exp((30-y(39))/50));
% taud = dss*(1-exp(-(y(39)+4.8)/6.2))/(0.035*(y(39)+4.8));
% %taud=Ad*Bd+Cd;
% ydot(4) = (dss-y(4))/taud;
% 
fss = 1/(1+exp((y(39)+35)/9))+0.6/(1+exp((50-y(39))/20));
% Af=21/(1+exp((-55-y(39))/13));
% Bf=20/(1+exp((y(39)+50)/5));
% %Cf=100/(1+exp((-20-y(39))/300));
% Cf=10.5;
% tauff = Af*Bf+Cf;
% ydot(5) = (fss-y(5))/tauff;
% 
% As=21*2/(1+exp((-70-y(39))/20));
% Bs=20*2/(1+exp((y(39)+50)/12));
% Cs=200/(1+exp((30-y(39))/20));
% taufs =As*Bs+Cs;
% ydot(41)= (fss-y(41))/taufs;
% 
% ydot(6) = 1.7*y(36)*(1-y(6))-11.9e-3*y(6); % fCa_junc   koff!!!!!!!!
% ydot(7) = 1.7*y(37)*(1-y(7))-11.9e-3*y(7); % fCa_sl
% fcaCaMSL= 0.1/(1+(0.01/y(37)));
% fcaCaj= 0.1/(1+(0.01/y(36)));
% fcaCaMSL=0;
% fcaCaj= 0;
% ibarca_j = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
% ibarca_sl = pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
% ibark = pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
% ibarna_j = pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
% ibarna_sl = pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
% I_Ca_junc = (Fjunc_CaL*ibarca_j*y(4)*(0.7 * y(5) + 0.3 * y(41))*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.5141;
% I_Ca_sl = (Fsl_CaL*ibarca_sl*y(4)*(0.7 * y(5) + 0.3 * y(41))*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.5141;
% I_Ca = I_Ca_junc+I_Ca_sl;
% I_CaK = (ibark*y(4)*(0.7 * y(5) + 0.3 * y(41))*(Fjunc_CaL*(fcaCaj+(1-y(6)))+Fsl_CaL*(fcaCaMSL+(1-y(7))))*Q10CaL^Qpow)*0.5141;
% I_CaNa_junc = (Fjunc_CaL*ibarna_j*y(4)*(0.7 * y(5) + 0.3 * y(41))*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.5141;
% I_CaNa_sl = (Fsl_CaL*ibarna_sl*y(4)*(0.7 * y(5) + 0.3 * y(41))*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.5141;
% I_CaNa = I_CaNa_junc+I_CaNa_sl;
% I_Catot = I_Ca+I_CaK+I_CaNa;

if celltype ==1
    ICaL_scale = 1.2;  %1
else
    ICaL_scale = 1; 
end

% I_Ca: L-type Calcium Current
dss = 1/(1+exp(-(y(39)+5+mutation_change_vCaL_LQT8*mutation_flag_LQT8)/6));
taud = dss*(1-exp(-(y(39)+5)/6.0))/(0.035*(y(39)+5));
% fss = 1/(1+exp((y(39)+35.06)/3.6))+0.6/(1+exp((50-y(39))/20));
tauf = 1/(0.0197*exp( -(0.0337*(y(39)+14.5))^2 )+0.02);
ydot(4) = (dss-y(4))/taud;
ydot(5) = (fss-y(5))/tauf;
ydot(6) = 1.7*y(36)*(1-y(6))-11.9e-3*y(6); % fCa_junc   koff!!!!!!!!
ydot(7) = 1.7*y(37)*(1-y(7))-11.9e-3*y(7); % fCa_sl
fcaCaMSL= 0.1/(1+(0.01/y(37)));
fcaCaj= 0.1/(1+(0.01/y(36)));
fcaCaMSL=0;
fcaCaj= 0;
%y(6)=0;
%y(7)=0;
ibarca_j = ICaL_scale*pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(36)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibarca_sl = ICaL_scale*pCa*4*(y(39)*Frdy*FoRT) * (0.341*y(37)*exp(2*y(39)*FoRT)-0.341*Cao) /(exp(2*y(39)*FoRT)-1);
ibark = ICaL_scale*pK*(y(39)*Frdy*FoRT)*(0.75*y(35)*exp(y(39)*FoRT)-0.75*Ko) /(exp(y(39)*FoRT)-1);
ibarna_j = ICaL_scale*pNa*(y(39)*Frdy*FoRT) *(0.75*y(32)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
ibarna_sl =ICaL_scale* pNa*(y(39)*Frdy*FoRT) *(0.75*y(33)*exp(y(39)*FoRT)-0.75*Nao)  /(exp(y(39)*FoRT)-1);
I_Ca_junc = (Fjunc_CaL*ibarca_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45*1;
I_Ca_sl = (Fsl_CaL*ibarca_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*0.45*1;
I_Ca = I_Ca_junc+I_Ca_sl;
I_CaK = (ibark*y(4)*y(5)*(Fjunc_CaL*(fcaCaj+(1-y(6)))+Fsl_CaL*(fcaCaMSL+(1-y(7))))*Q10CaL^Qpow)*0.45*1;

I_CaNa_junc = (Fjunc_CaL*ibarna_j*y(4)*y(5)*((1-y(6))+fcaCaj)*Q10CaL^Qpow)*0.45*1;
I_CaNa_sl = (Fsl_CaL*ibarna_sl*y(4)*y(5)*((1-y(7))+fcaCaMSL)*Q10CaL^Qpow)*.45*1;
I_CaNa = I_CaNa_junc+I_CaNa_sl;
I_Catot = I_Ca+I_CaK+I_CaNa;

%% I_ncx: Na/Ca Exchanger flux
if female_flag==0
    scaleNaCas=1;
elseif female_flag==1
    scaleNaCas=1.15;
end


Ka_junc = 1/(1+(Kdact/y(36))^2);
Ka_sl = 1/(1+(Kdact/y(37))^2);
s1_junc = exp(nu*y(39)*FoRT)*y(32)^3*Cao;
s1_sl = exp(nu*y(39)*FoRT)*y(33)^3*Cao;
s2_junc = exp((nu-1)*y(39)*FoRT)*Nao^3*y(36);
s2_sl = exp((nu-1)*y(39)*FoRT)*Nao^3*y(37);
% Shannon 2004
%s3_junc =KmCai*Nao^3*(1+(y(32)/KmNai)^3) + KmNao^3*y(36)*(1+y(36)/KmCai)               + KmCao*y(32)^3 + y(32)^3*Cao + Nao^3*y(36);
%s3_sl =  KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37)*(1+y(37)/KmCai)               + KmCao*y(33)^3 + y(33)^3*Cao + Nao^3*y(37);
% Weber 2001
s3_junc = KmCai*Nao^3*(1+(y(32)/KmNai)^3) + KmNao^3*y(36) + KmNai^3*Cao*(1+y(36)/KmCai) + KmCao*y(32)^3 + y(32)^3*Cao + Nao^3*y(36);
s3_sl =   KmCai*Nao^3*(1+(y(33)/KmNai)^3) + KmNao^3*y(37) + KmNai^3*Cao*(1+y(37)/KmCai) + KmCao*y(33)^3 + y(33)^3*Cao + Nao^3*y(37);

I_ncx_junc = Fjunc*IbarNCX*Q10NCX^Qpow*Ka_junc*(s1_junc-s2_junc)/s3_junc/(1+ksat*exp((nu-1)*y(39)*FoRT))*scaleNaCas;
I_ncx_sl = Fsl*IbarNCX*Q10NCX^Qpow*Ka_sl*(s1_sl-s2_sl)/s3_sl/(1+ksat*exp((nu-1)*y(39)*FoRT))*scaleNaCas;
I_ncx = I_ncx_junc+I_ncx_sl;

%% I_pca: Sarcolemmal Ca Pump Current
if female_flag == 0
    if celltype == 1; pCa_scale=0.88;
    elseif celltype == 0; pCa_scale=1; end
    
elseif female_flag == 1
    if celltype == 1; pCa_scale=1.6;
    elseif celltype == 0; pCa_scale=1.6; end
     
end
I_pca_junc = Fjunc*Q10SLCaP^Qpow*IbarSLCaP*y(36)^1.6/(KmPCa^1.6+y(36)^1.6)*pCa_scale;
I_pca_sl = Fsl*Q10SLCaP^Qpow*IbarSLCaP*y(37)^1.6/(KmPCa^1.6+y(37)^1.6)*pCa_scale;
I_pca = I_pca_junc+I_pca_sl;

%% I_cabk: Ca Background Current
if celltype ==1
    Gkb_scale = 0.6; %1
else
    Gkb_scale = 1;
end

I_cabk_junc = Gkb_scale*Fjunc*GCaB*(y(39)-eca_junc);
I_cabk_sl = Gkb_scale*Fsl*GCaB*(y(39)-eca_sl);
I_cabk = I_cabk_junc+I_cabk_sl;

%% SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
MaxSR = 15; MinSR = 1;
kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/y(31))^2.5);
koSRCa = koCa/kCaSR;
kiSRCa = kiCa*kCaSR;
RI = 1-y(14)-y(15)-y(16);
ydot(14) = (kim*RI-kiSRCa*y(36)*y(14))-(koSRCa*y(36)^2*y(14)-kom*y(15));   % R
ydot(15) = (koSRCa*y(36)^2*y(14)-kom*y(15))-(kiSRCa*y(36)*y(15)-kim*y(16));% O
ydot(16) = (kiSRCa*y(36)*y(15)-kim*y(16))-(kom*y(16)-koSRCa*y(36)^2*RI);   % I
J_SRCarel = ks*y(15)*(y(31)-y(36));          % [mM/ms]
% if t<2000
% J_SRCarel = ks*(y(31)-y(36));          % [mM/ms]
% end
if celltype ==1
    uptake_scale = 1.42;  %1
else
    uptake_scale = 1; 
end

J_serca = uptake_scale*1*Q10SRCaP^Qpow*Vmax_SRCaP*((y(38)/Kmf)^hillSRCaP-(y(31)/Kmr)^hillSRCaP)...
    /(1+(y(38)/Kmf)^hillSRCaP+(y(31)/Kmr)^hillSRCaP);
J_SRleak = Jleak_Multiplier*5.348e-6*(y(31)-y(36));           %   [mM/ms]

%% Sodium and Calcium Buffering
ydot(17) = kon_na*y(32)*(Bmax_Naj-y(17))-koff_na*y(17);        % NaBj      [mM/ms]
ydot(18) = kon_na*y(33)*(Bmax_Nasl-y(18))-koff_na*y(18);       % NaBsl     [mM/ms]

% Cytosolic Ca Buffers
ydot(19) = kon_tncl*y(38)*(Bmax_TnClow-y(19))-koff_tncl*y(19);            % TnCL      [mM/ms]
ydot(20) = kon_tnchca*y(38)*(Bmax_TnChigh-y(20)-y(21))-koff_tnchca*y(20); % TnCHc     [mM/ms]
ydot(21) = kon_tnchmg*Mgi*(Bmax_TnChigh-y(20)-y(21))-koff_tnchmg*y(21);   % TnCHm     [mM/ms]

if female_flag == 0
    if celltype == 1; CaM_multi=1.07;
    elseif celltype == 0; CaM_multi=1; end
    
elseif female_flag == 1
    if celltype == 1; CaM_multi=1.41;
    elseif celltype == 0; CaM_multi=1.21; end
end

ydot(22) = CaM_multi*kon_cam*y(38)*(Bmax_CaM-y(22))-koff_cam*y(22);                 % CaM       [mM/ms]
ydot(23) = kon_myoca*y(38)*(Bmax_myosin-y(23)-y(24))-koff_myoca*y(23);    % Myosin_ca [mM/ms]
ydot(24) = kon_myomg*Mgi*(Bmax_myosin-y(23)-y(24))-koff_myomg*y(24);      % Myosin_mg [mM/ms]
ydot(25) = kon_sr*y(38)*(Bmax_SR-y(25))-koff_sr*y(25);                    % SRB       [mM/ms]
%J_CaB_cytosol = sum(ydot(19:25)); % wrong formulation
J_CaB_cytosol = ydot(19)+ydot(20)+ydot(22)+ydot(23)+ydot(25);

% Junctional and SL Ca Buffers
ydot(26) = kon_sll*y(36)*(Bmax_SLlowj-y(26))-koff_sll*y(26);       % SLLj      [mM/ms]
ydot(27) = kon_sll*y(37)*(Bmax_SLlowsl-y(27))-koff_sll*y(27);      % SLLsl     [mM/ms]
ydot(28) = kon_slh*y(36)*(Bmax_SLhighj-y(28))-koff_slh*y(28);      % SLHj      [mM/ms]
ydot(29) = kon_slh*y(37)*(Bmax_SLhighsl-y(29))-koff_slh*y(29);     % SLHsl     [mM/ms]
J_CaB_junction = ydot(26)+ydot(28);
J_CaB_sl = ydot(27)+ydot(29);

%% Ion concentrations
% SR Ca Concentrations
ydot(30) = kon_csqn*y(31)*(Bmax_Csqn-y(30))-koff_csqn*y(30);       % Csqn      [mM/ms]
ydot(31) = J_serca-(J_SRleak*Vmyo/Vsr+J_SRCarel)-ydot(30);         % Ca_sr     [mM/ms] %Ratio 3 leak current
% ydot(31)=0;

% Sodium Concentrations
I_Na_tot_junc = I_nal_junc+I_Na_junc+I_nabk_junc+3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]
I_Na_tot_sl = I_nal_sl+I_Na_sl+I_nabk_sl+3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
I_Na_tot_sl2 = 3*I_ncx_sl+3*I_nak_sl+I_CaNa_sl;   % [uA/uF]
I_Na_tot_junc2 = 3*I_ncx_junc+3*I_nak_junc+I_CaNa_junc;   % [uA/uF]

ydot(32) = -I_Na_tot_junc*Cmem/(Vjunc*Frdy)+J_na_juncsl/Vjunc*(y(33)-y(32))-ydot(17);
ydot(33) = -I_Na_tot_sl*Cmem/(Vsl*Frdy)+J_na_juncsl/Vsl*(y(32)-y(33))...
   +J_na_slmyo/Vsl*(y(34)-y(33))-ydot(18);
%FluxNaSL=ydot(33);
% ydot(32) = 0;
% ydot(33) = 0;
ydot(34) =J_na_slmyo/Vmyo*(y(33)-y(34));             % [mM/msec] 
% ydot(34)=0;

% Potassium Concentration
I_K_tot = I_to+I_kr+I_ks+I_ki-2*I_nak+I_CaK+I_kp;     % [uA/uF]
% ydot(35) = -I_K_tot*Cmem/(Vmyo*Frdy);           % [mM/msec]
ydot(35) = 0; % -I_K_tot*Cmem/(Vmyo*Frdy);

% Calcium Concentrations
I_Ca_tot_junc = I_Ca_junc+I_cabk_junc+I_pca_junc-2*I_ncx_junc;                   % [uA/uF]
I_Ca_tot_sl = I_Ca_sl+I_cabk_sl+I_pca_sl-2*I_ncx_sl;            % [uA/uF]
ydot(36) = -I_Ca_tot_junc*Cmem/(Vjunc*2*Frdy)+J_ca_juncsl/Vjunc*(y(37)-y(36))...
    -J_CaB_junction+(J_SRCarel)*Vsr/Vjunc+J_SRleak*Vmyo/Vjunc;  % Ca_j
ydot(37) = -I_Ca_tot_sl*Cmem/(Vsl*2*Frdy)+J_ca_juncsl/Vsl*(y(36)-y(37))...
    + J_ca_slmyo/Vsl*(y(38)-y(37))-J_CaB_sl;   % Ca_sl
% ydot(36)=0;
% ydot(37)=0;
% ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol;%+J_ca_slmyo/Vmyo*(y(37)-y(38));    % [mM/msec]
ydot(38) = -J_serca*Vsr/Vmyo-J_CaB_cytosol +J_ca_slmyo/Vmyo*(y(37)-y(38));
% ydot(38)=0;
%if (t<15000)
%    ydot(41) = 0;
%    ydot(42) = 0;
%else
%ydot(41) = -I_Na*Cmem/(Vmyo*Frdy);
%ydot(42) = -I_Ca_sl*Cmem/(Vjunc*2*Frdy)*Vjunc/Vmyo;
%end


%% Simulation type
protocol = 'pace';

switch lower(protocol)
    case {'none',''}
        I_app = 0;
    case 'pace' % w/ current injection at rate 'rate'
		if mod(t, stimPeriod) <= 5
            I_app = 9.5;
        else
            I_app = 0.0;
        end            
end  

%% Membrane Potential
I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          % [uA/uF]
I_Cl_tot = I_ClCa+I_Clbk;                        % [uA/uF]
I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
%I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
%ydot(39) = -(I_Ca_tot+I_K_tot+I_Na_tot-I_app);
ydot(39) = -(I_tot-I_app);
vmax = ydot(39);
%% Outputs
% adjust output depending on the function call

if (nargin == 3)
    output = ydot;
elseif (nargin == 4) & strcmp(runType,'ydot')
    output = ydot;
elseif (nargin == 4) & strcmp(runType,'rates')
    output = r;
elseif (nargin == 4) & strcmp(runType,'currents')
    %currents = [I_Na I_nabk I_nak I_kr I_ks I_kp I_tos I_tof I_ki I_ClCa I_Clbk I_Catot I_ncx I_pca I_cabk J_serca*Vmyo/Vsr];
    %currents = [I_Na I_tof I_tos I_kr I_ks I_ClCa I_Catot J_SRCarel*Vsr/Vmyo J_SRleak RI I_ncx]; 
    %currents= [I_Catot J_SRCarel*2*Vsr*Frdy/Cmem];
    currents = [I_Catot I_Na I_NaL I_kr I_ks I_tos I_tof I_ki I_nak I_ncx];
    %currents=[];
    output = currents;
end