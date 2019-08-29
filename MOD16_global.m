%==============================%
% MOD16 algorithm Main Program %
%==============================%

% Code: Kun Zhang, Lanzhou University, Lanzhou, China.
% Date: First Version   -- 21/2/2017, Adelaide
%       Lastest Version -- 22/5/2018, Shenzhen

% References:
% Cleugh et al., 2007. Regional evaporation estimates from flux tower and MODIS satellite data
% Mu et al., 2007. Development of a global evapotranspiration algorithm based on MODIS and global meteorology data
% Mu et al., 2011. Improvements to a MODIS global terrestrial evapotranspiration algorithm
%--------------------------------------------------------------------------
function [ET,Ewet,Ttrans,Esoil] = MOD16_global(Input,G,flag,fc,LAI,BPLUT,conts)
% ET, MOD16 Evapotranspiration, mm
%----------------
% function input:
% Input     : 
%             1-Rn, net radiation, (W m-2)
%             2-Ta, air temperature, (C)
%             3-RH, relative humidity, (0.01)
%             4-Pa, air pressure, (kPa)
%             5-Tmin, minimum temperature, (C)
% G         : soil heat flux, (W m-2)
% flag      : 1-daytime, 0-nighttime
% LAI       : leaf area index
% BPLUT     : Biome Properties Look-Up Table
% conts     : the number of daytime hours in 24h
[row,clomn] = size(LAI);
%--------------------------------Adapted Parameters (Vegetation cover type)
Tmin_open   = BPLUT(:,:,1); % No inhibition to transpiration
Tmin_close  = BPLUT(:,:,2); % Nearly complete inhibition
VPD_close   = 0.001*BPLUT(:,:,3); % Nearly complete inhibition, Pa~kPa
VPD_open    = 0.001*BPLUT(:,:,4); % No inhibition to transpiration, Pa~kPa
gl_sh       = BPLUT(:,:,5); % Leaf conductance to sensible heat per unit LAI, m/s
gl_e_wv     = BPLUT(:,:,6); % Leaf conductance to evaporated water vapor per unit LAI, m/s
Cl          = BPLUT(:,:,7); % Mean potential stomatal conductance per unit leaf area,m/s
rblmin      = BPLUT(:,:,8); % Minimum value for rtotc, s/m
rblmax      = BPLUT(:,:,9); % Maximum value for rtotc, s/m
beta        = BPLUT(:,:,10); % beta, 
%--------------------------------------------------------------Main_Program
Rn   = Input(:,:,1);           % W m-2   
Ta   = Input(:,:,2)-273.16;    % C
RH   = 0.01.*Input(:,:,3);     % 0.01 
Pa   = 0.001.*Input(:,:,4);    % kPa
Tmin = Input(:,:,5)-273.16;    % C

[Ac,Asoil,VPD,delta,gamma,lambda,rcorr,rrc,fwet,rho] = Calv(Rn,Ta,RH,Pa,G,fc);
%--Module 1 Evaporation from wet canopy surface
[LEwetc]  = fEwetc(Pa,delta,lambda,rrc,VPD,fwet,LAI,Ac,fc,gl_sh,gl_e_wv,rho);
%--Module 2 Plant transpiration
[LEtrans] = fEtrans(delta,Ac,VPD,fc,rrc,Tmin,fwet,LAI,flag,Cl,gl_sh,Tmin_open,Tmin_close,VPD_open,VPD_close,gamma,rcorr,rho,row,clomn);
%--Module 3 Evaporation from soil surface
[LEsoil]  = fEsoil(delta,Asoil,fc,fwet,RH,rrc,VPD,VPD_open,VPD_close,rblmax,rblmin,gamma,rcorr,beta,rho,row,clomn);
%--Total daily evapotranspiration
Ewet   = Wm2mm(LEwetc, lambda, conts);
Ttrans = Wm2mm(LEtrans, lambda, conts);
Esoil  = Wm2mm(LEsoil, lambda, conts);
ET     = Ewet + Ttrans + Esoil;
end

%% W to mm 
function out = Wm2mm(inn, lambda, conts)
Aa = inn./lambda; % kg m-2 s-1
out = Aa.*3600.*conts; % kg m-2 day-1 ==> mm day-1
end

%% Precalculation for Model input
%--------------------------------%
% Precalculation for Model input %
%--------------------------------%
function [Ac,Asoil,VPD,delta,gamma,lambda,rcorr,rrc,fwet,rho]=Calv(Rn,Ta,RH,Pa,G,fc)
% Etrans, Plant transpiration, (W m-2)
%----------------
% function input:
% Rn        : net radiation, (W m-2)
% Ta        : air temperature, (C)
% RH        : relative humidity, (0.01)
% Pa        : air pressure, (kPa)
% G         : soil heat flux, (W m-2)
% fc        : vegetation cover fraction
%----------------
% function ouput:
% Ac        : the part of Net Radiation allocated to canopy, (W m-2)
% Asoil     : the part of Net Radiation allocated to soil surface (W m-2)
% VPD       : vapor pressure deficit, (kPa)
% delta     : Slope of saturation vapour pressure curve at Ta (kPa C-1)
% gamma     : psychrometric constant (kPa C-1)
% lambda    : Latent heat of water vaporization, (MJ kg-1)
% rcorr     : the correct patameter based on the Ta and Pa, (s/m)
% rrc       : resistance to radiative heat transfer through air, (s/m)
% fwet      : wet surface fraction
%----------------
Cp = 1013; % Specific heat (J kg-1 C-1)
eps = 0.622; % e (unitless) is the ratio of molecular weight of water to dry air (equal to 0.622)
% rho=1.292; % Density of air (kg m-3)
rho = air_density(Ta,RH,Pa); % Density of air (kg m-3)
sigma = 5.6703e-8; % Stefen-Boltzman's constant, (W m-2 K-4)
Tk = Ta+273.16; % C to K
lambda = 1e6.*(2.501-0.002361.*Ta); % Latent heat of water vaporization, (J kg-1)

es = 0.6108.*exp((17.27.*Ta)./(Ta+237.3)); % Saturation vapour pressure at Ta, (kPa)
ea = RH.*es; % Actual vapour pressure (kPa)
VPD = es-ea; % Vapor Pressure Deficit, (kPa)
delta = (4098.*es)./((Ta+237.3).^2); % Slope of saturation vapour pressure curve at Ta (kPa/degC)
gamma = (Cp.*Pa)./(eps*lambda); % Psychrometric constant (kPa/degC)
% refer to the wikupedia 0.622 Psychrometric constant

% rcorr=1./((101.3./Pa).*(Tk./293.15).^1.75);
a = (Tk./293.15);
b = sqrt(sqrt(a.*a.*a.*a.*a.*a.*a));
rcorr = 1./(101.3.*b./Pa);

% resistance to radiative heat transfer through air, s/m
rrc = rho.*Cp./(4.*sigma.*Tk.*Tk.*Tk); 

% Wet surface fraction
fwet = RH.*RH.*RH.*RH;
fwet(RH<0.7) = 0;
% Radiation partition
Ac = fc.*Rn; % The part of A allocated to the canopy
Asoil = (1-fc).*Rn-G; % The part of A partitioned on the soil surface
end
% air density module
function [ro] = air_density(t,hr,p)
% AIR_DENSITY calculates density of air
%  Usage :[ro] = air_density(t,hr,p)
%  Inputs:   t = ambient temperature (C)
%           hr = relative humidity [0.01]
%            p = ambient pressure [kPa]  
%  Output:  ro = air density [kg/m3]
%
%  Refs:
% 1)'Equation for the Determination of the Density of Moist Air' P. Giacomo  Metrologia 18, 33-40 (1982)
% 2)'Equation for the Determination of the Density of Moist Air' R. S. Davis Metrologia 29, 67-70 (1992)
%
% ver 1.0   06/10/2006    Jose Luis Prego Borges (Sensor & System Group, Universitat Politecnica de Catalunya)
% ver 1.1   05-Feb-2007   Richard Signell (rsignell@usgs.gov)  Vectorized

%-------------------------------------------------------------------------
T0 = 273.16;        % Triple point of water (aprox. 0C)
T = T0 + t;         % Ambient temperature in Kelvin
p = 1000 .* p;      % Kpa ==> Pa
hr = 100 .* hr;      % 0.01 ==> 100%
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% 1) Coefficients values
R  =  8.314510;           % Molar ideal gas constant   [J/(mol.K)]
Mv = 18.015*10^-3;        % Molar mass of water vapour [kg/mol]
Ma = 28.9635*10^-3;       % Molar mass of dry air      [kg/mol]

A =  1.2378847*10^-5;     % [K^-2]
B = -1.9121316*10^-2;     % [K^-1]
C = 33.93711047;          %
D = -6.3431645*10^3;      % [K]

a0 =  1.58123*10^-6;      % [K/Pa]
a1 = -2.9331*10^-8;       % [1/Pa]
a2 =  1.1043*10^-10;      % [1/(K.Pa)]
b0 =  5.707*10^-6;        % [K/Pa]
b1 = -2.051*10^-8;        % [1/Pa]
c0 =  1.9898*10^-4;       % [K/Pa]
c1 = -2.376*10^-6;        % [1/Pa]
d  =  1.83*10^-11;        % [K^2/Pa^2]
e  = -0.765*10^-8;        % [K^2/Pa^2]
%-------------------------------------------------------------------------
% 2) Calculation of the saturation vapour pressure at ambient temperature, in [Pa]
psv = exp(A.*(T.^2) + B.*T + C + D./T);   % [Pa]
%-------------------------------------------------------------------------
% 3) Calculation of the enhancement factor at ambient temperature and pressure
fpt = 1.00062 + (3.14*10^-8)*p + (5.6*10^-7)*(t.^2);
%-------------------------------------------------------------------------
% 4) Calculation of the mole fraction of water vapour
xv = hr.*fpt.*psv.*(1./p)*(10^-2);
%-------------------------------------------------------------------------
% 5) Calculation of the compressibility factor of air
Z = 1 - ((p./T).*(a0 + a1*t + a2*(t.^2) + (b0+b1*t).*xv + (c0+c1*t).*(xv.^2))) + ((p.^2 ./ T.^2).*(d + e.*(xv.^2)));
%-------------------------------------------------------------------------
% 6) Final calculation of the air density in [kg/m^3]
ro = (p.*Ma./(Z.*R.*T)).*(1 - xv.*(1-Mv./Ma));
end

%% Evaporation from wet canopy surface
%-------------------------------------%
% Evaporation from wet canopy surface %
%-------------------------------------%
function [ETwetc]= fEwetc(Pa,delta,lambda,rrc,VPD,fwet,LAI,Ac,fc,gl_sh,gl_e_wv,rho)
% ETwetc, Evaporation from wet canopy surface, (W m-2)
%----------------
% function input:
% Pa        : air pressure, (kPa)
% delta     : Slope of saturation vapour pressure curve at Ta (kPa C-1)
% lambda    : Latent heat of water vaporization, (MJ kg-1)
% Ta        : air temperature, (C)
% VPD       : vapor pressure deficit, (kPa)
% fwet
% LAI       : leaf area index
% Ac        : the part of Net Radiation allocated to the canopy, (W m-2)
% fc        : vegetation cover fraction
% gl_sh     : leaf conductance to sensible heat per unit LAI, (m/s)
% gl_e_wv   : leaf conductance to evaporated water vapor per unit LAI,(m/s)
%----------------
% rho=1.292; % Density of air (kg m-3)
Cp=1013; % Specific heat (J kg-1 K-1) = (J kg-1 C-1)

epsilon=0.622; % Ratio molecular weight of water vapour/dry air, (unitless)

% gl_sh, Leaf conductance to sensible heat per unit LAI, m/s
rhc=1./(gl_sh.*LAI.*fwet);  % wet canopy resistance to sensible heat, s/m

rhrc=(rhc.*rrc)./(rhc+rrc);  % aerodynamic resistance, s/m

rvc=1./(gl_e_wv.*LAI.*fwet); % wet canopy resistance, s/m

% delta, Slope of saturation vapour pressure curve at Ta (kPa/degC)
% gamma, Psychrometric constant (kPa/degC)
Q=((delta.*Ac)+(rho.*Cp.*VPD.*fc./rhrc)).*fwet;
W=delta+(Pa.*Cp.*rvc)./(lambda.*epsilon.*rhrc);
ETwetc=Q./W; % Evaporation from wet canopy surface, (W m-2)

ETwetc(isnan(ETwetc))=0; % delete the NaN to 0
end
%% Plant transpiration
%---------------------%
% Plant transpiration %
%---------------------%
function [Etrans]=fEtrans(delta,Ac,VPD,fc,rrc,Tmin,fwet,LAI,flag,Cl,gl_sh,Tmin_open,Tmin_close,VPD_open,VPD_close,gamma,rcorr,rho,row,clomn)
% Etrans, Plant transpiration, (W m-2)
%----------------
% function input:
% delta     : Slope of saturation vapour pressure curve at Ta (kPa C-1)
% Ac        : the part of Net Radiation allocated to the canopy, (W m-2)
% VPD       : vapor pressure deficit, (kPa)
% fc        : vegetation cover fraction
% Ta        : air temperature, (C)
% Tmin      : minimum air temperature, (C)
% fwet      :
% LAI       : leaf area index
% flag      : 1=day or 0=night
% Pa        : air pressure, (kPa)
% Cl        : the mean potential stomatal conductance per unit leaf area, (m/s)
% gl_sh     : leaf conductance to sensible heat per unit LAI, (m/s)
% Tmin_open : minmun temperature that has no inhibition to transpiration,(C)
% Tmin_close: minmun temperature that nearly complete inhibition to transpiration, (C)
% VPD_open  : VPD that has no inhibition to transpiration, (kPa)
% VPD_close : VPD that nearly complete inhibition, (kPa)
% gamma     : psychrometric constant (kPa C-1)
%----------------
% rho=1.292; % Density of air (kg m-3)
Cp=1013; % Specific heat (J kg-1 K-1)

% ra, aerodynamic resistance, (s/m)
ra=Ra(gl_sh,rrc); 

% rs, the dry canopy surface resistance to transpiration from the plant, (s/m)
rs = Rs(Tmin,VPD,fwet,LAI,flag,Cl,gl_sh,Tmin_open,Tmin_close,VPD_open,VPD_close,rcorr,row,clomn);

A=((delta.*Ac)+(rho.*Cp.*VPD.*fc./ra)).*(1-fwet);
B=delta+gamma.*(1+rs./ra);
Etrans=A./B;
Etrans(fc==0)=0;
end
function [ra]=Ra(gl_sh,rrc)
% ra, aerodynamic resistance, (s/m)
%----------------
% function input:
% Ta        : air temperature, (C)
% gl_sh     : leaf conductance to sensible heat per unit LAI, (m/s)
%----------------
gl_bl=gl_sh; % leaf-scale boundary layer conductance, m/s
rh=1./gl_bl; % resistance to convective heat transfer, s/m
% rr=(rho.*Cp)./(4.*sigma.*(Ta+273.16).^3); % s/m
rr=rrc;
ra=(rh.*rr)./(rh+rr);
end
function [rs]=Rs(Tmin,VPD,fwet,LAI,flag,Cl,gl_sh,Tmin_open,Tmin_close,VPD_open,VPD_close,rcorr,row,clomn)
% rs, the dry canopy surface resistance to transpiration from the plant, (s/m)
%----------------
% function input:
% Ta        : air temperature, (C)
% Tmin      : minimum air temperature, (C)
% VPD       : vapor pressure deficit, (kPa)
% fwet
% LAI
% flag      : 1=day or 0=night
% Pa        : Air pressure, (kPa)
% Cl        : the mean potential stomatal conductance per unit leaf area, (m/s)
% gl_sh     : leaf conductance to sensible heat per unit LAI, (m/s)
% Tmin_open : minmun temperature that has no inhibition to transpiration,(C)
% Tmin_close: Nearly complete inhibition, (C)
% VPD_open  : VPD that has no inhibition to transpiration, (kPa)
% VPD_close : VPD that nearly complete inhibition, (kPa)
%----------------
gcu = 0.00001; % Cuticular conductance per unit LAI (m/s)

% multiplier that limits potential stomatal conductance by minimum air temperature
mTmin = mT(Tmin,Tmin_open,Tmin_close,row,clomn);
% multiplier used to reduce the potential stomatal conductance when VPD is high enough to inhibit photosynthesis, (kPa)
mVPD = mV(VPD,VPD_open,VPD_close,row,clomn);

Gcu = gcu.*rcorr;  % m/s, Gcu is leaf cuticular conductane
Gs2 = gl_sh; % m/s, Gs2 is leaf boundary-layer conductance

Cc = zeros(row,clomn); % canopy conductance, m/s
if flag == 1 % Daytime
    Gsi = Cl.*mTmin.*mVPD.*rcorr; % m/s, daytime stomatal conductance
    xxx = LAI.*(1-fwet).*Gs2.*(Gsi+Gcu)./(Gsi+Gs2+Gcu);
    Cc(fwet<1 | LAI>0) = xxx(fwet<1 | LAI>0);
elseif flag == 0 % Nighttime
    Gsi = zeros(row,clomn); % m/s, nighttime stomatal conductance
    xxx = LAI.*(1-fwet).*Gs2.*(Gsi+Gcu)./(Gsi+Gs2+Gcu);
    Cc(fwet<1 | LAI>0) = xxx(fwet<1 | LAI>0);
end
rs = 1./Cc;
end
function [mTmin] = mT(Tmin,Tmin_open,Tmin_close,row,clomn)
% mTmin is a multiplier that limits potential stomatal conductance by
% minimum air temperature, (C)
mTmin = ones(row,clomn);
mTmin(Tmin <= Tmin_close) = 0.1;
ss = (Tmin-Tmin_close)./(Tmin_open-Tmin_close);
mTmin(Tmin>Tmin_close & Tmin<Tmin_open) = ss(Tmin>Tmin_close & Tmin<Tmin_open);
end
function [mVPD] = mV(VPD,VPD_open,VPD_close,row,clomn)
% mVPD is a multiplier used to reduce the potential stomatal conductance
% when VPD is high enough to inhibit photosynthesis, (kPa)
mVPD = ones(row,clomn);
mVPD(VPD>=VPD_close) = 0.1;
ss = (VPD_close-VPD)./(VPD_close-VPD_open);
mVPD(VPD>VPD_open & VPD<VPD_close) = ss(VPD>VPD_open & VPD<VPD_close);
end
%% Evaporation from soil surface
%-------------------------------%
% Evaporation from soil surface %
%-------------------------------%
function [ETsoil] = fEsoil(delta,Asoil,fc,fwet,RH,rrc,VPD,VPD_open,VPD_close,rblmax,rblmin,gamma,rcorr,beta,rho,row,clomn)
% ETsoil, Evaporation from soil surface, (W m-2)
%----------------
% function input:
% Pa        : air pressure, (kPa)
% delta     : Slope of saturation vapour pressure curve at Ta (kPa C-1)
% Asoil     : the part of Net Radiation allocated to soil surface (W m-2)
% fc        : vegetation cover fraction
% fwet
% RH        : Relative humidity, 0.01
% Ta        : air temperature, (C)
% VPD       : vapor pressure deficit, (kPa)
% VPD_open  : VPD that has no inhibition to transpiration, (kPa)
% VPD_close : VPD that nearly complete inhibition, (kPa)
% rblmax    : Maximum value for rtotc, (s/m)
% rblmin    : Minimum value for rtotc, (s/m)
% gamma     : psychrometric constant (kPa C-1)
%----------------
% rho    = 1.292; % Density of air (kg m-3)
Cp     = 1013; % Specific heat (J kg-1 K-1)
% beta   = 0.2; % (kPa)
% rtot, the total aerodynamic resistance to vapor transport, (s/m)
% ras, the aerodynamic resistance at the soil surface, (s/m)
[rtot,ras]=frasoil(VPD,VPD_open,VPD_close,rblmax,rblmin,rcorr,rrc,row,clomn);

A=delta.*Asoil+rho.*Cp.*(1-fc).*VPD./ras;
B=delta+gamma.*rtot./ras;
Etwetsoil=(A.*fwet)./B; % (W M-2)
Etsoilpot=(A.*(1-fwet))./B; % (W M-2)
ETsoil=Etwetsoil+Etsoilpot.*RH.^(VPD./beta); % (W M-2)
end
% rtot and ras
function [rtot,ras]=frasoil(VPD,VPD_open,VPD_close,rblmax,rblmin,rcorr,rrc,row,clomn)
% rtot, the total aerodynamic resistance to vapor transport, (s/m)
% ras, the aerodynamic resistance at the soil surface, (s/m)
%----------------
% function input:
% Ta        : air temperature, (C)
% Pa        : air pressure, (kPa)
% VPD       : vapor pressure deficit, (kPa)
% VPD_open  : VPD that has no inhibition to transpiration, (kPa)
% VPD_close : VPD that nearly complete inhibition, (kPa)
% rblmax    : Maximum value for rtotc, (s/m)
% rblmin    : Minimum value for rtotc, (s/m)
%----------------
rtotc=ftotc(VPD,VPD_open,VPD_close,rblmax,rblmin,row,clomn); %(s/m)
% rcorr=bsxfun(@rdivide,1,bsxfun(@times,bsxfun(@rdivide,101.3,Pa),(bsxfun(@rdivide,Tk,293.15)).^1.75));
% the total aerodynamic resistance to vapor transport, (s/m)
rtot=bsxfun(@times,rtotc,rcorr);
% rhs, the resistance to convective heat transfer (s/m)
rhs=rtot; 
% rrs, the resistance to radiative heat  transfer, (s/m)
rrs=rrc;
% rrs=bsxfun(@rdivide,(rho.*Cp),(4.*sigma.*Tk.^3)); 
ras=bsxfun(@rdivide,bsxfun(@times,rhs,rrs),bsxfun(@plus,rhs,rrs)); % (s/m)
end
function [rtotc]=ftotc(VPD,VPD_open,VPD_close,rblmax,rblmin,row,clomn)
% rtotc,
%----------------
% function input:
% VPD       : vapor pressure deficit, (kPa)
% VPD_open  : VPD that has no inhibition to transpiration, (kPa)
% VPD_close : VPD that nearly complete inhibition, (kPa)
% rblmax    : Maximum value for rtotc, (s/m)
% rblmin    : Minimum value for rtotc, (s/m)
%----------------
rtotc=zeros(row,clomn);
rtotc(VPD<=VPD_open)=rblmax(VPD<=VPD_open);
rtotc(VPD>=VPD_close)=rblmin(VPD>=VPD_close);
ss=rblmax-(((rblmax-rblmin).*(VPD_close-VPD))./(VPD_close-VPD_open));
rtotc(VPD>VPD_open & VPD<VPD_close)=ss(VPD>VPD_open & VPD<VPD_close);
end