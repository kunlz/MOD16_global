% The each day was divided into Daytime and Nighttime
% Forcing Data:
% -------------------------------------------------------------------------
% Daytime:
% Rn_D,  W m-2
% Ta_D,  C
% RH_D,  1
% Pa_D,  kPa
% Tmin_D, C

% Nighttime:
% Rn_N, W m-2
% Ta_N,  C
% RH_N,  1
% Pa_N,  kPa
% Tmin_N, C

% Fc,  vegetation cover fraction
% LAI,  leaf area index
% daycont, the number of daytime hours in 24h
% BPLUT, Biome Properties Look-Up Table
% -------------------------------------------------------------------------
% Get the Soil heat flux data
disp(' ******          Get the Soil heat flux data        ****** ')
[G_D, G_N] = Get_Soil_G(Ta_D, Ta_N, Rn_D, Rn_N, Fc, Tmin_close, Tann);

% Call MOD16
disp('        #####------ -------------------- ------##### ')
disp('        #####------ Call MOD16 algorithm ------##### ')
disp('        #####------ -------------------- ------##### ')

[ET_D, Ewet_D, Ttrans_D, Esoil_D] = MOD16_global(Rn_D, Ta_D, RH_D, Pa_D, Tmin_D, G_D, 1, Fc, LAI, BPLUT, daycont);
[ET_N, Ewet_N, Ttrans_N, Esoil_N] = MOD16_CMG(Rn_N, Ta_N, RH_N, Pa_N, Tmin_N, G_N, 0, Fc, LAI, BPLUT, 24 - daycont);

disp(' ****** Delete abnormal value according to MODIS Fc data ****** ')
% Total Ewet
Ewet = Ewet_D + Ewet_N;
Ewet(Fc > 1) = 0; Ewet = single(Ewet);
% Total Ttrans
Ttrans = Ttrans_D + Ttrans_N;
Ttrans(Fc > 1) = 0; Ttrans = single(Ttrans);
% Total Esoil
Esoil = Esoil_D + Esoil_N;
Esoil(Fc > 1) = 0; Esoil = single(Esoil);
% Total ET
ET = ET_D + ET_N;
ET(Fc > 1) = 0; 
ET = single(ET);

% write to
disp(' ****** Writting result to specific output folder ****** ')

save([fileout1 '\ET.Y' year '.DOY' DOY '.mat'], 'ET');
save([fileout2 '\Ewet.Y' year '.DOY' DOY '.mat'], 'Ewet');
save([fileout3 '\Trans.Y' year '.DOY' DOY '.mat'], 'Ttrans');
save([fileout4 '\Esoil.Y' year '.DOY' DOY '.mat'], 'Esoil');
nn = nn + 1;
disp(' ******              Writting finished            ****** ')
disp([' =======================  DOY : ' DOY ' Finished  ======================= '])
disp('                                 ')
