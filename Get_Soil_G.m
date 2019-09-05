function [G_day, G_night] = Get_Soil_G(Ta_day, Ta_night, Rn_day, Rn_night, fc, Tmin_close, Tann)
    % This function designed for calculating soil heat flux at daytime and nighttime,respectively.
    % ----------------
    % function input :

    % Ta_day   --  Daytime average air temperature,  C
    % Ta_night --  Nighttime average air temperature, C
    % Rn_day   --  Daytime available energy, w/m^2
    % Rn_night --  Nighttime available energy, w/m^2
    % Tmin_close -- The threshold value below which the stomata will close completely
    % Tann     --  Annual average aire temperature, C

    % ----------------

    Variation = Ta_day - Ta_night; % Variation between Tday and Tnight

    % Soil heat flux in Daytime
    G_day = SoilG(Rn_day, Ta_day, Tmin_close, Tann, Variation, fc);
    % Soil heat flux in Nighttime
    G_night = SoilG(Rn_night, Ta_night, Tmin_close, Tann, Variation, fc);

end

function [G] = SoilG(Ai, Ti, Tmin_close, Tann, Variation, fc)
    % note: Algorithm refer to the paper (Mu et al.,2011)
    Gsoil = ones(size(Ai, 1), size(Ai, 2));
    Gsoil(Tann >= 25) = 0;
    Gsoil(Tann < Tmin_close) = 0;
    Gsoil(Variation < 5) = 0;
    Gsoil(Tmin_close <= Tann & Tann < 25 & Variation >= 5) = 4.73 .* Ti(Tmin_close <= Tann & Tann < 25 & Variation >= 5) - 20.87;

    fl = abs(Gsoil) - 0.39 .* abs(Ai);

    Gsoil(fl > 0) = 0.39 .* abs(Ai(fl > 0));
    G = Gsoil .* (1 - fc);

end
