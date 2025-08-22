function q = get_equilibrium_loading(parm, P_vec, T)
% GET_EQUILIBRIUM_LOADING calculates multicomponent adsorption loading
% Inputs:
%   parm   - parameter struct with .isotherm (1 × N_gas struct array) and R
%   P_vec  - partial pressures [1 × N_gas] in kPa
%   T      - temperature in K
% Output:
%   q      - [1 × N_gas] equilibrium loading in mol/kg

n = parm.N_gas;
q = zeros(1, n);


% Calculate individual q_i
for i = 1:n
    gas_name = parm.gas_names{i};
    iso = parm.isotherm.(gas_name);
    Pi = P_vec(i);

    switch lower(parm.isotherm_type)
        case 'langmuir'
            bi = iso.b0 * exp(-iso.dH / (parm.R * T));
            q(i) = iso.qmax * (bi * Pi) / (1 + sum_bP);
        case 'linear'
            q(i) = iso.K * Pi;
        otherwise
            error('Unsupported isotherm model: %s', iso.iso_model);
    end
end

end
