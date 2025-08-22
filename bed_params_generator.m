
function parm = bed_params_generator(sim, gas, layer_config)
% BED_PARAMS_GENERATOR   Generate a bed parameter struct with layered adsorbents
%
% Inputs:
%   sim - simulation struct (needs sim.d_bed)
%   gas - struct with fields: gas_names, N_gas, MW, dH
%   layer_config - array of structs with fields:
%                  name, length (m), num_nodes, eps, rho_b, pellet_diameter

% === Basic Setup ===
parm.R = 8.314;               % J/mol/K
parm.d = sim.d_bed;           % common bed diameter
parm.cross_section = pi * (parm.d/2)^2;
parm.gas_names = gas.gas_names;
parm.N_gas = gas.N_gas;
parm.MW = gas.MW;
parm.dH = gas.dH;

% === Total Bed Layout ===
parm.NS = sum([layer_config.num_nodes]);          % total nodes
parm.L = sum([layer_config.length]);              % total bed length
parm.dz = [layer_config.length] ./ [layer_config.num_nodes];  % dz per layer (scalar per layer)
parm.dz = repelem(parm.dz, [layer_config.num_nodes]);         % expand to full profile

% === Initialize per-node properties ===
parm.rho_b = zeros(parm.NS, 1);  % bulk density [kg/m3]
parm.eps = zeros(parm.NS, 1);    % void fraction [-]
parm.dp = zeros(parm.NS, 1);     % pellet diameter [m]
parm.cp_s = zeros(parm.NS, 1);   % solid heat capacity [J/kg/K]
parm.isotherm = cell(parm.NS, 1);    % isotherm parameters per node
parm.k_LDF = 0.5 * ones(parm.NS, parm.N_gas);  % LDF constants per node

% === Assign Properties Layer by Layer ===
node_index = 1;
for i = 1:length(layer_config)
    lay = layer_config(i);
    ads_data = sim.adsorbent_db.(lay.name);

    idx = node_index : node_index + lay.num_nodes - 1;
    parm.rho_b(idx) = ads_data.rho_b;
    parm.eps(idx)   = ads_data.eps;
    parm.dp(idx)    = lay.pellet_diameter;
    parm.cp_s(idx)  = ads_data.cp_s;
    for j = 1:parm.N_gas
        parm.isotherm_type = ads_data.isotherm_type;
        parm.isotherm{idx,j} = ads_data.isotherm.(parm.gas_names{j});
    end
    node_index = node_index + lay.num_nodes;
end

% === Misc Constants ===
parm.cp_g = 29.0;              % J/mol/K (approximate, can be made gas-dependent)
parm.heat_source = 0;          % J/m3/s
parm.heat_sink = 0;            % J/m3/s
parm.Ergun_factor = 1e4;       % for velocity calc
parm.visc_func = @(T) 1e-5 * ones(size(T));  % constant viscosity (adjust if needed)

end