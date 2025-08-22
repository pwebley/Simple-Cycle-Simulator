
function dy = bed_ode(t, y, parm, inlet_bc, outlet_bc, bed_index)

% BED_ODE   Finite Volume ODE solver for one bed in a PSA system
%
%   y - vector of unknowns: [C_tot; C_i; T; q_i] for each spatial node
%   t - current simulation time
%   parm - struct of bed and gas/adsorbent properties
%   inlet_bc / outlet_bc - structs with .Q (flow) and .P as needed

%% Extract constants and geometry
NS = parm.NS;
N  = parm.N_gas;
dz = parm.L / NS;
eps = parm.eps;
rho_b = parm.rho_b;
A = parm.cross_section;

%% Index map
id_Ct  = 1:NS;
id_Ci  = reshape(NS+(1:NS*(N-1)), NS, N-1);
id_T   = NS*(N) + (1:NS);
id_qi  = reshape(NS*(N+1) + (1:NS*N), NS, N);

% Unpack state vector
Ct  = y(id_Ct);                % total gas conc
Ct = max(Ct, 1e-6);
Ci  = zeros(NS,N);
Ci(:,1:N-1) = y(id_Ci);        % species concentrations
Ci(:,N) = Ct - sum(Ci(:,1:N-1),2); % last species by difference
Ci = max(Ci, 1e-9);  % or add eps
T   = y(id_T);                 % temperature
T  = max(T, 1e-6);
qi  = y(id_qi);                % loading per species

% Compute mole fractions
yi = Ci ./ Ct;

% Gas pressure via ideal gas law
P = Ct .* parm.R .* T;

% Gas velocity via Ergun for interior nodes
mu = parm.visc_func(T);
M = yi * parm.MW(:);           % mixture MW
rho_g = P .* M / parm.R ./ T;  % gas density
v = zeros(NS+1,1);
for s = 2:NS
    dP = P(s-1) - P(s);
    v(s) = sign(dP)*sqrt(abs(dP) / parm.Ergun_factor / rho_g(s));
end
% Inlet boundary condition (molar flow)
v(1) = inlet_bc.Q / (A * Ct(1));  % at z = 0
% Outlet boundary condition from flow_network_update (not recomputed here)
v(NS+1) = outlet_bc.Q / (A * Ct(end));  % at z = L

% Upwind flux limiter
upwind = @(Y, v) ...
    (v >= 0) .* Y(1:end-1) + (v < 0) .* Y(2:end);

% Compute derivatives
dCt_dt = zeros(NS,1);
dCi_dt = zeros(NS,N);
dT_dt  = zeros(NS,1);
dqi_dt = zeros(NS,N);

% Compute solid loading terms
sum_term = zeros(NS,1);
q_eq = zeros(NS, N);
k_LDF = parm.k_LDF;  % assume this is a [1 x N] or scalar value in parm

for s = 1:NS
    P_vec = yi(s,:) .* P(s);  % partial pressures at node s
    q_eq(s,:) = get_equilibrium_loading(parm, P_vec, T(s));
end

dqi_dt = k_LDF .* (q_eq - qi);

% Species mass balances including solid-phase accumulation
for i = 1:(N-1)
    % --- Inlet & outlet mole fractions for species i ---
    yi_in  = inlet_bc.y(i);
    yi_out = outlet_bc.y(i);

    % --- Inlet and outlet temperatures ---
    T_in  = inlet_bc.T;
    T_out = outlet_bc.T;

    % --- Pressure at ends from bed state ---
    P_in  = P(1);
    P_out = P(end);

    % --- Total concentration at ends using ideal gas law ---
    Ct_in  = P_in  / (parm.R * T_in);
    Ct_out = P_out / (parm.R * T_out);

    % --- Build extended Ci vector for ghost nodes ---
    Ci_ext = zeros(NS+2,1);     % NS physical nodes + 2 ghost nodes
    Ci_ext(2:end-1) = Ci(:,i);  % Fill physical nodes

    % Left ghost node (inlet)
    Ci_ext(1) = yi_in * Ct_in;

    % Right ghost node (outlet)
    Ci_ext(end) = yi_out * Ct_out;

    % --- Compute species flux and balance ---
    flux_i = v .* upwind(Ci_ext, v);
    gas_term = -(flux_i(2:end) - flux_i(1:end-1)) / dz;

    % --- Species accumulation including solid phase ---
    dCi_dt(:,i) = (1/eps) * (gas_term - rho_b * dqi_dt(:,i));
end

% Build extended Ct vector for both inlet and outlet boundary
Ct_ext = zeros(NS+2,1);    % NS physical nodes + 2 ghost nodes
Ct_ext(2:end-1) = Ct;      % Fill physical nodes

Ct_ext(1) = Ct_in;   % ghost node from inlet
Ct_ext(end) = Ct_out;  % ghost node from outlet

% Upwinded total concentration at faces
Ct_face = upwind(Ct_ext,v);

% Compute total flux and spatial gradient
flux_Ct = v .* Ct_face;
gas_term_total = -(flux_Ct(2:end) - flux_Ct(1:end-1)) / dz;

% Final total mass balance
dCt_dt = (1/eps) * (gas_term_total - rho_b * sum(dqi_dt, 2));

% Heat capacities and enthalpies
cp_g = parm.cp_g;         % J/mol·K (gas-phase mixture)
cp_s = parm.cp_s;         % J/kg·K (solid)
dH   = parm.dH;           % -ΔH_i for each species (J/mol)
source = parm.heat_source;
sink   = parm.heat_sink;

% === Compute face temperature profile (upwinded) ===
T_ext = [T(1); T; T(end)];  % add ghost node on left and right
T_face = upwind(T_ext,v);

% === Compute convective energy flux ===
enthalpy_flux = v .* Ct_face .* cp_g .* T_face;  % W/m²

% === Spatial derivative (divergence) ===
convective_term = -(enthalpy_flux(2:end) - enthalpy_flux(1:end-1)) / dz;

% === Energy accumulation ===
denominator = eps * Ct .* cp_g + (1 - eps) * rho_b * cp_s;
sum_term = -rho_b * sum(dqi_dt .* repmat(dH(:)', NS, 1), 2);


% === Final energy balance ===
dT_dt = (1 ./ denominator) .* ...
        (convective_term - parm.R .* T .* dCt_dt + sum_term + source - sink);


% Return as vector
dy = [dCt_dt;
      reshape(dCi_dt(:,1:N-1), [], 1);
      dT_dt;
      reshape(dqi_dt, [], 1)];
end
