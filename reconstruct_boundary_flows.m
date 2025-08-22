function [Q_inlet, Q_outlet, T_inlet, T_outlet, Y_inlet, Y_outlet, P_inlet, P_outlet] = ...
    reconstruct_boundary_flows(t_all, y_all, parm_all, valves, nodes, sim)
% RECONSTRUCT_BOUNDARY_FLOWS
%   Recomputes and logs inlet/outlet boundary conditions.
%
% Outputs:
%   Q_inlet    - [Nt x num_beds] molar flow rate at bed inlets (mol/s)
%   Q_outlet   - [Nt x num_beds] molar flow rate at bed outlets (mol/s)
%   T_inlet    - [Nt x num_beds] inlet gas temperature (K)
%   T_outlet   - [Nt x num_beds] outlet gas temperature (K)
%   Y_inlet    - [Nt x num_beds x N_gas] mole fractions at inlet
%   Y_outlet   - [Nt x num_beds x N_gas] mole fractions at outlet
%   P_inlet    - [Nt x num_beds] gas pressure at inlet (Pa)
%   P_outlet   - [Nt x num_beds] gas pressure at outlet (Pa)

Nt = length(t_all);
n_beds = sim.num_beds;
N_gas = parm_all{1}.N_gas;
R = parm_all{1}.R;

% Preallocate arrays
Q_inlet  = zeros(Nt, n_beds);
Q_outlet = zeros(Nt, n_beds);
T_inlet  = zeros(Nt, n_beds);
T_outlet = zeros(Nt, n_beds);
Y_inlet  = cell(Nt, n_beds);
Y_outlet = cell(Nt, n_beds);
P_inlet  = zeros(Nt, n_beds);
P_outlet = zeros(Nt, n_beds);

for k = 1:Nt
    t = t_all(k);
    y = y_all(k, :)';

    [nodes, ~, inlet_bc, outlet_bc] = flow_network_update(t, y, parm_all, valves, nodes, sim);

    for b = 1:n_beds
        Q_inlet(k, b) = inlet_bc{b}.Q;
        Q_outlet(k, b) = outlet_bc{b}.Q;
        T_inlet(k, b) = inlet_bc{b}.T;
        T_outlet(k, b) = outlet_bc{b}.T;
        Y_inlet{k, b}  = inlet_bc{b}.y(:)';
        Y_outlet{k, b} = outlet_bc{b}.y(:)';

        Ct_in = sum(inlet_bc{b}.y(:)) * inlet_bc{b}.Ct;  % optional consistency check
        Ct_out = sum(outlet_bc{b}.y(:)) * outlet_bc{b}.Ct;

        % Pressure = Ct * R * T
        P_inlet(k, b)  = inlet_bc{b}.Ct  * R * inlet_bc{b}.T;
        P_outlet(k, b) = outlet_bc{b}.Ct * R * outlet_bc{b}.T;
    end
end

end
