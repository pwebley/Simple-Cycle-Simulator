
function dydt = ode_multibed_wrapper(t, y, parm_all, valves, nodes, sim)
% ODE_MULTIBED_WRAPPER   Computes dy/dt for all beds in the system

n_beds = numel(parm_all);
dydt = zeros(size(y));

% Compute boundary conditions for each bed
[nodes, Q_valves, inlet_bc, outlet_bc] = flow_network_update(t, y, parm_all, valves, nodes, sim);

% Compute derivatives for each bed
for i = 1:n_beds
    offset = (i - 1) * parm_all{i}.varlen;
    y_bed = y(offset + 1 : offset + parm_all{i}.varlen);
    dydt_bed = bed_ode(t, y_bed, parm_all{i}, inlet_bc{i}, outlet_bc{i}, i);
    dydt(offset + 1 : offset + parm_all{i}.varlen) = dydt_bed;
end

end
