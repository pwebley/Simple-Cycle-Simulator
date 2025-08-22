function cyclesim()
% CYCLESIM   Main simulation entry point for cyclic PSA/VSA simulation
tic;

%% === Load Centralized Simulation Configuration ===
sim = SimulationInput();  % replaces all hardcoded input

%% === Simulation Tolerances ===
sim.tol = 1e-5;            % Relative tolerance for ODE
abs_tol = 1e-7;            % Absolute tolerance

%% === Generate Parameters and Initial Global State ===
parm_all = cell(1, sim.num_beds);

if sim.restart_from_file
    fprintf('Loading previous state from %s...\n', sim.restart_file);
    load(sim.restart_file, 'y_all');  % assumes final bed state is last row
    y0 = y_all(end, :)';  % column vector

    % Generate parm_all (without computing initial_state)
    for i = 1:sim.num_beds
        parm = bed_params_generator(sim, sim.gas, sim.layers);
        parm.varlen = length(parm.initial_state);  % assumes varlen fixed
        parm_all{i} = parm;
    end

else
    % Build initial state and parameters from sim.init_conditions
    for i = 1:sim.num_beds
        parm = bed_params_generator(sim, sim.gas, sim.layers);
        NS = parm.NS;
        R = parm.R;
        init = sim.init_conditions{i};

        Ct = init.P0 ./ (R * init.T0);
        Ci = init.y0 .* Ct;
        Pp = init.y0 .* init.P0;
        q0 = zeros(NS, parm.N_gas);
        for j = 1:NS
            q0(j,:) = get_equilibrium_loading(parm, Pp(j,:), init.T0(j));
        end

        parm.initial_state = double([
            Ct;
            reshape(Ci(:,1:parm.N_gas-1), [], 1);
            init.T0;
            reshape(q0, [], 1)
        ]);
        parm.initial_state = parm.initial_state(:);
        parm.varlen = length(parm.initial_state);
        parm_all{i} = parm;
    end

    % Assemble full initial state vector
    y0 = zeros(sim.num_beds * parm_all{1}.varlen, 1);
    for i = 1:sim.num_beds
        y0((i-1)*parm_all{i}.varlen + 1 : i*parm_all{i}.varlen) = parm_all{i}.initial_state;
    end
end

%% === ODE Integration ===
opts = odeset('RelTol', sim.tol, 'AbsTol', abs_tol, 'MaxStep', 0.5);

[t_all, y_all] = ode15s(@(t, y) ...
    ode_multibed_wrapper(t, y, parm_all, sim.valves, sim.nodes, sim), ...
    [0, sim.n_cycles * sim.cycle_time], y0, opts);

%% === Postprocessing ===
[Q_inlet, Q_outlet, T_inlet, T_outlet, ...
 Y_inlet, Y_outlet, P_inlet, P_outlet] = ...
    reconstruct_boundary_flows(t_all, y_all, parm_all, sim.valves, sim.nodes, sim);

%% === Save Results ===
save('cycle_sim.mat', 't_all', 'y_all', ...
    'Q_inlet', 'Q_outlet', ...
    'T_inlet', 'T_outlet', ...
    'Y_inlet', 'Y_outlet', ...
    'P_inlet', 'P_outlet', ...
    'parm_all', 'sim');

%% === Save Valve Step Times ===
step_times = sim.step_times;
save('valve_step_times.mat', 'step_times');

%% === Optional Plot ===
if sim.output.plot_profiles
    plot_bed_profiles;  % Make sure this reads from sim if needed
end

fprintf('Simulation completed in %.2f seconds.\n', toc);
end
