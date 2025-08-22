
function sim = SimulationSettings()
% SIMULATIONSETTINGS defines solver tolerances, cycle timing, and output controls

%% ------------------TIME AND CYCLE CONTROL---------------------------
sim.t_cycle = 60;               % Duration of one full cycle [s]
sim.n_cycles = 5;               % Number of cycles to simulate
sim.t_total = sim.t_cycle * sim.n_cycles;

%% ------------------SOLVER OPTIONS-----------------------------------
sim.ode_solver = 'ode15s';      % Solver choice
sim.RelTol = 1e-6;              % Relative tolerance for ODE solver
sim.AbsTol = 1e-8;              % Absolute tolerance for ODE solver
sim.dt_max = 0.5;               % Max time step [s]
sim.dt_min = 1e-3;              % Min time step [s]

%% ------------------OUTPUT AND MONITORING----------------------------
sim.save_interval = 1.0;        % Interval for saving results [s]
sim.verbose = true;             % Print progress info to console
sim.use_mass_balance_check = true;

end
