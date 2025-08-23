function sim = NetworkSimulationInput()
% NETWORKSIMULATIONINPUT Defines a PSA cycle using the 9-state flow model.

%% === General Simulation Settings ===
sim.n_cycles = 1;
sim.cycle_time = 200; % Total cycle time (s)
sim.num_beds = 2;     % Number of beds in the system

% Species and gas properties (reuse from your existing setup)
sim.species_names = {'O2', 'CO2', 'N2'};
sim.n_species = numel(sim.species_names);
sim.species_index = containers.Map(sim.species_names, 1:sim.n_species);

% === Load Gas Properties ===
load('gas_database.mat', 'gas_database');
for i = 1:sim.n_species
    match = strcmp({gas_database.name}, sim.species_names{i});
    if any(match)
        sim.gas(i) = gas_database(find(match,1));
    else
        error(['Species ', sim.species_names{i}, ' not found in gas_database.mat']);
    end
end

%% === Bed & Adsorbent Configuration  ===
sim.bed_diameter = 0.041;
sim.n_layers = 1;          
sim.layers(1).length = 0.9;
sim.layers(1).num_nodes = 20;
sim.layers(1).adsorbent_name = 'ActivatedCarbon_1';
sim.bed_length = sum([sim.layers.length]);
sim.num_nodes = sum([sim.layers.num_nodes]);

% === Load Adsorbent Properties ===
load('adsorbent_database.mat', 'adsorbent_database');
for i = 1:sim.n_layers
    name = sim.layers(i).adsorbent_name;
    match = strcmp({adsorbent_database.name}, name);
    if any(match)
        sim.layers(i).properties = adsorbent_database(find(match,1));
    else
        error(['Adsorbent ', name, ' not found in adsorbent_database.mat']);
    end
end

%% === Tank Definitions ===
sim.tanks(1).name = 'Feed_Tank';
sim.tanks(1).type = 'infinite';
sim.tanks(1).P = 5e5;    % Pressure (Pa)
sim.tanks(1).T = 298.15; % Temperature (K)
sim.tanks(1).y = [0.21, 0.0004, 0.7896]; % Air composition

sim.tanks(2).name = 'Product_Tank';
sim.tanks(2).type = 'finite';
sim.tanks(2).volume = 10.0;
sim.tanks(2).P = 1.5e5;
sim.tanks(2).T = 298.15;
sim.tanks(2).y = [1.0, 0.0, 0.0]; % Assume pure O2 product

sim.tanks(3).name = 'Vent_Tank';
sim.tanks(3).type = 'infinite';
sim.tanks(3).P = 1e5;
sim.tanks(3).T = 298.15;
sim.tanks(3).y = [0.0, 1.0, 0.0]; % Assume pure CO2 waste

%% === CYCLE DEFINITION ===
% Define the duration of each step in the cycle
sim.step_times = [0, 60, 80, 85, 100, 160, 180, 185, 200]; % 8 steps

% Pre-allocate a structure array for steps
sim.step(1:8) = struct(); % For an 8-step cycle

%----------------------------------------------------------------------
% Step 1: Bed A in State 4 (Adsorption), Bed B in State 3 (Blowdown)
%----------------------------------------------------------------------
i = 1;
% --- Bed A Configuration (State 4: u_z0 > 0, u_zL > 0) ---
sim.step(i).BedA.state = 4;
% z=0: Flow IN. Source is Feed_Tank.
sim.step(i).BedA.z0.source = {"Feed_Tank"};
sim.step(i).BedA.z0.flow_law = {"valve"}; % Flow driven by P_Feed - P_BedA_z0
sim.step(i).BedA.z0.parameters = {struct('Cv', 1.5)};
% z=L: Flow OUT. Destination is Product_Tank.
sim.step(i).BedA.zL.destination = {"Product_Tank"};
sim.step(i).BedA.zL.flow_law = {"valve"}; % Flow driven by P_BedA_zL - P_Product
sim.step(i).BedA.zL.parameters = {struct('Cv', 1.2)};

% --- Bed B Configuration (State 3: u_z0 < 0, u_zL = 0) ---
sim.step(i).BedB.state = 3;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedB.z0.destination = {"Vent_Tank"};
sim.step(i).BedB.z0.flow_law = {"valve"}; % Flow driven by P_BedB_z0 - P_Vent
sim.step(i).BedB.z0.parameters = {struct('Cv', 2.0)};
% z=L: No flow.
sim.step(i).BedB.zL.flow_law = {"none"};

%----------------------------------------------------------------------
% Step 2: Bed A in State 4 (Adsorption), Bed B in State 2 (Purge)
%----------------------------------------------------------------------
i = 2;
% --- Bed A Configuration (State 4: u_z0 > 0, u_zL > 0) ---
sim.step(i).BedA.state = 4;
% z=0: Flow IN. Source is Feed_Tank.
sim.step(i).BedA.z0.source = {"Feed_Tank"};
sim.step(i).BedA.z0.flow_law = {"valve"};
sim.step(i).BedA.z0.parameters = {struct('Cv', 1.5)};
% z=L: Flow OUT. Destination is Product_Tank AND BedB_zL.
sim.step(i).BedA.zL.destination = {"Product_Tank", "BedB_zL"};
sim.step(i).BedA.zL.flow_law = {"valve_split"}; % ONE law for total flow
sim.step(i).BedA.zL.parameters = {struct('Cv', 1.2, 'split_frac', 0.2)};
% The solver will use P_BedA_zL and P_Product to get Q_total,
% then split it.

% --- Bed B Configuration (State 2: u_z0 < 0, u_zL < 0) ---
sim.step(i).BedB.state = 2;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedB.z0.destination = {"Vent_Tank"};
sim.step(i).BedB.z0.flow_law = {"valve"}; % Flow driven by P_BedB_z0 - P_Vent
sim.step(i).BedB.z0.parameters = {struct('Cv', 1.2)};
% z=L: Flow IN. Source is Bed A's zL.
sim.step(i).BedB.zL.source = {"BedA_zL"};
sim.step(i).BedB.zL.flow_law = {"linked"}; % Flow rate set by upstream split
sim.step(i).BedB.zL.parameters = {};

%----------------------------------------------------------------------
% Step 3: Pressure Equalization - Bed A to Bed B
% Bed A: State 6 (Depressurization), Bed B: State 8 (Repressurization)
%----------------------------------------------------------------------
i = 3;
% --- Bed A Configuration (State 6: u_z0 = 0, u_zL > 0) ---
sim.step(i).BedA.state = 6;
sim.step(i).BedA.z0.flow_law = {"none"};
% z=L: Flow OUT. Destination is Bed B's z=L.
sim.step(i).BedA.zL.destination = {"BedB_zL"};
sim.step(i).BedA.zL.flow_law = {"valve"}; % Flow driven by P_A_zL > P_B_zL
sim.step(i).BedA.zL.parameters = {struct('Cv', 2.5)}; % Valve coefficient for equalization line
% The solver will automatically use P at z=L for BedA and BedB.

% --- Bed B Configuration (State 8: u_z0 = 0, u_zL < 0) ---
sim.step(i).BedB.state = 8;
sim.step(i).BedB.z0.flow_law = {"none"};
% z=L: Flow IN. Source is Bed A's zL.
sim.step(i).BedB.zL.source = {"BedA_zL"};
sim.step(i).BedB.zL.flow_law = {"linked"}; % Flow rate is determined by the upstream valve
sim.step(i).BedB.zL.parameters = {};

%----------------------------------------------------------------------
% Step 4: Bed A in State 3 (Blowdown), Bed B in State 8 (Repressurization from Product)
%----------------------------------------------------------------------
i = 4;
% --- Bed A Configuration (State 3: u_z0 < 0, u_zL = 0) ---
sim.step(i).BedA.state = 3;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedA.z0.destination = {"Vent_Tank"};
sim.step(i).BedA.z0.flow_law = {"valve"}; % Flow driven by P_BedA_z0 - P_Vent
sim.step(i).BedA.z0.parameters = {struct('Cv', 2.0)};
% z=L: No flow.
sim.step(i).BedA.zL.flow_law = {"none"};

% --- Bed B Configuration (State 8: u_z0 = 0, u_zL < 0) ---
sim.step(i).BedB.state = 8;
% z=0: No flow.
sim.step(i).BedB.z0.flow_law = {"none"};
% z=L: Flow IN. Source is Product_Tank.
sim.step(i).BedB.zL.source = {"Product_Tank"};
sim.step(i).BedB.zL.flow_law = {"valve"}; % Flow driven by P_Product - P_BedB_zL
sim.step(i).BedB.zL.parameters = {struct('Cv', 1.5)};

%----------------------------------------------------------------------
% Step 5: Bed A in State 3 (Blowdown), Bed B in State 4 (Adsorption)
%----------------------------------------------------------------------
i = 5;
% --- Bed A Configuration (State 3: u_z0 < 0, u_zL = 0) ---
sim.step(i).BedA.state = 3;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedA.z0.destination = {"Vent_Tank"};
sim.step(i).BedA.z0.flow_law = {"valve"}; % Flow driven by P_BedA_z0 - P_Vent
sim.step(i).BedA.z0.parameters = {struct('Cv', 2.0)};
% z=L: No flow.
sim.step(i).BedA.zL.flow_law = {"none"};

% --- Bed B Configuration (State 4: u_z0 > 0, u_zL > 0) ---
sim.step(i).BedB.state = 4;
% z=0: Flow IN. Source is Feed_Tank.
sim.step(i).BedB.z0.source = {"Feed_Tank"};
sim.step(i).BedB.z0.flow_law = {"valve"}; % Flow driven by P_Feed - P_BedB_z0
sim.step(i).BedB.z0.parameters = {struct('Cv', 1.5)};
% z=L: Flow OUT. Destination is Product_Tank.
sim.step(i).BedB.zL.destination = {"Product_Tank"};
sim.step(i).BedB.zL.flow_law = {"valve"}; % Flow driven by P_BedB_zL - P_Product
sim.step(i).BedB.zL.parameters = {struct('Cv', 1.2)};

%----------------------------------------------------------------------
% Step 6: Bed B in State 4 (Adsorption), Bed A in State 2 (Purge)
%
%----------------------------------------------------------------------
i = 6;
% --- Bed B Configuration (State 4: u_z0 > 0, u_zL > 0) ---
sim.step(i).BedB.state = 4;
% z=0: Flow IN. Source is Feed_Tank.
sim.step(i).BedB.z0.source = {"Feed_Tank"};
sim.step(i).BedB.z0.flow_law = {"valve"};
sim.step(i).BedB.z0.parameters = {struct('Cv', 1.5)};
% z=L: Flow OUT. Destination is Product_Tank AND BedA_zL.
sim.step(i).BedB.zL.destination = {"Product_Tank", "BedA_zL"};
sim.step(i).BedB.zL.flow_law = {"valve_split"}; % ONE law for total flow
sim.step(i).BedB.zL.parameters = {struct('Cv', 1.2, 'split_frac', 0.2)};
% The solver will use P_BedB_zL and P_Product to get Q_total,
% then split it. 20% of flow is diverted to purge Bed A.

% --- Bed A Configuration (State 2: u_z0 < 0, u_zL < 0) ---
sim.step(i).BedA.state = 2;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedA.z0.destination = {"Vent_Tank"};
sim.step(i).BedA.z0.flow_law = {"valve"}; % Flow driven by P_BedA_z0 - P_Vent
sim.step(i).BedA.z0.parameters = {struct('Cv', 1.2)};
% z=L: Flow IN. Source is Bed B's zL.
sim.step(i).BedA.zL.source = {"BedB_zL"};
sim.step(i).BedA.zL.flow_law = {"linked"}; % Flow rate set by upstream split
sim.step(i).BedA.zL.parameters = {};

%----------------------------------------------------------------------
% Step 7: Pressure Equalization - Bed B to Bed A
% Bed B: State 6 (Depressurization), Bed A: State 8 (Repressurization)
% Mirror image of Step 3
%----------------------------------------------------------------------
i = 7;
% --- Bed B Configuration (State 6: u_z0 = 0, u_zL > 0) ---
sim.step(i).BedB.state = 6;
% z=0: No flow.
sim.step(i).BedB.z0.flow_law = {"none"};
% z=L: Flow OUT in +Z direction. Destination is Bed A's z=L.
sim.step(i).BedB.zL.destination = {"BedA_zL"};
sim.step(i).BedB.zL.flow_law = {"valve"}; % Flow driven by P_BedB_zL > P_BedA_zL
sim.step(i).BedB.zL.parameters = {struct('Cv', 2.5)}; % Valve coefficient for equalization line

% --- Bed A Configuration (State 8: u_z0 = 0, u_zL < 0) ---
sim.step(i).BedA.state = 8;
% z=0: No flow.
sim.step(i).BedA.z0.flow_law = {"none"};
% z=L: Flow IN in -Z direction. Source is Bed B's z=L.
sim.step(i).BedA.zL.source = {"BedB_zL"};
sim.step(i).BedA.zL.flow_law = {"linked"}; % Accepts the flow from Bed B
sim.step(i).BedA.zL.parameters = {};

%----------------------------------------------------------------------
% Step 8: Bed B in State 3 (Blowdown), Bed A in State 8 (Repressurization from Product)
% Mirror image of Step 4
%----------------------------------------------------------------------
i = 8;
% --- Bed B Configuration (State 3: u_z0 < 0, u_zL = 0) ---
sim.step(i).BedB.state = 3;
% z=0: Flow OUT. Destination is Vent_Tank.
sim.step(i).BedB.z0.destination = {"Vent_Tank"};
sim.step(i).BedB.z0.flow_law = {"valve"}; % Flow driven by P_BedB_z0 - P_Vent
sim.step(i).BedB.z0.parameters = {struct('Cv', 2.0)};
% z=L: No flow.
sim.step(i).BedB.zL.flow_law = {"none"};

% --- Bed A Configuration (State 8: u_z0 = 0, u_zL < 0) ---
sim.step(i).BedA.state = 8;
% z=0: No flow.
sim.step(i).BedA.z0.flow_law = {"none"};
% z=L: Flow IN. Source is Product_Tank.
sim.step(i).BedA.zL.source = {"Product_Tank"};
sim.step(i).BedA.zL.flow_law = {"valve"}; % Flow driven by P_Product - P_BedA_zL
sim.step(i).BedA.zL.parameters = {struct('Cv', 1.5)};

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

%% === Output Control ===
sim.output.save_profiles = true;
sim.output.plot_profiles = true;
sim.output.report_recovery = true;

%% === Restart Settings ===
sim.restart_from_file = false;
sim.restart_file = 'cycle_sim.mat';

if ~sim.restart_from_file
    for i = 1:sim.num_beds
        sim.init_conditions{i}.P0 = 1.01e5 * ones(sim.num_nodes, 1);
        sim.init_conditions{i}.T0 = 298 * ones(sim.num_nodes, 1);
        sim.init_conditions{i}.y0 = repmat([0.2, 0.00, 0.80], sim.num_nodes, 1);
    end
end

end