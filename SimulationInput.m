function sim = SimulationInput()

%% === General Simulation Settings ===
sim.n_cycles = 1;            % Number of cycles to run
sim.cycle_time = 840;         % Duration of one cycle (s)
sim.num_beds = 1;             % Number of beds in the system

%% === Species Setup ===
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

%% === Bed & Adsorbent Layer Configuration ===
sim.bed_diameter = 0.041;  % units, m
sim.n_layers = 1;          % number of adsorbent layers

sim.layers(1).length = 0.9;
sim.layers(1).num_nodes = 20;
sim.layers(1).adsorbent_name = 'ActivatedCarbon_1';

% sim.layers(2).length = 0.4;
% sim.layers(2).num_nodes = 16;
% sim.layers(2).adsorbent_name = 'ActivatedCarbon';

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


tanks(1).name = 'Tank1';
tanks(1).type = 'infinite';
tanks(1).isothermal = true;
tanks(1).temperature = 298.15;
tanks(1).pressure = 101325;
tanks(1).composition = [0.02, 0.15, 0.83, 0.00];

tanks(2).name = 'Tank2';
tanks(2).type = 'finite';
tanks(2).isothermal = true;
tanks(2).volume = 0.01;
tanks(2).temperature = 298.15;
tanks(2).pressure = 101325;
tanks(2).composition = [0, 1, 0, 0];

tanks(3).name = 'Tank3';
tanks(3).type = 'infinite';
tanks(3).isothermal = true;
tanks(3).volume = 0.01;
tanks(3).temperature = 298.15;
tanks(3).pressure = 101325;
tanks(3).composition = [0, 0, 1, 0];

sim.tanks = tanks;

%% === Node and Valve Setup ===
[sim.nodes, sim.valves, sim.step_times, sim.valve_schedule] = SetupValvesAndSchedule();
sim.num_valves = numel(sim.valves);
sim.node_names = unique([string({sim.valves.from}), string({sim.valves.to})]);

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
        sim.init_conditions{i}.y0 = repmat([0.02, 0.15, 0.83, 0.0], sim.num_nodes, 1);
    end
end

end

