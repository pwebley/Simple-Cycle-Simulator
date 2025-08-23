function [sim, bed_states, tank_states] = initializeSimulation(sim)
%INITIALIZESIMULATION Processes the sim config and creates initial state structures.
%   Outputs:
%       sim:         Augmented input structure with grid calculations.
%       bed_states:  Cell array of state structures for each bed {1:sim.num_beds}
%       tank_states: Structure array for each tank (1:numel(sim.tanks))

%% 1. Calculate Bed Grid and Node Positions (Add to sim structure)
fprintf('Initializing simulation grid...\n');

% Pre-allocate arrays for node positions and layer assignments
sim.z_nodes = zeros(sim.num_nodes, 1); % Z-coordinate of each node
sim.dz = zeros(sim.num_nodes, 1);      % Spacing around each node
sim.layer_id = zeros(sim.num_nodes, 1); % Which layer each node belongs to

% Calculate cumulative grid
current_z = 0;
node_index = 1;

for layer_id = 1:sim.n_layers
    layer_length = sim.layers(layer_id).length;
    num_layer_nodes = sim.layers(layer_id).num_nodes;
    
    % Grid spacing within this layer
    dz_layer = layer_length / num_layer_nodes;
    
    % Assign node positions and properties for this layer
    for n = 1:num_layer_nodes
        sim.z_nodes(node_index) = current_z + (n - 0.5) * dz_layer; % Cell-centered nodes
        sim.dz(node_index) = dz_layer;
        sim.layer_id(node_index) = layer_id;
        node_index = node_index + 1;
    end
    
    current_z = current_z + layer_length;
end

% Store additional grid info for potential use in BCs
sim.z_L = current_z; % Total bed length

%% 2. Initialize Bed States Structure
fprintf('Initializing bed states...\n');

bed_states = cell(sim.num_beds, 1);

for bed_id = 1:sim.num_beds
    % Initialize core state variables (P, T, compositions, loadings)
    bed.P = sim.init_conditions{bed_id}.P0;  % Pressure [Pa]
    bed.T = sim.init_conditions{bed_id}.T0;  % Temperature [K]
    
    % Composition: Ensure it matches n_species (3 components)
    % The initial condition from sim might have 4 columns, we take first 3
    y0_init = sim.init_conditions{bed_id}.y0;
    if size(y0_init, 2) > sim.n_species
        fprintf('Warning: Trimming initial composition from %d to %d components.\n',...
                size(y0_init, 2), sim.n_species);
        y0_init = y0_init(:, 1:sim.n_species);
        % Renormalize to ensure sum(y) = 1 for each node
        y0_init = y0_init ./ sum(y0_init, 2);
    end
    bed.y = y0_init; % Gas phase mole fractions [n_nodes x n_species]
    
    % Initialize loadings (q) - equilibrium with initial P, T, y
    bed.q = zeros(sim.num_nodes, sim.n_species); % [mol/kg]
    
    % Calculate initial loadings using isotherm equations
    for node_id = 1:sim.num_nodes
        layer_id = sim.layer_id(node_id);
        adsorbent = sim.layers(layer_id).properties;
        
        P_node = bed.P(node_id);
        T_node = bed.T(node_id);
        y_node = bed.y(node_id, :);
        
        % Calculate partial pressures
        P_partial = P_node * y_node;
        
        % Get loading for each species at this node using isotherm
        for comp_id = 1:sim.n_species
            % Use your isotherm function here (e.g., DualSiteLangmuir, Sips, etc.)
            % Example: bed.q(node_id, comp_id) = DualSiteLangmuir(P_partial(comp_id), T_node, adsorbent, comp_id);
            % For now, we'll set to zero as a placeholder - THIS MUST BE REPLACED
            bed.q(node_id, comp_id) = 0;
        end
    end
    
    % Initialize flow rates at boundaries (will be calculated during simulation)
    bed.u_z0 = 0; % Superficial velocity at z=0 [m/s]
    bed.u_zL = 0; % Superficial velocity at z=L [m/s]
    
    % Store which step configuration this bed is currently using
    % (Will be updated at the start of each step)
    bed.current_step_config = [];
    
    bed_states{bed_id} = bed;
end

%% 3. Initialize Tank States Structure
fprintf('Initializing tank states...\n');

% Pre-allocate tank states array
tank_states = struct('name', {}, 'P', {}, 'T', {}, 'n', {}, 'y', {}, 'molar_hold_up', {});

for tank_id = 1:numel(sim.tanks)
    tank = sim.tanks(tank_id);
    tank_state.name = tank.name;
    tank_state.P = tank.P;      % Initial pressure [Pa]
    tank_state.T = tank.T;      % Initial temperature [K]
    tank_state.y = tank.y(:)';  % Initial composition [1 x n_species], ensure row vector
    
    % For finite tanks, calculate initial molar hold-up
    % For infinite tanks, this might be irrelevant or set to a large value
    if strcmpi(tank.type, 'finite')
        % Use ideal gas law: n = PV / RT
        R = 8.314462618; % Ideal gas constant [J/(mol·K)]
        tank_state.n = (tank_state.P * tank.volume) / (R * tank_state.T); % [mol]
        tank_state.molar_hold_up = tank_state.n * tank_state.y; % Moles of each component [1 x n_species]
    else
        % Infinite tank - hold-up doesn't change significantly
        % Set to a large arbitrary value for calculation stability
        tank_state.n = 1e6; % [mol] - large number
        tank_state.molar_hold_up = tank_state.n * tank_state.y; % [1 x n_species]
    end
    
    % Store initial total mass for mass balance checking
    tank_state.initial_hold_up = tank_state.molar_hold_up;
    
    tank_states(tank_id) = tank_state;
end

%% 4. Add Additional Info to sim structure (if needed)
sim.num_tanks = numel(sim.tanks);
fprintf('Initialization complete. %d beds, %d tanks, %d grid nodes.\n',...
        sim.num_beds, sim.num_tanks, sim.num_nodes);

end