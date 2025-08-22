% === Bed Profile Summary Plotter with Boundary Analysis and Auto Export ===

% Parameters
dt = 1;  % seconds between points
output_dir = 'plots';  % Create if needed

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% === Load required files ===
load('valve_step_times.mat', 'step_times');
load('gas_database.mat');
load('adsorbent_database.mat');
load('cycle_sim.mat', 't_all', 'y_all', ...
     'T_inlet', 'T_outlet', ...
     'P_inlet', 'P_outlet', ...
     'Y_inlet', 'Y_outlet', ...
     'Q_inlet', 'Q_outlet');


% === Simulation setup ===
sim.cycle_time = 134;
sim.gas_names = {'H2O', 'O2', 'N2'};
sim.adsorbent_name = 'Zeolite5A';
sim.num_beds = 2;

% === Reconstruct parameters ===
parm_all = cell(1, sim.num_beds);
for i = 1:sim.num_beds
    parm_all{i} = bed_params_generator(sim.gas_names, sim.adsorbent_name, gas_db, adsorbent_db);
end

% === Mass and Energy Balance ===
dt_all = diff(t_all);  % Typically constant, but check if needed
% Convert Y_inlet and Y_outlet into numeric arrays:
n_times = length(t_all);
n_beds = size(Q_inlet, 2);
n_species = length(sim.gas_names);  % or 3 in your case

Y_inlet_array  = zeros(n_times, n_beds, n_species);
Y_outlet_array = zeros(n_times, n_beds, n_species);

for i = 1:n_times
    for b = 1:n_beds
        Y_inlet_array(i, b, :)  = Y_inlet{i, b};
        Y_outlet_array(i, b, :) = Y_outlet{i, b};
    end
end

% Initialize accumulators
n_species = length(Y_inlet{1,1});
mol_in_total = zeros(sim.num_beds,1);
mol_out_total = zeros(sim.num_beds,1);
mol_in_species = zeros(sim.num_beds, n_species);
mol_out_species = zeros(sim.num_beds, n_species);

for bed = 1:sim.num_beds
    mol_in_total(bed)  = trapz(t_all, Q_inlet(:,bed));
    mol_out_total(bed) = trapz(t_all, Q_outlet(:,bed));
    
    for g = 1:n_species
        mol_in_species(bed, g) = trapz(t_all, Q_inlet(:,bed)  .* Y_inlet_array(:,bed,g));
        mol_out_species(bed,g) = trapz(t_all, Q_outlet(:,bed) .* Y_outlet_array(:,bed,g));    end
end

% display mass balance
fprintf('\n=== Mass Balance Summary ===\n');
for bed = 1:sim.num_beds
    fprintf('Bed %d: Total In = %.3f mol, Total Out = %.3f mol\n', ...
        bed, mol_in_total(bed), mol_out_total(bed));
    for g = 1:n_species
        fprintf('  Species %s: In = %.3f mol, Out = %.3f mol\n', ...
            sim.gas_names{g}, mol_in_species(bed,g), mol_out_species(bed,g));
    end
end

% === Step-Resolved Mass Balance ===
fprintf('\n=== Mass Balance Summary by Step ===\n');
n_steps = length(step_times) - 1;
for step = 1:n_steps
    t_start = step_times(step);
    t_end   = step_times(step+1);
    step_idx = t_all >= t_start & t_all <= t_end;
    t_seg = t_all(step_idx);
    
    fprintf('--- Step %d: Time = %.2f to %.2f s ---\n', step, t_start, t_end);
    for bed = 1:sim.num_beds
        Q_in_seg  = Q_inlet(step_idx, bed);
        Q_out_seg = Q_outlet(step_idx, bed);

        mol_in_step  = trapz(t_seg, Q_in_seg);
        mol_out_step = trapz(t_seg, Q_out_seg);

        fprintf('  Bed %d: Total In = %.3f mol, Total Out = %.3f mol\n', bed, mol_in_step, mol_out_step);

        for g = 1:n_species
            Y_in_g  = cellfun(@(x) x(g), Y_inlet(step_idx, bed));
            Y_out_g = cellfun(@(x) x(g), Y_outlet(step_idx, bed));
            
            Q_in_g  = Q_in_seg .* Y_in_g;
            Q_out_g = Q_out_seg .* Y_out_g;

            mol_in_g  = trapz(t_seg, Q_in_g);
            mol_out_g = trapz(t_seg, Q_out_g);

            fprintf('    Species %s: In = %.3f mol, Out = %.3f mol\n', sim.gas_names{g}, mol_in_g, mol_out_g);
        end
    end
end



% === Interpolation of data for plotting ===

t_uniform = 0:dt:sim.cycle_time;
y_interp = interp1(t_all, y_all, t_uniform);
% Interpolate boundary conditions to match t_uniform
Q_inlet    = interp1(t_all, Q_inlet,    t_uniform);
Q_outlet   = interp1(t_all, Q_outlet,   t_uniform);
T_inlet    = interp1(t_all, T_inlet,    t_uniform);
T_outlet   = interp1(t_all, T_outlet,   t_uniform);
P_inlet    = interp1(t_all, P_inlet,    t_uniform);
P_outlet   = interp1(t_all, P_outlet,   t_uniform);
Y_inlet_interp  = zeros(length(t_uniform), n_beds, n_species);
Y_outlet_interp = zeros(length(t_uniform), n_beds, n_species);

for b = 1:n_beds
    for g = 1:n_species
        Y_inlet_interp(:, b, g)  = interp1(t_all, Y_inlet_array(:,b,g),  t_uniform);
        Y_outlet_interp(:, b, g) = interp1(t_all, Y_outlet_array(:,b,g), t_uniform);
    end
end

step_labels = strcat("Step ", string(1:length(step_times)));

% disp(fieldnames(parm_all{1}));
NS = parm_all{1}.NS;
N_gas = parm_all{1}.N_gas;
A = parm_all{1}.cross_section;
R = parm_all{1}.R;
L = parm_all{1}.L;           % Total bed length (m)
dz = L / NS;                      % Spatial step size
z_centers = dz * (0.5:NS-0.5);    % Node center locations

z_nodes = [1, round(NS/2), NS];   % Node indices
z_locs  = z_centers(z_nodes);     % Actual physical locations

% Generate labels with true locations
z_labels = arrayfun(@(z) sprintf('z = %.2f m', z), z_locs, 'UniformOutput', false);
gas_names = sim.gas_names;
colors = lines(N_gas);

% === Index map ===
id_Ct  = 1:NS;
id_Ci  = reshape(NS+(1:NS*(N_gas-1)), NS, N_gas-1);
id_T   = NS*(N_gas) + (1:NS);
id_qi  = reshape(NS*(N_gas+1) + (1:NS*N_gas), NS, N_gas);

% === Internal bed profiles ===
for bed = 1:sim.num_beds
    varlen = NS * (1 + (N_gas - 1) + 1 + N_gas);
    y_bed = y_interp(:, (bed-1)*varlen + (1:varlen));
    Ct_t = y_bed(:, id_Ct);
    T_t = y_bed(:, id_T);
    P_t = Ct_t .* R .* T_t;

    v_t = zeros(length(t_uniform), NS+1);
    for k = 1:length(t_uniform)
        Ct = y_bed(k, id_Ct)';
        T = y_bed(k, id_T)';
        P = Ct .* R .* T;
        Ci = zeros(NS, N_gas);
        Ci(:,1:N_gas-1) = reshape(y_bed(k, id_Ci), NS, N_gas-1);
        Ci(:,N_gas) = max(0, Ct - sum(Ci(:,1:N_gas-1), 2));
        Ct = sum(Ci, 2);
        yi = Ci ./ Ct;
        MW_mix = yi * parm_all{bed}.MW(:);
        rho_g = P(:) .* MW_mix / R ./ T(:);
        for s = 2:NS
            dP = P(s-1) - P(s);
            v_t(k,s) = sign(dP) * sqrt(abs(dP) / parm_all{bed}.Ergun_factor / rho_g(s));
        end
        v_t(k,1) = Q_inlet(k,bed) / (A * Ct(1));
        v_t(k,NS+1) = Q_outlet(k,bed) / (A * Ct(end));
    end

    mol_flow = v_t(:,2:NS+1) .* Ct_t * A;  % Use face velocities at internal nodes (s=2:NS+1)

    % === Plot results ===
    for idx = 1:length(z_nodes)
        z = z_nodes(idx);
        fig = figure('Name', sprintf('Bed %d — %s', bed, z_labels{idx}));
        tiledlayout(3,2, 'TileSpacing', 'compact');
        sgtitle(sprintf('Bed %d — %s', bed, z_labels{idx}));

        nexttile;
        plot(t_uniform, P_t(:,z)/1000); ylabel('Pressure (kPa)'); title('Pressure'); hold on; grid on;
grid on;

        nexttile;
        % Extract total concentration and compute mole fractions properly
        Ct = y_bed(:, id_Ct);  % Nt x NS
        yi_all = zeros(length(t_uniform), N_gas);  % Will store mole fractions at z

        for i = 1:N_gas
            if i < N_gas
                Ci = reshape(y_bed(:, id_Ci(:,i)), [], NS);  % Nt x NS
            else
                % Compute C_last = Ct - sum(Ci)
                Ci_sum = zeros(size(Ct));
                for j = 1:N_gas-1
                    Ci_sum = Ci_sum + reshape(y_bed(:, id_Ci(:,j)), [], NS);
                end
                Ci = max(0, Ct - Ci_sum);  % Nt x NS
            end
            yi_all(:, i) = Ci(:,z) ./ Ct(:,z);  % Extract mole fraction at node z
        end

        % Plot all species
        for i = 1:N_gas
            plot(t_uniform, yi_all(:, i), 'DisplayName', gas_names{i}); hold on; grid on;
        end
        ylabel('Mole fraction'); title('Composition'); legend(gas_names, 'Location','northeastoutside');

        nexttile;
        for i = 1:N_gas
            qi = reshape(y_bed(:, id_qi(:,i)), [], NS);
            plot(t_uniform, qi(:,z)); hold on; grid on;
        end
        ylabel('q_i (mol/kg)'); title('Adsorbed loading');
        legend(gas_names, 'Location','northeastoutside');

        nexttile;
        plot(t_uniform, y_bed(:, id_T(z))); ylabel('T (K)'); title('Temperature');

        nexttile;
        plot(t_uniform, v_t(:,z)); ylabel('v (m/s)'); title('Velocity');

        nexttile;
        plot(t_uniform, mol_flow(:,z)); ylabel('mol/s'); title('Molar Flow Rate');

        
        % Add step markers
        ax_list = findall(fig, 'Type', 'axes');
        for ax = ax_list'
            axes(ax); hold on; grid on;
            for s = 1:length(step_times)
                xline(step_times(s), '--k', 'HandleVisibility','off');
                if s < length(step_labels)
                    yl = ylim;
                    x_shift = 0.5;  % seconds (adjust as needed)
                    text(step_times(s) + x_shift, yl(2), step_labels(s), ...
                        'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize', 7);
                end
            end
        end

        end

        % === Save Figure ===
        fname_base = sprintf('bed%d_z%d', bed, z);
        saveas(fig, fullfile(output_dir, [fname_base '.png']));
        savefig(fig, fullfile(output_dir, [fname_base '.fig']));
        % exportgraphics(fig, fullfile(output_dir, [fname_base '.pdf'])); % optional high-quality
        % close(fig);
end

% === Boundary plots ===
for bed = 1:sim.num_beds
    fig = figure('Name', sprintf("Bed %d Boundary Conditions", bed), 'Position', [100 100 1000 800]);
    tiledlayout(4,1, 'TileSpacing', 'compact');
    sgtitle(sprintf('Boundary Conditions — Bed %d', bed));

    nexttile;
    plot(t_uniform, Q_inlet(:,bed), 'b', ...
         t_uniform, Q_outlet(:,bed), 'r--');
    grid on;
    ylabel('Flow (mol/s)'); title('Flow Rate'); legend('Inlet','Outlet');

    nexttile;
    plot(t_uniform, T_inlet(:,bed), 'b', ...
         t_uniform, T_outlet(:,bed), 'r--');
    grid on;
    ylabel('Temp (K)'); title('Temperature'); legend('Inlet','Outlet');

    nexttile; hold on; grid on;
    for g = 1:N_gas
        plot(t_uniform, Y_inlet_interp(:,bed,g), '-',  ...
            'Color', colors(g,:), 'DisplayName', [gas_names{g} ' In']);
        plot(t_uniform, Y_outlet_interp(:,bed,g), '--', ...
            'Color', colors(g,:), 'DisplayName', [gas_names{g} ' Out']);
    end
    ylabel('Mole Fraction'); title('Composition'); legend('show');

    nexttile;
    plot(t_uniform, P_inlet(:,bed), 'b', ...
         t_uniform, P_outlet(:,bed), 'r--');
    grid on;
    ylabel('Pressure (Pa)'); xlabel('Time (s)'); title('Pressure'); legend('Inlet','Outlet');

    % Add step markers without adding to legend
    ax_list = findall(fig, 'Type', 'axes');
    for ax = ax_list'
        axes(ax); hold on; grid on;
        for s = 1:length(step_times)
            xline(step_times(s), '--k', 'HandleVisibility','off');
            if s < length(step_times)
                yl = ylim;
                text(step_times(s), yl(2), step_labels(s), ...
                    'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize', 6);
            end
        end
    end

    % === Save Figure ===
    fname_base = sprintf('bed%d_boundary', bed);
    saveas(fig, fullfile(output_dir, [fname_base '.png']));
    savefig(fig, fullfile(output_dir, [fname_base '.fig']));
    % exportgraphics(fig, fullfile(output_dir, [fname_base '.pdf']));
    close(fig);
end
