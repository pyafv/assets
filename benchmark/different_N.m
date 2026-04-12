%% Benchmark finite Voronoi force build/update vs N
% Mirrors the structure of your Python timing script.

clear; clc;

%% Reproducibility (matches intent of np.random.seed(42))
rng(42);

%% N sweep: (10**linspace(1,3,10)).astype(int)
Ns = floor(10.^linspace(1, 3, 10));   % 10 values, cast to int-like
Ns = unique(Ns, 'stable');            % ensure no accidental duplicates

%% Pre-generate initial points for each N (like all_pts in Python)
all_pts = cell(numel(Ns), 1);
for i = 1:numel(Ns)
    N = Ns(i);
    all_pts{i} = rand(N, 2);          % in [0,1]
end

%% Simulation parameters (use same physics settings across N)
radius = 1;                 % l in your paper/code
phi = 0.5;                  % packing fraction
P0 = 4.8;
A0 = pi;
tension_difference = 0.2;   % Lambda
delta_t = 0.01;

T = 10.0;                   % match your Python T=10
num_steps = round(T / delta_t);

n_reps = 3;                 % run 3 times per N

%% Timing results: rows correspond to Ns, cols to repetitions
all_times = zeros(numel(Ns), n_reps);

%% Main sweep
for i = 1:numel(Ns)
    N_cell = Ns(i);
    fprintf('Running simulation for N=%d\n', N_cell);

    % Per-N box size (your demo formula)
    box_size = sqrt(N_cell*pi/phi);

    % Per-N parameter lists
    K_P    = ones(N_cell, 1);
    A0_list = A0 * ones(N_cell, 1);
    P0_list = P0 * ones(N_cell, 1);

    % Initial points (reused across repetitions, like Python)
    pts0 = all_pts{i};

    for r = 1:n_reps
        % Reset state for each repetition (important for fair timing)
        center_xy = pts0 * box_size;

        tStart = tic;
        for step = 1:num_steps
            % Build finite Voronoi configuration
            [cell_chain, edgelist, vertex_position, area_list, perimeter_list] = ...
                make_finite_voronoi_pbc(center_xy, radius, box_size);

            % Compute forces
            center_force = get_finite_voronoi_force(center_xy, cell_chain, edgelist, ...
                vertex_position, radius, ones(N_cell, 1), ...
                A0_list, K_P, P0_list, tension_difference, ...
                area_list, perimeter_list, box_size);

            % Update positions (Euler), enforce PBC
            center_xy = mod(center_xy + center_force * delta_t, box_size);
        end
        all_times(i, r) = toc(tStart);
    end
end

%% Save times (same output shape idea as Python: (#Ns) x 3)
writematrix(all_times, 'different_N_times_matlab_new.csv');

%% (Optional) Save Ns too (so the CSV is self-describing)
% writematrix(Ns(:), 'different_N_Ns.csv');

disp('Done. Wrote different_N_times_matlab.csv');
