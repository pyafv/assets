%% set parameters
N_cell = 100; % number of cells in the system
radius = 1; % the radius of cell, i.e. l in the paper
phi = 0.5; % packing fraction
P0 = 4.8;
A0 = pi;
tension_difference = 0.2; % \Lambda in the paper
box_size = sqrt(N_cell*pi/phi); % box_size could be a scalar L, or a 1x2 vector [Lx, Ly]
delta_t = 0.01;
% Parameters for tissue energy, see Eq. (3) in the paper
K_P = 1*ones(N_cell, 1);
A0_list = A0*ones(N_cell, 1);
P0_list = P0*ones(N_cell, 1);
% set random seed for this demo
rng(3);
% initial the cell center positions in the center of the box
center_xy = ((rand(N_cell, 2)-0.5)*0.3+0.5)*box_size;

writematrix(center_xy, 'init_pts.csv');

[cell_chain, edgelist, vertex_position, area_list, perimeter_list] = ...
make_finite_voronoi_pbc(center_xy, radius, box_size);
% calculate the interaction force
center_force = get_finite_voronoi_force(center_xy, cell_chain, edgelist, ...
vertex_position, radius, ones(N_cell, 1), ...
A0_list, K_P, P0_list, tension_difference, ...
area_list, perimeter_list, box_size);

writematrix(center_force, 'init_forces.csv');


tic;
for time = 0:delta_t:30
% generate finite Voronoi configuration
[cell_chain, edgelist, vertex_position, area_list, perimeter_list] = ...
make_finite_voronoi_pbc(center_xy, radius, box_size);
% calculate the interaction force
center_force = get_finite_voronoi_force(center_xy, cell_chain, edgelist, ...
vertex_position, radius, ones(N_cell, 1), ...
A0_list, K_P, P0_list, tension_difference, ...
area_list, perimeter_list, box_size);
% update cell center positions
center_xy = mod(center_xy + center_force*delta_t, box_size);
end
toc

% save final pts and force
writematrix(center_xy, 'final_pts.csv');

[cell_chain, edgelist, vertex_position, area_list, perimeter_list] = ...
make_finite_voronoi_pbc(center_xy, radius, box_size);
% calculate the interaction force
center_force = get_finite_voronoi_force(center_xy, cell_chain, edgelist, ...
vertex_position, radius, ones(N_cell, 1), ...
A0_list, K_P, P0_list, tension_difference, ...
area_list, perimeter_list, box_size);

writematrix(center_force, 'final_forces.csv');

% figure, hold on, axis off
% draw_finite_voronoi(cell_chain, edgelist, vertex_position, box_size, center_xy, radius)

