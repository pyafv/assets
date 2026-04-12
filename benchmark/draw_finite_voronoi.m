function draw_finite_voronoi(cell_chain, edgelist, ...
                                                vertex_position, box_size, ...
                                                center_xy, radius)

% <input>
% box_size=[Lx, Ly]: the x and y periodic boundary condition distance
% center_xy: input [] if vertex model

N_cell = numel(cell_chain);

axis equal
hold on

if ~isempty(center_xy)
    plot(center_xy(:, 1), center_xy(:, 2), 'r.')
end

% periodic boundary condition
Lx = box_size(1);
Ly = box_size(end);
rectangle('position', [0 0 Lx Ly])
    
for i_c = 1:N_cell

    chain = cell_chain{i_c};
    chain_len = numel(chain);

    for i_e = 1:chain_len
        i_edge = abs(chain(i_e));
        if ~edgelist(i_edge, 1)
            Ximg = vertex_position(edgelist(i_edge, 2:3), 1:2);
            Ximg = pbc_relocate(center_xy(i_c, :), Ximg, box_size);
            plot(Ximg(:, 1), Ximg(:, 2), 'color', 'b', 'linewidth', 1);
        else
            arc_angle = edgelist(i_edge, 6):0.01:edgelist(i_edge, 7);
            plot(center_xy(i_c, 1)+radius*cos(arc_angle), ...
                    center_xy(i_c, 2)+radius*sin(arc_angle), 'm');
        end
    end
end
    
    
    
    
    

    
    
    
    
    