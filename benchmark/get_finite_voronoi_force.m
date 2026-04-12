function center_force = ...
        get_finite_voronoi_force(center_xy, cell_chain, edgelist, vertex_position, ...
        radius, K_A, A0_list, K_P, P0_list, tension_difference, ...
        area_list, perimeter_list, box_size)

% rely on: pbc_relocate, 
	
N_cell = size(center_xy, 1);
N_edge = size(edgelist, 1);
N_vertex = size(vertex_position, 1);
if ~exist('box_size', 'var')
    box_size = sqrt(N_cell);
end
Lx = box_size(1);
Ly = box_size(end);


% calculate cell-edge adjacency matrix
counter_cell_edge_adj = zeros(N_cell, N_edge);
clock_cell_edge_adj = zeros(N_cell, N_edge);
for i_c = 1:N_cell
    chain = cell_chain{i_c};
    counter_cell_edge_adj(i_c, chain(chain>0)) = 1;
    clock_cell_edge_adj(i_c, -chain(chain<0)) = 1;
end
cell_edge_adj = counter_cell_edge_adj + clock_cell_edge_adj;


% calculate dl/dh & da/dh (straight line) and 
dl_dh = zeros(N_edge, 2*N_vertex); % edge vs. v_x, v_y
counter_da_dh = zeros(N_edge, 2*N_vertex);
clock_da_dh = zeros(N_edge, 2*N_vertex);
% calculate dl/dtheta & da/dtheta (arc)
dl_dtheta1 = zeros(N_edge, N_edge);
dl_dtheta2 = zeros(N_edge, N_edge);
da_dtheta1 = zeros(N_edge, N_edge);
da_dtheta2 = zeros(N_edge, N_edge);
for i_e = 1:N_edge % loop over edges to get dl / dh (sl) and dl / dz (para)
    
    % straight line or arc edge
    counter_c_xy = center_xy(counter_cell_edge_adj(:, i_e)>0, 1:2);
    v1_id = edgelist(i_e, 2);
    v2_id = edgelist(i_e, 3);
    if v1_id==0, continue, end % isolated cell
    
    counter_v1_xy = pbc_relocate(counter_c_xy, vertex_position(v1_id, 1:2), box_size);
    counter_v2_xy = pbc_relocate(counter_c_xy, vertex_position(v2_id, 1:2), box_size);
    
    % straight line or arc edge
    % da / dhx
    da_dh_1x = 1/2*( counter_v2_xy(2)-counter_v1_xy(2) );
    counter_da_dh(i_e, v1_id) = counter_da_dh(i_e, v1_id) + da_dh_1x;
    counter_da_dh(i_e, v2_id) = counter_da_dh(i_e, v2_id) + da_dh_1x;
    % da / dhy
    da_dh_1y = -1/2*( counter_v1_xy(1)+counter_v2_xy(1) );
    counter_da_dh(i_e, v1_id+N_vertex) = counter_da_dh(i_e, v1_id+N_vertex) + da_dh_1y;
    counter_da_dh(i_e, v2_id+N_vertex) = counter_da_dh(i_e, v2_id+N_vertex) - da_dh_1y;
    
    
    if ~edgelist(i_e, 1) % straight line
        % dl / dhx
        dl_dh_1x = ( counter_v1_xy(1)-counter_v2_xy(1) ) / norm(counter_v1_xy-counter_v2_xy);
        dl_dh(i_e, v1_id) = dl_dh(i_e, v1_id) + dl_dh_1x;
        dl_dh(i_e, v2_id) = dl_dh(i_e, v2_id) - dl_dh_1x;
        % dl / dhy
        dl_dh_1y = ( counter_v1_xy(2)-counter_v2_xy(2) ) / norm(counter_v1_xy-counter_v2_xy);
        dl_dh(i_e, v1_id+N_vertex) = dl_dh(i_e, v1_id+N_vertex) + dl_dh_1y;
        dl_dh(i_e, v2_id+N_vertex) = dl_dh(i_e, v2_id+N_vertex) - dl_dh_1y;
        
        clock_c_xy = center_xy(clock_cell_edge_adj(:, i_e)>0, 1:2);
        clock_v1_xy = pbc_relocate(clock_c_xy, vertex_position(v1_id, 1:2), box_size);
        clock_v2_xy = pbc_relocate(clock_c_xy, vertex_position(v2_id, 1:2), box_size);
        % da / dhx
        clock_da_dh(i_e, v1_id) = clock_da_dh(i_e, v1_id) + da_dh_1x;
        clock_da_dh(i_e, v2_id) = clock_da_dh(i_e, v2_id) + da_dh_1x;
        % da / dhy
        da_dh_1y = -1/2*( clock_v1_xy(1)+clock_v2_xy(1) );
        clock_da_dh(i_e, v1_id+N_vertex) = clock_da_dh(i_e, v1_id+N_vertex) + da_dh_1y;
        clock_da_dh(i_e, v2_id+N_vertex) = clock_da_dh(i_e, v2_id+N_vertex) - da_dh_1y;
        
    else % arc edge
        
        theta_diff = edgelist(i_e, 7) - edgelist(i_e, 6);
        dl_dtheta1(i_e, i_e) = -radius;
        dl_dtheta2(i_e, i_e) = radius;
        % sector minus triangle
        da_dtheta1(i_e, i_e) = -0.5*radius^2*(1-cos(theta_diff));
        da_dtheta2(i_e, i_e) = 0.5*radius^2*(1-cos(theta_diff));

        
    end
    
end

    
    
%% calculate dtheta / dr for arc edge
% everything is relocated around [xe ye]
dtheta1_dr = zeros(N_edge, 2*N_cell);
dtheta2_dr = zeros(N_edge, 2*N_cell);    
for i_e = 1:N_edge
    if any(edgelist(i_e, 1:2)==0), continue, end % straight line or isolated edge

    for theta_order = 1:2

        theta_vid = edgelist(i_e, 1+theta_order);
        theta = edgelist(i_e, 5+theta_order);
        main_cid = edgelist(i_e, 4);
        neighbor_cid = sum(edgelist(~edgelist(:, 1)&any(edgelist(:, 2:3)==theta_vid, 2), 4:5)) - main_cid;
        neighbor_cxy = pbc_relocate(center_xy(main_cid, :), center_xy(neighbor_cid, :), box_size);
        xx = neighbor_cxy(1) - center_xy(main_cid, 1);
        yy = neighbor_cxy(2) - center_xy(main_cid, 2);
        whole_bunch = sign(yy)*sqrt((4*radius^2-xx^2-yy^2)/(xx^2+yy^2));
		if ~isreal(whole_bunch)
			continue % 211014 quick update, not sure is correct
		end
        option1 = get_angle([0 0], [xx-yy*whole_bunch, yy+xx*whole_bunch]);
        option2 = get_angle([0 0], [xx+yy*whole_bunch, yy-xx*whole_bunch]);

        theta_diff = abs([option1 option2]-mod(theta, 2*pi));
        if theta_diff(1)<theta_diff(2) % option1 is correct
            dx_value = -((yy*(radius*xx^3 + radius*xx*yy^2 + sqrt(-(radius^2*yy^2*(xx^2 + yy^2)*(-4*radius^2 + xx^2 + yy^2)))))/((xx^2 + yy^2)*sqrt(-(radius^2*yy^2*(xx^2 + yy^2)*(-4*radius^2 + xx^2 + yy^2)))));
            dy_value = xx/(xx^2 + yy^2) - (radius*yy^2)/sqrt(-(radius^2*yy^2*(xx^2 + yy^2)*(-4*radius^2 + xx^2 + yy^2)));
        else
            dx_value = (yy*(radius*xx^3 + radius*xx*yy^2 - sqrt(-(radius^2*yy^2*(xx^2 + yy^2)*(-4*radius^2 + xx^2 + yy^2)))))/((xx^2 + yy^2)*sqrt(-(radius^2*yy^2*(xx^2 + yy^2)*(-4*radius^2 + xx^2 + yy^2))));
            dy_value = xx/(xx^2 + yy^2) + (radius*yy^2)/sqrt(-(radius^2*yy^2*(xx^2 + yy^2)*(-4*radius^2 + xx^2 + yy^2)));
        end
        if theta_order==1
            dtheta1_dr(i_e, neighbor_cid) = dx_value;
            dtheta1_dr(i_e, main_cid) = -dx_value;
            dtheta1_dr(i_e, neighbor_cid+N_cell) = dy_value;
            dtheta1_dr(i_e, main_cid+N_cell) = -dy_value;
        else
            dtheta2_dr(i_e, neighbor_cid) = dx_value;
            dtheta2_dr(i_e, main_cid) = -dx_value;
            dtheta2_dr(i_e, neighbor_cid+N_cell) = dy_value;
            dtheta2_dr(i_e, main_cid+N_cell) = -dy_value;
        end

    end
end
    

%% calculate straight line dh_dr
dh_dr = zeros(2*N_vertex, 2*N_cell);
for i_v = 1:N_vertex
    vertex_info = vertex_position(i_v, :);
    joint_v_x = vertex_info(1);
    joint_v_y = vertex_info(2);
    if vertex_info(3)==3 % eee, inner vertex
        e1_id = vertex_info(4);
        e2_id = vertex_info(5);
        e3_id = vertex_info(6);
        
        
        x1 = pbc_relocate(joint_v_x, center_xy(e1_id, 1), Lx);
        y1 = pbc_relocate(joint_v_y, center_xy(e1_id, 2), Ly);
        x2 = pbc_relocate(joint_v_x, center_xy(e2_id, 1), Lx);
        y2 = pbc_relocate(joint_v_y, center_xy(e2_id, 2), Ly);
        x3 = pbc_relocate(joint_v_x, center_xy(e3_id, 1), Lx);
        y3 = pbc_relocate(joint_v_y, center_xy(e3_id, 2), Ly);

        % d hx / d e1x
        dh_dr(i_v, e1_id) = -((y2-y3)*(2*x1*(x3*(y2-y1)+x2*(y1-y3))+(y1-y2)*(x3^2+(y1-y3)*(y2-y3))+x2^2*(y3-y1)+ x1^2*(y3-y2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
        % d hx / d e2x
        dh_dr(i_v, e2_id) = ((x3^2*(y1-y2)-2*x2*(x3*(y1-y2)+x1*(y2-y3))+x2^2*(y1-y3)+(x1^2+(y1-y2)*(y1-y3))*(y2-y3))*(y1-y3))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
        % d hx / d e3x
        dh_dr(i_v, e3_id) = -((y1-y2)*(x3^2*(y2-y1)+2*x2*x3*(y1-y3)+(x1^2+(y1-y2)*(y1-y3))*(y2-y3)+x2^2*(y3-y1)+ 2*x1*x3*(y3-y2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);

        % d hx / d e1y
        dh_dr(i_v, e1_id + N_cell) = ((y2-y3)*(x1^2*(x2-x3)+x2^2*x3+x3*(y1-y2)^2-x2*(x3^2+(y1-y3)^2)+ x1*(-x2^2+x3^2+2*y1*y2-y2^2-2*y1*y3+y3^2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
        % d hx / d e2y
        dh_dr(i_v, e2_id + N_cell) = -((y1-y3)*(x1^2*(x2-x3)+x2^2*x3-x2*x3^2-x3*(y1-y2)^2+x1*(-x2^2+x3^2+(y2-y3)^2)+ x2*(y1-y3)*(y1-2*y2+y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
        % d hx / d e3y
        dh_dr(i_v, e3_id + N_cell) = ((y1-y2)*(x1^2*(x2-x3)+x2^2*x3+x2*(-x3^2+(y1-y3)^2)-x1*(x2^2-x3^2+(y2-y3)^2)- x3*(y1-y2)*(y1+y2-2*y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);

        % d hy / d e1x
        dh_dr(i_v + N_vertex, e1_id) = -((x2-x3)*(-((y1-y2)*(x3^2+(y1-y3)*(y2-y3)))+x2^2*(y1-y3)+x1^2*(y2-y3)+ 2*x1*(x3*(y1-y2)+x2*(y3-y1))))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
        % d hy / d e2x
        dh_dr(i_v + N_vertex, e2_id) = -((x1-x3)*(x3^2*(y1-y2)-2*x2*(x3*(y1-y2)+x1*(y2-y3))+x2^2*(y1-y3)+ (x1^2+(y1-y2)*(y1-y3))*(y2-y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
        % d hy / d e3x
        dh_dr(i_v + N_vertex, e3_id) = ((x1-x2)*(x3^2*(y2-y1)+2*x2*x3*(y1-y3)+(x1^2+(y1-y2)*(y1-y3))*(y2-y3)+x2^2*(y3-y1)+ 2*x1*x3*(y3-y2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);

        % d hy / d e1y
        dh_dr(i_v + N_vertex, e1_id + N_cell) = -((x2-x3)*(x1^2*(x2-x3)+x2^2*x3+x3*(y1-y2)^2-x2*(x3^2+(y1-y3)^2)+ x1*(-x2^2+x3^2+2*y1*y2-y2^2-2*y1*y3+y3^2)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
        % d hy / d e2y
        dh_dr(i_v + N_vertex, e2_id + N_cell) = ((x1-x3)*(x1^2*(x2-x3)+x2^2*x3-x2*x3^2-x3*(y1-y2)^2+x1*(-x2^2+x3^2+(y2-y3)^2)+ x2*(y1-y3)*(y1-2*y2+y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);
        % d hy / d e3y
        dh_dr(i_v + N_vertex, e3_id + N_cell) = -((x1-x2)*(x1^2*(x2-x3)+x2^2*x3+x2*(-x3^2+(y1-y3)^2)-x1*(x2^2-x3^2+(y2-y3)^2)- x3*(y1-y2)*(y1+y2-2*y3)))/(2*(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))^2);

    else % ee, circle vertex
        e1_id = vertex_info(4);
        e2_id = vertex_info(5);

        x1 = pbc_relocate(joint_v_x, center_xy(e1_id, 1), Lx);
        y1 = pbc_relocate(joint_v_y, center_xy(e1_id, 2), Ly);
        x2 = pbc_relocate(joint_v_x, center_xy(e2_id, 1), Lx);
        y2 = pbc_relocate(joint_v_y, center_xy(e2_id, 2), Ly);
        
        vx_option1 = (x1^3 - x1^2*x2 + x2^3 + x1*(-x2^2 + (y1 - y2)^2) + x2*(y1 - y2)^2 - sqrt( (-x1^2 + 2*x1*x2 - x2^2 - (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2) );
        vx_option2 = (x1^3 - x1^2*x2 + x2^3 + x1*(-x2^2 + (y1 - y2)^2) + x2*(y1 - y2)^2 + sqrt( (-x1^2 + 2*x1*x2 - x2^2 - (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2)))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2) );
        if abs(joint_v_x-vx_option1)<abs(joint_v_x-vx_option2) % use vx_option1
            % d hx / d e1x
            dh_dr(i_v, e1_id) = 0.5 + (2*radius^2*(x1 - x2)*(y1 - y2)^2)/((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));
            % d hx / d e2x
            dh_dr(i_v, e2_id) = 0.5 - (2*radius^2*(x1 - x2)*(y1 - y2)^2)/((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));

            % d hx / d e1y
            dh_dr(i_v, e1_id + N_cell) = ((-4*radius^2*(x1 - x2)^2 + (x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)^2)*(y1 - y2))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));
            % d hx / d e2y
            dh_dr(i_v, e2_id + N_cell) = -((-4*radius^2*(x1 - x2)^2 + (x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)^2)*(y1 - y2))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));

            % d hy / d e1x
            dh_dr(i_v + N_vertex, e1_id) = -((x1^4 - 4*x1^3*x2 + x2^4 - 4*x1*x2*(x2^2 + (y1 - y2)^2) + 2*x1^2*(3*x2^2 + (y1 - y2)^2) - 4*radius^2*(y1 - y2)^2 + 2*x2^2*(y1 - y2)^2 + (y1 - y2)^4)*(y1 - y2))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));
            % d hy / d e2x
            dh_dr(i_v + N_vertex, e2_id) = ((x1^4 - 4*x1^3*x2 + x2^4 - 4*x1*x2*(x2^2 + (y1 - y2)^2) + 2*x1^2*(3*x2^2 + (y1 - y2)^2) - 4*radius^2*(y1 - y2)^2 + 2*x2^2*(y1 - y2)^2 + (y1 - y2)^4)*(y1 - y2))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));

            % d hy / d e1y
            dh_dr(i_v + N_vertex, e1_id + N_cell) = 0.5 - (2*radius^2*(x1 - x2)*(y1 - y2)^2)/((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));
            % d hy / d e2y
            dh_dr(i_v + N_vertex, e2_id + N_cell) = 0.5 + (2*radius^2*(x1 - x2)*(y1 - y2)^2)/((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));

        else % use vx_option2
            % d hx / d e1x
            dh_dr(i_v, e1_id) = 0.5 - (2*radius^2*(x1 - x2)*(y1 - y2)^2)/((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));
            % d hx / d e2x
            dh_dr(i_v, e2_id) = 0.5 + (2*radius^2*(x1 - x2)*(y1 - y2)^2)/((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));

            % d hx / d e1y
            dh_dr(i_v, e1_id + N_cell) = -((-4*radius^2*(x1 - x2)^2 + (x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)^2)*(y1 - y2))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));
            % d hx / d e2y
            dh_dr(i_v, e2_id + N_cell) = ((-4*radius^2*(x1 - x2)^2 + (x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)^2)*(y1 - y2))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));

            % d hy / d e1x
            dh_dr(i_v + N_vertex, e1_id) = ((x1^4 - 4*x1^3*x2 + x2^4 - 4*x1*x2*(x2^2 + (y1 - y2)^2) + 2*x1^2*(3*x2^2 + (y1 - y2)^2) - 4*radius^2*(y1 - y2)^2 + 2*x2^2*(y1 - y2)^2 + (y1 - y2)^4)*(y1 - y2))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));
            % d hy / d e2x
            dh_dr(i_v + N_vertex, e2_id) = -((x1^4 - 4*x1^3*x2 + x2^4 - 4*x1*x2*(x2^2 + (y1 - y2)^2) + 2*x1^2*(3*x2^2 + (y1 - y2)^2) - 4*radius^2*(y1 - y2)^2 + 2*x2^2*(y1 - y2)^2 + (y1 - y2)^4)*(y1 - y2))/(2*(x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));

            % d hy / d e1y
            dh_dr(i_v + N_vertex, e1_id + N_cell) = 0.5 + (2*radius^2*(x1 - x2)*(y1 - y2)^2)/((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));
            % d hy / d e2y
            dh_dr(i_v + N_vertex, e2_id + N_cell) = 0.5 - (2*radius^2*(x1 - x2)*(y1 - y2)^2)/((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*sqrt(-((x1^2 - 2*x1*x2 + x2^2 + (y1 - y2)^2)*(y1 - y2)^2*(-4*radius^2 + x1^2 - 2*x1*x2 + x2^2 + y1^2 - 2*y1*y2 + y2^2))));
        
        end
        
    end
        
end


%% get dl/dr and da/dr and thus de/dr
dl_dtheta1 = sparse(dl_dtheta1);
dtheta1_dr = sparse(dtheta1_dr);
dl_dtheta2 = sparse(dl_dtheta2);
dtheta2_dr = sparse(dtheta2_dr);
dl_dh = sparse(dl_dh);
dh_dr = sparse(dh_dr);
clock_da_dh = sparse(clock_da_dh);
counter_da_dh = sparse(counter_da_dh);
da_dtheta1 = sparse(da_dtheta1);
da_dtheta2 = sparse(da_dtheta2);
cell_edge_adj = sparse(cell_edge_adj);
clock_cell_edge_adj = sparse(clock_cell_edge_adj);
counter_cell_edge_adj = sparse(counter_cell_edge_adj);

dln_dr = dl_dtheta1*dtheta1_dr + dl_dtheta2*dtheta2_dr;
dl_dr = dl_dh*dh_dr + dln_dr;

counter_da_dr = counter_da_dh*dh_dr + da_dtheta1*dtheta1_dr + da_dtheta2*dtheta2_dr;
clock_da_dr = clock_da_dh*dh_dr;

% (C x 2C) = (C x 2C) .* ((C x E) * (E x 2C)) + (C x E) * (E x 2C)
d_EP_dr = 2*repmat(K_P.*(perimeter_list-P0_list), 1, 2*N_cell).*(cell_edge_adj*dl_dr) + ...
			tension_difference*cell_edge_adj*dln_dr;

% (C x 2C) = (C x 2C) .* ((C x E) * (E x 2C) - (C x E) * (E x 2C))
d_EA_dr = 2*repmat(K_A.*(area_list-A0_list), 1, 2*N_cell).* ...
			(counter_cell_edge_adj*counter_da_dr-clock_cell_edge_adj*clock_da_dr);

center_force = -reshape(full(sum(d_EP_dr+d_EA_dr, 1)), [N_cell 2]);

% clear nan force
nan_cid = find(isnan(center_force(:, 1)));
for i_id = 1:numel(nan_cid)
    i_c = nan_cid(i_id);
    if numel(cell_chain{i_c})==1 % isolated
        center_force(i_c, :) = [0, 0];
    end
end


end % end of the whole function