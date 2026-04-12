function [cell_chain, edgelist, vertex_position, area_list, perimeter_list] = ...
                    make_finite_voronoi_pbc(center_xy, radius, box_size)
% rely on: pbc_relocate

N_cell = size(center_xy, 1);
if ~exist('box_size', 'var')
    box_size = sqrt(N_cell);
end


center_xy = [center_xy; ...
        center_xy + box_size*[1 0]; ...
        center_xy - box_size*[1 0]];

center_xy = [center_xy; ...
        center_xy + box_size*[0 1]; ...
        center_xy - box_size*[0 1]];
			


DT = delaunayTriangulation(center_xy);
[big_vertex_position, big_cell_chain] = voronoiDiagram(DT);



cell_chain = big_cell_chain(1:N_cell);

central_to_big_vid_dict = find(all(big_vertex_position>=0, 2)&all(big_vertex_position<box_size, 2));

big2central_vid = zeros(1, size(big_vertex_position, 1));
big2central_vid(central_to_big_vid_dict) = 1:numel(central_to_big_vid_dict);

vertex_position = mod(big_vertex_position(central_to_big_vid_dict, :), box_size);

outside_v_candidate = zeros(20*N_cell, 1);
for i_c = 1:N_cell
    outside_v_candidate(cell_chain{i_c}) = 1;
end
outside_v_candidate(central_to_big_vid_dict) = 0;
outside_vid = find(outside_v_candidate);

shifted_ov_xy = mod(big_vertex_position(outside_vid, :), box_size);
for i_o = 1:numel(outside_vid)
    i_o_xy = shifted_ov_xy(i_o, :);
    [~, min_pos] = min(sum(abs(i_o_xy - vertex_position), 2));
    big2central_vid(outside_vid(i_o)) = min_pos;
end




cell_chain = cellfun(@(c) big2central_vid(c), cell_chain, 'uniformoutput', 0);
% now we have a regular Voronoi tessallation
N_origin_vertex = size(vertex_position, 1);
% cell_vertex_adj = zeros(N_cell, N_origin_vertex);
vertex_position = [vertex_position, 3*ones(N_origin_vertex, 1), zeros(N_origin_vertex, 3)];
for i_c = 1:N_cell
    i_chain = cell_chain{i_c};
    %remove duplicated vertices
    if numel(i_chain)>numel(unique(i_chain))
        for pos = numel(i_chain):-1:1
            if (pos>1 && i_chain(pos) == i_chain(pos-1)) || (pos==1 && i_chain(pos) == i_chain(end))
                i_chain(pos) = [];
            end
        end
        cell_chain{i_c} = i_chain;
    end
    i_chain_len = numel(i_chain);
    for i_v_idx = 1:i_chain_len
        i_v = i_chain(i_v_idx);
        if ~vertex_position(i_v, 4)
            vertex_position(i_v, 4) = i_c;
        elseif ~vertex_position(i_v, 5)
            vertex_position(i_v, 5) = i_c;
        else
            vertex_position(i_v, 6) = i_c;
        end
    end
end



% % next construct bound Voronoi
% split far vertex
% [out_vid, in_vid, new_vid, new_v_xy, c1, c2]
new_vertex_mat = zeros(10*N_cell, 7);
new_vertex_count = 0;
% edgelist format:
% (straight edge, shared by two cells): [1/2/3, v1id, v2id, cell1, cell 2, old v1, old v2]
% (arc edge, owned by one cell): [4, v1id, v2id, cell1, 0, v1 angle, v2 angle]
edgelist = zeros(20*N_cell, 7);
inner_vertex_bool = ones(size(vertex_position, 1), 1);
edge_count = 0;
vid2row = zeros(N_origin_vertex^2+N_origin_vertex, 1);
old_vid2row = zeros(N_origin_vertex^2+N_origin_vertex, 1);
% if stay 1 after loop, then it is a real inner vertex
for i_c = 1:N_cell
    xc = center_xy(i_c, 1);
    yc = center_xy(i_c, 2);
    chain = cell_chain{i_c}; % counterclockwise vertex chain
    
    Ximg = pbc_relocate([xc yc], vertex_position(chain, 1:2), box_size);
    v_dist = vecnorm([xc yc]-Ximg, 2, 2);
    across_check = 2*(v_dist>radius)+([v_dist(2:end); v_dist(1)]>radius);
    % in-in: 0; in-out: 1; out-in: 2; out-out: 3.
    chain_len = numel(chain);
    
    last_vertex_info = [1 0]; % in: 1; circle: 0.
    head_vertex_info = [];
    new_chain = zeros(1, 50);
    new_chain_pointer = 1;
    for i_v = 1:chain_len
        i_v_id = chain(i_v);
        i_next_pos = mod(i_v, chain_len)+1;
        i_next_id = chain(i_next_pos);
            
        switch across_check(i_v)

            case 0 % in-in, result in one straight edge
                if i_v_id<i_next_id
                    edge_sign = 1;
                else
                    edge_sign = -1;
                end
                edge_ends = sort([i_v_id, i_next_id]);
                vid_key = edge_ends(1)*N_origin_vertex+edge_ends(2);
                edge_row_pos = vid2row(vid_key);
                is_old_edge = edge_row_pos>0;
                if is_old_edge
                    new_chain(new_chain_pointer) = edge_sign*edge_row_pos;
                    new_chain_pointer = new_chain_pointer + 1;
                    edgelist(edge_row_pos, 5) = i_c;
                else
                    % add straight
                    edgelist(edge_count+1, :) = [1, edge_ends, i_c, 0, edge_ends]; % leave cell2 id empty
                    vid2row(vid_key) = edge_count+1;
                    edge_count = edge_count+1;
                    new_chain(new_chain_pointer) = edge_sign*edge_count;
                    new_chain_pointer = new_chain_pointer + 1;
                end
                last_vertex_info = [1 i_next_id];
                if isempty(head_vertex_info)
                    head_vertex_info = [1 i_v_id];
                end

            case {1, 2} % 1: in-out; 2: out-in, result in one straight edge and one arc edge
                if across_check(i_v)==1
                    in_vid = i_v_id;
                    out_vid = i_next_id;
                    in_pos = i_v;
                    out_pos = i_next_pos;
                    edge_sign = 1;
                    head_sign = 1;
                else
                    in_vid = i_next_id;
                    out_vid = i_v_id;
                    in_pos = i_next_pos;
                    out_pos = i_v;
                    edge_sign = -1;
                    head_sign = 0;
                end
                inner_vertex_bool(out_vid) = 0;

                old_vid_key = in_vid*N_origin_vertex+out_vid; % [in_vid out_vid];
                edge_row_pos = old_vid2row(old_vid_key);
                is_old_edge = edge_row_pos>0;
                if is_old_edge
                    new_v_id = edgelist(edge_row_pos, 3);
                    new_v_xy = pbc_relocate([xc yc], new_vertex_mat(new_v_id-N_origin_vertex, 4:5), box_size);
                    new_vertex_mat(new_v_id-N_origin_vertex, 7) = i_c;
                else % solve new vertex xy and add a row to the edgelist
                    out_x = Ximg(out_pos, 1);
                    in_x = Ximg(in_pos, 1);
                    out_y = Ximg(out_pos, 2);
                    in_y = Ximg(in_pos, 2);
                    A = in_y-out_y;
                    B = out_x-in_x;
                    C = in_x*out_y - out_x*in_y;
                    d = (A*xc+B*yc+C)/sqrt(A^2+B^2);
                    foot_x = (B*B*xc - A*B*yc - A*C) / (A*A + B*B);
                    foot_y = (-A*B*xc + A*A*yc - B*C) / (A*A + B*B);
                    foot_to_new_len = sqrt(radius^2-d^2);
                    foot_to_out_len = sqrt((out_x-foot_x)^2 + (out_y-foot_y)^2);
                    new_v_xy = [foot_x foot_y] + foot_to_new_len/foot_to_out_len ...
                                        *[out_x-foot_x, out_y-foot_y];

                    new_v_id = N_origin_vertex+new_vertex_count+1;
                    new_vertex_mat(new_vertex_count+1, :) = [in_vid, out_vid, new_v_id, new_v_xy, i_c, -1];
                    new_vertex_count = new_vertex_count + 1;
                end

                if across_check(i_v)==2 && last_vertex_info(1)==0 % circle - circle -> arc edge
                    last_v_id = last_vertex_info(2);
                    last_v_xy = new_vertex_mat(last_v_id - N_origin_vertex, 4:5);
                    last_v_xy = pbc_relocate([xc yc], last_v_xy, box_size);
                    last_v_angle = get_angle([xc yc], last_v_xy);
                    new_v_angle = get_angle([xc yc], new_v_xy);
                    if new_v_angle<last_v_angle
                        new_v_angle = new_v_angle + 2*pi;
                    end
                    % add arc
                    edgelist(edge_count+1, :) = [4, last_v_id, new_v_id, i_c, -1, last_v_angle, new_v_angle];
                    edge_count = edge_count+1;
                    new_chain(new_chain_pointer) = edge_count;
                    new_chain_pointer = new_chain_pointer + 1;
                elseif isempty(head_vertex_info)
                    head_vertex_info = [head_sign new_v_id];
                end

                last_vertex_info = [0 new_v_id];
                if is_old_edge
                    edgelist(edge_row_pos, 5) = i_c;
                    new_chain(new_chain_pointer) = edge_sign*edge_row_pos;
                    new_chain_pointer = new_chain_pointer + 1;
                else
                    % add straight
                    edgelist(edge_count+1, :) = [2, in_vid, new_v_id, i_c, 0, in_vid, out_vid];
                    old_vid2row(old_vid_key) = edge_count+1;
                    edge_count = edge_count+1;
                    new_chain(new_chain_pointer) = edge_sign*edge_count;
                    new_chain_pointer = new_chain_pointer + 1;
                end


            case 3 % out-out, if last vertex is a circle, and a squeeze new edge, then another arc edge

                inner_vertex_bool([i_v_id i_next_id]) = 0;
                if i_v_id<i_next_id
                    edge_sign = 1;
                else
                    edge_sign = -1;
                end
                old_ends_pair = sort([i_v_id, i_next_id]);
                old_vid_key = old_ends_pair(1)*N_origin_vertex+old_ends_pair(2);
                edge_row_pos = old_vid2row(old_vid_key);
                is_old_edge = edge_row_pos>0;
                if is_old_edge
                    new_v1_id = edgelist(edge_row_pos, 2+(i_v_id>i_next_id));
                    new_v1_xy = new_vertex_mat(new_v1_id-N_origin_vertex, 4:5);
                    new_v1_xy = pbc_relocate([xc yc], new_v1_xy, box_size);
                    new_v2_id = edgelist(edge_row_pos, 2+(i_v_id<=i_next_id));
                    new_vertex_mat(new_v1_id-N_origin_vertex, 7) = i_c;
                    new_vertex_mat(new_v2_id-N_origin_vertex, 7) = i_c;
                else % solve new vertex xy and add a row to new_vertex_mat
                    v1_x = Ximg(i_v, 1);
                    v2_x = Ximg(i_next_pos, 1);
                    v1_y = Ximg(i_v, 2);
                    v2_y = Ximg(i_next_pos, 2);
                    A = v2_y-v1_y;
                    B = v1_x-v2_x;
                    C = v2_x*v1_y - v1_x*v2_y;
                    d = (A*xc+B*yc+C)/sqrt(A^2+B^2);
                    if abs(d)>radius, continue, end
                    foot_x = (B*B*xc - A*B*yc - A*C) / (A*A + B*B);
                    foot_y = (-A*B*xc + A*A*yc - B*C) / (A*A + B*B);
                    if (foot_x-v1_x)*(foot_x-v2_x)>0, continue, end
					% check whether the mirror location has a cell
                    mirror_c_xy = [2*foot_x-xc, 2*foot_y-yc];
                    shifted_all_cxy = pbc_relocate([xc yc], center_xy(1:N_cell, 1:2), box_size);
                    if min(vecnorm(mirror_c_xy-shifted_all_cxy, 2, 2))>1e-5
                        continue
                    end
					
                    foot_to_new_len = sqrt(radius^2-d^2);
                    foot_to_v1_len = sqrt((v1_x-foot_x)^2 + (v1_y-foot_y)^2);
                    new_v1_xy = [foot_x foot_y] + foot_to_new_len/foot_to_v1_len ...
                                        *[v1_x-foot_x, v1_y-foot_y];
                    new_v2_xy = [foot_x foot_y] - foot_to_new_len/foot_to_v1_len ...
                                        *[v1_x-foot_x, v1_y-foot_y];

                    new_v1_id = N_origin_vertex+new_vertex_count+1;
                    new_v2_id = N_origin_vertex+new_vertex_count+2;
                    new_vertex_mat(new_vertex_count+1:new_vertex_count+2, :) = ...
                        [i_v_id, i_next_id, new_v1_id, new_v1_xy, i_c, 0; ...
                        i_next_id, i_v_id, new_v2_id, new_v2_xy, i_c, 0];
                    new_vertex_count = new_vertex_count + 2;
                end

                if last_vertex_info(1)==0 % circle - circle -> arc edge
                    last_v_id = last_vertex_info(2);
                    last_v_xy = new_vertex_mat(last_v_id - N_origin_vertex, 4:5);
                    last_v_xy = pbc_relocate([xc yc], last_v_xy, box_size);
                    last_v_angle = get_angle([xc yc], last_v_xy);
                    new_v_angle = get_angle([xc yc], new_v1_xy);
                    if new_v_angle<last_v_angle
                        new_v_angle = new_v_angle + 2*pi;
                    end
                    % add arc
                    edgelist(edge_count+1, :) = [4, last_v_id, new_v1_id, i_c, -1, last_v_angle, new_v_angle];
                    edge_count = edge_count+1;
                    new_chain(new_chain_pointer) = edge_count;
                    new_chain_pointer = new_chain_pointer + 1;
                elseif isempty(head_vertex_info)
                    head_vertex_info = [0 new_v1_id];
                end

                last_vertex_info = [0 new_v2_id];
                if is_old_edge
                    new_chain(new_chain_pointer) = edge_sign*edge_row_pos;
                    new_chain_pointer = new_chain_pointer + 1;
                    edgelist(edge_row_pos, 5) = i_c;
                else
                    if i_v_id<i_next_id
                        correspond_new_v_pair = [new_v1_id new_v2_id];
                    else
                        correspond_new_v_pair = [new_v2_id new_v1_id];
                    end
                    % add straight
                    edgelist(edge_count+1, :) = [3, correspond_new_v_pair, i_c, 0, old_ends_pair];
                    old_vid2row(old_vid_key) = edge_count+1;
                    edge_count = edge_count+1;
                    new_chain(new_chain_pointer) = edge_sign*edge_count;
                    new_chain_pointer = new_chain_pointer + 1;
                end


        end   

    end



    % deal with possible tail - head arc edge
    if ~isempty(head_vertex_info) && head_vertex_info(1)==0 ...
                    && last_vertex_info(1)==0
        last_v_id = last_vertex_info(2);
        last_v_xy = new_vertex_mat(last_v_id - N_origin_vertex, 4:5);
        last_v_xy = pbc_relocate([xc yc], last_v_xy, box_size);
        last_v_angle = get_angle([xc yc], last_v_xy);
        head_v_id = head_vertex_info(2);

        head_v_xy = new_vertex_mat(head_v_id - N_origin_vertex, 4:5);
        head_v_xy = pbc_relocate([xc yc], head_v_xy, box_size);
        head_v_angle = get_angle([xc yc], head_v_xy);
        if head_v_angle<last_v_angle
            head_v_angle = head_v_angle + 2*pi;
        end
        % add arc
        edgelist(edge_count+1, :) = [4, last_v_id, head_v_id, i_c, -1, last_v_angle, head_v_angle];
        edge_count = edge_count+1;
        new_chain(new_chain_pointer) = edge_count;
        new_chain_pointer = new_chain_pointer + 1;
    elseif isempty(head_vertex_info)
        % add circle
        edgelist(edge_count+1, :) = [4, 0, 0, i_c, -1, 0, 2*pi];
        edge_count = edge_count+1;
        new_chain(new_chain_pointer) = edge_count;
        new_chain_pointer = new_chain_pointer + 1;
    end
    
    cell_chain{i_c} = new_chain(1:new_chain_pointer-1);

end
edgelist(edge_count+1:end, :) = [];
new_vertex_mat(new_vertex_count+1:end, :) = [];
% edgelist(edgelist(:, 4)==98, :)
% edgelist(edgelist(:, 5)==338, :)
% vertex_position([414 102], :)
%% build a new vid dict
vertex_position = [vertex_position; new_vertex_mat(:, 4:5), ...
    2*ones(new_vertex_count, 1), new_vertex_mat(:, 6:7), zeros(new_vertex_count, 1)];

bound_vid_dict = [0; (1:N_origin_vertex).'; ...
    new_vertex_mat(:, 3)];
sort_out_pool = unique([new_vertex_mat(:, 2); find(inner_vertex_bool==0)]);
for i_o = numel(sort_out_pool):-1:1
    i_out = sort_out_pool(i_o);
    bound_vid_dict(i_out+2:end) = bound_vid_dict(i_out+2:end) - 1;
end
vertex_position(sort_out_pool, :) = [];
edgelist(:, 2:3) = bound_vid_dict(edgelist(:, 2:3)+1);
edgelist(:, 1) = round(edgelist(:, 1), 3)==4; % 0: straight; 1: arc.



if nargout <= 3
    return
end

%% calculate the perimeter and area lists
perimeter_list = zeros(N_cell,1);
area_list = zeros(N_cell,1);
for i_e = 1:size(edgelist, 1)
    
    i_edge_info = edgelist(i_e, :);
    i_c = i_edge_info(4);
        
    if ~i_edge_info(1) % straight edge
        i_c_pair = i_edge_info(4:5);
        i_center = center_xy(i_c, :);
        Ximg = pbc_relocate(i_center, vertex_position(i_edge_info(2:3), 1:2), box_size);
        % add perimeter
        perimeter_list(i_c_pair) = perimeter_list(i_c_pair) + norm(Ximg(1, :)-Ximg(2, :));

        % add area
        area_list(i_c_pair) = area_list(i_c_pair) + 0.5*abs( ...
            i_center(1)*(Ximg(1, 2)-Ximg(2, 2)) + ...
            Ximg(1, 1)*(Ximg(2, 2)-i_center(2)) + ...
            Ximg(2, 1)*(i_center(2)-Ximg(1, 2)) );
    else % arc edge
        perimeter_list(i_c) = perimeter_list(i_c) + (i_edge_info(7)-i_edge_info(6))*radius;
        area_list(i_c) = area_list(i_c) + 0.5*(i_edge_info(7)-i_edge_info(6))*radius^2;
    end

end

end

