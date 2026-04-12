function points = pbc_relocate(center, points, box_size)

if isempty(points), return, end

if size(points, 2) > size(center, 2)
    n_column = size(points, 2);
    n_dim = size(center, 2);
    origin_points = points;
    origin_point_mean = zeros(size(points, 1), n_dim);
    for i_c = 1:size(center, 2)
        origin_point_mean(:, i_c) = mean(points(:, i_c:n_dim:end), 2);
    end
    points = origin_point_mean;
end



for i = 1:size(points, 1)
    i_point = points(i, :);
    difference = (i_point - center)./box_size + 0.5;
     
    % we want this difference to be between 0 and 1
    % because (difference - 0.5) should be between -0.5 and 0.5
    points(i, :) = i_point - floor(difference).*box_size;
end
 
 
% average_center = mean(points, 1);
% for i = 1:size(points, 1)
%     i_point = points(i, :);
%     difference = (i_point - average_center)./box_size + 0.5;
%      
%     % we want this difference to be between 0 and 1
%     % because (difference - 0.5) should be between -0.5 and 0.5
%     points(i, :) = i_point - floor(difference).*box_size;
% end
 
 


if exist('origin_points', 'var')
    deviation = repmat(points - origin_point_mean, 1, n_column/n_dim);
    points = origin_points + deviation;
end

end % end of function