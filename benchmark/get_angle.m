function angle_list = get_angle(center_xy, end_xy)

N = size(end_xy, 1);
angle_list = zeros(N, 1);

if N>1 && numel(center_xy)==2
    center_xy = repmat(center_xy, N, 1);
end


for i = 1:N
    if norm(end_xy(i, :)-center_xy(i, :))==0
        angle_list(i) = 0;
        continue
    end

    angle = atan( (end_xy(i, 2)-center_xy(i, 2))/(end_xy(i, 1)-center_xy(i, 1)) );
    if end_xy(i, 1)<center_xy(i, 1)
        angle = mod(angle+pi, 2*pi);
    else
        angle = mod(angle, 2*pi);
    end
    
    angle_list(i) = angle;
end



end % end of function