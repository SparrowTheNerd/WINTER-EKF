clear; clc;

v1 = [9.64792-1.78739 0.18566 -1.78739];
v2 = [0 0 -9.80665];

v1_norm = v1 / norm(v1);
v2_norm = v2 / norm(v2);

axis = cross(v1_norm, v2_norm);
angle = acos(dot(v1_norm, v2_norm));

if norm(axis) == 0
    if dot(v1_norm, v2_norm) > 0
        q = quaternion(1, 0, 0, 0)
    else
        q = quaternion(0, 1, 0, 0)
    end
else
    axis_norm = axis / norm(axis);
    q = [cos(angle/2), sin(angle/2)*axis_norm(1), sin(angle/2)*axis_norm(2), sin(angle/2)*axis_norm(3)]
end

