function h = magMeasurementFcn(x,mag0)
    q_conj = [x(1) x(2) x(3) x(4)];
    
    globalquat = [0, mag0];

    temp = quatmultiply(q_conj,globalquat);

    rotQuat = quatmultiply(temp,[x(1) -x(2) -x(3) -x(4)]);

    h = rotQuat(2:4);
end