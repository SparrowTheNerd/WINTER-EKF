function h = magMeasurementFcn(x,mag0)
    q_state = [x(1) x(2) x(3) x(4)];
    h = quatrotate(q_state,mag0);
end