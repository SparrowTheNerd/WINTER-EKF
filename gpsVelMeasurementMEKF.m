function [dv,hv] = gpsVelMeasurementMEKF(v, gpsV)

    dv = gpsV - v;
    hv = zeros(2,18);
    hv(1,4) = 1;
    hv(2,5) = 1;

end