function [db,hb] = baroMeasurementMEKF(pZ,baro)

    db = baro-pZ;
    hb = zeros(1,18);
    hb(9) = 1;
end