function [dg, hg] = gpsMeasurementMEKF(pX,pY,lat,lon,lat0,lon0)

% returns [E N U] in meters
denu = dllh2denu([lat0,lon0,0],[lat,lon,0]);
dg = [denu(1:2)] - [pX pY];

hg = zeros(2,18);
hg(1,7) = 1; hg(2,8) = 1;

end