function [dm, hm] = magMeasurementMEKF(q,m,m0)

mPred = quatrotIB(q,m0);
dm = (m - mPred);
hm = zeros(3,18);
hm(1:3,1:3) = skew(quatDCM(quatconj(q'))*m0');
hm(1:3,16:18) = -eye(3);

end