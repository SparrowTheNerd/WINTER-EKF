% creates state transition matrix F for the MEKF implementation
% p suffix is previous iterations a-posteriori estimate and IMU outputs

function F = stateTransitionMEKF(q,w,f,qp,wp,fp)

I3 = eye(3);
F = zeros(18,18);

qAvg = quatAverage(qp,q);
wx = -skew((wp+w)/2);
Cbif = ( (-quatDCM(q)*skew(fp)) + (-quatDCM(qp)*skew(f)) ) / 2;
Cbi = -quatDCM(qAvg);

F(1:3,1:3) = wx; 
F(4:6,1:3) = Cbif; 
F(7:9,4:6) = I3;
F(1:3,10:12) = -I3; 
F(4:6,13:15) = Cbi;

end