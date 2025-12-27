% creates state transition matrix F for the MEKF implementation

function F = stateTransitionMEKF(q,w,f)
Cbi = quatDCM(q);
I3 = eye(3);

F = zeros(18,18);

F(1:3,1:3) = skew(w); F(4:6,1:3) = -Cbi*skew(f); F(7:9,4:6) = I3;
F(1:3,10:12) = -I3; F(4:6,13:15) = -Cbi;

end