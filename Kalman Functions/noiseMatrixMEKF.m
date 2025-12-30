function G = noiseMatrixMEKF(q,sigw,sigf,sigbw,sigbf,sigbm)

Cbi = quatDCM(q);
I3 = eye(3);

G = zeros(18,18);

G(1:3,1:3) = -sigw*I3;
G(4:6,4:6) = -Cbi*(sigf*I3);
G(10:12,10:12) = sigbw*I3;
G(13:15,13:15) = sigbf*I3;
G(16:18,16:18) = sigbm*I3;

end