function Qd = noiseCovarianceMEKF(dT,w,f,Bw,Bf,Bm)

Qd = zeros(18,18);
I3 = eye(3);

Iw = w^2*I3; If = f^2*I3; IBw = Bw^2*I3; IBf = Bf^2*I3; IBm = Bm^2*I3;

Qd(1:3,1:3) = Iw*dT + IBw*dT^3/3;
Qd(10:12,1:3) = -IBw*dT^2/2;
Qd(4:6,4:6) = If*dT + IBf*dT^3/3;
Qd(7:9,4:6) = If*dT^2/2 + IBf*dT^4/8;
Qd(13:15,4:6) = -IBf*dT^2/2;
Qd(4:6,7:9) = IBf*dT^4/8 + If*dT^2/2;
Qd(7:9,7:9) = If*dT^3/3 + IBf*dT^5/20;
Qd(13:15,7:9) = -IBf*dT^3/6;
Qd(1:3,10:12) = -IBw*dT^2/2;
Qd(10:12,10:12) = IBw*dT^2/2;
Qd(4:6,13:15) = -IBf*dT^2/2;
Qd(7:9,13:15) = -IBf*dT^3/6;
Qd(13:15,13:15) = IBf*dT;
Qd(16:18,16:18) = IBm*dT;

end