function Qd = noiseCovarianceMEKF(dT,w,a,Bw,Ba,Bm)

Qd = zeros(18,18);
I3 = eye(3);

Iw = w^2*I3; Ia = a^2*I3; IBw = Bw^2*I3; IBa = Ba^2*I3; IBm = Bm^2*I3;

Qd(1:3,1:3) = Iw*dT + IBw*dT^3/3;
Qd(10:12,1:3) = -IBw*dT^2/2;
Qd(4:6,4:6) = Ia*dT + IBa*dT^3/3;
Qd(7:9,4:6) = Ia*dT^2/2 + IBa*dT^4/8;
Qd(13:15,4:6) = -IBa*dT^2/2;
Qd(4:6,7:9) = IBa*dT^4/8 + Ia*dT^2/2;
Qd(7:9,7:9) = Ia*dT^3/3 + IBa*dT^5/20;
Qd(13:15,7:9) = -IBa*dT^3/6;
Qd(1:3,10:12) = -IBw*dT^2/2;
Qd(10:12,10:12) = IBw*dT^2/2;
Qd(4:6,13:15) = -IBa*dT^2/2;
Qd(7:9,13:15) = -IBa*dT^3/6;
Qd(13:15,13:15) = IBa*dT;
Qd(16:18,16:18) = IBm*dT;

end