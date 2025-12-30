function x = integrationMEKF(x,u,bias,dT)
% Quaternions here represent a body-to-global rotation. Accelerations will
% be rotated into the global frame for position and velocity

qw = x(1); qx = x(2); qy = x(3); qz = x(4);
px = x(5); py = x(6); pz = x(7);
vx = x(8); vy = x(9); vz = x(10);

axb = bias(1); ayb = bias(2); azb = bias(3);
wxb = bias(4); wyb = bias(5); wzb = bias(6);

ax = u(1)-axb; ay = u(2)-ayb; az = u(3)-azb;
wx = u(4)-wxb; wy = u(5)-wyb; wz = u(6)-wzb;

% Quaternion Integration:
% https://ahrs.readthedocs.io/en/latest/filters/ekf.html#prediction-step 
qw1 = qw - (dT/2)*wx*qx - (dT/2)*wy*qy - (dT/2)*wz*qz;
qx1 = qx + (dT/2)*wx*qw - (dT/2)*wy*qz + (dT/2)*wz*qy;
qy1 = qy + (dT/2)*wx*qz + (dT/2)*wy*qw - (dT/2)*wz*qx;
qz1 = qz - (dT/2)*wx*qy + (dT/2)*wy*qx + (dT/2)*wz*qw;

qNorm = norm([qw1 qx1 qy1 qz1]);
qw1 = qw1/qNorm;
qx1 = qx1/qNorm;
qy1 = qy1/qNorm;
qz1 = qz1/qNorm;

% Rotate accelerations into the global frame
% for some reason, the aerospace toolbox quatrotate expects a world-to-body
% quaternion, so we must do the conjugate
ag = quatrotBI([qw1 qx1 qy1 qz1],[ax ay az]);

ag(3) = ag(3)+9.80665;

% Integrate acceleration for position and velocity
vx1 = vx + ag(1)*dT; 
vy1 = vy + ag(2)*dT;
vz1 = vz + ag(3)*dT;
px1 = px + vx1*dT;
py1 = py + vy1*dT;
pz1 = pz + vz1*dT;

x = [qw1;
     qx1;
     qy1;
     qz1;
     px1;
     py1;
     pz1;
     vx1;
     vy1;
     vz1];
end