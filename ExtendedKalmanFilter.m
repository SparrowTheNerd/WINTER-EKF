clear; clc; close;

% States are wxyz quaternion, xyz position, xyz velocity
% Control inputs are xyz accelerometer, xyz gyroscope
% Measurements are xyz GPS, x barometer

tS = 0.01; %sample time (100hz)

sigAccel = (0.002*9.81)^2; % m/s^2 rms
sigGyro = 0.1^2; % dps rms
sigMag = 0.11^2; % degrees rms
sigBaro = 1^2; % meters rms
sigGPS = 0.5^2; % meters rms
sigGPSVel = 0.05^2; % m/s rms


%0.1 for quats, 10m position, 10m/s velocity
Q = [0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 10.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 10.0 0.0;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 10.0];

R_baro = [sigBaro];
R_gps = [sigGPS 0 0;
         0 sigGPS 0;
         0 0 sigGPS];
R_gpsVel = [sigGPSVel 0 0;
            0 sigGPSVel 0;
            0 0 sigGPSVel];

dat = readtable('SampleFlightData.csv');

Xnn = [1,0,0,0,0,0,0,0,0,0];
Pnn = [10 0 0 0 0 0 0 0 0 0;
       0 10 0 0 0 0 0 0 0 0;
       0 0 10 0 0 0 0 0 0 0;
       0 0 0 10 0 0 0 0 0 0;
       0 0 0 0 50 0 0 0 0 0;
       0 0 0 0 0 50 0 0 0 0;
       0 0 0 0 0 0 50 0 0 0;
       0 0 0 0 0 0 0 50 0 0;
       0 0 0 0 0 0 0 0 50 0;
       0 0 0 0 0 0 0 0 0 50];

statesX = zeros(10,height(dat));
statesX(:,1) = Xnn;
for i = 1:height(dat)-1
    % Prediction Step
    dT = (dat.Time(i+1)-dat.Time(i)); %timestep
    u = [dat.aX(i) dat.aY(i) dat.aZ(i) dat.gX(i)*pi/180 dat.gY(i)*pi/180 dat.gZ(i)*pi/180]; %control vector
    Xn1n = stateTransitionFcn(Xnn,u,dT); %Xn+1,n
    J = @(x) stateTransitionFcn(x,u,dT); %function definition for taking jacobian numerically
    dF = numericalJacobian(J,Xnn); %jacobian dFdX
    
    Pn1n = dF*Pnn*dF' + Q;
    
    %Correction Step
    if(~isnan(dat.Prs(i)))
        zn = 44330*(1-(dat.Prs(i)/101325)^(0.190284))-188;
        h = Xn1n(5);
        dh = [0 0 0 0 1 0 0 0 0 0];

        Kn = Pn1n * dh' / (dh*Pn1n*dh'+R_baro);
        Pnn = (eye(10)-Kn*dh)*Pn1n*((eye(10)-Kn*dh)')+Kn*R_baro*Kn';
        Xnn = Xn1n + Kn*(zn-h);
    else
        Xnn = Xn1n;
    end
    statesX(:,i+1) = Xn1n;
end

posList = statesX(5:7,:);
quatList = statesX(1:4,:);

% v = VideoWriter('ekf_animation_justlaunch.mp4','MPEG-4');
% v.FrameRate = 30;
% open(v);
% for i = 1:height(dat)-2100
%     q = quaternion(quatList(:,i)');
%     quatRot = quaternion(sqrt(2),0,sqrt(2),0); %rotate so that X is up in animation
%     qRot = quatRot*q;
%     pos = quatrotate(quatRot,posList(:,i)');
%     poseplot(qRot,pos,'ScaleFactor',0.2);
% 
%     timeStr = sprintf('Time: %.2f s', dat.Time(i));
%     title(timeStr, 'FontSize', 14);
% 
%     drawnow;
% 
%     % frame = getframe(gcf);
%     % writeVideo(v, frame);
% end
% close(v);

figure(1);
plot(dat.Time,statesX(5,:));
xlabel("Time (s)");
ylabel("Altitude (m)");
title("Altitude vs Time");

figure(2);
plot(dat.Time,statesX(6,:))
xlabel("Time (s)");
ylabel("Crossrange X (m)");
title("Crossrange X vs Time");

figure(3);
plot(dat.Time,statesX(7,:))
xlabel("Time (s)");
ylabel("Crossrange Y (m)");
title("Crossrange Y vs Time");

quatStates = statesX(1:4,:);
eulAngs = zeros(height(dat),3);
for i = 1:height(dat)
    eulAngs(i,:) = quat2eul(quatStates(:,i)',"XYZ")*180/pi;
end

figure(4);
plot(dat.Time(1:500),eulAngs((1:500),1))
xlabel("Time (s)");
ylabel("Roll (deg)");
title("Roll Angle vs Time");
figure(5);
plot(dat.Time(1:500),eulAngs((1:500),2))
xlabel("Time (s)");
ylabel("Pitch (deg)");
title("Pitch Angle vs Time");
figure(6);
plot(dat.Time(1:500),eulAngs((1:500),3))
xlabel("Time (s)");
ylabel("Yaw (deg)");
title("Yaw Angle vs Time");