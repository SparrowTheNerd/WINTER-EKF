% all calculation details in https://apps.dtic.mil/sti/tr/pdf/ADA588831.pdf

clear; clc; close;

SetupEnv();

% Kalman error states are quat angle, velocity, position, gyro bias, accel bias, mag bias
% inertial states are quaternion, position, velocity
% quaternion here represents a body to inertial rotation

tS = 0.01; %sample time (100hz)
[ORDat,dat] = ORDataImport(tS);

% sigX^2 is noise covariance
sigAccel = (0.005*9.81); % m/s^2 rms
sigGyro = deg2rad(0.1); % dps rms
sigMag = 0.0004; % Gauss rms
sigBaro = 3; % meters rms
sigGPS = 5; % meters rms

Rmag = eye(3)*sigMag^2;
Rbaro = sigBaro^2;
Rgps = eye(2)*sigGPS^2;
% sigGPSVel = 0.05; % m/s rms

% pnX is process noise
pnAccel = 0.05;
pnGyro = deg2rad(1);
pnBf = 0.05;
pnBw = 0.01;
pnBm = 0.0005;

G = 9.80665;

alt0 = atmospalt(dat.Prs(1)); %initial pressure altitude
qPos = [ORDat.qW(1) ORDat.qX(1) ORDat.qY(1) ORDat.qZ(1)]; % a-posteriori
wPos = [0 0 0]; aPos = [dat.aX dat.aY dat.aZ]; % previous IMU measurements
wBias = [0; 0; 0]; aBias = [0; 0; 0]; mBias = [0; 0; 0]; % IMU bias vectors

mag0 = [dat.mX(1) dat.mY(1) dat.mZ(1)]; %initial magnetic vector
mag0 = quatrotBI(qPos,mag0);

LL0 = [ORDat.lat(1) ORDat.lon(1)];

inertialState = [qPos'; 0; 0; alt0; 0; 0; 0;]; % initial rocket state
inertialStateList = zeros(10,size(ORDat,1));
inertialStateList(:,1) = inertialState;

errStateList = zeros(18,size(ORDat,1));

Pnn = zeros(18,18);
Qd = noiseCovarianceMEKF(tS,pnGyro,pnAccel,pnBw,pnBf,pnBm); % time-invariant noise covariance matrix

for i = 1:size(ORDat,1)-1

    % 1) measurement / storage
    w = [dat.gX(i) dat.gY(i) dat.gZ(i)];
    a = [dat.aX(i) dat.aY(i) dat.aZ(i)];
    % a-priori inertial state
    inertialStateN = integrationMEKF(inertialStateList(:,i),[a w],[aBias wBias],tS); 
    
    
    % 2) covariance prediction
    F = stateTransitionMEKF(inertialStateN(1:4),w,a,qPos,wPos,aPos);
    Phi = eye(18) + F*tS + 1/2*F^2*tS^2;
    Pn1n = Phi*Pnn*Phi'+Qd;
    % cond(Pn1n)

    if (mod(i,2) == 0)
        % 3 & 4) residual mappings & kalman gain
        % P = (I - K*H)*P*(I - K*H)' + K*R*K';
        I = eye(18);
        [dm, hm] = magMeasurementMEKF(inertialStateN(1:4)',[dat.mX(i) dat.mY(i) dat.mZ(i)],mag0);
        Km = Pn1n*hm' / (hm*Pn1n*hm' + Rmag);
        dXm = Km*dm';
        % Pm = (I-Km*hm)*Pn1n;
        Pm = (I - Km*hm)*Pn1n*(I-Km*hm)' + Km*Rmag*Km';

        alt = atmospalt(dat.Prs(i));
        [db, hb] = baroMeasurementMEKF(inertialStateN(7),alt-alt0);
        Kb = Pm*hb' / (hb*Pm*hb' + Rbaro);
        dXb = dXm + Kb*db;
        % Pb = (I-Kb*hb)*Pm;
        Pb = (I - Kb*hb)*Pm*(I-Kb*hb)' + Kb*Rbaro*Kb';


        if(mod(i,1) == 0)
            [dg, hg] = gpsMeasurementMEKF(inertialStateN(5),inertialStateN(6),dat.lat(i),dat.lon(i),LL0(1),LL0(2));
            Kg = Pb*hg' / (hg*Pb*hg' + Rgps);
            dXg = dXb + Kg*dg';
            % Pg = (I-Kg*hg)*Pb;
            Pg = (I - Kg*hg)*Pb*(I-Kg*hg)' + Kg*Rgps*Kg';

            Xnn = dXg;
            Pnn = Pg;
        else
            Xnn = dXb;
            Pnn = Pb;
        end

        % 5) update full states and calibrations
        inertialState = zeros(10,1);
        inertialState(1:4) = quatmultiply(inertialStateN(1:4)',[1 (Xnn(1:3)/2)'])';
        inertialState(1:4) = inertialState(1:4)/norm(inertialState(1:4));
        inertialState(5:7) = inertialStateN(5:7) + Xnn(7:9);
        inertialState(8:10) = inertialStateN(8:10) + Xnn(4:6);
        wBias = wBias + Xnn(10:12);
        aBias = aBias + Xnn(13:15);
        mBias = mBias + Xnn(16:18);

        errStateList(:,i) = Xnn;
        inertialStateList(:,i+1) = inertialState;

        % dtheta = Xnn(1:3);
        % G = eye(18);
        % G(1:3,1:3) = eye(3) - skew(dtheta);
        % Pnn = G * Pnn * G.';
        Xnn = 0;
        wPos = w; aPos = a; qPos = inertialState(1:4);

    else
        inertialState = inertialStateN;
        inertialStateList(:,i+1) = inertialState;
        Pnn = Pn1n;
    end


end
wBias
aBias
mBias
quatList = inertialStateList(1:4,:);
posList = inertialStateList(5:7,:);
%% Figure Plotting
figure(1);
subplot(3,1,1);
plot(dat.Time(1:end),inertialStateList(5,1:end));
xlabel("Time (s)");
ylabel("Crossrange X (m)");
title("Crossrange X vs Time");
hold on;
plot(ORDat.t,ORDat.relPosX);
legend("Simulated","Truth")
grid on;
subplot(3,1,2);
plot(dat.Time(1:end),inertialStateList(6,1:end))
xlabel("Time (s)");
ylabel("Crossrange Y (m)");
title("Crossrange Y vs Time");
hold on;
plot(ORDat.t,ORDat.relPosY);
legend("Simulated","Truth")
grid on;
subplot(3,1,3);
plot(dat.Time(1:end),inertialStateList(7,1:end))
xlabel("Time (s)");
ylabel("Altitude (m)");
title("Altitude vs Time");
hold on;
plot(ORDat.t,ORDat.relPosZ);
legend("Simulated","Truth")
grid on;

eulAngs = zeros(height(dat),3);
eulAngsTrue = zeros(height(dat),3);
for i = 1:height(dat)
    eulAngs(i,:) = quat2eul(quatList(:,i)',"ZYX")*180/pi;
    eulAngsTrue(i,:) = quat2eul([ORDat.qW(i) ORDat.qX(i) ORDat.qY(i) ORDat.qZ(i)],"ZYX")*180/pi;
end

figure(2);
axP = subplot(3,1,1);
plot(dat.Time,eulAngs(:,1));
xlabel("Time (s)");
ylabel("Z-Axis (deg)");
title("Z-Axis Angle vs Time");
hold on;
plot(ORDat.t,eulAngsTrue(:,1));
legend("Simulated","Truth")
grid on;
axY = subplot(3,1,2);
plot(dat.Time,eulAngs(:,2));
xlabel("Time (s)");
ylabel("Y-Axis (deg)");
title("Y-Axis Angle vs Time");
hold on;
plot(ORDat.t,eulAngsTrue(:,2));
legend("Simulated","Truth")
grid on;
axR = subplot(3,1,3);
plot(dat.Time,eulAngs(:,3));
% hold on;
% plot(dat.Time,dat.mY*100);
xlabel("Time (s)");
ylabel("X-Axis (deg)");
title("X-Axis vs Time");
hold on;
plot(ORDat.t,eulAngsTrue(:,3));
legend("Simulated","Truth")
% linkaxes([axP axY axR],'y')
% axP.YLim = [-180 180];
grid on;

figure(3)
qwp = subplot(4,1,1);
plot(dat.Time,inertialStateList(1,1:end));
xlabel("Time (s)");
title("Quat W vs Time");
hold on;
plot(ORDat.t,ORDat.qW);
legend("Simulated","Truth")
grid on;
qxp = subplot(4,1,2);
plot(dat.Time,inertialStateList(2,1:end));
xlabel("Time (s)");
title("Quat X vs Time");
hold on;
plot(ORDat.t,ORDat.qX);
legend("Simulated","Truth")
grid on;
qyp = subplot(4,1,3);
plot(dat.Time,inertialStateList(3,1:end));
xlabel("Time (s)");
title("Quat Y vs Time");
hold on;
plot(ORDat.t,ORDat.qY);
legend("Simulated","Truth")
grid on;
qzp = subplot(4,1,4);
plot(dat.Time,inertialStateList(4,1:end));
xlabel("Time (s)");
title("Quat Z vs Time");
hold on;
plot(ORDat.t,ORDat.qZ);
legend("Simulated","Truth")
grid on;
linkaxes([qwp qxp qyp qzp],'y')
qwp.YLim = [-1 1];

%%
run_animation(quatList, posList, dat) ;

function run_animation(quatList, posList, dat)   
    
% v = VideoWriter('flight_of_fuckin_bumblebee.mp4','MPEG-4');
% v.FrameRate = 60;
% open(v);

    q = quaternion(quatList(:,1)');
    patch = poseplot(q);
    
    patch.ScaleFactor = 150;
    xlabel("X")
    ylabel("Y")
    zlabel("Z")
    
    for i = 1:5:height(dat)
        q = quatList(:,i);
        q = [q(1) q(2) q(3) q(4)];
        q = quaternion(q);
        pos = posList(:,i);
        set(patch, Orientation=q, Position=pos); hold on
        plot3(pos(1), pos(2), pos(3), '.b', 'MarkerSize', 2); hold on
    
        xlim([-100, 600]);
        ylim([-100, 600]);
        zlim([0, 4500]);
    
        set(gca,'ZDir','normal')  
        title(sprintf("t = %0.2f", dat.Time(i)))
        drawnow
        % frame = getframe(gcf);
        % writeVideo(v, frame);
    end
% close(v);
end