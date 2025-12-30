clear; clc; close;

SetupEnv();

% States are wxyz quaternion, xyz position, xyz velocity
% Control inputs are xyz accelerometer, xyz gyroscope
% Measurements are xyz GPS, x barometer

tS = 0.01; %sample time (100hz)
[ORDat,dat] = ORDataImport(tS);

sigAccel = (0.005*9.81)^2; % m/s^2 rms
sigGyro = deg2rad(0.1)^2; % dps rms
sigMag = 0.0004^2; % Gauss rms
sigBaro = 3^2; % meters rms
sigGPS = 0.5^2; % meters rms
sigGPSVel = 0.05^2; % m/s rms

G = 9.80665;

%0.1 for quats, 10m position, 10m/s velocity
Q = [0.0075 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0075 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0075 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0075 0.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.0 0.0 0.0;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 5.0 0.0;
     0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 10.0];

R_baro = [sigBaro];
R_gps = [sigGPS 0 0;
         0 sigGPS 0;
         0 0 sigGPS];
R_gpsVel = [sigGPSVel 0 0;
            0 sigGPSVel 0;
            0 0 sigGPSVel];
R_mag = [sigMag 0 0;
         0 sigMag 0;
         0 0 sigMag];

alt0 = atmospalt(dat.Prs(1)); %initial pressure altitude

% ====================

%this uses the last measured (landed) gravity vector along with the raw
%quaternion output from the EKF to determine approximately the true launch
%angle of the rocket.

% gf = [dat.aX(1) dat.aY(1) dat.aZ(1)]; %gravity vector
% 
% GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;     %%GG, FFi, UU create rotm
%               norm(cross(A,B)) dot(A,B)  0;
%               0              0           1];
% FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
% UU = @(Fi,G) Fi*G/(Fi);
% 
% Rfb = (quat2rotm([1 0 0 0]));
% v = Rfb*gf';
% Rft = UU(FFi(v,[0; 0; 1]),GG(v,[0; 0; 1]));
% vec0 = Rft*[1;0;0];
% 
% q0 = quatvec2([1;0;0],vec0);
q0 = [ORDat.qW(1) ORDat.qX(1) ORDat.qY(1) ORDat.qZ(1)];
mag0 = [dat.mX(1) dat.mY(1) dat.mZ(1)]; %initial magnetic vector
% mag0 = mag0/norm(mag0);

mag0 = quatrotate(quatconj(quaternion(q0)),mag0);

% =======================

Xnn = [q0(1),q0(2),q0(3),q0(4),0,0,0,0,0,0];

Pnn = [0.2 0 0 0 0 0 0 0 0 0;
       0 0.2 0 0 0 0 0 0 0 0;
       0 0 0.2 0 0 0 0 0 0 0;
       0 0 0 0.2 0 0 0 0 0 0;
       0 0 0 0 5 0 0 0 0 0;
       0 0 0 0 0 5 0 0 0 0;
       0 0 0 0 0 0 5 0 0 0;
       0 0 0 0 0 0 0 5 0 0;
       0 0 0 0 0 0 0 0 5 0;
       0 0 0 0 0 0 0 0 0 5];

statesX = zeros(10,height(dat));

covariance = zeros(10,10,height(dat));

statesX(:,1) = Xnn;
accelVecs = zeros(3,height(dat));

for i = 1:height(dat)-1
    % Prediction Step
    dT = (dat.Time(i+1)-dat.Time(i)); %timestep
    u = [dat.aX(i) dat.aY(i) dat.aZ(i) dat.gX(i) dat.gY(i) dat.gZ(i)]; %control vector
    Xn1n = stateTransitionFcn(Xnn,u,dT); %Xn+1,n
    J = @(x) stateTransitionFcn(x,u,dT); %function definition for taking jacobian numerically
    dF = numericalJacobian(J,Xnn); %jacobian dFdX

    Pn1n = dF*Pnn*dF' + Q;

    % Correction Step
    % barometer correction
    if(~isnan(dat.Prs(i)))
        zn = atmospalt(dat.Prs(i))-alt0; %atmospheric pressure altitude
        h = Xn1n(7);
        dh = [0 0 0 0 0 0 1 0 0 0];

        Kn = Pn1n * dh' / (dh*Pn1n*dh'+R_baro);
        PnnB = (eye(10)-Kn*dh)*Pn1n*((eye(10)-Kn*dh)')+Kn*R_baro*Kn';
        XnnB = Xn1n + Kn*(zn-h);
        q = XnnB(1:4);
        q = q/norm(q);
        XnnB(1:4) = q;
    else
        PnnB = Pn1n;
        XnnB = Xn1n;
    end
    % magnetometer correction
    zn = [dat.mX(i) dat.mY(i) dat.mZ(i)];
    % % % zn = zn/norm(zn);
    J = @(x) magMeasurementFcn(x,mag0);
    h = J(XnnB);
    dh = numericalJacobian(J,XnnB);

    Kn = PnnB * dh' / (dh*PnnB*dh'+R_mag);
    Pnn = (eye(10)-Kn*dh)*PnnB*((eye(10)-Kn*dh)')+Kn*R_mag*Kn';
    Xnn = XnnB + Kn*((zn-h)');
    q = Xnn(1:4);
    q = q/norm(q);
    Xnn(1:4) = q;

    % Xnn = XnnB;
    % Pnn = PnnB;
    statesX(:,i+1) = Xnn;
    covariance(:,:,i+1) = Pnn;
    

    accelVecs(:,i+1) = quatrotate([Xnn(1) -Xnn(2:4)'],[dat.aX(i),dat.aY(i),dat.aZ(i)]);
    % accelVecs(3,i+1) = accelVecs(3,i+1) - 1;
end

posList = statesX(5:7,:);
quatList = statesX(1:4,:);

%% Figure Plotting
figure(1);
subplot(3,1,1);
plot(dat.Time,statesX(5,:));
xlabel("Time (s)");
ylabel("Crossrange X (m)");
title("Crossrange X vs Time");
hold on;
plot(ORDat.t,ORDat.relPosX);
legend("Simulated","Truth")
grid on;
subplot(3,1,2);
plot(dat.Time,statesX(6,:))
xlabel("Time (s)");
ylabel("Crossrange Y (m)");
title("Crossrange Y vs Time");
hold on;
plot(ORDat.t,ORDat.relPosY);
legend("Simulated","Truth")
grid on;
subplot(3,1,3);
plot(dat.Time,statesX(7,:))
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
plot(dat.Time,statesX(1,:));
xlabel("Time (s)");
title("Quat W vs Time");
hold on;
plot(ORDat.t,ORDat.qW);
legend("Simulated","Truth")
grid on;
qxp = subplot(4,1,2);
plot(dat.Time,statesX(2,:));
xlabel("Time (s)");
title("Quat X vs Time");
hold on;
plot(ORDat.t,ORDat.qX);
legend("Simulated","Truth")
grid on;
qyp = subplot(4,1,3);
plot(dat.Time,statesX(3,:));
xlabel("Time (s)");
title("Quat Y vs Time");
hold on;
plot(ORDat.t,ORDat.qY);
legend("Simulated","Truth")
grid on;
qzp = subplot(4,1,4);
plot(dat.Time,statesX(4,:));
xlabel("Time (s)");
title("Quat Z vs Time");
hold on;
plot(ORDat.t,ORDat.qZ);
legend("Simulated","Truth")
grid on;
linkaxes([qwp qxp qyp qzp],'y')
qwp.YLim = [-1 1];

normQ = zeros(1,length(quatList));
for i = 1:length(quatList)
    normQ(i) = norm(quatList(:,i));
end
%%
numState = 4;
truthDat = ORDat.qZ;
plot(dat.Time,sqrt(squeeze(covariance(numState,numState,:))),"Color",'b');
hold on;
plot(dat.Time,-sqrt(squeeze(covariance(numState,numState,:))),"Color",'b');
hold on;
plot(dat.Time,statesX(numState,:)'-truthDat);
%% 6DoF Animation
% v = VideoWriter('ekf_animation_justlaunch.mp4','MPEG-4');
% v.FrameRate = 30;
% open(v);
for i = 1:2:height(dat)
    q = quaternion(quatList(:,i)');
    pos = posList(:,i);
    poseplot(q,pos,'ScaleFactor',0.2);

    timeStr = sprintf('Time: %.2f s', dat.Time(i));
    title(timeStr, 'FontSize', 14); 

    drawnow;

    % frame = getframe(gcf);
    % writeVideo(v, frame);
end
% close(v);

%%
stepsize = 20;
endpos = size(ORDat,1);
quiver3(posList(1,1:stepsize:endpos), posList(2,1:stepsize:endpos), posList(3,1:stepsize:endpos), accelVecs(1,1:stepsize:endpos),accelVecs(2,1:stepsize:endpos),accelVecs(3,1:stepsize:endpos),2); axis equal;

%%
function q = quatvec2(v1,v2)
    v1_norm = v1 / norm(v1);
    v2_norm = v2 / norm(v2);
    
    axis = cross(v1_norm, v2_norm);
    angle = acos(dot(v1_norm, v2_norm));
    
    if norm(axis) == 0
        if dot(v1_norm, v2_norm) > 0
            q = quaternion(1, 0, 0, 0);
        else
            q = quaternion(0, 1, 0, 0);
        end
    else
        axis_norm = axis / norm(axis);
        q = [cos(angle/2), sin(angle/2)*axis_norm(1), sin(angle/2)*axis_norm(2), sin(angle/2)*axis_norm(3)];
    end
end

%%
run_animation(quatList, posList, dat) ;

function run_animation(quatList, posList, dat)   
    
% v = VideoWriter('flight_of_fuckin_bumblebee.mp4','MPEG-4');
% v.FrameRate = 60;
% open(v);

    q = quaternion(quatList(:,1)');
    patch = poseplot(q);
    
    patch.ScaleFactor = 15;
    xlabel("X")
    ylabel("Y")
    zlabel("Z")
    
    for i = 1:5:height(dat)
        q = quaternion(quatList(:,i)');
        pos = posList(:,i);
        set(patch, Orientation=q, Position=pos); hold on
        plot3(pos(1), pos(2), pos(3), '.b', 'MarkerSize', 2); hold on
    
        xlim([-200, 200]);
        ylim([-200, 200]);
        zlim([0, 460]);
        
    
        set(gca,'ZDir','normal')  
        title(sprintf("t = %0.2f", dat.Time(i)))
        drawnow
        % frame = getframe(gcf);
        % writeVideo(v, frame);
    end
% close(v);
end