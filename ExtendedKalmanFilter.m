clear; clc; close;

SetupEnv();
dat = readtable("rocket_sd_data_fixed.csv");

% States are wxyz quaternion, xyz position, xyz velocity
% Control inputs are xyz accelerometer, xyz gyroscope
% Measurements are xyz GPS, x barometer

tS = 0.01; %sample time (100hz)

sigAccel = (0.002*9.81)^2; % m/s^2 rms
sigGyro = 0.1^2; % dps rms
sigMag = 20^2; % degrees rms
sigBaro = 5^2; % meters rms
sigGPS = 0.5^2; % meters rms
sigGPSVel = 0.05^2; % m/s rms

G = 9.80665;

magCoords = [38.82914122590703 -77.80889365885524];

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
R_mag = [sigMag 0 0;
         0 sigMag 0;
         0 0 sigMag];

alt0 = atmospalt(dat.Prs(1)); %initial pressure altitude

% ====================

%this uses the last measured (landed) gravity vector along with the raw
%quaternion output from the EKF to determine approximately the true launch
%angle of the rocket.

gf = [dat.aX(1) dat.aY(1) dat.aZ(1)]; %gravity vector

GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;     %%GG, FFi, UU create rotm
              norm(cross(A,B)) dot(A,B)  0;
              0              0           1];
FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
UU = @(Fi,G) Fi*G/(Fi);

Rfb = (quat2rotm([1 0 0 0]));
v = Rfb*gf';
Rft = UU(FFi(v,[0; 0; 1]),GG(v,[0; 0; 1]));
vec0 = Rft*[1;0;0];

q0 = quatvec2([1;0;0],vec0);

mag0 = [dat.mZ(1) dat.mY(1) dat.mX(1)]; %initial magnetic vector
mag0 = mag0/norm(mag0);

mag0 = quatrotate(quaternion(q0),mag0);

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
    u = [dat.aX(i)*G dat.aY(i)*G dat.aZ(i)*G dat.gX(i)*pi/180 dat.gY(i)*pi/180 dat.gZ(i)*pi/180]; %control vector
    Xn1n = stateTransitionFcn(Xnn,u,dT); %Xn+1,n
    J = @(x) stateTransitionFcn(x,u,dT); %function definition for taking jacobian numerically
    dF = numericalJacobian(J,Xnn); %jacobian dFdX

    Pn1n = dF*Pnn*dF' + Q;

    %Correction Step
    %barometer correction
    % if(~isnan(dat.Prs(i)))
    %     zn = atmospalt(dat.Prs(i))-alt0; %atmospheric pressure altitude
    %     h = Xn1n(7);
    %     dh = [0 0 0 0 0 0 1 0 0 0];
    % 
    %     Kn = Pn1n * dh' / (dh*Pn1n*dh'+R_baro);
    %     PnnB = (eye(10)-Kn*dh)*Pn1n*((eye(10)-Kn*dh)')+Kn*R_baro*Kn';
    %     XnnB = Xn1n + Kn*(zn-h);
    %     q = XnnB(1:4);
    %     q = q/norm(q);
    %     XnnB(1:4) = q;
    % else
        PnnB = Pn1n;
        XnnB = Xn1n;
    % end
    %magnetometer correction
    zn = [dat.mZ(i) dat.mY(i) dat.mX(i)];
    zn = zn/norm(zn);
    J = @(x) magMeasurementFcn(x,mag0);
    h = J(XnnB);
    dh = numericalJacobian(J,XnnB);

    Kn = PnnB * dh' / (dh*PnnB*dh'+R_mag);
    Pnn = (eye(10)-Kn*dh)*PnnB*((eye(10)-Kn*dh)')+Kn*R_mag*Kn';
    Xnn = XnnB + Kn*((zn-h)');
    q = Xnn(1:4);
    q = q/norm(q);
    Xnn(1:4) = q;

    statesX(:,i+1) = Xnn;
    covariance(:,:,i+1) = Pnn;
    % Xnn = XnnB;
    % Pnn = PnnB;

    accelVecs(:,i+1) = quatrotate(Xnn(1:4)',[dat.aX(i),dat.aY(i),dat.aZ(i)]);
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
subplot(3,1,2);
plot(dat.Time,statesX(6,:))
xlabel("Time (s)");
ylabel("Crossrange Y (m)");
title("Crossrange Y vs Time");
subplot(3,1,3);
plot(dat.Time,statesX(7,:))
xlabel("Time (s)");
ylabel("Altitude (m)");
title("Altitude vs Time");

quatStates = statesX(1:4,:);
eulAngs = zeros(height(dat),3);
for i = 1:height(dat)
    eulAngs(i,:) = quat2eul(quatStates(:,i)',"ZYX")*180/pi;
end

figure(2);
axP = subplot(3,1,1);
plot(dat.Time,eulAngs(:,1));
xlabel("Time (s)");
ylabel("Pitch (deg)");
title("Pitch Angle vs Time");
axY = subplot(3,1,2);
plot(dat.Time,eulAngs(:,2));
xlabel("Time (s)");
ylabel("Yaw (deg)");
title("Yaw Angle vs Time");
axR = subplot(3,1,3);
plot(dat.Time,eulAngs(:,3));
hold on;
plot(dat.Time,dat.mY*100);
xlabel("Time (s)");
ylabel("Roll (deg)");
title("Roll Angle vs Time");
linkaxes([axP axY axR],'y')
axP.YLim = [-180 180];

% figure(3);
% plot(dat.Time,statesX(8,:));

normQ = zeros(1,length(quatList));
for i = 1:length(quatList)
    normQ(i) = norm(quatList(:,i));
end

%%

plot(sqrt(squeeze(covariance(1,1,:))));
hold on;
plot(-sqrt(squeeze(covariance(1,1,:))));

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
endpos = 2500;
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
    
    patch.ScaleFactor = 60;
    xlabel("X")
    ylabel("Y")
    zlabel("Z")
    
    for i = 1:5:height(dat)
        q = quaternion(quatList(:,i)');
        pos = posList(:,i);
        set(patch, Orientation=q, Position=pos); hold on
        plot3(pos(1), pos(2), pos(3), '.b', 'MarkerSize', 2); hold on
    
        xlim([-1000, 1000]);
        ylim([-1000, 1000]);
        zlim([0, 1000]);
    
        set(gca,'ZDir','normal')  
        title(sprintf("t = %0.2f", dat.Time(i)))
        drawnow
        % frame = getframe(gcf);
        % writeVideo(v, frame);
    end
% close(v);
end