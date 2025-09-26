clear; clc; close;

SetupEnv();
dat = readtable("rocket_sd_data_fixed.csv");

% States are wxyz quaternion, xyz position, xyz velocity
% Control inputs are xyz accelerometer, xyz gyroscope
% Measurements are xyz GPS, x barometer

tS = 0.01; %sample time (100hz)
G = 9.80665;

alt0 = atmospalt(dat.Prs(1)); %initial pressure altitude

% ====================

%this uses the first measured gravity vector along with the raw
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

% =======================

Xnn = [q0(1),q0(2),q0(3),q0(4),0,0,0,0,0,0];

statesX = zeros(10,height(dat));

statesX(:,1) = Xnn;
for i = 1:height(dat)-1
    % Prediction Step
    dT = (dat.Time(i+1)-dat.Time(i)); %timestep
    u = [dat.aX(i)*G dat.aY(i)*G dat.aZ(i)*G dat.gX(i)*pi/180 dat.gY(i)*pi/180 dat.gZ(i)*pi/180]; %control vector
    Xn1n = stateTransitionFcn(Xnn,u,dT); %Xn+1,n
    J = @(x) stateTransitionFcn(x,u,dT); %function definition for taking jacobian numerically
    dF = numericalJacobian(J,Xnn); %jacobian dFdX
    statesX(:,i+1) = Xn1n;
    Xnn = Xn1n;
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
