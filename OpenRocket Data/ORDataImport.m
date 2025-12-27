%  Processes CSV data input from InertialExport,
%  an OpenRocket simulation extension

function [ORDat,inertialDat] = ORDataImport(dT)
% dT = 0.01;

ORraw = readtable("OpenRocket Data/simulation-003.csv");

TT = table2timetable(ORraw, 'RowTimes', seconds(ORraw.t));
% retime resamples time data and resolves irregular timing, like from OR
TT_100Hz = retime(TT, seconds(0:dT:max(ORraw.t)), 'linear');

ORDat = timetable2table(TT_100Hz, 'ConvertRowTimes', true);
ORDat = removevars(ORDat,"Time");
clear TT TT_100Hz ORraw

[~,~,P,~,~,~] = atmosisa(ORDat.relPosZ);

% create magnetometer measurements
mXYZ = zeros(size(ORDat,1),3);
for i=1:size(ORDat,1)
    mXYZ(i,:) = wrldmagm(ORDat.relPosZ(i),ORDat.lat(i),ORDat.lon(i),decyear(2025,12,25),'2025');
    % rotate magnetometer measurements into body frame
    mXYZrot = quatrotate([ORDat.qW ORDat.qX ORDat.qY ORDat.qZ],mXYZ)*1e-5;
end

inertialDatArr = [ORDat.t P ORDat.bodyAccX ORDat.bodyAccY ORDat.bodyAccZ ORDat.bodyGyrX ORDat.bodyGyrY ORDat.bodyGyrZ mXYZrot(:,1) mXYZrot(:,2) mXYZrot(:,3)];
inertialDat = array2table(inertialDatArr,"VariableNames",{'Time','Prs','aX','aY','aZ','gX','gY','gZ','mX','mY','mZ'});

accel = quatrotate(quatconj([ORDat.qW ORDat.qX ORDat.qY ORDat.qZ]),[inertialDat.aX inertialDat.aY inertialDat.aZ]);
accel(:,3) = accel(:,3) - 9.80665;
accelBody = quatrotate(([ORDat.qW ORDat.qX ORDat.qY ORDat.qZ]),accel);
inertialDat.aX = accelBody(:,1); inertialDat.aY = accelBody(:,2); inertialDat.aZ = accelBody(:,3);

clear accel accelBody

% add noise and bias to measurements

sigAccel = (0.005*9.81); % m/s^2 rms
sigGyro = 0.5; % dps rms
sigMag = 0.0004; % Gauss rms
sigBaro = 3; % meters rms
% sigGPS = 0.5^2; % meters rms
% sigGPSVel = 0.05^2; % m/s rms

inertialDat.Prs = inertialDat.Prs + sigBaro * randn(size(ORDat,1),1);
inertialDat.aX = inertialDat.aX   + sigAccel * randn(size(ORDat,1),1);
inertialDat.aY = inertialDat.aY   + sigAccel * randn(size(ORDat,1),1);
inertialDat.aZ = inertialDat.aZ   + sigAccel * randn(size(ORDat,1),1);
inertialDat.gX = inertialDat.gX   + sigGyro * randn(size(ORDat,1),1)+0.1;
inertialDat.gY = inertialDat.gY   + sigGyro * randn(size(ORDat,1),1);
inertialDat.gZ = inertialDat.gZ   + sigGyro * randn(size(ORDat,1),1);
inertialDat.mX = inertialDat.mX   + sigMag * randn(size(ORDat,1),1);
inertialDat.mY = inertialDat.mY   + sigMag * randn(size(ORDat,1),1);
inertialDat.mZ = inertialDat.mZ   + sigMag * randn(size(ORDat,1),1);

inertialDat.Prs(2:2:end) = nan;
end
