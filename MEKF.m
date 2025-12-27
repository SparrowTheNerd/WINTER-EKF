clear; clc; close;

SetupEnv();

% 18 states are 3 ea. of quat, vel, pos, gyro bias, accel bias, mag bias errors
% quaternion here represents a body to inertial rotation

tS = 0.01; %sample time (100hz)
[ORDat,dat] = ORDataImport(tS);
