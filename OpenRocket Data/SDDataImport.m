function [ORDat, inertialDat] = SDDataImport(dT, csvFile)
%SDDATAIMPORT Import an sdread.py telemetry CSV for EKF replay testing.
%
%   [ORDat,inertialDat] = SDDATAIMPORT(dT,csvFile) resamples a logged
%   telemetry CSV at dT seconds. inertialDat contains the actual logged IMU,
%   magnetometer, barometric-altitude-derived pressure, and GPS measurements.
%
%   ORDat intentionally contains no simulated truth: all non-position fields
%   are zero. GPS supplies its latitude, longitude, and local N/E/U position.
%   relPosX is north (m), relPosY is east (m), and relPosZ is up (m), each
%   relative to the first valid GPS fix. This is a local convention chosen for
%   replay; it is not inferred from the EKF state vector.
%
%   csvFile is optional and defaults to "FLIGHTDATA17.csv".

arguments
    dT (1,1) double {mustBePositive}
    csvFile (1,1) string = "FLIGHTDATA02.csv"
end

raw = readtable(csvFile, "VariableNamingRule", "preserve");
required = ["time", "accel_x", "accel_y", "accel_z", "highG_x", "highG_y", "highG_z", ...
    "gyro_x", "gyro_y", "gyro_z", "mag_x", "mag_y", "mag_z", ...
    "lat", "lon", "gpsAlt", "gpsNVel", "gpsEVel", "baroAlt"];
missing = required(~ismember(required, string(raw.Properties.VariableNames)));
if ~isempty(missing)
    error("SDDataImport:MissingColumns", ...
        "CSV is not an sdread.py export. Missing: %s", strjoin(missing, ", "));
end

signalNames = ["accel_x", "accel_y", "accel_z", "highG_x", "highG_y", "highG_z", "gyro_x", "gyro_y", "gyro_z", ...
    "mag_x", "mag_y", "mag_z", "baroAlt", "lat", "lon", "gpsAlt", "gpsNVel", "gpsEVel"];
signals = raw(:, cellstr(signalNames));
TT = table2timetable(signals, "RowTimes", seconds(raw.time));
sampleTimes = seconds(raw.time(1):dT:raw.time(end));
TT = retime(TT, sampleTimes, "linear");
% Colon indexing normalizes both timetable row times and variables to columns
% across MATLAB releases.
time = seconds(TT.Time(:));

lat = TT.lat(:);
lon = TT.lon(:);
gpsAlt = TT.gpsAlt(:);
accelX = -TT.accel_z(:); accelY = TT.accel_y(:); accelZ = TT.accel_x(:);
gyroX = -TT.gyro_z(:); gyroY = TT.gyro_y(:); gyroZ = TT.gyro_x(:);
highGX = -TT.highG_z(:); highGY = -TT.highG_y(:); highGZ = -TT.highG_x(:); 
magX = TT.mag_x(:); magY = TT.mag_y(:); magZ = TT.mag_z(:);
baroAlt = TT.baroAlt(:);

enu = dllh2denu([lat(1) lon(1) gpsAlt(1)], [lat lon gpsAlt]);
gpsVelNorth = TT.gpsNVel; gpsVelEast = TT.gpsEVel;

% ORDat mirrors common OpenRocket truth names used by existing EKF harnesses.
% Only GPS-derived position is populated; every other truth field is zero.
n = height(TT);
z = zeros(n, 1);
ORDatArr = [time, z, z, z, z, enu(:,1), enu(:,2), enu(:,3), z, z, z, ...
    z, z, z, z, z, z, lat, lon, gpsAlt];
ORDatArr = ORDatArr(7732:end,:);
ORDat = array2table(ORDatArr, "VariableNames", {'t', 'qW', 'qX', 'qY', 'qZ', 'relPosX', 'relPosY', 'relPosZ', 'worldVelX', 'worldVelY', 'worldVelZ', 'bodyAccX', 'bodyAccY', 'bodyAccZ','bodyGyrX', 'bodyGyrY', 'bodyGyrZ', 'lat', 'lon', 'gpsAlt'});

% sdread.py logs barometric altitude, not raw pressure. Convert it to Pa to
% preserve the Prs interface of the original OpenRocket importer.
[~, ~, pressure] = atmosisa(baroAlt);
pressure = pressure(:);
inertialDatArr = [time, pressure, accelX, accelY, accelZ, highGX, highGY, highGZ, gpsVelEast, gpsVelNorth, deg2rad(gyroX), deg2rad(gyroY), deg2rad(gyroZ), magX, magY, magZ, lat, lon];
inertialDatArr = inertialDatArr(7732:end,:);
inertialDat = array2table(inertialDatArr, "VariableNames", {'Time', 'Prs', 'aX', 'aY', 'aZ', 'hGX', 'hGY', 'hGZ', 'vX', 'vY', ...
    'gX', 'gY', 'gZ', 'mX', 'mY', 'mZ', 'lat', 'lon'});

end