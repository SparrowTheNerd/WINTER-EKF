clear; clc;
device = serialport("COM6",115200);
device.configureTerminator("CR/LF");
i = 0;
% r = rateControl(200);
while(true)
    i = i+1;
    s = readline(device);
    d(i,:) = textscan(s,'%f%f%f', 'Delimiter', ',', 'CollectOutput',1);
    % waitfor(r);
end
%%
dat = cell2mat(d);
save("gyro.mat","dat","-append");
r.statistics