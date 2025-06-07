function SetupEnv()

simPaths = {
    genpath(fullfile('Kalman Functions'))
    genpath(fullfile('Sensor Measurements'))
};

simPaths = strjoin(simPaths,';');
addpath(simPaths);

end