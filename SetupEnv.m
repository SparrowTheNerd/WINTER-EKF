function SetupEnv()

simPaths = {
    genpath(fullfile('Kalman Functions'))
    genpath(fullfile('Sensor Measurements'))
    genpath(fullfile('OpenRocket Data'))
};

simPaths = strjoin(simPaths,';');
addpath(simPaths);

end