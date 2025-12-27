function SetupEnv()

simPaths = {
    genpath(fullfile('Kalman Functions'))
    genpath(fullfile('Sensor Measurements'))
    genpath(fullfile('OpenRocket Data'))
    genpath(fullfile('Helper Functions'))
    genpath(fullfile('Quaternion Logic'))
};

simPaths = strjoin(simPaths,';');
addpath(simPaths);

end