clear;
clc;
load('SensorCalDataLong.mat');
varnames = sensorDat.Properties.VariableNames;

maxNumM = 100;

% dat = [sensorDat.ICMaX,sensorDat.ICMaY,sensorDat.ICMaZ];
% dat = sensorDat.Baro(~isnan(sensorDat.Baro));

for j=1:1:size(sensorDat,2)
    
    if(j == size(sensorDat,2))
        Fs = 200/3;
        omega = sensorDat.Baro(~isnan(sensorDat.Baro));
    else 
        Fs = 200;
        omega = table2array(sensorDat(:,j));
    end
    t0 = 1/Fs;
    
    L = size(omega,1);
    maxM = 2.^floor(log2(L/2));
    m = logspace(log10(1), log10(maxM), maxNumM).';
    m = floor(m); % m must be an integer.
    m = unique(m); % Remove duplicates.

    [avar, tau] = allanvar(omega, m, Fs);
    adev = sqrt(avar);

    % Find the index where the slope of the log-scaled Allan deviation is equal
    % to the slope specified.
    slope = -0.5;
    logtau = log10(tau);
    logadev = log10(adev);
    dlogadev = diff(logadev) ./ diff(logtau);
    [~, i] = min(abs(dlogadev - slope));

    % Find the y-intercept of the line.
    b = logadev(i) - slope*logtau(i);

    % Determine the angle random walk coefficient from the line.
    logN = slope*log(1) + b;
    N = 10^logN;
    
    slope = 0.5;
    [~, i] = min(abs(dlogadev - slope));
    
    % Find the y-intercept of the line.
    b = logadev(i) - slope*logtau(i);
    % Determine the rate random walk coefficient from the line.
    logK = slope*log10(3) + b;
    K = 10^logK;

    slope = 0;
    [~, i] = min(abs(dlogadev - slope));
    
    % Find the y-intercept of the line.
    b = logadev(i) - slope*logtau(i);
    % Determine the bias instability coefficient from the line.
    scfB = sqrt(2*log(2)/pi);
    logB = b - log10(scfB);
    B = 10^logB;

    fprintf("%s:\t N=%.6f\tK=%.6f\tB=%.6f\n",varnames{j},N,K,B);
    tauN = 1;
    lineN = N ./ sqrt(tau);
    tauK = 3;
    lineK = K .* sqrt(tau/3);
    tauB = tau(i);
    lineB = B * scfB * ones(size(tau));
    tauParams = [tauN, tauK, tauB];
    params = [N, K, scfB*B];
    figure
    loglog(tau, adev, tau, [lineN, lineK, lineB], '--', ...
        tauParams, params, 'o')
    title(sprintf('Allan Deviation of %s with Noise Parameters',varnames{j}))
    xlabel('\tau')
    ylabel('\sigma(\tau)')
    legend('$\sigma$', '$\sigma_N$', ...
        '$\sigma_K$', '$\sigma_B$', 'Interpreter', 'latex')
    text(tauParams, params, {'N', 'K', '0.664B'})
    grid on
    axis equal
end