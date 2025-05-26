function J = numericalJacobian(f, x)
    epsilon = 1e-6;
    n = length(x);
    fx = f(x);
    m = length(fx);
    J = zeros(m, n);
    
    for i = 1:n
        dx = zeros(n,1);
        dx(i) = epsilon;
        J(:,i) = (f(x + dx) - f(x - dx)) / (2 * epsilon);
    end
end
