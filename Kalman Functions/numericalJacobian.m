function J = numericalJacobian(f, x)
    epsilon = 1e-6;
    n = length(x);
    fx = f(x);
    m = length(fx);
    J = zeros(m, n);
    
    for i = 1:n
        dx = zeros(n,1);
        dx(i) = epsilon;
        x1 = x + dx;
        x2 = x - dx;
        x1(1:4) = x1(1:4) / norm(x1(1:4));
        x2(1:4) = x2(1:4) / norm(x2(1:4));
        J(:,i) = (f(x1) - f(x2)) / (2 * epsilon);
    end
end
