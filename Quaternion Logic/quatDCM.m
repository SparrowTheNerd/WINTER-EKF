% (DCM*V is equivalent action to quatrotBI(q,V) )

function DCM = quatDCM(q)

Q = [q(1); q(2); q(3); q(4)]; % make sure it is vertical

DCM = 2*Q(2:4)*Q(2:4)'  +  eye(3) * ( Q(1)^2 - Q(2:4)'*Q(2:4) ) - (2*Q(1) * skew(-Q(2:4)));

end