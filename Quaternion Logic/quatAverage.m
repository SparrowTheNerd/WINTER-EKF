% old quat q1, new quat q2

function qAvg = quatAverage(Q1,Q2)
% make sure quat vector orientation is proper
q1 = [Q1(1) Q1(2) Q1(3) Q1(4)]; q2 = [Q2(1) Q2(2) Q2(3) Q2(4)];

r = quatmultiply(quatconj(q1),q2);
uMag = abs(2*acos(r(1)));

u = r(2:4) * (uMag/(sin(uMag/2)));

rn = [cos(uMag/4) u/uMag*sin(uMag/4)];

qAvg = quatmultiply(q1,rn);
qAvg = qAvg/norm(qAvg);

end