% rotates vector from body to inertial frame

function vecprod = quatrotBI(Q,v)

q = [Q(1) Q(2) Q(3) Q(4)]; % ensure it is horizontal
qconj = [q(1) -q(2) -q(3) -q(4)];

temp = quatmultiply(q,[0 v(1) v(2) v(3)]);

quatprod = quatmultiply(temp,qconj);

vecprod = quatprod(2:4);

end