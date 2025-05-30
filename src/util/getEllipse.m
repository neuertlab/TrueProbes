<<<<<<< HEAD
<<<<<<< HEAD
function [x,y] = getEllipse(r1,r2,C,phi)
beta = linspace(0.1,2*pi,100);
x = r1*cos(phi)*cos(beta) - r2*sin(phi)*sin(beta);
y = r1*sin(phi)*cos(beta) + r2*cos(phi)*sin(beta);
x = x + C(1,1);
y = y + C(1,2);
end
=======
function [x,y] = getEllipse(r1,r2,C,phi)
beta = linspace(0.1,2*pi,100);
x = r1*cos(phi)*cos(beta) - r2*sin(phi)*sin(beta);
y = r1*sin(phi)*cos(beta) + r2*cos(phi)*sin(beta);
x = x + C(1,1);
y = y + C(1,2);
end
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
=======
function [x,y] = getEllipse(r1,r2,C,phi)
beta = linspace(0.1,2*pi,100);
x = r1*cos(phi)*cos(beta) - r2*sin(phi)*sin(beta);
y = r1*sin(phi)*cos(beta) + r2*cos(phi)*sin(beta);
x = x + C(1,1);
y = y + C(1,2);
end
>>>>>>> 08410c48414cbfd1141b5d6a99035e1f365fbe06
