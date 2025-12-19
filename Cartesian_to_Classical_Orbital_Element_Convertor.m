%{
Caden N. Dennis
Cartesian to Classical Elements Convertor
%} 

format long g

'Input System Parameters';

%Spacecraft orbiting earth with cartesian orbital elements at some epoch

mu = 398600;             % Gravitational Parameter of Earth (km^3/s^2)
x = 1.8655e+05;          % x comp of position vector
y = 1.1690e+05;          % y comp of position vector
z = 2.6114e+05;          % z comp of position vector
xdot = 0.0254;           % x comp of velocity vector
ydot = -0.0987;          % y comp of velocity vector
zdot = -0.5750;          % z comp of velocity vector

r = sqrt(((x^2)+(y^2)+(z^2))) %radial distance (cross referenced --> good)

v = sqrt((xdot^2)+(ydot^2)+(zdot^2)) %velocity scalar (also cr --> good)

c1 = ((v^2)/2) - (mu/r);
a = -mu/(2*c1)

rbar = [x;y;z];          %position vector matrix
vbar = [xdot;ydot;zdot]; %velocity vector matrix
hbar = cross(rbar,vbar);

h = rssq(hbar)                % angular momentum scalar

e = sqrt((-(h^2)/(mu*a))+1)    %eccintricity 

x1 = dot(rbar,vbar);         %const term for dot product in ebar claculation

ebar = (1/mu)*((((v^2)-(mu/r))*rbar)-(x1)*vbar); %eccintricity vector
 
echeck = rssq(ebar)       % checking that e = echeck 

zbar = [0;0;1];           % define the zhat vector

x2 = dot(hbar,zbar);       % constant term for dot product in i calculation

i = acos(x2/h) *(180/pi)  % inclination of orbit (degrees)


x3 = cross(zbar,hbar);     %constant term for dot cross product in RAAN calc

c1 = rssq(x3);            % calculation term

abar = x3/c1;              % Direction of Ascending node

ynegcheck = abar(2,1);       %Selecting y component of ahat

xbar = [1;0;0]; %define the xhat matrix

x4 = dot(xbar,abar);     %Calculation term for finding RAAN


% Quad check for RAAN
if ynegcheck <0
    RAAN = 360 - (acos(x4) *(180/pi))
else 
    RAAN = acos(x4) *(180/pi)
end




znegcheck = ebar(3,1);       %Selecting z component of ehat

x5 = dot(abar,ebar)/(e);     %Calculation term for finding omega 

%Quad check for omega
if znegcheck <0
    omega = 360 - (acos(x5) *(180/pi))
else 
    omega = acos(x5) *(180/pi)
end

nunegcheck = dot(rbar,vbar);

x6 = dot(ebar,rbar)/(e*r);          %Calculation term for nu


%Quad check for true anomaly
if nunegcheck <0
    nu = 360 - (acos(x6)*(180/pi))
else 
    nu = acos(x6)*(180/pi)
end

%nu = nu*(pi/180)                   %Degree to radian for E calculation

x7 = (((1+e)/(1-e))^(1/2));
x8 = atan(nu/2);

E = (2*atan(x8/x7))+pi         %Eccintric anomoly (rad) used to get M
M = E-(e*sin(E));

M = M*(180/pi)                 %Mean Anomoly (degrees)



