%{
Caden N. Dennis
%}

close all
clc
format long g

function [x,y,z,xdot,ydot,zdot] = Classical2Cart(a,e,i,RAAN,w,M);

'Input System Parameters';
%a,e, i in degrees, raan in degrees, w in degrees, M in RADIANS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END SYSTEM INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = sqrt((mu)/(a^3));       % Mean Motion
h = sqrt(mu*(a*(1-(e^2)))); % Specific Angular Momentum 


%{
Note incliniation is 100 degrees, retrograde orbit

First, need to get Eccentric anomaly from M numerically: M = E-esinE
Took this function directly from "MATLAB Scripts"
%} 



E = kepler_E(e, M)



nu = 2*atan(((((1+e)/(1-e))^(1/2))*tan(E/2)));

nu = nu/pi*180           %Comes out at negative when E is 4 rad, quad check
                         %Looks right, as -162degrees matches that, 
                         %Spacecraft on approach 
theta = w + nu;          %Theta for Direction Cosine Matrix               
                    

%{
Now to define the Direction Cosine Matrix relating the Inertial Vernal 
Equinox frame to polar rotating frame. Definine matrix as mcr for matrix
column row, so m13 would be column 1 row 3 of the matrix
%}

m11 = cosd(RAAN)*cosd(theta)-(sind(RAAN)*cosd(i)*sind(theta));
m12 = sind(RAAN)*cosd(theta)+(cosd(RAAN)*cosd(i)*sind(theta));
m13 = sind(i)*sind(theta);
m21 = -sind(theta)*cosd(RAAN)-(sind(RAAN)*cosd(i)*cosd(theta));
m22 = -sind(RAAN)*sind(theta)+(cosd(RAAN)*cosd(i)*cosd(theta));
m23 = sind(i)*cosd(theta);
m31 = sind(RAAN)*sind(i);
m32 = -cosd(RAAN)*sind(i);
m33 = cosd(i);

DCM = [m11,m21,m31;m12,m22,m32;m13,m23,m33]


r = (((a*(1-(e^2)))/(1+(e*cosd(nu))))) % radial distance (km)

rmatrix = [r;0;0]; %r vector rhat, thetahat, what;

rbar = DCM*rmatrix

rcheck = rssq(rbar) %Root sum square to check rbar correct (r=rbar then T)

gamma = atand((e*sind(nu))/(1+(e*cosd(nu)))) %Flight Path Angle

% Quad check, negative Flight path angle, SC on approach, good.

v = sqrt((2*(((mu/r)-((mu)/(2*a)))))) %Velocity at point on orbit

vmatrix = [v*sind(gamma);v*cosd(gamma);0]; %velocity vector matrix;

vbar = DCM*vmatrix

vcheck = rssq(vbar) %Root sum square checking vbar (vbar = v, true)

hcheck = cross(rmatrix,vmatrix) %Double checking v and r, to get same h



% Found this function in Appendix D "MATLAB Scripts" 


function E = kepler_E(e, M)
% 
%{
This function uses Newton’s method to solve Kepler’s
equation E - e*sin(E) = M for the eccentric anomaly,
given the eccentricity and the mean anomaly.
E - eccentric anomaly (radians)
e - eccentricity, passed from the calling program
M - mean anomaly (radians), passed from the calling program
pi - 3.1415926...
User m-functions required: none
%}
% ––––––––––––––––––––––––––––––––––––––––––––––
%...Set an error tolerance:
error = 1.e-8;
%...Select a starting value for E:
if M < pi
E = M + e/2;
else
E = M - e/2;
end
%...Iterate on Equation 3.17 until E is determined to within
%...the error tolerance:
ratio = 1;
while abs(ratio) > error
ratio = (M-(E - e*sin(E)))/(-1 + e*cos(E));
E = E - ratio;
end
end %kepler_E

end


