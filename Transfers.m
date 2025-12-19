%{
Caden N. Dennis
%}

clear all
close all
clc
format long g

'Input System Parameters';

JD1_4 = 2459883;
JD2_4 = 2460069;
JD1_Mat=[2459883:1:(2459883+1555)];
JD2_Mat=[2459883:1:(2459883+1555)];
JJ = 1;
n = 10;
tol = 10^-14;
kmax = 100;
mu = 132700000000; % Gravitational Parameter of the Sun (km^3/s^2)
mumoon = 398600; %Gravitational Parameter of Earth (km^3/s^2)

%Setup Jualian Date Matrices 
[JD1_Mat3,JD2_Mat3]=meshgrid(JD1_Mat,JD2_Mat);

for ii = 1:778*2    %~ 778 days in one synodic period (2.13 years)
   
[JD1_Matrix(ii,:)] = ephem(JD1_Mat(ii),"Earth");
[JD2_Matrix(ii,:)] = ephem(JD2_Mat(ii),"Mars");
[JD3_Matrix(ii,:)] = ephem(JD1_Mat(ii),"Moon");


[JD1_Matrix2(ii,:)] = Classical2Cart(JD1_Matrix(ii,1),JD1_Matrix(ii,2),JD1_Matrix(ii,3),JD1_Matrix(ii,4),JD1_Matrix(ii,5),JD1_Matrix(ii,6),mu);
[JD2_Matrix2(ii,:)] = Classical2Cart(JD2_Matrix(ii,1),JD2_Matrix(ii,2),JD2_Matrix(ii,3),JD2_Matrix(ii,4),JD2_Matrix(ii,5),JD2_Matrix(ii,6),mu);
[JD3_Matrix2(ii,:)] = Classical2Cart(JD3_Matrix(ii,1),JD3_Matrix(ii,2),JD3_Matrix(ii,3),JD3_Matrix(ii,4),JD3_Matrix(ii,5),JD3_Matrix(ii,6),mumoon);

end
%Direction Cosine Matrix
MoonDCM=[1 0 0;0 cosd(23.4144) sind(23.4144);0 -sind(23.4144) cosd(23.4144)];

for i=1:778*2
    JD3_Matrix2(i,1:3)=MoonDCM*transpose(JD3_Matrix2(i,1:3));
    JD3_Matrix2(i,4:6)=MoonDCM*transpose(JD3_Matrix2(i,4:6));
end
JD3_Matrix2=JD1_Matrix2+JD3_Matrix2;

for    ii = 1:778*2
JD1=JD1_Mat(ii);
%JD2=JD2_Mat(ii);

for jj = 1:778*2
%JD1=JD1_Mat(jj);
JD2=JD2_Mat(jj);    
    
TOF = ((JD2-JD1)*86400);
[TOF_Matrix(ii,jj)] = TOF;

% function F = CelsiustoFar(X)
 
[A,P,V1,V2,conv] = Lambert(JD3_Matrix2(ii,1:3),JD2_Matrix2(jj,1:3),TOF,mu,JJ,n,tol,kmax);

%patched conic method below : Vinf+ Moon = Vinf+earth -Vmoon
amoon = JD3_Matrix(ii,1);
emoon = JD3_Matrix(ii,2);
Mmoon = JD3_Matrix(ii,6);
Emoon = kepler_E(emoon, Mmoon);
rmoon = amoon*(1-emoon*cos(Emoon));
Vmoon = sqrt(2*((-mumoon/(2*amoon))+(mumoon/rmoon)));

%Vinfmoon = V1 - Vmoon;

%Vinfmoon = norm(transpose(V1)-JD3_Matrix2(ii,4:6));

%C3 = (norm(transpose(Vinfmoon)-JD1_Matrix2(ii,4:6)))^2;
C3 = (norm(transpose(V1)-JD3_Matrix2(ii,4:6)))^2;
V_inf_minus=norm(transpose(V2)-JD2_Matrix2(jj,4:6));
V100km2=sqrt(2*(((V_inf_minus^2)/2)+((4.2828*10^4)/(100+3390))));

%{ 
if  C3 < 80
    C3_Mat(jj,ii)= (norm(transpose(V1)-JD3_Matrix2(ii,4:6)))^2;
    
    
else C3 > 80 
    C3_Mat(jj,ii) =NaN;
  
end
%}

if  (C3 > 40) || (JD1_Mat3(jj,ii)>JD2_Mat3(jj,ii)) %% checking for absurdities 
    
    C3_Mat(jj,ii) = NaN;

%{    
else C3 < 80 
    C3_Mat(jj,ii)= (norm(transpose(V1)-JD3_Matrix2(ii,4:6)))^2;
%} 
else
    C3_Mat(jj,ii)= (norm(transpose(V1)-JD3_Matrix2(ii,4:6)))^2;
    
end

%sc orbits moon circularly at 100km alt:
vsc1 = sqrt(2*((4903/1837.5)-(4903/3675)));

%delta V = V new - V old
[DeltaVMat(jj,ii)] = sqrt(C3_Mat(jj,ii)) - vsc1;

%{
if (JD1_Mat3(jj,ii) >= JD2_Mat3(jj,ii))%+(186)
    C3_Mat(jj,ii) = NaN;
end
%} 

if V100km2 <= 20
    V100km(jj,ii)=V100km2;
end

if JD1_Mat3(jj,ii) > JD2_Mat3(jj,ii)
    V100km(jj,ii) =NaN;
end
    
if V100km2 > 20  
    V100km(jj,ii) =NaN;
    
end

%C3_Mat(jj,ii)= (norm(transpose(V1)-JD1_Matrix2(ii,4:6)))^2;
end

%C3  = V_Inf_Plus_Mag^2

ii %Just a counting term displayed to make sure code running

end

%Launch_prd=C3_Mat(:,1:30);
Min_C3=min(C3_Mat(:))
Launch_prd2=V100km(1:30,:);

%MaxEntry=max(Launch_prd2(:))
%MinEntry=min(Launch_prd2(1:30))


figure(1)
contour(JD1_Mat3,JD2_Mat3,C3_Mat,'ShowText','on')
xlabel(['Julian Departure Date'])
ylabel(['Julian Arrival Date'])
title('Moon to Mars C3 Over Two Synodic Periods')

figure(2)
contour(JD1_Mat3,JD2_Mat3,V100km,'ShowText','on')
xlabel(['Julian Departure Date'])
ylabel(['Julian Arrival Date'])
title('Mars 100km Atmospheric Entry Velocity (km/s) Over Two Synodic Periods')

figure(3)
contour(JD1_Mat3,JD2_Mat3,DeltaVMat,'ShowText','on')
xlabel(['Julian Departure Date'])
ylabel(['Julian Arrival Date'])
title('Moon to Transfer Delta V Over Two Synodic Periods')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1 = [129140000,70740000,61.6167]
V1 = [-14.383,26.604,0.0000265]

R2 = [112750000,-185710000,-105480]
V2 = [20.484,14.8256,0.0082]

%Cartesian to Classical gross non function for now cause I left my
%function file on a campus computer sad

cartesian(1) = 129140000;
cartesian(2) = 70740000;
cartesian(3) = 61.616;
cartesian(4) = -14.383;
cartesian(5) = 26.604;
cartesian(6) = 0.0000265;

r_new=sqrt(cartesian(1)^2+cartesian(2)^2+cartesian(3)^2);
V_new=sqrt(cartesian(4)^2+cartesian(5)^2+cartesian(6)^2);
a_new=((-mu)/(((V_new^2)/2)-(mu/r_new)))/2;
rbar_new=[cartesian(1) cartesian(2) cartesian(3)];
Vbar_new=[cartesian(4) cartesian(5) cartesian(6)];
hbar_new=cross(rbar_new,Vbar_new);
h_new=sqrt(hbar_new(1)^2+hbar_new(2)^2+hbar_new(3)^2);
e_new=sqrt(-((h_new^2/(mu*a_new))-1));
ebar_new=(1/mu).*((V_new^2-(mu/r_new)).*rbar_new-(dot(rbar_new,Vbar_new)).*Vbar_new);
ebar_new_check=sqrt(ebar_new(1)^2+ebar_new(2)^2+ebar_new(3)^2);
in_new=(acos((dot(hbar_new,[0;0;1]))/h_new))*(180/pi);
node_line_part=cross([0;0;1],hbar_new);
node_line_part2=sqrt(node_line_part(1)^2+node_line_part(2)^2+node_line_part(3)^2);
node_line=node_line_part/node_line_part2;
omega_new=360-((acos((dot(node_line,ebar_new))/e_new))*(180/pi));
if node_line(2) < 0
    Omega_new=360-((acos(dot([1;0;0],node_line)))*(180/pi));
    omega_new=360-((acos((dot(node_line,ebar_new))/e_new))*(180/pi));
else
    Omega_new=(acos(dot([1;0;0],node_line)))*(180/pi);
    omega_new=(acos((dot(node_line,ebar_new))/e_new))*(180/pi);
end
v_new=(acos((dot(ebar_new,rbar_new))/(e_new*r_new)));
if (dot(rbar_new,Vbar_new)) < 0
    v_new=(2*pi)-v_new;
end
E_new=(2*atan((tan(v_new/2))/(((1+e_new)/(1-e_new))^(1/2))));
if E_new < 0
    E_new=(2*pi)+(2*atan((tan(v_new/2))/(((1+e_new)/(1-e_new))^(1/2))));
else
    E_new=(2*atan((tan(v_new/2))/(((1+e_new)/(1-e_new))^(1/2))));
end
M_new=(E_new-e_new*sin(E_new))*(180/pi);
classical_new=[a_new,e_new,in_new,Omega_new,omega_new,M_new];
disp(classical_new)


%[sca,sce,sci,scRANN,scomega,scM] = classical_new;

sca = classical_new(1)
sce = classical_new(2)
sci = classical_new(3)
scRANN = classical_new(4)
scomega = classical_new(5)
scM = classical_new(6)

%[sc.a,sc.e,sc.i,sc.RANN,sc.omega,sc.M] = cart2classic(R1,V1,mu);
scp = sca*(1-sce^2);

%%JDtrasnfer orbit aka lambert stuff

nu1 = 2*pi-acos(scp/(norm(R1)*sce)-1/sce);
nu2 = acos(scp/(norm(R2)*sce)-1/sce);
nu = linspace(nu1,nu2,650);


JDVec = linspace(2459883,2459883+700,650);

for i = 1:length(JDVec)

    orbitEarthCL = ephem(JDVec(i),'Earth');
    orbitMarsCL = ephem(JDVec(i),'Mars');

    [orbitEarthCA(:,i)] = Classical2Cart(orbitEarthCL(1),orbitEarthCL(2), ...
        orbitEarthCL(3),orbitEarthCL(4),orbitEarthCL(5),orbitEarthCL(6),mu);


    [orbitMarsCA(:,i)] = Classical2Cart(orbitMarsCL(1),orbitMarsCL(2), ...
        orbitMarsCL(3),orbitMarsCL(4),orbitMarsCL(5),orbitMarsCL(6),mu);

    scr = scp/(1+sce*cos(nu(i)));
    f = 1-(scr/scp)*(1-cos(nu(i)-nu1));
    g = scr*norm(R1)/sqrt(mu*scp)*sin(nu(i)-nu1);
    fdot = sqrt(mu/scp)*tan((nu(i)-nu1)/2)*((1-cos(nu(i)-nu1)/scp)-1/scr-1/norm(R1));
    gdot = 1-norm(R1)/scp*(1-cos(nu(i)-nu1));

    rbar(:,i) = f*R1+g*V1;
    vbar(:,i) = fdot*R1+gdot*V1;
    AU = 1.496e8; 

%     figure(3)
%     plot3(orbit.EarthCA(1,:),orbit.EarthCA(2,:),orbit.EarthCA(3,:),'color','#00FFFF','LineWidth',1)
%     hold on
%     plot3(orbit.MarsCA(1,:),orbit.MarsCA(2,:),orbit.MarsCA(3,:),'color','#FF0000','LineWidth',1)
%     plot3(rbar(1,:),rbar(2,:),rbar(3,:),'color','#000000','LineWidth',1)
%     plot3(0,0,0,'o','color','#EDB120','MarkerFaceColor','#EDB120')
%     xlim([-3*AU 3*AU])
%     ylim([-3*AU 3*AU])
%     zlim([-3*AU 3*AU])

end

AU = 1.496e8; 

figure(4)
plot3(orbitEarthCA(1,:),orbitEarthCA(2,:),orbitEarthCA(3,:),'color','#00FFFF','LineWidth',1)
hold on
plot3(orbitMarsCA(1,:),orbitMarsCA(2,:),orbitMarsCA(3,:),'color','#FF0000','LineWidth',1)
plot3(rbar(1,:),rbar(2,:),rbar(3,:),'color','#000000','LineWidth',1)
plot3(0,0,0,'o','color','#EDB120','MarkerFaceColor','#EDB120')
xlim([-3*AU 3*AU])
ylim([-3*AU 3*AU])
zlim([-3*AU 3*AU])
title('Moon to Mars Trajectory')

% Found this function in Appendix D "MATLAB Scripts" 

function E = kepler_E(e, M)
%                      
%{
This function uses Newton’s method to solve Kepler’s
equation E - e*sin(E) = M for the eccentric anomaly,
2/14/22 2:26 PM S:\Orbits II\Classical_to_Cartesian_Orb... 3 of 3
given the eccentricity and the mean anomaly.
E - eccentric anomaly (radians)
e - eccentricity, passed from the calling program
M - mean anomaly (radians), passed from the calling program
pi - 3.1415926...
User m-functions required: none
%}
% 
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

function [A,P,V1,V2,conv] = Lambert(R1,R2,TOF,mu,JJ,n,tol,kmax)
%{
    Programmer: Grant Hecht
    Date:       3/11/2019
    File:       Lambert.m
    Purpose:    This function solves Lambert's Problem using Battin's
                method.

    Inputs:
    R1:     Vector to Pt#1
    R2:     Vector to Pt#2
    TOF:    Time of Flight
    mu:     Gravitational parameter for central body
    JJ:     Integer that determines initial guess for x
                (set JJ = 1 for an ellipse)
                (set JJ = 0 for a parabola or hyperbola)
    n:      Number of continued fracton levels
                (set n = 0 for default of 10)
                (recomend setting n >= 100 for most accurate results)
    tol:    Tolerance to exit iterations
    kmax:   Maximum Iterations
    

    Outputs:
    V1:     Velocity Vector at Pt#1
    V2:     Velocity Vector at Pt#2
    conv:   Boolean to indicate convergence.

%}
% Defalt Value for n
n_default = 10;

% If does not converge, set as false
conv = true;

% Sets n to defalt value if n = 0 is passed
if n == 0
    n = n_default;
end

% Find Transfer Angle
ta = acos(dot(R1,R2)/(norm(R1)*norm(R2)));
if R1(1)*R2(2)-R1(2)*R2(1) < 0
    ta = 2*pi - ta;
end

% Find Chord
c = sqrt(norm(R1)^2 + norm(R2)^2 - 2*norm(R1)*norm(R2)*cos(ta));

% Find Semi-Perimeter
s = (norm(R1) + norm(R2) + c)/2;

% Find Lambda
lambda = sqrt(norm(R1)*norm(R2))*cos(ta/2)/s;

% Find w
w = atan((norm(R2)/norm(R1))^0.25)-(pi/4);

% Find l
if (0 < ta) && (ta < pi)
    l = (sin(ta/4)^2+tan(2*w)^2)/(sin(ta/4)^2+tan(2*w)^2+cos(ta/2));
elseif (pi <= ta) && (ta < 2*pi)
    l = (cos(ta/4)^2+tan(2*w)^2-cos(ta/2))/(cos(ta/4)^2+tan(2*w)^2);
else
    fprintf('Cannot Compute for Transfer Angle of 0 or 360 degrees.');
    return
end

% Find m
m = (8*mu*TOF^2)/(s^3*(1+lambda)^6);

% For Eliptical Transfer Orbit Use x = l for Initial Guess
% For Hyperbolic or Parabolic Transfer Orbit Use x = 0.
if JJ == 0
    x = 0;
else
    x = l;
end

%Define Velocity Vectors
V1 = zeros(3,1);
V2 = zeros(3,1);

%Define delta x and counter
DX = 100;
k  = 0;

%Iterate to solve
while abs(DX) > tol 
    
    % Breaks Itteration and sets conv = false if k = kmax 
    if k >= kmax
        conv = false;
        break
    end
    
    % Calculates Continued Fraction PHI for 'n' Levels
    eta = x/(sqrt(1+x)+1)^2;
    f   = 1;
    % Itterate to Calculate levels 4 -> n
    for j = n:-1:4
        ceta = j^2/((2*j)^2-1);
        f = 1 + ceta*eta/f;
    end
    % Finishes Calculation of PHI with levels 1 -> 3
    PHI = 8*(sqrt(1+x)+1)/(3+1/(5+eta+(9/7)*eta/f));

    h1 = (l+x)^2*(1+3*x+PHI)/((1+2*x+l)*(4*x+PHI*(3+x)));
    h2 = m*(x-l+PHI)/((1+2*x+l)*(4*x+PHI*(3+x)));
    B  = 27*h2/(4*(1+h1)^3);
    u  = B/(2*(sqrt(1+B)+1));

    % Calculates Continued Fraction K(u) for 'n' Levels
    f = 1;
    r = n/2-1;
    % Itterate to Calculate Levels 3 -> n
    for j = r:-1:1
        g2n  = 2*(3*j+1)*(6*j-1)/(9*(4*j-1)*(4*j+1));
        g2n1 = 2*(3*j+2)*(6*j+1)/(9*(4*j+1)*(4*j+3));
        f = 1 + g2n*u/(1 + g2n1*u/f);   
    end
    % Finishes Calculation of K(u) with levels 1 -> 2
    K = (1/3)/(1+(4/27)*u/f);
    
    % Calculates New Values for y and x
    yNew = ((1+h1)/3)*(2+sqrt(1+B)/(1+2*u*K^2));
    xNew = sqrt(((1-l)/2)^2+m/yNew^2)-(1+l)/2;
    
    % Compares xNew with x
    DX = abs(xNew - x);
    
    % Sets x and y to xNew and yNew
    x = xNew;
    y = yNew;
    
    k = k + 1;
end

%Computes Orbit Parameters and Initial and Final Velocity
A  = m*s*(1+lambda)^2/(8*x*y^2);
P0 = c^2*(1+x)^2/(16*A*x);
P  = 4*norm(R1)*norm(R2)*P0*sin(ta/2)^2/c^2;
F  = 1-(norm(R2)/P)*(1-cos(ta));
G  = norm(R1)*norm(R2)*sin(ta)/sqrt(mu*P);
FDOT = sqrt(mu/P)*tan(ta/2)*((1-cos(ta))/P-1/norm(R2)-1/norm(R1));
GDOT = 1-(norm(R1)/P)*(1-cos(ta));
for K=1:3
    V1(K)=(1/G)*(R2(K)-(F*R1(K)));
    V2(K)=FDOT*R1(K)+GDOT*V1(K);
end

end







