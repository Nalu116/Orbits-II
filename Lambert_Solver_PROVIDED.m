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
