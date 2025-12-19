function[X] = ephem(JDi,BODY)
 %{
    Programmer: Dylan Clay
    Date:       3/3/2022
    File:       ephem.m V2.1
    Purpose:    This function pulls and propagates spice data from a .mat file

    Inputs:
    JDi = user specified julian date (Version 1 valid for JD 2459595 -- > 2463245, inclusive)
    (In calendar dates (1/17/2022 -->  1/15/2032)
    BODY = planet of interest, ex: "Earth"

    Outputs:
    X:    Classical Elements
    X(1) - a (km)
    X(2) - e
    X(3) - i (rads)
    X(4) - RAAN (rads)
    X(5) - arg. periapsis (w) (rads)
    X(6) - M (rads)
      
      BODY numbers
      1 - Sun
      2 - Mercury
      3 - Venus
      4 - Earth
      5 - Moon
      6 - Mars
      7 - Phobos
      8 - Deimos
      9 - Jupiter
      10 - Saturn
      11 - Neptune
      12 - Uranus
       
      Note: for the Moon, Phobos, and Deimos the returned vector will be in the parent bodies' inertial frame
      Note: if you'd like to increase speed, reduce the 'tol' variable to 1*10^(-12)

%}
states = load('ephimeredes.mat');

%Internal Variables
AUtoKM = 1.49597870691*10^8;
JDS = 2459594;
JDM = 2463244;
body = 0;
cont = 0;
tol = 1*10^(-15);

if (JDS < JDi) && (JDi <=JDM)
  cont = 1;
  JDe = floor(JDi - JDS); %Whole number component
  JDd = JDi - floor(JDi); %Decimal component
else
  disp('Julian Date not in accepted range')
end

if strcmp(BODY,"Sun")
  body = 1;
  au = states.Positions(JDe,:,body);
elseif strcmp(BODY,"Mercury")
    body = 2;
    au = states.Positions(JDe,:,body);
elseif strcmp(BODY,"Venus")
    body = 3;
    au = states.Positions(JDe,:,body);
elseif strcmp(BODY,"Earth")
    body = 4;
        au = states.Positions(JDe,:,body);
elseif strcmp(BODY,"Moon")
    body = 5;
    au = states.Positions(JDe,:,body) - states.Positions(JDe,:,4);
elseif strcmp(BODY,"Mars")
    body = 6;
    au = states.Positions(JDe,:,body) ;
elseif strcmp(BODY,"Phobos")
    body = 7;
    au = states.Positions(JDe,:,body) - states.Positions(JDe,:,6);
elseif strcmp(BODY,"Deimos")
    body = 8;
    au = states.Positions(JDe,:,body) - states.Positions(JDe,:,6);
elseif strcmp(BODY,"Jupiter")
    body = 9;
    au = states.Positions(JDe,:,body)
elseif strcmp(BODY,"Saturn")
    body = 10;
    au = states.Positions(JDe,:,body)
elseif strcmp(BODY,"Neptune")
    body = 11;
    au = states.Positions(JDe,:,body)
elseif strcmp(BODY,"Uranus")
    body = 12;
    au = states.Positions(JDe,:,body)
else
  disp('Not a valid body')
end

if cont == 1
  %definition
  dT = JDd * 86400;
  mu = 1.32712440018*10^11;
  if (body == 5)
    mu = 3.98604418*10^5;
  else if (body == 7) || (body == 8)
    mu = 4.282837*10^4;
  end
  end
  km = au*AUtoKM;
  km(4) = km(4)/86400;
  km(5) = km(5)/86400;
  km(6) = km(6)/86400;
  xhat = [1;0;0;];
  yhat = [0;1;0;];
  zhat = [0;0;1;];
  if strcmp(BODY,"Moon")
      r = [km(1);km(2);km(3);];
      V = [km(4);km(5);km(6);];
      toeq = [1,0,0;0,cosd(23.4144),-sind(23.4144);0,sind(23.4144),cosd(23.4144);];
      r = toeq*r;
      V = toeq*V;
      km = [r;V;];
  end
  if strcmp(BODY,"Phobos") || strcmp(BODY,"Deimos")
      r = [km(1);km(2);km(3);];
      V = [km(4);km(5);km(6);];
      toeq = [1,0,0;0,cosd(25.2),-sind(25.2);0,sind(25.2),cosd(25.2);];
      r = toeq*r;
      V = toeq*V;
      km = [r;V;];
  end
    r = sqrt(km(1)^2 + km(2)^2 + km(3)^2);
    V = sqrt(km(4)^2 + km(5)^2 + km(6)^2);  
    a = 1/((2/r)-(V^2)/mu);
    n = sqrt(mu/a^3);
    rbar = [km(1),km(2),km(3)];
    Vbar = [km(4),km(5),km(6)];
    rbaru = rbar/r;
    Vbaru = Vbar/V;
    h = cross(rbar,Vbar);
    ebar = cross(Vbar,h)/mu - rbar/r;
    e = norm(ebar);
    i = acosd((dot(h,zhat))/(norm(h)));
    abar = (cross(zhat',h));
    ahat = abar/norm(abar);
    if(ahat(2)>=0)
      o = acosd(dot(xhat,ahat)); 
    else
      o = 360 - acosd(dot(xhat,ahat));
    end
    if(ebar(3)>=0)
      w = acosd((dot(ahat,ebar))/(e*norm(ahat)));
    else
       w = 360 - acosd((dot(ahat,ebar))/(e*norm(ahat)));
    end
    if(dot(rbar,Vbar)>=0)
      E = acos((a-r)/(a*e));
    else
      E = 2*pi - acos((a-r)/(a*e));
    end 
    M = E - e*sin(E);
  
 %propagation
  if (0 < JDd)
    Mf = M + dT*n;
    Eg = Mf + e*sin(Mf) + (e^2/2)*sin(2*Mf);
    
    %Newton's Method
    f=@(Q) Q - e*sin(Q);
   df=@(Q) 1 - e*cos(Q);
   x1 = Mf;                           
   x2 = Eg; 
   dif = x2 - x1; 
   k = 1;
   kmax = 15;
   while( abs(dif) > tol ) 
     x2 = x1 - (f(x1)-Mf)/df(x1);
     dif = x2-x1;
     x1 = x2;
       if k >= kmax
        conv = false;
      break
    end
     k = k + 1;
  end
   Ef = x2;
    dE = Ef - E;
    rf = a*(1-e*cos(Ef));
   
   %Building final vectors
   f = 1 - (a/r)*(1-cos(dE));
   g = dT - sqrt(a^3/mu)*(dE-sin(dE));
   fdot = (-sqrt(mu*a)*sin(dE))/(rf*r);
   gdot = 1 - (a/rf)*(1-cos(dE));
   
   rfinal = f*rbar + g*Vbar;
   vfinal = fdot*rbar + gdot*Vbar;
  
else
    rfinal = [km(1),km(2),km(3)];
    vfinal = [km(4),km(5),km(6)];
end
end
    r = norm(rfinal);
    V =norm(vfinal);
    rbar = rfinal;
    Vbar = vfinal;
    if(dot(rbar,Vbar)>=0)
      E = acos((a-r)/(a*e));
    else
      E = 2*pi - acos((a-r)/(a*e));
    end  
    M = E - e*sin(E);
    i = i*pi/180;
    w = w*pi/180;
    o = o*pi/180;
X = [a;e;i;o;w;M;];
end