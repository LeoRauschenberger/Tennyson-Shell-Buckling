% Analytical Solution by Tennyson for Cylinder
%
% by Leo Rauschenberger
% Version 3.0
%
clear; format;
format short g;

%%%%%%%%%% Definitions %%%%%%%%%%%
%Dimensions Cylinder
L = 510;
R = 250;
C = 2*pi*R;
t_ply  = 0.125; 

% Number of half-waves
m = 9
m_0 = 2*m
n = 13
% Initial imperfection to thickness of shell mu = w0/t 
% is 0 for perfect
mu = 0

% Ply engineering properties
% Unit: Mpa
E1   = 123551;
E2   = 8708;
nu12 = .32; 
G12  = 5695; 

% Laminate definition
Nplies = 10;
% ply angles in degrees, from top to bottom(!)
thetadt = [-53 53 -38 38 -22 22 90 90 -30 30 ];     % for Z17
%thetadt = [-38 38 90 90 90 90 -68 68 -38 38];       % for FL5
% ply angles in degrees, from bottom
thetadb = fliplr(thetadt);         
t      = Nplies * t_ply ;

for i = 1:Nplies;
  zbar(i) = - (t + t_ply)/2 + i*t_ply;
end;

% -----------Q matrix--------------------
denom = 1 - nu12^2 * E2 / E1 ;
Q11 = E1 / denom;
Q12 = nu12 * E2 / denom;
Q22 = E2 / denom;
Q66 = G12;

Q = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66]; 
S = inv(Q); 

% ------------ABD matrices---------------
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

for i = 1:Nplies;
  theta  = thetadb(i) * pi / 180; % ply i angle in radians, from bottom
  cs = cos(theta) ;
  si = sin(theta) ;
  T = [ cs^2 si^2 2*cs*si; si^2 cs^2 -2*cs*si; -cs*si cs*si (cs^2 - si^2)];
  Qbar = inv(T) * Q * (inv(T))' ;
  Sbar = inv(Qbar);
  
  A = A + Qbar * t_ply;
  B = B + Qbar * t_ply * zbar(i); 
  D = D + Qbar * (t_ply * zbar(i)^2  + t_ply^3 / 12);
end;

% ------------Outputs---------------------
%Qbar
%Sbar
ABD = [A B; B D] 
% ------------Partial inverse-------------
Ac = inv(A);
Bc = -Ac*B;
Dc = D-B*Ac*B;
ABDc = [Ac Bc; Bc Dc]
% ------------Replacement values----------
p       = pi*m*R/L;
alpha   = Dc(1,1) / sqrt(Ac(2,2)*Dc(1,1));
beta    = Bc(2,1) / sqrt(Ac(2,2)*Dc(1,1));
gamma   = 1       / sqrt(Ac(2,2)*Dc(1,1));
roh     = sqrt(2/(R*gamma))*p;
tau     = sqrt(2/(R*gamma))*n;            
%%%-----------a11,a12,...-----------------
for j = 0:2;
    a1j = Ac(2,2)*(2*j-1)^4*roh^4       + (2*Ac(1,2)+Ac(3,3))*(2*j-1)^2*roh^2*tau^2         + Ac(1,1)*tau^4; 
    a2j = 2*Ac(2,3)*(2*j-1)^3*roh^3+tau + 2*Ac(1,3)*(2*j-1)*roh*tau^3;
    
    b1j = Bc(2,1)*(2*j-1)^4*roh^4       + (Bc(1,1)+Bc(2,2)-2*Bc(3,3))*(2*j-1)^2*roh^2*tau^2 + Bc(1,2)*tau^4; 
    b2j = (Bc(3,1)-2*Bc(2,3))*(2*j-1)^3*roh^3*tau + (Bc(3,2)-2*Bc(1,3))*(2*j-1)*roh*tau^3;
    b1jS = b1j -(2*(2*j-1)^2*roh^2)/gamma;
    
    d1j = Dc(1,1)*(2*j-1)^4*roh^4        + 2*(Dc(1,2)+2*Dc(3,3))*(2*j-1)^2*roh^2*tau^2        + Dc(2,2)*tau^4;
    deltaj = a1j^2-a2j^2;
    
    if j == 0;
        a10 = a1j;
        a20 = a2j;
        b10 = b1j;
        b20 = b2j;
        b10S = b1jS;
        d10 = d1j;
        delta0 = deltaj;
    elseif j == 1;
        a11 = a1j;
        a21 = a2j;
        b11 = b1j;
        b21 = b2j;
        b11S = b1jS;
        d11 = d1j;
        delta1 = deltaj;
    elseif j == 2;
        a12 = a1j;
        a22 = a2j;
        b12 = b1j;
        b22 = b2j;
        d12 = d1j;
        delta2 = deltaj;
    end
end

%%%
C1 = roh^2 + (1 - 2*roh^2*beta)^2/(4*roh^2);

A1 = d11 + ((a11*b11S - a21*b21)*b11S + (a11*b21 - a21*b11S)*b21)/delta1;
A2 = (4*alpha*roh^2)/gamma;
A3 = 2*C1*mu*t*roh^2*tau^2*((a10*b10S - a20*b20)/delta0 + (a11*b11S - a21*b21)/delta1);
A4 = mu*t*tau^2*alpha*(1-2*roh^2*beta);
A5 = 4*roh^4*tau^4*C1^2*mu^2*t^2*(a10/delta0 + a12/delta2);

%%%
o1 = A2;
o2 = -(A1 + 2*A2*C1 + A4);
o3 = 2*A1*C1 + A2*C1^2 + A4*C1 + A3;
o4 = -(A1*C1^2 + A3*C1 + A5);
o5 = 27*o1^2*o4^2 - 18*o1*o2*o3*o4 + 4*o1*o3^3 + 4*o2^3*o4 - o2^2*o3^2;

%%% 3rd order equation
% o1*lambda^3 + o2*lambda^2 + o3*lambda + o4 = 0 
p = [o1 o2 o3 o4];
lambda1 = roots(p);

%%% finding the smallest lambda
lambdamin = 10000
for m = 1:20;
    for n = 0:25;
        temp_lambda = lambda1(3)
        if temp_lambda < lambdamin
            lambdamin = temp_lambda;
        end
    end
end

%%%%%%% RESULTS %%%%%%%
%%% (1) ------------------Perfect Shell------------------
lambda_crit_perf = sqrt(1+beta^2)-beta
F_krit_perf_2 = 4*pi*alpha*lambda_crit_perf

%%% (2) ------------------Non-Axisymmetric---------------
% finding the perfect shell again
F_kr_lamdamin = 4*pi*alpha*lambdamin


