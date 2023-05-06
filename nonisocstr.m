function f = nonisocstr(t, x)

Ca = x(1);
T = x(2);

% Data
ko = 460; % pre-exponential factot
E = 1380; % Activation engergy
Cao = 0.4; % initial concentration of compoent A
tau = 0.18; % residence time
Tc = 298.15;
k = 79;
Cp = 32;

Hr = -151080 + 2*(T - 298.15);
ra = -ko*exp(-E/T);

% the modelling equations are
dCadt = (Cao - Ca)/tau + (ra * Ca);
dTdt = (-Hr/Cp) * (-(ra * Ca)/Cao) - (((1 + k)*(T - Tc))/tau);
f = [dCadt; dTdt];