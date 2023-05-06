function xdot = vdv_ode(t, x)
%
% Solventhe two differential equations modeling
% the van de vusse reaction
%scheme in an isothermal CSTR. The are the concentration
% of A and B in the reactor.
% 
%[t, x] = ode45(vdv_ode, [0 5], x0]
% integrate from t = 0 to t = 5 min, with initial conditions
% ca0 = x0(1) and cb0 = x0(2), and x0 is a column vector 
%b.w. bequette

%since the states arepassed to this routine in the x vector, 
%convert to natrual notation

ca = x(1);
cb = x(2);

% the parameters are:
k1 = 5/6; %rate constant for A --> B (min^-1)
k2 = 5/3; %rate constant for B --> C (min^-1)
k3 = 1/6; %rate constant for 2A --> D (mol/(min^-1))
% the input values are:

fov = 4/7; % dilution rate (min^-1)
caf = 10; % mol/1

% the modeling equations are;
dcadt = fov*(caf - ca) - k1*ca - k3*ca^2;
dcbdt = -fov*cb + k1*ca - k2*cb;

xdot = [dcadt; dcbdt];
end

