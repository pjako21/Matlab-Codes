function dxdt = contcstr(t, C)
%
% Solves the two differential equations modeling
%am series of CSTR. The are the concentration
% of A, B and C in the reactor. 
% A + B ----> C

Ca1 = C(1);
Ca2 = C(3);
Ca3 = C(5);
Cb1 = C(2);
Cb2 = C(4);
Cb3 = C(6);

%Data
v0 = 6; % inlet flowrate to the first reactor
v = 12; % volume of the reactor vessels
V = 200;
Ca0 = 2; % inital concentration A
Cb0 = 2; 
k = 0.5; 

% the modeling equations are;
dadt = (v0*Ca0-v*Ca1-k*V*Ca1*Ca2)/V;
dbdt = (v0*Cb0-v*Cb1-k*V*Ca1*Cb1)/V;
dcdt = (v*Ca1-v*Ca2-k*V*Ca2*Cb2)/V;
dddt = (v*Cb1-v*Cb2-k*V*Ca2*Cb2)/V;
dedt =(v*Ca2-v*Ca3-k*V*Ca3*Cb3)/V; 
dfdt = (v*Cb2-v*Cb3-k*V*Ca3*Cb3)/V;

dxdt = [dadt; dbdt; dcdt; dddt; dedt; dfdt];
end


