function [u,v,p,rho] = get_primitives(U)
% *************************************************************************
% Extract primitve values from U
% *************************************************************************
%% Extract primitive values
gamma = 1.4;
rho = U(1);
u = U(2) / U(1);
v = U(3) / U(1);
V = sqrt(u*u+v*v);
p = (U(4)/U(1) - 0.5*V) * rho * (gamma-1);