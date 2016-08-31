%% Computes the derivative of the fluxes wrt the u vector
%% Setting up the matrices
syms u v p rho;
gamma = 1.4;
e = p / ( rho * (gamma-1) );
VEL = 0.5*(u^2+v^2);
E = e + VEL;
U = [ rho, rho*u, rho*v, rho*E];
F = [rho*u,rho*u*u+p,rho*u*v,rho*u*E+p*u];
G = [rho*v,rho*u*v,rho*v*v+p,rho*v*E+p*v];
X = [p,u,v,rho];

%% Compute the jacobians

dFdX = jacobian(F,X);
dUdX = jacobian(U,X);
dGdX = jacobian(G,X);

dFdU = dFdX/dUdX;
dGdU = dGdX/dUdX;
