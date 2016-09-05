function [F,H] = construct_fluxes(U,J)
% *************************************************************************
% Set intermediate fluxes from U
% *************************************************************************
%% Set up the fluxes
gamma = 1.4;
F = zeros(size(U));
G = F;
H = F;
%% Assign values
for i = 1:size(U,3)
    for j = 1:size(U,2)
        [u,v,p,rho] = get_primitives(U(:,j,i));
        e = p / ( rho * (gamma-1) );
        VEL = sqrt(u^2 + v^2);
        E = e + 0.5 * VEL^2;
        
        F(:,j,i) = [rho*u    rho*u^2+p    rho*u*v     rho*u*E*+p*u];
        G(:,j,i) = [rho*v    rho*u*v      rho*v^2+p   rho*v*E+p*v];
        H(:,j,i) = J(3,j,i) .* F(:,j,i,1) + J(4,j,i) .* G(:,j,i,1);
    end
end