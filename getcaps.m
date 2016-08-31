%% Generate reconstructed fluxes from weno5
function [fcap_r,fcap_l] = getcaps(U,F,G,dx,dy,t)
%% Split the flux
nX = size(F,3);
nY = size(F,2);
gamma = 1.4;
dFdU = @(u,v,p,rho) [                                            0,                                     1,          0,       0;
                                                 v^2/5 - (4*u^2)/5,                               (8*u)/5,   -(2*v)/5,     2/5;
                                                              -u*v,                                     v,          u,       0;
                      -(3*rho*u^3 + 3*rho*u*v^2 + 35*p*u)/(10*rho), (rho*u^2 + 5*rho*v^2 + 35*p)/(10*rho), -(2*u*v)/5, (7*u)/5];
dGdU = @(u,v,p,rho) [                                            0,          0,                                     1,       0;
                                                              -u*v,          v,                                     u,       0;
                                                 u^2/5 - (4*v^2)/5,   -(2*u)/5,                               (8*v)/5,     2/5;
                      -(3*rho*u^2*v + 3*rho*v^3 + 35*p*v)/(10*rho), -(2*u*v)/5, (5*rho*u^2 + rho*v^2 + 35*p)/(10*rho), (7*v)/5] ;                 
maxf = zeros(4,1);
maxg = maxf;
for i = 1:nX
    for j = 1:nY
        %  extract primitive variables
        rho = U(1,j,i,t);
        u = U(2,j,i,t) / U(1,j,i,t);
        v = U(3,j,i,t) / U(1,j,i,t);
        V = sqrt(u*u+v*v);
        p = (U(4,j,i,t)/U(1,j,i,t) - 0.5*V) * rho * (gamma-1);
        A1 = max(dFdU(u,v,p,rho),[],2);
        A2 = max(dGdU(u,v,p,rho),[],2);
        
        maxf = 0.5 * (A1+maxf+ abs(A1-maxf));
        maxg = 0.5*( A2+maxg+abs(A2-maxg));
    end
end

%% iteration to reconstruct the parts
for i = 4:nX-3
    for j = 4:nY-3
        ;
    end
end
fcap_l=0;
fcap_r = 0;