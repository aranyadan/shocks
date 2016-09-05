%% Generate reconstructed fluxes from weno5
function [resx,resy] = get_derivs(U,F,H,dX,dY,J)
%% compute alpha for lax Friedrich's splitting
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
maxh = maxf;
for i = 1:nX
    for j = 1:nY
        %  extract primitive variables
        [u,v,p,rho] = get_primitives(U(:,j,i));
        A1 = max(dFdU(u,v,p,rho),[],2);
        A2 = dFdU(u,v,p,rho).* J(3,j,i) + dGdU(u,v,p,rho) .* J(4,j,i);
        A3 = max(A2,[],2);
        
        maxf = 0.5 * (A1+maxf+ abs(A1-maxf));
        maxh = 0.5 * (A3+maxh+ abs(A3-maxh));
    end
end
F_p = zeros(4,nY,nX);
F_m = F_p;
H_p = F_p;
H_m = F_p;

%% split the fluxes
for i = 1:nX
    for j = 1:nY
        for k = 1:4
            F_p(k,j,i) = 0.5*(F(k,j,i) + maxf(k).*U(k,j,i));
            F_m(k,j,i) = 0.5*(F(k,j,i) - maxf(k).*U(k,j,i));
            H_p(k,j,i) = 0.5*(H(k,j,i) + maxh(k).*U(k,j,i));
            H_m(k,j,i) = 0.5*(H(k,j,i) - maxh(k).*U(k,j,i));
        end
    end
end

%% compute derivatives
resx = WENO5(F_p,F_m,dX,1);
resy = WENO5(H_p,H_m,dY,2);

%% pad the derivatives with 0th order interpolation
for i = 1:nX
    for j = 1:nY
        if(i<4)
            resx(:,j,i) = resx(:,j,4);
        end
        if(j<4)
            resy(:,j,i) = resy(:,4,i);
        end
        if(i>nX-3)
            resx(:,j,i) = resx(:,j,nX-3);
        end
        if(j>nY-3)
            resy(:,j,i) = resy(:,nY-3,i);
        end
    end
end