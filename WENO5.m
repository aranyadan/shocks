    %% Calculation of fcap
function res = WENO5(w,f,dx,i,j)
%% Right Flux
% Choose the positive fluxes, 'v', to compute the left cell boundary flux:
% $u_{i+1/2}^{-}$

        vmm = f(i-2,j);
        vm  = f(i-1,j);
        vp  = f(i+1,j);
        vpp = f(i+2,j);
        v   = f(i,j);
% Polynomials
p0n = (2*vmm - 7*vm + 11*v)/6;
p1n = ( -vm  + 5*v  + 2*vp)/6;
p2n = (2*v   + 5*vp - vpp )/6;

% Smooth Indicators (Beta factors)
B0n = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2; 
B1n = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2;
B2n = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% Constants
d0n = 1/10; d1n = 6/10; d2n = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alpha2n = d2n./(epsilon + B2n).^2;
alphasumn = alpha0n + alpha1n + alpha2n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w2n = alpha2n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n + w2n.*p2n;

%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i+1/2}^{+}$ 
   umm = f(i-1,j);
   um  = f(i+1,j);
   up  = f(i+2,j);
   upp = f(i+3,j);
   u   = f(i,j);
  % umm = f(i,j-1);
  % um  = f(i,j+1);
  % up  = f(i,j+2);
  % upp = f(i,j+3);
% Polynomials
p0p = ( -umm + 5*u + 2*um )/6;
p1p = ( 2*u + 5*um  - up   )/6;
p2p = (11*u  - 7*up + 2*upp)/6;

% Smooth Indicators (Beta factors)
B0p = 13/12*(umm-2*u+um  ).^2 + 1/4*(umm-4*u+3*um).^2; 
B1p = 13/12*(u -2*um +up ).^2 + 1/4*(u-up).^2;
B2p = 13/12*(um  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2;

% Constants
d0p = 3/10; d1p = 6/10; d2p = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alpha2p = d2p./(epsilon + B2p).^2;
alphasump = alpha0p + alpha1p + alpha2p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;
w2p = alpha2p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p + w2p.*p2p;

tt=w(i+1,j);
pp=w(i,j);
abar=(vp-v)./(tt-pp);
if(abar>=0)
    fc=hp/dx;
else
    fc = hn/dx;
end
res = fc;
