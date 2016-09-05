tic;
clc;
clear all;
%% grid generation
% zeta = x
% eta = (y - ys(x)) / h(x)

wedges = 1;
thetad = 45;
nt = 5;
DISPLAYMESH = 0;
A = 2; % wedge start x
H = 2.5; % height in the beginning
dz = 0.1; finalz = 3.9;
nz = finalz/dz + 1;
zeta = 0:dz:(nz-1)*dz;
nx = nz;

det = 0.02;
eta = 0:det:1;
neta = (1/det)+1;
ny = neta;

if DISPLAYMESH == 1
    x = zeros(length(zeta),neta);
    y = zeros(length(zeta),neta);
    for i = 1:length(zeta)
        if(zeta(i)<=A)
            h = H;
            ys = 0;
        elseif(zeta(i)>A)
            h = H - (zeta(i)-A) * tand(thetad);
            ys = (zeta(i)-A) * tand(thetad);
        end
        y(i,:) = eta(:).*h + ys;

        xtemp = zeros(1,neta);
        xtemp(:) = zeta(i);
        x(i,:) = xtemp;
    end
    hold on;
    for i = 1:size(x,2)
        scatter(x(:,i),y(:,i),5,'filled');
        plot(x(:,i),y(:,i));
    end
end

% From now, X = zeta, Y = eta
X = zeta;
Y = eta;
nX = nz;
nY = neta;
dX = dz;
dY = det;

C = 0.5;
dt = C * dY;


%% Jacobian matrix
% [dzeta , deta] = J [dx , dy]
J = zeros(4,nY,nX);
for i = 1:nX
    for j = 1:nY
        J(1,j,i) = 1; % dzeta/dx 
        J(2,j,i) = 0; % dzeta/dy
        if ( X(i) <= A )
            h = H;
            J(3,j,i) = 0;  % deta/dx
            J(4,j,i) = 1/h;  % deta/dy
        elseif ( X(i) > A )
            h = H - ( X(i) - A ) * tand(thetad);
            J(3,j,i) = ( (Y(j) - 1)/h )*tand(thetad);  % deta/dx
            J(4,j,i) = 1/h; % deta/dy
        end
    end
end
%% Get initial conditions
Minf = 4;            % initial Mach no
rhoinf = 1.225;      % Density of air
pinf = 1.01325e5;    % Atmospheric Pressure
cv = .718;           % Specific Heat at Constant Volume
R = 8.3144598;       % Gas constant
gamma = 1.4;         % for air
ainf = sqrt(gamma * pinf / rhoinf);

U=zeros(4,nY,nX,nt);
F = U;
G = U;
H = U;

for i = 1:nX
    for j=1:nY
        p = pinf;
        rho = rhoinf;
        e = p / ( rho * (gamma-1) );
        v = 0;
        if(j == 0)
            u = 0;
        else
            u = Minf * ainf;
        end
        VEL = sqrt(u^2 + v^2);
        E = e + 0.5 * VEL^2;
        U(:,j,i,1) = [rho      rho*u        rho*v       rho*E];
        F(:,j,i,1) = [rho*u    rho*u^2+p    rho*u*v     rho*u*E*+p*u];
        G(:,j,i,1) = [rho*v    rho*u*v      rho*v^2+p   rho*v*E+p*v];
        H(:,j,i,1) = J(3,j,i) .* F(:,j,i,1) + J(4,j,i) .* G(:,j,i,1);
    end
end

%% Solver Loop
for i = 2:nt
    U_0 = U(:,:,:,i-1);
    
    [F_0,H_0] = construct_fluxes(U_0,J);
    [resx,resy] = get_derivs(U_0,F_0,H_0,dX,dY,J);
    
    U_1 = U_0 - dt * (resx+resy);   % RK3 step 1
    
    [F_1,H_1] = construct_fluxes(U_1,J);
    [resx,resy] = get_derivs(U_1,F_1,H_1,dX,dY,J);
    
    U_2 = (3/4) * U_0 + (1/4) * U_1 - (1/4) * dt * (resx+resy);  % RK3 Step 2
    
    [F_2,H_2] = construct_fluxes(U_2,J);
    [resx,resy] = get_derivs(U_2,F_2,H_2,dX,dY,J);
    
    U_3 = (1/3) * U_0 + (2/3) * U_2 - (2/3) * dt * (resx+resy);  % RK3 Step 3
    
    for j = 1:nX
        for k = 1:nY
            U(:,k,j,i) = U_3(:,k,j);
        end
    end
    
    percent = i*100/nt;
    display(['Calculated for time step ' num2str(i) ' of ' num2str(nt) ' (' num2str(percent) '%)']);
end
toc;