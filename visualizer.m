function [] = visualizer(U,X,Y)
% *************************************************************************
% Function to visualize any flux/derivative in 2d space
% *************************************************************************
%% Check number of plots
num = size(U,1);
nr = max([floor(num/2) 1]);
nc = ceil(num/2);
%% reshape the vectors
x = reshape(X,[size(X,1)*size(X,2),1]);
y = reshape(Y,size(x));
%% Plots
for i = 1:num
    subplot(nr,nc,i);
    z = reshape(U(i,:,:),size(x));
    % Delaunay triangulation
    tri = delaunay(x,y);
    trisurf(tri,x,y,z);
    view(2);
    title(['i=' num2str(i)]);
    xlabel('x');
    ylabel('y');
    grid on; box on;
    shading interp;
end

