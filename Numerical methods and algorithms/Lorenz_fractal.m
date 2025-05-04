clear all; close all; clc;

r = 28;
sigma = 10;
beta = 8./3;
y0 = 3*sqrt(2);
eps = 0.0001;
iter = 1001;
maxIter = 1000;

root1 = [0, 0, 0]; root2 = [sqrt(beta*(r-1)), sqrt(beta*(r-1)), r-1]; root3 = [-sqrt(beta*(r-1)), -sqrt(beta*(r-1)), r-1];
xmin = -100; xmax = 100; zmin = -100; zmax = 100;

Z = zeros(iter, iter);
xlist = linspace(xmin, xmax, iter);
zlist = linspace(zmin, zmax, iter);

for i = 1:iter
    for j = 1:iter
        x = xlist(i);
        z = zlist(j);
        [x1, y1, z1] = Newton(r, sigma, beta, x, y0, z, maxIter);
        if max(abs([x1, y1, z1] - root1)) < eps
            pixel = 1;
        elseif max(abs([x1, y1, z1] - root2)) < eps
            pixel = 2;
        elseif max(abs([x1, y1, z1] - root3)) < eps
            pixel = 3;
        else
            pixel = 4;
        end
        %i = round(1 + (iter - 1) * (x - xmin)./(xmax - xmin));
        %j = round(1 + (iter - 1) * (z - zmin)./(zmax - zmin));
        Z(i, j) = pixel;
    end
end



figure;
map = [1 0 0; 0 1 0; 0 0 1; 0 0 0]; colormap(map);

image([xmin xmax], [zmin zmax], Z);
set(gca, 'YDir', 'normal');


function [x, y, z] = Newton(r, sigma, beta, x0, y0, z0, maxNewtoniter)

threshold = 1e-3;
converge = false;
iter_count = 0;

x = x0;
y = y0;
z = z0;

while ~converge && iter_count < maxNewtoniter
    J = [-sigma sigma 0; r-z -1 -x; y x -beta];
    b = [sigma*(y-x); x*(r-z)-y; x*y-beta*z];
    X = J \ (-b);
    x = x + X(1);
    y = y + X(2);
    z = z + X(3);
    iter_count = iter_count + 1;
    if abs(max(X)) < threshold 
        converge = true;
    end
end

end
