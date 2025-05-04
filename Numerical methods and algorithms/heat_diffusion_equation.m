clear all; close all; clc;

%%%%% Define the square and grid parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 1; %square is 2L x 2L 
N = 100; %# of intervals in x and y directions 
n = N+1; %# of gridpoints in x,y directions including boundaries 
h = 2*L/N; %grid size in x,y directions 
x = -L + (0:N)*h; %x values on the grid 
y = -L + (0:N)*h; %y values on the grid 
[X,Y] = meshgrid(x,y); 

%%%%% Definition of the indices associated with the boundaries %%%%%%%%%%%%%%%%%%% 

% boundary_index = [bottom, left, top, right] 
boundary_index = [ 1:n, 1:n:1+(n-1)*n, ... 
    1+(n-1)*n:n*n, n:n:n*n ];

%%%%% Diffusion constant and time-step parameters %%%%%%%%%%%%%%%%%%%

D = 1; 
dt = h^2/(2*D); %borderline stability of FTCS scheme 
alpha = dt*D/h^2; %equation parameter 
nsteps = 1000; %number of time steps

%%%%% CONSTRUCT THE MATRIX AND COMPUTE LU DECOMPOSITION %%%%%%%%%%%%%%%%%%%%

nx = n;
ny = n;
Ntot = nx * ny;

Ix = speye(nx);
Iy = speye(ny);

% Tridiagonal matrices construction
ex = alpha * ones(nx,1);
Tx = spdiags([ex -2*ex ex], [-1 0 1], nx, nx);

ey = alpha * ones(ny,1);
Ty = spdiags([ey -2*ey ey], [-1 0 1], ny, ny);

% Kronecker product
Lmat = kron(Iy, Tx) + kron(Ty, Ix);

% System matrix for the implicit equation
A = speye(Ntot) - alpha * Lmat;

% Dirichlet conditions
for j = 1:n
    for i = 1:n
        if i == 1 || i == n || j == 1 || j == n
            idx = (j-1)*n + i;
            A(idx,:) = 0;
            A(idx,idx) = 1;
        end
    end
end

% LU decomposition for the linear system
[Lf, Uf] = lu(A);


%%%%% Define initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = zeros(n,n,nsteps); 
sigma = L/4; 
u(:,:,1) = 1/(2*pi*sigma^2)*exp(-0.5*(X.^2+Y.^2)/sigma^2); 
u(1,:,1) = 0; u(n,:,1) = 0; u(:,1,1) = 0; u(:,n,1) = 0; %b.c. 

%%%%% ADVANCE SOLUTION u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for m = 2:nsteps
    % Column vector reshaping for the previous solution
    u_vec_prev = reshape(u(:,:,m-1), [], 1);

    % Implicit system evaluation
    u_vec_next = Uf \ (Lf \ u_vec_prev);

    % 2D matrix reconstruction with boundary conditions
    u_next = reshape(u_vec_next, n, n);
    u_next(1,:) = 0; u_next(end,:) = 0;
    u_next(:,1) = 0; u_next(:,end) = 0;

    % Storing of the result
    u(:,:,m) = u_next;
end

%%%% Plot with animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

figure('units','normalized','outerposition',[0 0 1 1]) 
s=surf(X,Y,u(:,:,1)); zlim([0, 2.6]); 
xlabel('$x$','Interpreter','latex','FontSize',14); 
ylabel('$y$','Interpreter','latex','FontSize',14); 
zlabel('$u(x,y,t)$','Interpreter','latex','FontSize',14); 
title('Solution of the 2D diffusion equation','Interpreter','latex','FontSize',16); 
pause(1) 

% Animation step
for j = 2:nsteps 
    s.ZData = u(:,:,j); 
    zlim([0, max(u(:,:,j), [], 'all')]);
    pause(0.01); 
end
