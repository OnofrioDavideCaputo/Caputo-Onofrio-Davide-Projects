clear all; close all; clc;

e = 0.7; m1 = 1; m2 = 4; 
T = 2*pi./(1-e).^1.5; tspan = linspace(0,T,1000); 
x0 = -1; u0 = 0;
y0 = 0; v0 = sqrt(1 + e);
options = odeset('RelTol',1.e-6);
x1 = zeros(1, length(tspan) + 1); y1 = zeros(1, length(tspan) + 1);
x2 = zeros(1, length(tspan) + 1); y2 = zeros(1, length(tspan) + 1);

%%%%% Solve differential equations for x and y using ode45 
%%%%% with arguments tspan and options. 
%%%%% Determine x1, y1 and x2, y2 

[t,xy] = ode45(@(t,xy) orbit(xy), tspan, [x0, u0, y0, v0], options);
x = xy(:,1); y = xy(:,3);

for i = 1:length(t)
    x1(i) = -m2/(m1 + m2) * x(i);
    y1(i) = -m2/(m1 + m2) * y(i);
    x2(i) =  m1/(m1 + m2) * x(i);
    y2(i) =  m1/(m1 + m2) * y(i);
end


k = 0.1; 
R1 = k*(m1)^(1/3); R2 = k*(m2)^(1/3); %radius of masses 
theta = linspace(0,2*pi); 
figure; axis equal; hold on; set(gcf,'color','w'); 
axis off; 
xlim([-5,5]); ylim([-5,5]); 
planet = fill(R1*cos(theta)+x1(1), R1*sin(theta)+y1(1),'b'); 
sun = fill(R2*cos(theta)+x2(1), R2*sin(theta)+y2(1),'r'); 
pause(1); 
nperiods = 5;

for j = 1:nperiods 
    for i = 1:length(t) 
        planet.XData = R1*cos(theta)+x1(i); planet.YData = R1*sin(theta)+y1(i); 
        sun.XData = R2*cos(theta)+x2(i); sun.YData = R2*sin(theta)+y2(i); 
        drawnow; 
    end 
end


function dxy_dt = orbit(xy)
x = xy(1); u = xy(2); y = xy(3); v = xy(4);
dxy_dt = [u; (-x)./((x^2 + y^2)^(3/2)); v; (-y)./((x^2 + y^2)^(3/2))];
end
