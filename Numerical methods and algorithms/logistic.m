clear all; close all; clc;

mu = linspace(0, 4, 1000);
fixed = [];
muvals = [];
transient = 400;
iterations = 500;

for el = mu
    x = 0.5;
    for i=1:iterations
        x = el.*x.*(1 - x);
        if i > transient
            fixed = [fixed, x];
            muvals = [muvals, el];
        end
    end
end


figure;
plot(muvals, fixed, '.', 'MarkerSize', 1);
axis([-1,5,-0.2,1.2])
xlabel('\mu');
ylabel('Punti fissi stabili x^*');
title('Punti fissi stabili della mappa logistica');
grid on;



function df = derivate(x, mu)
    df = mu - 2 .* mu .* x; 
end