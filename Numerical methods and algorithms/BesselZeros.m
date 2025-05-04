clear all; close all; clc;

J = @(n, x) integral(@(theta) (1/pi).*cos(x.*sin(theta) - n.*theta), 0, pi);
xmax = 25;
max_n = 4;
max_zeros = 5;
X = linspace(0, xmax, 1000);
x0 = 0;
ZERO = zeros(max_n+1, max_zeros);

zeros_guess=[2.4,5.5,8.6,11.8,15;
             3.8,7,10,13,16.4;
             5.1,8.4,11.6,15,18;
             6,9.7,13,16,19.4;
             7.5,11,14,18,21;
             8.7,12,16,19,22];



figure; hold on; grid on; 

for n = 0:max_n
    Jn = @(x) J(n, x);
    for i = 1:max_zeros
        ZERO(n + 1, i) = fzero(Jn, zeros_guess(n+1, i));
    end
    Jn_values = arrayfun(@(x) Jn(x), X);
    plot(X, Jn_values, 'MarkerSize', 2);
end

axis([-1,30,-3,3])
xlabel('x');
ylabel('y');
title('Bessel functions');