% per dei valori differenti di n  mi ricavo il valore della delta man mano
% che andiamo avanti e questo ci permette di valutare il valore del limite 

%all'interno di questo limite devo andare a valutare il valore di m_n che
%non conosciamo. Lo posso fare partendo da m_0 ed m_1 che possiamo ricavare
%analiticamente, mentre il valore di m_2 lo ricaviamo facendo un'assunzione
%iniziale sul valore di delta_1. Il valore di m_n lo ricaviamo mediante il
%metodo di Newton che dovrebbe convergere rapidamente dopo qualche
%iterazione.

%i valori da inserire all'interno del metodo di Newton devono essere
%calcolati mediante la mappa logistica che deve essere iterata un tot
%numero di volte per farle passare il transiente.

clear all; clc;

x0 = 0.5;
dx0 = 0;
iterations = 500;
nmax = 13;
tol = 0.0001;

m = zeros(1, nmax);
delta = zeros(1, nmax);
r = roots([1,-4,0,8]);
m(1) = r(2);
m(2) = r(1);
delta(2) = 5;

convergenza = false;

for n = 3:nmax
    m(n) = m(n-1) + (m(n-1) - m(n-2))./(delta(n-1));
    mu = m(n);
    while ~convergenza
        x = x0;
        dx = dx0;
        for i = 1:2^n
            dx = x*(1 - x) + mu*(1 - 2*x)*dx;
            x = mu*x*(1 - x);
        end
        mu1 = mu - (x - 0.5)./(dx);
        if abs(mu1 - mu) < tol
            convergenza = true;
            m(n) = mu1;
        else
            mu = mu1;
        end
    end
    delta(n) = (m(n-1) - m(n-2))./(m(n) - m(n-1));
    convergenza = false;
end


