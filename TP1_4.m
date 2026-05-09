clc, clear all, close all;

La = 3.66e-4;
J = 5e-9;
Ra = 55.6;
B = 0;
Ki = 6.49e-3;
Km = 6.54e-3;
dt = 1e-7;
t = 0:dt:0.1;
[ia, wr, theta] = deal(zeros(size(t)));
Va = 12;
TL = 0;

for k = 1: length(t)-1
    dia = -(Ra/La)*ia(k) -(Km/La)*wr(k) +(1/La)*Va;
    dwr = (Ki/J)*ia(k) -(B/J)*wr(k) -(1/J)*TL;
    dtheta = wr(k);
    %x(k+1) = x(k) + dt*dx Integracion de Euler
    ia(k+1) = ia(k) + dt*dia;
    wr(k+1) = wr(k) + dt*dwr;
    theta(k+1) = theta(k) + dt*dtheta;
end

Tmax = Ki * max(ia)

subplot(2,1,1), plot(t, ia), xlabel('t'), ylabel('i_a(t)');
subplot(2,1,2), plot(t, wr), xlabel('t'), ylabel('w_r(t)');