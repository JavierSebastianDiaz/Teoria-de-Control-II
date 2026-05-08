clc, clear all, close all;

excel = xlsread('Curvas_Medidas_RLC_2026.xls');
t = excel(:,1);
i = excel(:,2);
Vc = excel(:,3);
u = excel(:,4);
Vo = excel(:,5);

t0 = find(Vc>0 ,1);
x = 200;
t1 = t0+x;
t2 = t0+x*2;
t3 = t0+x*3;
% plot(Vc(1:5000)),hold on,scatter(t1,Vc(t1)),hold on,scatter(t2,Vc(t2)),hold on,scatter(t3,Vc(t3)),grid on

K = max(Vc);
k1 = Vc(t1)/K - 1;
k2 = Vc(t2)/K - 1;
k3 = Vc(t3)/K - 1;
b = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3;
a1 = (k1*k2 + k3 - sqrt(b))/(2*(k1^2 + k2));
a2 = (k1*k2 + k3 + sqrt(b))/(2*(k1^2 + k2));
%beta = (2*k1^3 + 3*k1*k2 + k3 - sqrt(b))/sqrt(b)
T1 = -(t(t1)-t(t0))/log(a1);
T2 = -(t(t1)-t(t0))/log(a2);
%T3 = beta*(T1 - T2) + T1

%num = [K*T3 K];
num = [1];
den = conv([T1 1],[T2 1]);
sys = tf(num, den);

[Vc_aprox, x2] = lsim(sys, u, t);
Vc_p_aprox = diff(Vc_aprox)/(t(1)-t(2));
C = max(i)/max(Vc_p_aprox)
%          1
%   ---------------
%   s^2LC + sRC + 1
L = den(1,1)/C
R = den(1,2)/C