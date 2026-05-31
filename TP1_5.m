clc, clear all, close all; 

excel = xlsread('Curvas_Medidas_Motor_2026.xls');
t  = excel(:,1);
wr = excel(:,2);
ia = excel(:,3);
Va = excel(:,4);
TL = excel(:,5);

t0 = find(Va>0 ,1);
x = fix(find(wr==max(wr), 1, 'first')/5);
t1 = t0+x;
t2 = t0+x*2;
t3 = t0+x*3;

step = max(Va);
M = [wr ia];
[K T1 T2 T3] = deal(zeros(1, 2));
for k = 1:2
    K(k) = M(end,k)/step;
    k1 = (1/step)*M(t1,k)/K(k) - 1;
    k2 = (1/step)*M(t2,k)/K(k) - 1;
    k3 = (1/step)*M(t3,k)/K(k) - 1;
    b = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3;
    a1 = (k1*k2 + k3 - sqrt(b))/(2*(k1^2 + k2));
    a2 = (k1*k2 + k3 + sqrt(b))/(2*(k1^2 + k2));
    % beta = (2*k1^3 + 3*k1*k2 + k3 - sqrt(b))/sqrt(b);
    beta = (k1 + a2)/(a1-a2);
    T1(k) = -(t(t1)-t(t0))/log(a1);
    T2(k) = -(t(t1)-t(t0))/log(a2);
    T3(k) = beta*(T1(k) - T2(k)) + T1(k);
end
    
sys_wr = tf(K(1)*[T3(1) 1], conv([T1(1) 1],[T2(1) 1]));
sys_ia = tf(K(2)*[T3(2) 1], conv([T1(2) 1],[T2(2) 1]));
y_wr = lsim(sys_wr, Va, t, [0;0]);
y_ia = lsim(sys_ia, Va, t, [0;0]);
subplot(2,1,1), plot(t, wr, 'k', 'LineWidth', 3), hold on, plot(t, y_wr, '--r','LineWidth', 3)
plot(t(t1),wr(t1),'ob', 'LineWidth', 5),plot(t(t2),wr(t2),'ob', 'LineWidth', 5),plot(t(t3),wr(t3),'ob', 'LineWidth', 5);
title('Modelado de Motor CC', 'Fontsize', 28);
legend('Datos','Modelo', 'Puntos Chen', 'Fontsize', 16);
ylabel('\omega_r(t) [rad/s]', 'Fontsize', 28);
grid on, grid minor;
subplot(2,1,2), plot(t, ia, 'k', 'LineWidth', 3), hold on, plot(t, y_ia, '--r','LineWidth', 3)
plot(t(t1),ia(t1),'ob', 'LineWidth', 5),plot(t(t2),ia(t2),'ob', 'LineWidth', 5),plot(t(t3),ia(t3),'ob', 'LineWidth', 5);
ylabel('i_a(t) [A]', 'Fontsize', 28), xlabel('t [s]', 'Fontsize', 28);
grid on, grid minor;

inicio_TL = find(TL>0, 1);
final_TL = find(TL>0, 1, 'last');
K_TL = (wr(inicio_TL) - wr(final_TL)) / max(TL);

Ra = (-T1(1)*T2(1) + T1(1)*T3(2) + T2(1)*T3(2))/(K(2)*T3(2)^2);
La = T1(1)*T2(1)/(K(2)*T3(2));
Km = (T1(1)*T2(1) + T3(2)^2 - T3(2)*(T1(1) + T2(1)))/(K(1)*T3(2)^2);
Ki = K(1)*Ra/K_TL;
J  = K(2)*Ki*T3(2)/K(1);
B  = K(2)*Ki/K(1);
fprintf('Ra = %f\nLa = %f\nKm = %f\nKi = %f\nJ  = %f\nB  = %f\n',Ra,La,Km,Ki,J,B)