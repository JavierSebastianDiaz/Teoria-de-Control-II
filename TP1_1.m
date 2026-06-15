clc, clear all, close all;

R = 2200;
L = 0.5;
C = 1e-5;

A = [-R/L -1/L; 1/C 0];
B = [1/L; 0];
C = [R 0];
D = 0;

t = 0:1e-4:0.5;
u = 12*square(2*pi*4*t); %se cambio la frecuencia de la señal cuadrada para obtener curvas similares a las de la figura 1-2

sys_ss = ss(A, B, C, D);

[Vo, t, x] = lsim(sys_ss, u, t);
i = x(:, 1);
Vc = x(:, 2);

figure(1);
plot(t, u, '--r', 'LineWidth', 3), hold on, plot(t, Vo, 'k', 'LineWidth', 3), xlabel('t'), ylabel('V'), legend('V_e', 'V_o'), grid on;

figure(2);
subplot(4,1,1), plot(t, u, 'k', 'LineWidth', 3), ylabel('V_e(t)'), grid on;
subplot(4,1,2), plot(t, Vo, 'k', 'LineWidth', 3), ylabel('V_o(t)'), grid on;
subplot(4,1,3), plot(t, i, 'k', 'LineWidth', 3), ylabel('i(t)'), grid on;
subplot(4,1,4), plot(t, Vc, 'k', 'LineWidth', 3), ylabel('V_c(t)'), xlabel('t'), grid on;