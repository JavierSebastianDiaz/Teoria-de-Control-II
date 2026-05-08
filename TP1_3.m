clc, clear all, close all;

excel = xlsread('Curvas_Medidas_RLC_2026.xls');
t = excel(:,1);
i = excel(:,2);
Vc = excel(:,3);
u = excel(:,4);
Vo = excel(:,5);

R = 219.9786;
L = 1.0006;
C = 2.2002e-04;

A = [-R/L -1/L; 1/C 0];
B = [1/L; 0];
C = [R 0];
D = 0;

sys_ss = ss(A, B, C, D);

[Vo_aprox, t, x] = lsim(sys_ss, u, t);
i_aprox = x(:, 1);
Vc_aprox = x(:, 2);

plot(t, i,'LineWidth', 3), hold on, plot(t, i_aprox, '--','LineWidth', 3);
xlabel('t'),ylabel('V');
legend('V_e', 'V_o');
grid on