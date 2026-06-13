clc, clear all, close all;

excel = xlsread('Curvas_Medidas_Motor_2026.xls');
t = excel(:,1);
TL = excel(:,5);

Ra = 0.999976;
La = 1.023575;
Km = 0.099994;
Ki = 10.010591;
J = 2.048650;
B = 0.500588; 

dt = t(2) - t(1);
[ia, wr, theta, Va, error] = deal(zeros(size(t)));
theta_ref = 1;

Kp = 1;
Ki_control = 1; 
Kd = 8;
integral = 0;
error_prev = 0;
X=[0,0,0]

for k = 1: length(t)-1
    error(k) = theta_ref - theta(k);
    
    integral = integral + error(k)*dt;
    derivada = (error(k) - error_prev)/dt;
    
    Va(k) = Kp*error(k) + Ki_control*integral + Kd*derivada;
    
    error_prev = error(k);
    X=modmotor(dt, X, [Va(k), TL(k)]);
    
    ia(k+1)=X(3);
    wr(k+1)=X(2);
    theta(k+1)=X(1);
end

subplot(4,1,1), plot(t, ia, 'k', 'LineWidth', 3); 
xlabel('t'), ylabel('i_a(t)'),grid on;
legend('Corriente del Motor');
subplot(4,1,2), plot(t, wr, 'k', 'LineWidth', 3); 
xlabel('t'), ylabel('\omega_r(t)'),grid on;
legend('Velocidad del Motor');
subplot(4,1,3), plot(t, theta, 'k', 'LineWidth', 3);
xlabel('t'), ylabel('\theta(t)'),hold on;
plot(t, theta_ref*ones(size(t)), 'r--', 'LineWidth', 3),grid on;
legend('Angulo del motor', 'Angulo de Referencia');
subplot(4,1,4), plot(t(2:end), Va(2:end), 'k', 'LineWidth', 3); 
xlabel('t'), ylabel('V_a(t)  T_L'),hold on;
plot(t(2:end), TL(2:end), 'r', 'LineWidth', 3),grid on;
legend('Accion de Control', 'Torque');