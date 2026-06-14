clc, clear all,close all;

Ra = 0.999976;
La = 1.023575;
Km = 0.099994;
Ki = 10.010591;
J = 2.048650;
B = 0.500588; 

t=(0:0.01:100)';
dt = t(2) - t(1);
[ia, wr, theta, Va, error] = deal(zeros(size(t)));
theta_ref = [(pi/2)*ones(2500,1);-(pi/2)*ones(2500,1);(pi/2)*ones(2500,1);-(pi/2)*ones(2501,1)];
TL = [10*ones(2500,1);zeros(2500,1);10*ones(2500,1);zeros(2501,1)];

integral = 0;
error_prev = 0;
X=[0;0;0];
mat_A=[0, 1, 0;0, -B/J, Ki/J; 0, -Km/La, -Ra/La];
mat_B=[0; 0; 1/La];
mat_C=[1,0,0];
mat_Co=[1, 0, 0; 0, 1, 0];
mat_Aa=[mat_A,zeros(3,1);-mat_C, 0];
mat_Ba=[mat_B;0];
q=diag([ 1, 1/8^2, 1/3^2, 15000]);r=1;
Ka=lqr(mat_Aa, mat_Ba, q, r);
eig(mat_Aa-mat_Ba*Ka)

for k = 1: length(t)-1
    error(k) = theta_ref(k) - theta(k);
    
    integral = integral + error(k)*dt;
    Va(k)=-Ka*[X ; integral];

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
plot(t, theta_ref, 'r--', 'LineWidth', 3),grid on;
legend('Angulo del Motor', 'Angulo de Referencia');
subplot(4,1,4), plot(t(2:end), Va(2:end), 'k', 'LineWidth', 3); 
xlabel('t'), ylabel('V_a(t)  T_L'),hold on;
plot(t(2:end), TL(2:end), 'r', 'LineWidth', 3),grid on;
legend('Accion de Control', 'Torque');