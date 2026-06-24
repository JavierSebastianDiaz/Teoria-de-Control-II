clc, clear all, close all;

excel = xlsread('Curvas_Medidas_Motor_2026.xls');
t = excel(:,1);
TL = excel(:,5);
[ia, wr, theta, Va, error] = deal(zeros(size(t)));
theta_ref = [(pi/2)*ones(1500,1);-(pi/2)*ones(1500,1);(pi/2)*ones(1000,1)];

%Parámetros calculados
Ra = 0.999976;
La = 1.023575;
Km = 0.099994;
Ki = 10.010591;
J = 2.048650;
B = 0.500588; 

%Cálculo del Controlador
mat_A=[0, 1, 0;0, -B/J, Ki/J; 0, -Km/La, -Ra/La];
mat_B=[0; 0; 1/La];
mat_C=[1,0,0];
mat_Co=[1, 0, 0; 0, 1, 0];
mat_Aa=[mat_A,zeros(3,1);-mat_C, 0];
mat_Ba=[mat_B;0];
% q=diag([ 1, 1/8^2, 1/3^2, 15000]);r=1; %sin restricciones en la accion de control
q=diag([ 1, 1/8^2, 1/3^2, 15000]);r=1e4; %con restricciones
Ka=lqr(mat_Aa, mat_Ba, q, r);
eig(mat_Aa-mat_Ba*Ka)

%Cálculo del Observador
mat_Ao= mat_A';
mat_Bo= mat_Co';
q_o=diag([1,1,1]);
r_o=diag([1e1,1e1]);
K_o=lqr(mat_Ao, mat_Bo, q_o, r_o)';
eig(mat_A-K_o*mat_Co)

dt = 0.01;
integral = 0;
X=[0;0;0];
X_ang=[0;0;0];
NL=5;
for k = 1: length(t)-1
    error(k) = theta_ref(k) - theta(k);
    Y=[theta(k); wr(k)];
    integral = integral + error(k)*dt;
    
    Va(k)=-Ka*[X_ang ; integral]; %con observador
%     Va(k)=-Ka*[X ; integral]; %sin observador
    
    %No linealidad de la acción de control
    if Va(k)>NL
        Va(k)=Va(k)-NL;
    end
    if Va(k)<-NL
        Va(k)=Va(k)+NL;
    end
    if Va(k)<NL && Va(k)>-NL
        Va(k)=0;
    end
    
    X=modmotor(dt, X, [Va(k), TL(k)]);
    X_ang_p=mat_A*X_ang+mat_B*Va(k)+K_o*(Y-mat_Co*X_ang);
    X_ang= X_ang+X_ang_p*dt;

    ia(k+1)=X(3);
    wr(k+1)=X(2);
    theta(k+1)=X(1);
end

subplot(4,1,1), plot(t, ia, 'k', 'LineWidth', 3); 
ylabel('i_a(t)', 'Fontsize', 28),grid on;
legend('Corriente del Motor');
% title('Control sin restricciones ni observador de Motor CC', 'Fontsize', 28);
% title('Control restringiendo la tensión a 5 V', 'Fontsize', 28); %a
% title('Control con Observador (\omega y \theta)', 'Fontsize', 28); %b
title('Control con no linealidad en la acción de control', 'Fontsize', 28); %c
subplot(4,1,2), plot(t, wr, 'k', 'LineWidth', 3); 
ylabel('\omega_r(t)', 'Fontsize', 28),grid on;
legend('Velocidad del Motor');
subplot(4,1,3), plot(t, theta, 'k', 'LineWidth', 3);
ylabel('\theta(t)', 'Fontsize', 28),hold on;
plot(t, theta_ref, 'r--', 'LineWidth', 3),grid on;
legend('Angulo del motor', 'Angulo de Referencia');
subplot(4,1,4), plot(t(2:end), Va(2:end), 'k', 'LineWidth', 3); 
xlabel('t [s]', 'Fontsize', 28), ylabel('V_a(t)  T_L', 'Fontsize', 28),hold on;
plot(t(2:end), TL(2:end), 'r', 'LineWidth', 3),grid on;
legend('Accion de Control', 'Torque');