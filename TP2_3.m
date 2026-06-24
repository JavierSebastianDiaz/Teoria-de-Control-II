clc;clear all, close all;
Ts=0.1;KMAX=1000;T=Ts*KMAX;
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
TamanioFuente=20;
%pkg load signal; pkg load control; %Solo una vez
%Condiciones iniciales
alfa(1)=.5; colorc='.g';colorl='g';
% alfa(1)=.1; colorc='.r';colorl='r';
% alfa(1)=.8; colorc='.b';colorl='b';
alfa(1)=2.8;Ci=2.8; colorc='.k';colorl='k';
%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
% estado=[p(i); p_p(i); alfa(i); omega(i)]
% Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0];
% Mat_Bc=[0; 1/M; 0; -1/(long*M)];
Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(long*M) -g*(m+M)/(long*M) 0];
Mat_Bc=[0; 1/M; 0; 1/(long*M)];
Mat_C=[1 0 0 0; 0 0 1 0];
sys_c=ss(Mat_Ac,Mat_Bc,Mat_C,[0]);
sys_d=c2d(sys_c,Ts,'zoh');
Mat_A=sys_d.a;
Mat_B=sys_d.b;
H=[0;0;0;0];d_tao=Ts/100;tao=0;
for hh=1:100
    dH=expm(Mat_Ac*tao)*Mat_Bc*d_tao;
    H=H+dH;
    tao=tao+d_tao;
end
Mat_B=H;
Mat_A=expm(Mat_Ac*Ts);
Mat_A1= Mat_A;
Mat_B1=Mat_B;
m=1;
Mat_Ac=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(long*M) -g*(m+M)/(long*M) 0];
Mat_Bc=[0; 1/M; 0; 1/(long*M)];
sys_c2=ss(Mat_Ac,Mat_Bc,Mat_C,[0]);
sys_d2=c2d(sys_c,Ts,'zoh');
Mat_A2=sys_d2.a;
Mat_B2=sys_d2.b;

Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B];%Matriz Controlabilidad
rango=rank(Mat_M);
%Cálculo del controlador por asignación de polos
auto_val=eig(Mat_A);
c_ai=poly(auto_val);
Mat_W=[c_ai(4) c_ai(3) c_ai(2) 1 ;c_ai(3) c_ai(2) 1 0 ;c_ai(2) 1 0 0 ;1 0 0 0 ];
Mat_T=Mat_M*Mat_W;
A_controlable=inv(Mat_T)*Mat_A*Mat_T %Verificación de que T esté bien
%Ubicación de los polos de lazo cerrado en mui:
mui(1)=0.9997; mui(2)=0.994568; mui(3)= conj(mui(2)); mui(4)=0.99; %Para Ts=0.001
% mui=mui*.98;
mui(1)=0.8590123457; mui(2)=0.978; mui(3)= conj(mui(2)); mui(4)=0.9019;%Para Ts=0.01
% mui=mui*0;
% mui=mui*.9998;
mui(1)=0.67987; mui(2)=0.6798; mui(3)= conj(mui(2)); mui(4)=0.79;%Para Ts=0.1
% mui=mui*0.999;
alfa_i=poly(mui);
Mat_Q=diag([1,1/9^2,1/.5^2,1/.5^2]);
Mat_R=1e2;
Mat_Qo=diag([1,1/9^2,1/.5^2,1/.5^2]);
Mat_Ro=diag([1e1,1e0]);
% estado=[p(i); p_p(i); alfa(i); omega(i)]
K=fliplr(alfa_i(2:5)-c_ai(2:5))*inv(Mat_T);
K=dlqr(Mat_A, Mat_B, diag([1,1,1e1,1e1]),1e0);
Ko1=dlqr(Mat_A', Mat_C', diag([1,1,1,1/.5^2]), diag([1e-2,1e-2]))'; %calculo de observador
Ko2=dlqr(Mat_A2', Mat_C', diag([1,1,1e2,1e2]), diag([1e-5,1e-3]))';
K1=K;
eig(Mat_A-Mat_B*K)
K2=dlqr(Mat_A2, Mat_B2, diag([1,1,1,1]), Mat_R*1e-1);
G=inv(Mat_C(1,:)*inv(eye(4)-Mat_A+Mat_B*K)*Mat_B);
G1=G;
G2=inv(Mat_C(1,:)*inv(eye(4)-Mat_A2+Mat_B2*K2)*Mat_B2);
abs(eig(Mat_A-Mat_B*K))
t=0; x=[0;0;alfa(1);0];
p(1)=x(1); p_p(1)=x(2); alfa(1)=x(3); omega(1)=x(4);
%Implementación del controlador en el modelo no lineal en tiempo continuo.
x=[0;0;alfa(1);0];
pl=x(1); p_pl=x(2); alfal=x(3); omegal=x(4);
p=x(1); p_p=x(2); alfa=x(3); 
omega=0; tita_pp(1)=0;h=Ts/20; u=[];i=1;
u_k(1)=0;
x_ang=[0; 0; 0; 0];
Ko=Ko1;
Xop=[0 0 pi 0]';%x=[0 0 alfa(1) 0]';
ref=10;
for ki=1:KMAX
    ul=-K*(x_ang-Xop)+G*ref;
    Y=Mat_C*(x-Xop);
%     x=Mat_A*x+Mat_B*ul;
    x=Mat_A*(x-Xop)+Mat_B*ul+Xop;
    y_ang=Mat_C*(x_ang-Xop);
    x_ang=Mat_A*(x_ang-Xop)+Mat_B*ul+Xop+Ko*(Y-y_ang);
    pl(ki)=x(1);
    p_pl(ki)=x(2);
    alfal(ki)=x(3);
    omegal(ki)=x(4);
    u_kl(ki)=ul;
    if x(1)>9.9 && Ts*ki>3
        ref=0;
        Mat_A=Mat_A2;
        Mat_B=Mat_B2;
        m=1;
        Ko=Ko2;
        G=G2;
        K=K2;
    end
end
% pause
K=K1;
G=G1;
Mat_A=Mat_A1;
Mat_B=Mat_B1;
tl=(0:KMAX-1)*Ts;
x=[0;0;Ci;0];
ref = 10;
x_ang=[0; 0; 0; 0];
Ko=Ko1;
m=0.1;
i=1;
flag=0;
for ki=1:KMAX
    u1(ki)=-K*(x_ang-Xop)+G*ref;
    Y=Mat_C*(x-Xop);
    for kii=1:Ts/h
        u(i)=u1(ki);
        p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))- Fricc*p_p(i));
        tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
        p_p(i+1)=p_p(i)+h*p_pp;
        p(i+1)=p(i)+h*p_p(i);
        omega(i+1)=omega(i)+h*tita_pp;
        alfa(i+1)=alfa(i)+h*omega(i);
        x=[p(i+1); p_p(i+1); alfa(i+1); omega(i+1)];
        i=i+1;
    end
    y_ang=Mat_C*(x_ang-Xop);
    x_ang=Mat_A*(x_ang-Xop)+Mat_B*u1(ki)+Xop+Ko*(Y-y_ang);
    if x(1)>9.9 && flag==0
        flag=1;
        ref=0;
        x_ang=[0;0;0;0];
        m=1;
        Mat_A=Mat_A2;
        Mat_B=Mat_B2;
%         Ko=Ko2;
%         G=G2;
%         K=K2;
    end
end
u(i)=u1(ki);t=0:h:T;
figure(1);
subplot(3,2,1);plot(tl,omegal,colorl,t,omega,colorc);grid on;
title('Velocidad ángulo','FontSize',TamanioFuente);hold on;
subplot(3,2,2);plot(tl,alfal,colorl,t,alfa,colorc);grid on;
title('Ángulo','FontSize',TamanioFuente);hold on;
legend('Continuo', 'Discreto');
subplot(3,2,3);plot(tl,pl,colorl,t,p,colorc);grid on;title('Posición carro','FontSize',TamanioFuente);hold on;
subplot(3,2,4);plot(tl,p_pl,colorl,t,p_p,colorc);grid on;title('Velocidad carro','FontSize',TamanioFuente);hold on;
subplot(3,1,3);plot(tl,u_kl,colorl,t,u,colorc);grid on;title('Acción de control','FontSize',TamanioFuente);xlabel('Tiempo en Seg.','FontSize',TamanioFuente);hold on;
figure(2);
subplot(2,1,1);plot(alfal,omegal,colorl,alfa,omega,colorc);grid on;xlabel('Ángulo','FontSize',TamanioFuente);ylabel('Velocidad angular','FontSize',TamanioFuente);hold on;
legend('Continuo', 'Discreto');
subplot(2,1,2);plot(pl,p_pl,colorl,p,p_p,colorc);grid on;xlabel('Posición carro','FontSize',TamanioFuente);ylabel('Velocidad carro','FontSize',TamanioFuente);hold on;
