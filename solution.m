%% Dados do problema
clear all; clc;

h = 0.001;
t0 = 0.0;
tf = 0.5;

t = t0:h:tf;
n = length(t); 



%% Parâmetros 
a1 = 1;      a2 = 1;
l1 = 0.5;    l2 = 0.5;
m_l1 = 50;   m_l2 = 50;
I_l1 = 10;   I_l2 = 10;
k_r1 = 100;  k_r2 = 100; 
m_m1 = 5;    m_m2 = 5;
I_m1 = 0.01; I_m2 = 0.01;
g = 9.81;


%% Início dos cálculos

%olhar no outro código como calcular a posição (x,y) do ef

% Condição inicial
u = [0.0;1.0];
y0 = [0.1; 0.1; 0.1; 0.1];
gravidade = [100.0;0.0];
%posição inicial é de (0.2,0)
t = t0:h:tf;
n = length(t); 
Y = y0;
y = y0;

for i = 1:n-1

    % transformando thetas em y's
    theta1 = y(1, i);
    theta2 = y(2, i);
    % aqui deve ser calulada a posição (da mesma forma em que foi feito no
    % código anterior)
    %plotar também o comportamento do robô


    theta1_ponto = y(3, i);
    theta2_ponto = y(4, i);

    b11 = I_l1 + m_l1*(l1^2) + (k_r1^2)*I_m1 + I_l2 + m_l2*(a1^2 + l2^2 + 2*a1*l2*cos(theta2)) + I_m2 + m_m2*(a1^2);
    b12 = I_l2 + m_l2*(l2^2 + a1*l2*cos(theta2)) + k_r2*I_m2;
    b21 = b12;
    b22 = I_l2 + m_l2*(l2^2) + (k_r2^2)*I_m2;

    h_c = - m_l2*a1*l2*sin(theta2);

    c11 = h_c*theta2_ponto;
    c12 = h_c*(theta1_ponto + theta2_ponto);
    c21 = -h_c*theta1_ponto;
    c22 = 0;

    g1 = (m_l1*l1 + m_m2*a1 + m_l2*a1)*g*cos(theta1) + m_l2*l2*g*cos(theta1+theta2);
    g2 = m_l2*l2*g*cos(theta1+theta2);
    
    
    Mi = [b11, b12; b21, b22];
    C = [c11, c12; c21, c22];
    g_vet = [g1;g2];
    gravidade(1:2,i+1) = g_vet;
    
 
    A = [zeros(2,2), eye(2,2); inv(Mi)*0, -inv(Mi)*C];
    B = [zeros(2,2);inv(Mi)];

    Y = rungekutta(t(i),Y,h,A,B,g_vet,u);
    % x vira t arrumar depois
    t(i+1) = t(i) + h;
    y(1:4,i+1) = Y;

end

%% Plot dos yRunge-Kutta e da tabela de valores

%table_data = [ y'];  
%table_data = num2cell(table_data);

%figure;
%uitable('Data',table_data,'ColumnName',{'y'},'Position',[140 100 400 300])
%title('Valores de y')

figure;
plot(t,y(1,:),'k',t,y(2,:),'g',t,y(3,:),'r',t,y(4,:),'b')
title('Aproximações para Runge-Kutta')
legend('y_1','y_2','y_3','y_4')
xlabel('x')
ylabel('f(x,y)')

figure;
plot(t,gravidade(1,:),'k',t,gravidade(2,:),'b')
