clear all;clc

%% Parâmetros 
a1 = 1;      a2 = 1;
l1 = 0.5;    l2 = 0.5;
m_l1 = 50;   m_l2 = 50;
I_l1 = 10;   I_l2 = 10;
k_r1 = 100;  k_r2 = 100; 
m_m1 = 5;    m_m2 = 5;
I_m1 = 0.01; I_m2 = 0.01;
g = 9.81;

%% Condição inicial

x_0 = 0.2;
y_0 = 0;

[t1,t2] = calcAngulos(x_0,y_0,a1,a2);

y0 = [t1; t2; 0.1; 0.1];
gravidade = [100.0; 0.10];

h = 0.005;
t0 = 0.0;
tf = 0.5;
t = t0:h:tf;

n = length(t); 
Y = y0;
y = y0;

% tau1 = zeros(1,101);
% tau2 = zeros(1,101);

theta_ponto_ponto = [25;25];

% 2 coisas: passo do theta e passo do tempo(ja defini --> fazer por int)
% supor de 0 a 1 grau (theta) em meio ao 0 a 1, integro o sistema de 0 a 1 
% no proximo passo ir de 1 a 2 --> integro o sistema de 1 a 2
% outra saida: descobrir o tempo necessário para ir até pi/2 --> aí não
% preciso de 2 for

for i = 1:n-1

    %provavelmente colocar um for aqui dentro 
    %um para a mudança do ângulo(por fora) e outro para integração(por
    %dentro)

    b11 = I_l1 + m_l1*(l1^2) + (k_r1^2)*I_m1 + I_l2 + m_l2*(a1^2 + l2^2 + 2*a1*l2*cos(y(2,i))) + I_m2 + m_m2*(a1^2);
    b12 = I_l2 + m_l2*(l2^2 + a1*l2*cos(y(2,i))) + k_r2*I_m2;
    b21 = b12;
    b22 = I_l2 + m_l2*(l2^2) + (k_r2^2)*I_m2;

    h_c = - m_l2*a1*l2*sin(y(2,i));

    c11 = h_c*y(4,i);
    c12 = h_c*(y(3,i) + y(4,i));
    c21 = -h_c*y(3,i);
    c22 = 0;

    g1 = (m_l1*l1 + m_m2*a1 + m_l2*a1)*g*cos(y(1,i)) + m_l2*l2*g*cos(y(1,i)+y(2,i));
    g2 = m_l2*l2*g*cos(y(1,i)+y(2,i));
    
    Mi = [b11, b12; b21, b22];
    C = [c11, c12; c21, c22];

    g_vet = [g1;g2];
    gravidade(1:2,i+1) = g_vet;

    u1 = (I_l1 + m_l1*(l1^2) + (k_r1^2)*I_m1 + I_l2 + m_l2*(a1^2 + l2^2 + 2*a1*l2*cos(y(2,i))) +I_m2 + m_m2*(a1^2))*theta_ponto_ponto(1) + (I_l2 + m_l2*(l2^2 + a1*l2*cos(y(2,i))) + k_r2*I_m2)*theta_ponto_ponto(2)... 
        -2*m_l2*a1*l2*sin(y(2,i))*y(3,i)*y(4,i) - m_l2*a1*l2*sin(y(2,i))*y(4,i)^2 + (m_l1*l1 + m_m2*a1 + m_l2*a1)*g*cos(y(1,i)) + m_l2*l2*g*cos(y(1,i)+y(2,i));

    u2 = (I_l2 + m_l2*(l2^2 + a1*l2*cos(y(2,i))) + k_r2*I_m2)*theta_ponto_ponto(1) + (I_l2 + m_l2*l2^2 + k_r2^2*I_m2)*theta_ponto_ponto(2) + m_l2*a1*l2*sin(y(2,i))*y(3,i)^2 ...
        +m_l2*l2*g*cos(y(1,i)+y(2,i));

    tau1(i) = u1;
    tau2(i) = u2;

    u = [u1;u2];

    theta_ponto_ponto = inv(Mi)*(u - C*y(3:4,i) - g_vet);  
    aceleracao(:,i) = theta_ponto_ponto;
    %aqui? porque não funciona 

    A = [zeros(2,2), eye(2,2); inv(Mi)*0, -inv(Mi)*C];
    B = [zeros(2,2);inv(Mi)];

    Y = rungekutta(t(i),Y,h,A,B,g_vet,u);
    y(1:4,i+1) = Y;

    % if y(1,i+1) >= pi/2 || y(2,i+1) >= pi/2
    %     break;  % Sai do loop se o ângulo atingiu pi/2
    % end

    t(i+1) = t(i) + h;

    %sai da equação da dinâmica --> B*q_dot_dot + C*q_dot + g = tau
    % theta1_ponto_ponto = inv(Mi)*(u - C*y(3,i) - g_vet);
    % theta2_ponto_ponto = inv(Mi)*(u - C*y(4,i) - g_vet);
    
end
% 
 %% Plot dos resultados

% Cinemática direta

theta1_f = t1 + pi/2;
theta2_f = t2 + pi/2;

theta1 = linspace(t1, theta1_f, length(t));
theta2 = linspace(t2, theta2_f, length(t));

posx = a1*cos(theta1) + a2*cos(theta1 + theta2);
posy = a1*sin(theta1) + a2*sin(theta1 + theta2);

figure(1);

for i = 1:length(t)
    plot([0, a1*cos(theta1(i))], [0, a1*sin(theta1(i))], 'r', 'LineWidth', 2); % Primeiro elo
    hold on;
    plot([a1*cos(theta1(i)), posx(i)], [a1*sin(theta1(i)), posy(i)], 'b', 'LineWidth', 2); % Segundo elo
    plot(posx(1:i), posy(1:i), 'g', 'LineWidth', 1); % Trajetória
    plot(x_0, y_0, 'rs', 'markersize', 10); % Posição inicial
    axis equal;
    axis([-(a1+a2)-0.1 (a1+a2)+0.1 -(a1+a2)-0.1 (a1+a2)+0.1]);
    title('Movimento do manipulador');
    pause(0.05);
    hold off;
end

 figure(2);
 plot(t,y(1,:),'k',t,y(2,:),'g',t,y(3,:),'r',t,y(4,:),'b')
 title('Aproximações para Runge-Kutta')
 legend('y_1','y_2','y_3','y_4')
 xlabel('Tempo [s]')
 ylabel('Valor')
 
 figure(3);
 plot(t,gravidade(1,:),'k',t,gravidade(2,:),'b')
 title('Componente Gravitacional')
 xlabel('[s]')
 ylabel('[Nm]')

 figure(4);
 plot(t,y(1,:),'k',t,y(2,:),'b')
 title('Posição das juntas')
 legend('junta_1','junta_2')
 xlabel('[s]')
 ylabel('[rad]')

 figure(5);
 plot(t,y(3,:),'r',t,y(4,:),'g')
 title('Velocidade das juntas')
 legend('junta_1','junta_2')
 xlabel('[s]')
 ylabel('[rad/s]')

 figure(6);
 plot(t,,'r',t,theta_ponto_ponto(2,:),'g')
 title('Aceleração das juntas')
 legend('junta_1','junta_2')
 xlabel('[s]')
 ylabel('[rad/s^2]')

 figure(7);
 plot(t,tau1,'r',t,tau2,'g')
 title('Torque das juntas')
 legend('junta_1','junta_2')
 xlabel('[s]')
 ylabel('[Nm]')
