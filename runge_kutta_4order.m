%% Dados do problema 
f_xy = @(r,t) t;  % função f(x,y)
h = 0.1; % passo
% intervalo 
inicio = 0.0;
fim = 1.0;
%% Início dos cálculos
x = inicio:h:fim;
y = zeros(1,length(x)); 
% condição inicial
y(1) = 1.0;

for i=1:(length(x)-1)  % precisa do -1 para que o ylen seja igual a xlen       
    % para 2 funções 
    k_1 = f_xy(x(i),y(i));
    k_2 = f_xy(x(i)+0.5*h,y(i)+0.5*h*k_1);
    k_3 = f_xy((x(i)+0.5*h),(y(i)+0.5*h*k_2));
    k_4 = f_xy((x(i)+h),(y(i)+k_3*h));
    y(i+1) = y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  
  
end

[w,s] = ode45(f_xy,[inicio:h:fim],1); %ode45 é o próprio runge-kutta de 4 ordem 5 passos 

solucao_edo = exp(x); %solução exata

%% Plot dos yRunge-Kutta e da tabela de valores

table_data = [x', y'];  
%table_data = num2cell(table_data);

figure;
uitable('Data',table_data,'ColumnName',{'x','y'},'Position',[140 100 400 300])
%title('Valores de x e y')

figure;
plot(x,y,'b',x,solucao_edo,'or')
title('Aproximações para Runge-Kutta')
xlabel('x')
ylabel('f(x,y)')





