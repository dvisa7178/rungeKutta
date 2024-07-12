%% Dados do problema
inicio = 0.0;
fim = 1.0;
n = 10;
h = (fim-inicio)/n; % passo
%% Início dos cálculos
x = inicio:h:fim;
y = zeros(1,length(x)); 
% valor inicial
y0 = 1;
y(1,1) = y0;
Y = y0;

for i=1:n  % precisa do -1 para que o ylen seja igual a xlen   
    Y = rungekutta(x(i),Y,h);
    x(i+1) = x(i)+h;
    y(1:2,i+1) = Y;
end

%% Plot dos yRunge-Kutta e da tabela de valores

table_data = [x', y'];  
%table_data = num2cell(table_data);

figure;
uitable('Data',table_data,'ColumnName',{'x','y'},'Position',[140 100 400 300])
%title('Valores de x e y')

figure;
plot(x,y,'b')
%plot(x,y,'b',x,solucao_edo,'or')
title('Aproximações para Runge-Kutta')
xlabel('x')
ylabel('f(x,y)')





