% Очищаем рабочее пространство
clear all;
clc;

% Задаем диапазон значений x и t
x = linspace(0, 1, 50);
t = linspace(0, 1, 50);

% Создаем сетку для x и t
[X, T] = meshgrid(x, t);

% Вычисляем значения y = e^t * sinh(x)
Y = exp(T) .* sinh(X);

% Строим 3D график
figure;
surf(X, T, Y);
title('График аналитического решения y = e^t \cdot sinh(x)', 'Fontsize', 25);
xlabel('x');
ylabel('t');
zlabel('y');
colormap('summer');
colorbar; % Добавляем цветовую шкалу
grid on;

