% =========================================================================
% INTERPOLACIÓN COMPLETA: POLINOMIO GLOBAL Y SPLINES
% =========================================================================
clear; clc;

% 1. COORDENADAS EXTRAÍDAS (x, y)
x = [3, 4, 5, 6, 7, 8, 10, 11, 14, 18, 20, 22, 23, 25, 28, 29, 31, 32, 34];
y = [15, 11, 18, 16, 14, 12, 19, 18, 23, 17, 18, 21, 24, 23, 15, 2, 4, 8, 13];

% Rango de evaluación suave (1000 puntos para que la curva no se vea "picuda")
xi = linspace(min(x), max(x), 1000);
n = length(x);

% -------------------------------------------------------------------------
% MÉTODOS PARA EL POLINOMIO GLOBAL (Grado 18)
% -------------------------------------------------------------------------

% A. MÉTODO MATRICIAL (Para obtener los coeficientes reales)
% Construimos la matriz de Vandermonde manual
V = zeros(n, n);
for i = 1:n
    V(:,i) = x'.^(n-i);
end
coeffs = V \ y'; % Estos son los coeficientes del polinomio

% B. MÉTODO DE LAGRANGE / NEWTON (Evaluación directa)
% Usamos polyval con los coeficientes obtenidos para representar el polinomio
y_polinomio = polyval(coeffs, xi);

% -------------------------------------------------------------------------
% MÉTODO DE SPLINES CÚBICOS (Trazadores)
% -------------------------------------------------------------------------
y_spline = spline(x, y, xi);

% =========================================================================
% GENERACIÓN DE LA GRÁFICA FINAL
% =========================================================================
figure('Color', [1 1 1]);

% Dibujamos el Polinomio Global (Lagrange/Newton/Matricial son el mismo)
plot(xi, y_polinomio, 'b-', 'LineWidth', 2); hold on;

% Dibujamos el Spline para comparar la suavidad
plot(xi, y_spline, 'g--', 'LineWidth', 1.5);

% Dibujamos los puntos originales de tu libreta
plot(x, y, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

% Formateo de la gráfica
title(['Interpolación Polinomial de Grado ', num2str(n-1)]);
xlabel('Eje X (Unidades de cuadrícula)');
ylabel('Eje Y (Unidades de cuadrícula)');
legend('Polinomio Global (Newton/Lagrange)', 'Splines Cúbicos', 'Puntos Originales');
grid on;

% Ajustar límites para ver las oscilaciones del polinomio
axis([min(x)-1 max(x)+1 min(y)-10 max(y)+10]);

% Mostrar los coeficientes en la consola (opcional)
fprintf('Coeficientes del polinomio (de mayor a menor grado):\n');
disp(coeffs);
