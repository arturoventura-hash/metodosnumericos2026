% =========================================================================
% SCRIPT DE INTERPOLACIÓN: MATRICIAL, NEWTON, LAGRANGE Y SPLINES
% =========================================================================
clear; clc;

% 1. ENTRADA DE DATOS (Extraídos de la imagen)
x = [3, 4, 5, 6, 7, 8, 10, 11, 14, 18, 20, 22, 23, 25, 28, 29, 31, 32, 34];
y = [15, 11, 18, 16, 14, 12, 19, 18, 23, 17, 18, 21, 24, 23, 15, 2, 4, 8, 13];

% Rango de evaluación para las gráficas (más puntos para suavidad)
xi = linspace(min(x), max(x), 1000);
n = length(x);

% -------------------------------------------------------------------------
% A. MÉTODO MATRICIAL (Vandermonde)
% -------------------------------------------------------------------------
% El polinomio es y = a1*x^(n-1) + a2*x^(n-2) + ... + an
V = zeros(n, n);
for i = 1:n
    V(:,i) = x'.^(n-i);
end
a_mat = V \ y'; % Resolver el sistema de ecuaciones lineales
y_mat = polyval(a_mat, xi);

% -------------------------------------------------------------------------
% B. MÉTODO DE NEWTON (Diferencias Divididas)
% -------------------------------------------------------------------------
D = zeros(n, n);
D(:,1) = y';
for j = 2:n
    for i = j:n
        D(i,j) = (D(i,j-1) - D(i-1,j-1)) / (x(i) - x(i-j+1));
    end
end

% Evaluación del polinomio de Newton en el rango xi
y_newton = D(1,1);
for j = 2:n
    prod_term = 1;
    for k = 1:j-1
        prod_term = prod_term .* (xi - x(k));
    end
    y_newton = y_newton + D(j,j) * prod_term;
end

% -------------------------------------------------------------------------
% C. MÉTODO DE LAGRANGE
% -------------------------------------------------------------------------
y_lagrange = zeros(size(xi));
for i = 1:n
    Li = ones(size(xi));
    for j = 1:n
        if i ~= j
            Li = Li .* (xi - x(j)) / (x(i) - x(j));
        end
    end
    y_lagrange = y_lagrange + y(i) * Li;
end

% -------------------------------------------------------------------------
% D. MÉTODO SPLINES CÚBICOS
% -------------------------------------------------------------------------
% Octave utiliza splines cúbicos naturales por defecto con la función spline
y_spline = spline(x, y, xi);

% =========================================================================
% VISUALIZACIÓN DE RESULTADOS
% =========================================================================

figure('Name', 'Comparativa de Metodos Numéricos', 'NumberTitle', 'off');

% Subplot 1: Polinomios Globales (Todos coinciden matemáticamente)
subplot(2,1,1);
plot(x, y, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Puntos (Datos)'); hold on;
plot(xi, y_mat, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Matricial');
plot(xi, y_newton, 'b--', 'LineWidth', 1, 'DisplayName', 'Newton');
plot(xi, y_lagrange, 'g:', 'LineWidth', 1, 'DisplayName', 'Lagrange');
title('Interpolación Polinomial Global (Grado 18)');
legend('show'); grid on;
ylabel('y');

% Subplot 2: Splines Cúbicos (El más estable)
subplot(2,1,2);
plot(x, y, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Puntos (Datos)'); hold on;
plot(xi, y_spline, 'm-', 'LineWidth', 2, 'DisplayName', 'Spline Cúbico');
title('Interpolación por Splines Cúbicos (Trazadores)');
legend('show'); grid on;
xlabel('x'); ylabel('y');

fprintf('--- Informe de Ejecución ---\n');
fprintf('Número de puntos procesados: %d\n', n);
fprintf('Grado del polinomio global: %d\n', n-1);
fprintf('Nota: Los métodos Matricial, Newton y Lagrange generan el mismo polinomio.\n');
