function constitutiveGauss(a, b, f, w, m, N)
fprintf('Составная Квадратурная формула Гаусса с %d узлами с чилом промежутков деления [%f, %f] равным %d: \n', N, a, b, m);
syms t;
phi = w*f;
h = (b - a)/m;
Z = a : h : b;
J = 0;
P = t.^(0 : 1);
for i = 3 : N + 1
    P(1, i) = (2*(i - 1) - 1)/(i - 1) * t * P(1, i - 1) - ((i - 1) - 1)/(i - 1) * P(1, i - 2);
end
fprintf('Многочлен Лежандра степени %d:\n', N);
disp(expand(P(1, N + 1)));
L = P(1, N + 1);
disp('Корни многочлена Лежандра:');
T = solve(L);
disp(T);
C = 2./((1 - T.*T).*(subs(diff(L, t), t, T)).^2);
A = h/2 * C;
disp('Коэффициенты КФ:');
disp(C);
for i = 1 : m 
    X = (Z(1, i + 1) - Z(1, i))/2 .* T + (Z(1, i + 1) + Z(1, i))/2;
    Y = phi(X);
    I = A'*Y;
    J = J + I;
end
fprintf('Приближенное значение интеграла %6.8f\n', J);
end 

