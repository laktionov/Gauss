function pointsGauss(a, b, f, w, N)
fprintf('Квадратурная формула Гаусса с %d узлами: \n', N);
syms x;
M = zeros(1, 2*N);
B = zeros(N);
d = zeros(N, 1);
for i = 1 : 2*N
    M(1, i) = int(w*x^(i-1), a, b);
end
disp('Моменты весовой функции: ');
disp(M);
for i = 1 : N
    B(:, i) = M(1, N - i + 1 : 2*N - i)';
end
d = -1 * M(1, N + 1 : 2*N)';
K = B \ d;
disp('Ортогональный многочлен: ');
X = x.^(N : -1 : 0);
omega = X * [1; K];
disp(omega);
S = roots([1; K]);

V = zeros(N);
for i = 1 : N
    V(:, i) = S(i).^(0 : N - 1); 
end
disp('Узлы :');
disp(S');
A = V \ M(1, 1 : N)';
disp('Коэффициенты КФ:');
disp(A);
Y = f(S);
J = A' * Y;
fprintf('Приближенное значение интеграла %6.8f\n', J);
fprintf('\n')
end

