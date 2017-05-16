function pointsGauss(a, b, f, w, N)
fprintf('������������ ������� ������ � %d ������: \n', N);
syms x;
M = zeros(1, 2*N);
B = zeros(N);
d = zeros(N, 1);
for i = 1 : 2*N
    M(1, i) = int(w*x^(i-1), a, b);
end
disp('������� ������� �������: ');
disp(M);
for i = 1 : N
    B(:, i) = M(1, N - i + 1 : 2*N - i)';
end
d = -1 * M(1, N + 1 : 2*N)';
K = B \ d;
disp('������������� ���������: ');
X = x.^(N : -1 : 0);
omega = X * [1; K];
disp(omega);
S = roots([1; K]);

V = zeros(N);
for i = 1 : N
    V(:, i) = S(i).^(0 : N - 1); 
end
disp('���� :');
disp(S');
A = V \ M(1, 1 : N)';
disp('������������ ��:');
disp(A);
Y = f(S);
J = A' * Y;
fprintf('������������ �������� ��������� %6.8f\n', J);
fprintf('\n')
end

