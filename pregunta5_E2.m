
function pregunta5_E2()
  clc; clear;

  mvalues = [16; 32; 64; 128; 256];

  display('Metodo 1: HSS');
  for i=1:5
    m = mvalues(i);
    h = 1 / (m + 1);
    [W,T] = calc_W_T(m,h);
    [f,g] = calc_f_g(m,h);

    disp(["Caso", num2str(i), ": m=", num2str(m)]);
    pregunta1(W, T, f, g);
  endfor

end

function [W,T]=calc_W_T(m,h)
  A = tridiag(-1,2,-1,m);
  Vm = (h^-2)*A;
  I1 = eye(m);
  I2 = eye(m^2);
  cita1 = -10;
  cita2 = 1;
  K =  kron(I,Vm) + kron(Vm,I);
  W = K+(cita1*I);
  T = cita2*I;
end

function [f,g]=calc_f_g(m,h)
  A = tridiag(-1,2,-1,m);
  Im = eye(m,1);
  f = A*Im;
  g = A*Im;

  end

function KronAB=kron(A,B)
  [m, n] = size(A);
  [p, q] = size(B);
  KronAB = zeros(m * p, n * q);
  for i = 1:m
    for j = 1:n
      KronAB((i - 1) * p + 1:i * p, (j - 1) * q + 1:j * q) = A(i, j) * B;
    end
  end
end

function A = tridiag(a, b, c, n)
  A = b * eye(n);
  A = A + a * diag(ones(n-1, 1), -1);
  A = A + c * diag(ones(n-1, 1), 1);
end
