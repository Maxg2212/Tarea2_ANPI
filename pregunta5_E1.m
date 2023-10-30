% La funcion pregunta5_E1 ejecuta todos los metodos para distintos valores de m.
% Sintaxis de la funcion: pregunta5_E1()
function pregunta5_E1()
  clc; clear;

  mvalues=[16;32;64;128;256];

  display('Metodo 1: HSS');
  for i=1:2
    m=mvalues(i);
    h=1/(m+1);
    [f,g]=calc_f_g(m,h);
    [W,T]=calc_W_T(m,h);

    disp(["Caso", num2str(i), ": m=", num2str(m)]);
    pregunta1(W, T, f, g);
  endfor


  display('Método 3: MHSS');
  for i = 1 : 2
    m = mvalues(i);
    h = 1/(m+1);
    [f,g]=calc_f_g(m,h);
    [W,T]=calc_W_T(m,h);

    iterMax = 5000;
    tol = 1e-12;

    disp(["Caso", num2str(i), ": m = ", num2str(m)]);
    metodo_3(W,T,f,g,iterMax,tol)

  endfor

  display('Metodos 4: QR y Eliminacion Gausseana');
  for i=1:2
    m=mvalues(i);
    h = 1/(m+1);
    [f,g]=calc_f_g(m,h);
    [W,T]=calc_W_T(m,h);

    disp(["Caso", num2str(i), ": m=", num2str(m)]);
    pregunta4(W,T,f,g);
  endfor


end

% La funcion calc_W_T calcula los valores de las matrices W y T.
% Sintaxis de la funcion: [W,T]=calc_W_T(m,h)
% Parametros de entrada:
%         m = Tamaño de matriz
%         h = Constante obtenido de forma previa
% Parametros de salida:
%         W = Matriz de medida m x m
%         T = Matriz de medida m x m
function [W,T]=calc_W_T(m,h)
  bm=(1/((h)^2))*tridiag(-1,2,-1,m);
  I1=eye(m);
  I2=eye(m^2);
  tau=h;
  K=(kron(I,bm))+(kron(bm,I));
  W=(K+((3-sqrt(3))/tau))*I;
  T=(K+((3+sqrt(3))/tau))*I;
  %A = tridiag(-1, 2, -1, m);
  %Vm = (1/h^2)*A;
  %cita1 = -10;
  %cita2 = 1;
  %K = kron(I,Vm) + kron(Vm,I);
  %W = K + (cita1*I);
  %T = cita2 * I;

end

% La funcion calc_f_g calcula los valores de las matrices f y g.
% Sintaxis de la funcion: [f,g]=calc_f_g(m,h)
% Parametros de entrada:
%         m = Tamaño de matriz
%         h = Constante obtenido de forma previa
% Parametros de salida:
%         f = Matriz de medida m x 1
%         g = Matriz de medida m x 1
function [f,g]=calc_f_g(m,h)
  f=zeros(m,1);
  g=zeros(m,1);
  for j=1:m
    fj=(1-i)*j/(h*(j+1)^2);
    f(j)=real(fj);
    g(j)=imag(fj);
  end
end

% La funcion kron realiza el calculo del producto de Kronecker.
% Sintaxis de la funcion: KronAB=kron(A,B)
% Parametros de entrada:
%         A = Matriz de medida m x n
%         B = Matriz de medida p x q
% Parametros de salida:
%         KronAB = Matriz de medida m^2 x m^2
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

% La funcion tridiag realiza una matriz tridiagonal con lo valores mostrados.
% Sintaxis de la funcion: A = tridiag(a, b, c, n)
% Parametros de entrada:
%         a = Primer valor de la tridiagonal de la matriz
%         b = Segundo valor de la tridiagonal de la matriz
%         c = Tercer valor de la tridiagonal de la matriz
%         n = Tamaño de la matriz
% Parametros de salida:
%         A = Matriz tridiagonal producida
function A = tridiag(a, b, c, n)
  A = b * eye(n);
  A = A + a * diag(ones(n-1, 1), -1);
  A = A + c * diag(ones(n-1, 1), 1);
end

function metodo_3(W,T,p,q,iterMax,tol)
  tic

  # Variable compleja
  i = sqrt(-1);
  # Se definen las matrices iniciales

  # matrices complejas
  A = W + i*T;
  b = p + i*q;
  m = length (W);

  # valor inicial x0
  x_k = zeros(m, 1);

  # Determina el valor propio más pequeño y más
  # grande de la matriz
  lambda_min = min(eig(W));
  lambda_max = max(eig(W));

  # Calcula el valor de alpha*
  alpha = sqrt(lambda_min * lambda_max);

  # Calculo de M_alpha
  # Matrices intermedias para el calculo
  Im = eye(size(W));
  alpha_Im = alpha*Im;

  alpha_Im_T = alpha_Im + T;
  alpha_Im_iW = alpha_Im + i*W;
  alpha_Im_W = alpha_Im + W;
  alpha_Im_iT = alpha_Im - i*T;

  # Calculo de matrices inversas para el calculo de M_alpha
  alpha_Im_T_Inv = mldivide(alpha_Im_T, eye(size(alpha_Im_T)));
  alpha_Im_W_Inv = mldivide(alpha_Im_W, eye(size(alpha_Im_W)));

  M_alpha = (alpha_Im_T_Inv)*(alpha_Im_iW)*(alpha_Im_W_Inv)*(alpha_Im_iT);

  # Calculo de N_alpha
  N_alpha = (1 - i)*alpha*(alpha_Im_T_Inv)*(alpha_Im_W_Inv);

  # Cálculo de la sucesión desde k = 0 hasta iterMax

  for k = 0: iterMax
    x_k1 = M_alpha*x_k + N_alpha*b;

    err = norm(A*x_k - b);

    x_k = x_k1;
    if err <= tol*norm(b);
      break
    end

  end

  # Resultados

  display(['Iteraciones: k = ',num2str(k)]);
  display(['Error:= ',num2str(err)]);
  display('Tiempo de ejecucion');
  toc
end




