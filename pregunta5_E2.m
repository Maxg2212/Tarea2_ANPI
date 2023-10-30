%Funcion pregunta5_E2() implementa la ecuacion compleja de Helmholtz para obtener matrices de diferentes tamaños.
% que nos permiten obtener diferentes soluciones del problema de sistemas de ecuaciones con valores complejos.
% ai como ver la eficiencia de cada uno de los metodos implementados en las demas preguntas.
% Sintaxis de la funcion: pregunta5_E2().
%
% Parametros de entrada:
%
% Parametros de salida:
%
function pregunta5_E2()
  clc; clear;

  mvalues = [16; 32; 64; 128; 256];

  display('Metodo 1: HSS');
  for i=1:5
    m = mvalues(i);
    h = 1 / (m + 1);
    [W,T] = calc_W_T(m,h);
    A = W+i*T;
    [f,g] = calc_f_g(m,A,h);

    disp(["Caso", num2str(i), ": m=", num2str(m)]);
    pregunta1(W, T, f, g);
  endfor

  fprintf('\n');
  display('Metodo 2: PNHSS & PS*HSS');
  for i=1:5
    m = mvalues(i);
    h = 1 / (m + 1);
    [W, T] = calc_W_T(m, h);
    A = W+i*T;
    [f,g] = calc_f_g(m,A,h);

    fprintf('\n');
    disp(["Caso", num2str(i), ": m=", num2str(m)]);
    pregunta2(W, T, f, g);
  endfor

  display('Metodo 3: MHSS');
  for i=1:5
    m = mvalues(i);
    h = 1 / (m + 1);
    [W,T] = calc_W_T(m,h);
    A = W+i*T;
    [f,g] = calc_f_g(m,A,h);

    disp(["Caso", num2str(i), ": m=", num2str(m)]);
    pregunta3(W, T, f, g, 5000, 1e-6);
  endfor

  fprintf('\n');
  display('Metodos 4: QR y Eliminacion Gausseana');
  for i=1:5
    m = mvalues(i);
    h = 1 / (m + 1);
    [W,T] = calc_W_T(m,h);
    A = W+i*T;
    [f,g] = calc_f_g(m,A,h);

    disp(["Caso", num2str(i), ": m=", num2str(m)]);
    pregunta4(W, T, f, g);
  endfor
end

%Funcion [W,T]=calc_W_T(m,h) implementa un algoritmo que nos permite crear los valores para las matrices de tamano m*m.
%
% Sintaxis de la funcion: [W,T]=calc_W_T(m,h)
%
% Parametros de entrada:
%             m = el tamano de la matriz que se desea crear.
%             h = es el tamano de la malla que se desea crear.
%
% Parametros de salida:
%             W = la primera matriz del sistema que sumada a T equivalen a A.
%             T = la segunda matriz del sistema que multiplicada por i y sumada con W equivale a A.
%
function [W,T]=calc_W_T(m,h)
  A = tridiag(-1,2,-1,m);
  Vm = (1/h^2)*A;
  cita1 = -1;
  cita2 = 1;
  I1 = eye(m);
  I2= eye(m^2);
  K =  kron(I,Vm) + kron(Vm,I);
  W = K+(cita1*I1);
  T = cita2*I1;
  W = W*h^2;
  T = T*h^2;

end


%Funcion [f,g]=calc_f_g(m) implementa un algoritmo que nos permite crear los valores para las matrices de tamano m*m.
%
% Sintaxis de la funcion: [f,g]=calc_f_g(m)
%
% Parametros de entrada:
%             m = el tamano de la matriz que se desea crear.
%
% Parametros de salida:
%             f = la primera matriz del sistema que sumada a g equivalen a b.
%             g = la segunda matriz del sistema que multiplicada por i y sumada con f equivale a b.
%
function [f,g]=calc_f_g(m,A,h)
  Im = ones(m,1);
  f = A*Im;
  g = A*Im;
  f=f*h^2;
  g=g*h^2;
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
