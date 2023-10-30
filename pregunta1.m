%Funcion pregunta1 implementa el metodo de simple Hermitian and skew Hermitian splitting method.
% para obtener la solucion del problema de sistemas de ecuaciones con valores complejos.
% Sintaxis de la funcion: pregunta1(W, T, p, q)
% Parametros de entrada:
%         W = Matriz de medida m x m
%         T = Matriz de medida m x m
%         p = Matriz de medida m x 1
%         q = Matriz de medida m x 1
% Parametros de salida:
function pregunta1(W, T, p, q)
  %clc; clear;
  A = W + i*T;
  b = p + i*q;
  m = length(b);
  x0 = zeros(m,1);
  iterMax = 5000;
  tol = 1e-6;

  %Medir el tiempo usando tic y toc.
  tic;
  [Sa,err,iter]=HSS(A, b, x0, iterMax, tol);
  tiempo = toc;

  fprintf('Error = %.4f \n', err);
  fprintf('Tiempo de ejecución = %.4f segs\n', tiempo);
  fprintf('Iteraciones = %.4f \n', iter);
  display('-------------------------------------------------------------');


end

% La funcion HSS aproxima la solucion de un sistema de ecuaciones usando el simple Hermitian and skew-Hermitian splitting methods.
%Sintaxis de la funcion: [Sa,err,iter]=HSS(A, b, x0, iterMax, tol).
%Parametros de entrada:
%               A = Matriz de la combinacion W + iT.
%               b = Matriz de la combinacion p + iq.
%               x0 = La solucion cero del sistema de ecuaciones.
%               iterMax = Numero que representa el maximo de iteraciones.
%               tol = Numero que representa la tolerancia minima del problema.
%Parametros de salida:
%               Sa = la matriz con la aproximacion de los resultados.
%               err = El error dado por la formula de parada.
%               iter = La cantidad de iteraciones que se llevaron a cabo.
function [Sa,err,iter]=HSS(A, b, x0, iterMax, tol)
    %Obtener el tamaño de la matriz
    n = size(A);
    %Matrices

    W = real(A);
    T = imag(A);
    I = eye(n, n);

    %Solucion inicial
    x=x0;

    %Iterar hasta que se cumpla el criterio de parada o el maximo de iteraciones.
    for k=1:iterMax
      %z = (inv(I + W) * (I - i*T) * x) + (inv(I + W) * b);
      %x = ((inv(I + i*T) * (I - W) * z) + (inv(I + i*T) * b));

      inv_I_W = (I+W)\I;
      inv_I_iT = (I+i*T)\I;
      z = (inv_I_W * (I - i*T) * x) + (inv_I_W * b);
      x = ((inv_I_iT * (I - W) * z) + (inv_I_iT * b));

      %Calculo del error
      error=norm(A*x-b);
      tolerance = tol * norm(b);
      if abs(error)<=abs(tolerance);
        break;
      end
    end
    err=error;
    iter=k;
    Sa=x;
 end


