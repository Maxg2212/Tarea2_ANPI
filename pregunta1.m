function pregunta1()
  clc; clear;
  A = [12 -2 6 -2; -2 5 2 1; 6 2 9 -2; -2 1 -2 1] + (i*[6 2 7 2; 2 7 1 1; 7 1 9 0; 2 1 0 10]);
  b = [9;-7;-5;7] + (i*[12;-4;17;-2]);
  x0 = [0; 0; 0; 0];
  iterMax = 1000;
  tol = 1e-12;

  %Medir el tiempo usando tic y toc.
  tic;
  [Sa,err,iter]=HSS(A, b, x0, iterMax, tol);
  tiempo = toc;

  fprintf('Error = %.4f \n', err);
  fprintf('Tiempo de ejecución = %.4f segs\n', tiempo);
  fprintf('Iteraciones = %.4f \n', iter);


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
      if abs(error)<=tol;
        err=error;
        iter=k;
        Sa=x;
        break;
      end
      err=error;
      iter=k;
      Sa=x;
    end
 end


