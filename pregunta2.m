%Funcion pregunta2 implementa el metodo de parameterized HSS (PNHSS) y
% parameterized single-step HSS (PS*HSS) para obtener la solucion del problema de
% sistemas de ecuaciones con valores complejos.
% Sintaxis de la funcion: pregunta2(W, T, p, q)
% Parametros de entrada:
%         W = Matriz de medida m x m
%         T = Matriz de medida m x m
%         p = Matriz de medida m x 1
%         q = Matriz de medida m x 1


function pregunta2(W, T, p, q)
  %clc; clear;

  A = W + i*T;
  b = p + i*q;

  iterMax = 5000;
  tol = 1e-6;

  %Medicion del tiempo y muestra de error y tolerancia.
  fprintf('\nMetodo PNHSS\n');
  tic;
  [Ans,err,iter] = metodo_PNHSS(A, b, iterMax, tol);
  tiempo = toc;

  fprintf('Error = %.4f \n', err);
  fprintf('Tiempo de ejecución = %.4f segs\n', tiempo);
  fprintf('Iteraciones = %.4f \n', iter);

  fprintf('\nMetodo PS*HSS\n');
  tic;
  [Ans,err,iter] = metodo_PSHSS(A, b, iterMax, tol);
  tiempo = toc;

  fprintf('Error = %.4f \n', err);
  fprintf('Tiempo de ejecución = %.4f segs\n', tiempo);
  fprintf('Iteraciones = %.4f \n', iter);

endfunction

% La funcion PNHSS aproxima la solucion de un sistema de ecuaciones usando el
%parameterized HSS
%Sintaxis de la funcion: [Ans, err, iter] = metodo_PNHSS(A, b, iterMax, tol)
%Parametros de entrada:
%               A = Matriz de la combinacion W + iT.
%               b = Matriz de la combinacion p + iq.
%               iterMax = Numero que representa el maximo de iteraciones.
%               tol = Numero que representa la tolerancia minima del problema.
%Parametros de salida:
%               Ans = la matriz con la aproximacion de los resultados.
%               err = El error dado por la formula de parada.
%               iter = La cantidad de iteraciones que se llevaron a cabo.
function [Ans, err, iter] = metodo_PNHSS(A, b, iterMax, tol)

  %Valores iniciales
  n = length(b);
  I = eye(n, n);
  x = zeros(n, 1);

  %Matrices
  W = real(A);
  T = imag(A);

  %Constantes omega y alpha
  w = 1;
  a = 1;

  %Itera el metodo hasta que se cumpla el criterio de parada o el maximo de
  %iteraciones.
  for k = [1:iterMax]

    inv_wW_T = inv(w*W + T);
    x = (inv_wW_T * ((-i)*(w*T - W)) * x) + (inv_wW_T * ((w - i) * b));
    inv_aI_wW_T = inv(a*I + w*W + T);
    x = (inv_aI_wW_T * (((a*I) - (i*((w*T) - W))) * x)) + (inv_aI_wW_T * ((w - i) * b));

    %Calculo de error y almacenamiento de resultados.
    error=norm(A*x-b);
    Ans = x;
    err = error;
    iter = k;
    if (error <= tol);
      break
    endif
  endfor

endfunction

% La funcion PS*HSS aproxima la solucion de un sistema de ecuaciones usando el
%parameterized single-step HSS
%Sintaxis de la funcion: [Ans, err, iter] = metodo_PSHSS(A, b, iterMax, tol)
%Parametros de entrada:
%               A = Matriz de la combinacion W + iT.
%               b = Matriz de la combinacion p + iq.
%               iterMax = Numero que representa el maximo de iteraciones.
%               tol = Numero que representa la tolerancia minima del problema.
%Parametros de salida:
%               Ans = la matriz con la aproximacion de los resultados.
%               err = El error dado por la formula de parada.
%               iter = La cantidad de iteraciones que se llevaron a cabo.
function [Ans, err, iter] = metodo_PSHSS(A, b, iterMax, tol)

  %Valores iniciales
  n = length(A);
  I = eye(n);
  x = zeros(n, 1);

  %Matrices
  W = real(A);
  T = imag(A);

  %Constantes omega y alpha
  w = 1;
  a = 1;

  %Itera el metodo hasta que se cumpla el criterio de parada o el maximo de
  %iteraciones.
  for k = [1:iterMax]

    inv_wW_T = inv(w*W + T);
    x = (inv_wW_T * ((-i)*(w*T - W)) * x) + (inv_wW_T * ((w - i) * b));

    %Calculo de error y almacenamiento de resultados.
    error=norm(A*x-b);
    Ans = x;
    err = error;
    iter = k;
    if (error <= tol);
      break
    endif
  endfor

endfunction
