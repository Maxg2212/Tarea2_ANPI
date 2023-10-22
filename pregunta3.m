function pregunta3 ()
  clc; clear;

    # Matrices iniciales
    W = [12, -2, 6, -2;
         -2, 5, 2, 1;
          6, 2, 9, -2;
         -2, 1, -2, 1];
    T = [6, 2, 7, 2;
         2, 7, 1, 1;
         7, 1, 9, 0;
         2, 1, 0, 10];

    # Vectores iniciales
     p = [9;
          -7;
          -5;
           7];
     q = [12;
          -4;
          17;
          -2];
      # Iteracion maxima
      iterMax = 1000;

      # Tolerancia
      tol = 1e-12;

      #Llamada a funcion MHSS
      MHSS (W, T, p, q, iterMax, tol)

endfunction


function  MHSS (W, T, p, q, iterMax, tol)

  # Variable compleja
  j = sqrt(-1);
  n = length(p);
  Im = eye(size(W));

  # Valor inicial de xk (es decir, x0)
  xk = zeros(4,1);

  # Matrices iniciales en cero
  A = zeros(n,n);
  b = zeros(n,1);
  M = zeros(n,n);
  N = zeros(n,n);

  # Matrices Complejas
  A = W + j * T;
  b = p + j * q;

  # Determinación de el valor más pequeño y el valor más grande de una matriz,
  # respectivamente
  lambda_min = min(eig(W)); #Se calcula el valor más pequeño de la matriz
  lambda_max = max(eig(W)); #Se calcula el valor más grande de la matriz

  # Calcula el valor de alfa*
  alfa_ast = sqrt(lambda_min * lambda_max);

  #Calculo de matrices inversas
  Inv_alfaxImplusT = mldivide(alfa_ast * Im + T, eye(size(alfa_ast * Im + T)));
  Inv_alfaxImplusW = mldivide(alfa_ast * Im + W, eye(size(alfa_ast * Im + W)));

  # Calcula el valor M(alfa)
  M = (Inv_alfaxImplusT) * ((alfa_ast * Im + j * W)) * (Inv_alfaxImplusW) * (alfa_ast * Im - j * T);

  # Calcula el valor de N(alfa)
  N = (1 - j) * alfa_ast * (Inv_alfaxImplusT) * (Inv_alfaxImplusW);

  # Calcular sucesion desde i = 0 hasta iterMax
  for i = 0 : iterMax
      xk_next = (M * xk + N * b);

      error = norm(xk - b); #Critero de parada

      xk = xk_next; #Actualizacion de valor xk
      if error <= tol
        break
      endif

endfor

  # Resultados
  display('Resultados del metodo MHSS')
  display('Aproximacion de xk = [')
  display(xk')
  display(']')
  display(['Valor de alfa = ', num2str(alfa_ast)])
  display(['Iteraciones: i = ', num2str(i)])
  display(['Error: ||A*xk - b||_2 = ', num2str(error)])


endfunction
