function pregunta2()
  clc; clear;

  W = [12, -2, 6, -2;
       -2, 5, 2, 1;
        6, 2, 9, -2;
       -2, 1, -2, 1];

  T = [6, 2, 7, 2;
       2, 7, 1, 1;
       7, 1, 9, 0;
       2, 1, 0, 10];

  A = W + i*T;

  p = [9;
      -7;
      -5;
       7];

  q = [12;
       -4;
       17;
       -2];

  b = p + i*q;

  iterMax = 1000;
  tol = 1e-12;
  metodo_PNHSS(A, b, iterMax, tol)

endfunction

function Ans = metodo_PNHSS(A, b, iterMax, tol)

  n = size(A);
  I = eye(n);
  W = real(A);
  T = imag(A);
  x = zeros(4, 1);
  w = 1;
  a = 1;

  for k = [1:iterMax]

    inv_wW_T = inv(w*W + T);
    x = (inv_wW_T * ((-i)*(w*T - W)) * x) + (inv_wW_T * ((w - i) * b));
    inv_aI_wW_T = inv(a*I + w*W + T);
    x = (inv_aI_wW_T * (((a*I) - (i*((w*T) - W))) * x)) + (inv_aI_wW_T * ((w - i) * b));

    error=norm(A*x-b);
    if (error <= tol);
      Ans=x;
      break
    endif
  endfor


endfunction

