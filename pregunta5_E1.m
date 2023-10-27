

function pregunta5_E1()
  clc; clear;

  mvalues=[16;32;64;128;256];

  display('Metodo 1: HSS');
  for i=1:5
    m=mvalues(i);
    h=1/(m+1);
    [f,g]=calc_f_g(m,h);
    [W,T]=calc_W_T(m,h);

    disp(["Caso", num2str(i), ": m=", num2str(m)]);
    pregunta1(W, T, f, g);
  endfor

  display('Metodos 4: QR y Eliminacion Gausseana');
  for i=1:5
    m=mvalues(i);
    h = 1/(m+1);
    [f,g]=calc_f_g(m,h);
    [W,T]=calc_W_T(m,h);

    disp(["Caso", num2str(i), ": m=", num2str(m)]);
    pregunta4(W,T,f,g);
  endfor


end

function [W,T]=calc_W_T(m,h)
  bm=(1/((h)^2))*tridiag(-1,2,-1,m);
  I1=eye(m);
  I2=eye(m^2);
  tau=h;
  K=(kron(I,bm))+(kron(bm,I));
  W=(K+((3-sqrt(3))/tau))*I;
  T=(K+((3+sqrt(3))/tau))*I;

end

function [f,g]=calc_f_g(m,h)
  f=zeros(m,1);
  g=zeros(m,1);
  for j=1:m
    fj=(1-i)*j/(h*(j+1)^2);
    f(j)=real(fj);
    g(j)=imag(fj);
  end
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


