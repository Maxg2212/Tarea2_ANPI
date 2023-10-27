% La funcion pregunta4 implementa metodos de eliminacion gausseana y factorizacion QR
% para obtener la solucion del problema de sistemas de ecuaciones con valores complejos.
% Sintaxis de la funcion: pregunta4(W, T, p, q)
% Parametros de entrada:
%         W = Matriz de medida m x m
%         T = Matriz de medida m x m
%         p = Matriz de medida m x 1
%         q = Matriz de medida m x 1
% Parametros de salida:
function pregunta4(W, T, p, q)
  %W=zeros(4,4);
  %W=[12 -2 6 -2; -2 5 2 1; 6 2 9 -2; -2 1 -2 1];

  %T=zeros(4,4);
  %T=[6 2 7 2; 2 7 1 1; 7 1 9 0; 2 1 0 10];


  %p=zeros(4,1);
  %p=[9;-7;-5;7];

  %q=zeros(4,1);
  %q=[12;-4;17;-2];

  s1=elim_gauss(W,T,p,q);
  s2=qr_method(W,T,p,q);

end


% La funcion sust_atras aplica el metodo de sustitucion hacia atras para resolver sistemas de ecuaciones
% Sintaxis de la funcion: x=sust_atras(A,b)
% Parametros de entrada:
%         A = Matriz de medida 2m x 2m
%         b = Matriz de medida 2m x 1
% Parametros de salida:
%         x = Matriz de medida 2m x 2m
function x=sust_atras(A,b)
  m=length(b);
  x=zeros(m,1);
  for i=m:-1:1
    aux=0;
    for j=i+1:m
      aux+=A(i,j)*x(j);
    end
    x(i)=(1/A(i,i))*(b(i)-aux);
  end
end


% La funcion triang_sup convierte una matriz de entrada en triangular superior
% Sintaxis de la funcion: [Ar,br]=triang_sup(A,b)
% Parametros de entrada:
%         A = Matriz de medida 2m x 2m
%         b = Matriz de medida 2m x 1
% Parametros de salida:
%         Ar = Matriz de medida 2m x 2m
%         br = Matriz de medida 2m x 1
function [Ar,br]=triang_sup(A,b)
  m=size(A,1);
  At=[A b];
  for k=1:m-1
    for i=k+1:m
      p=At(i,k)/At(k,k);
      for j=k:m+1
        At(i,j)=At(i,j)-p*At(k,j);
      end
    end
  end
  Ar=At(:,1:m);
  br=At(:,m+1);
end


% La funcion elim_gauss procesa las matrices de entrada para encontrar una solucion para x
% mediante el metodo de eliminacion Gausseana
% Sintaxis de la funcion: x=elim_gauss(W,T,p,q)
% Parametros de entrada:
%         W = Matriz de medida m x m
%         T = Matriz de medida m x m
%         p = Matriz de medida m x 1
%         q = Matriz de medida m x 1
% Parametros de salida:
%         x = Matriz de medida 2m x 1
function x=elim_gauss(W,T,p,q)
  display('---------------------------------------------------------------------')

  display('Metodo Eliminacion Gausseana:');
  tic;

  M=[W -T; T W];

  d=[p;q];

  [At,bt]=triang_sup(M,d);
  z=sust_atras(At,bt);

  %u=[z(1);z(2);z(3);z(4)];
  %v=[z(5);z(6);z(7);z(8)];
  %x=u+(i*v);
  %display(x);

  %-
  m=size(z,1);
  maux = m/2;
  u=zeros(maux,1);
  v=zeros(maux,1);

  for j=1:maux
    u(j)=z(j);
  endfor

  iaux=1;
  for j=maux+1:m
    v(iaux)=z(j);
    iaux=iaux+1;
  endfor


  x=u+(i*v);
  %display(x);

  %-

  A=W+(i*T);
  b=p+(i*q);
  err=norm((A*x)-b);
  time=toc;

  display('Error:');
  display(err);
  fprintf('\n');
  display('Tiempo de ejecucion:');
  display(time);
  fprintf('\n');
  fprintf('\n');
  display('---------------------------------------------------------------------')


end


% La funcion qr_method procesa las matrices de entrada para encontrar una solucion para x
% mediante el metodo QR
% Sintaxis de la funcion: x=qr_method(W,T,p,q)
% Parametros de entrada:
%         W = Matriz de medida m x m
%         T = Matriz de medida m x m
%         p = Matriz de medida m x 1
%         q = Matriz de medida m x 1
% Parametros de salida:
%         x = Matriz de medida 2m x 1
function x=qr_method(W,T,p,q)
  display('---------------------------------------------------------------------')


  display('Metodo QR:');
  tic;

  M=[W -T; T W];

  d=[p;q];

  [Q, R] = qr(M);
  c=Q'*d;
  x1=sust_atras(R,c);

  %u=[x1(1);x1(2);x1(3);x1(4)];
  %v=[x1(5);x1(6);x1(7);x1(8)];
  %x=u+(i*v);
  %display(x);

  %-
  m=size(x1,1);
  maux = m/2;
  u=zeros(maux,1);
  v=zeros(maux,1);

  for j=1:maux
    u(j)=x1(j);
  endfor

  iaux=1;
  for j=maux+1:m
    v(iaux)=x1(j);
    iaux=iaux+1;
  endfor


  x=u+(i*v);
  %display(x);

  %-

  A=W+(i*T);
  b=p+(i*q);
  err=norm((A*x)-b);
  time=toc;

  display('Error:');
  display(err);
  fprintf('\n');
  display('Tiempo de ejecucion:');
  display(time);
  fprintf('\n');
  fprintf('\n');
  display('---------------------------------------------------------------------')

end
