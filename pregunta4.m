function pregunta4()
  W=zeros(4,4);
  W=[12 -2 6 -2; -2 5 2 1; 6 2 9 -2; -2 1 -2 1];

  T=zeros(4,4);
  T=[6 2 7 2; 2 7 1 1; 7 1 9 0; 2 1 0 10];


  p=zeros(4,1);
  p=[9;-7;-5;7];

  q=zeros(4,1);
  q=[12;-4;17;-2];

  z1=sol1(W,T,p,q);
  z2=sol2(W,T,p,q);

end

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

function x=elim_gauss_sust_atras(A,b)
  [At,bt]=triang_sup(A,b);
  x=sust_atras(At,bt);
end

function x=resQR(A,b)
  [Q, R] = qr(A);
  c=Q'*b;
  x=sust_atras(R,c);
end

function z=sol1(W,T,p,q)
  display('---------------------------------------------------------------------')

  display('Metodo Eliminacion Gausseana:');
  tic;

  M=[W -T; T W];

  d=[p;q];

  z=elim_gauss_sust_atras(M,d);

  u=[z(1);z(2);z(3);z(4)];
  v=[z(5);z(6);z(7);z(8)];
  x=u+(i*v);
  %display(x);

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

function z=sol2(W,T,p,q)
  display('Metodo QR:');
  tic;

  M=[W -T; T W];

  d=[p;q];

  z=resQR(M,d);

  u=[z(1);z(2);z(3);z(4)];
  v=[z(5);z(6);z(7);z(8)];
  x=u+(i*v);
  %display(x);

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
