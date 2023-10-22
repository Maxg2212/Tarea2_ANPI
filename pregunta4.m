function pregunta4()
  W=zeros(4,4);
  W(1,1)=12; W(1,2)=-2; W(1,3)=6; W(1,4)=-2;
  W(2,1)=-2; W(2,2)=5; W(2,3)=2; W(2,4)=1;
  W(3,1)=6; W(3,2)=2; W(3,3)=9; W(3,4)=-2;
  W(4,1)=-2; W(4,2)=1; W(4,3)=-2; W(4,4)=1;

  T=zeros(4,4);
  T(1,1)=6; T(1,2)=2; T(1,3)=7; T(1,4)=2;
  T(2,1)=-2; T(2,2)=7; T(2,3)=1; T(2,4)=1;
  T(3,1)=7; T(3,2)=1; T(3,3)=9; T(3,4)=-0;
  T(4,1)=2; T(4,2)=1; T(4,3)=0; T(4,4)=10;

  p=zeros(4,1);
  p(1)=9;
  p(2)=-7;
  p(3)=-5;
  p(4)=7;

  q=zeros(4,1);
  q(1)=12;
  q(2)=-4;
  q(3)=17;
  q(4)=-2;

  z=sol(W,T,p,q);
  display(sol);

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
  x=elim_gauss_sust_atras(R,c);
end

function z=sol(W,T,p,q)
  M=zeros(2,2);
  M(1,1)=W;
  M(1,2)=-T;
  M(2,1)=T;
  M(2,2)=W;

  d=zeros(2,1);
  d(1)=p;
  d(2)=q;

  z=zeros(2,1);
  z=resQR(M,d);
end
