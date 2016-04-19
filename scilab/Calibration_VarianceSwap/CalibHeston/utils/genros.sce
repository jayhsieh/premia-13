function [f,g,ind]=genros(x,ind)
  n=size(x,'*');
  dzs(2)=100.0d+0;
  f=1.0d+0
  for i=2:n
    f=f + dzs(2)*(x(i)-x(i-1)**2)**2 + (1.0d+0-x(i))**2
  end
  g=ones(n,1);
  g(1)=-4.0d+0*dzs(2)*(x(2)-x(1)**2)*x(1)
  nm1=n-1
  for i=2:n-1
    ip1=i+1
    g(i)=2.0d+0*dzs(2)*(x(i)-x(i-1)**2)
    g(i)=g(i) -4.0d+0*dzs(2)*(x(ip1)-x(i)**2)*x(i) - 2.0d+0*(1.0d+0-x(i))
  end
  g(n)=2.0d+0*dzs(2)*(x(n)-x(nm1)**2) - 2.0d+0*(1.0d+0-x(n))
endfunction


