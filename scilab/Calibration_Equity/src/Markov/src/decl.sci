function [c]=Call_BS(s,k,t,r,q,sigma)
  d=(log(s./k)+(r+sigma*sigma/2)*t)/(sigma*sqrt(t));
  d2=d-sigma*sqrt(t);
  c=s.*cdfnor("PQ",d,0*ones(d),ones(d))-exp(-r*t)*k.*cdfnor("PQ",d2,0*ones(d2),ones(d2))
endfunction
function [c]=Put_BS(s,k,t,r,q,sigma)
  d=(log(s./k)+(r+sigma*sigma/2)*t)/(sigma*sqrt(t));
  d2=d-sigma*sqrt(t);
  c=-s.*cdfnor("PQ",-d,0*ones(d),ones(d))+exp(-r*t)*k.*cdfnor("PQ",-d2,0*ones(d2),ones(d2))
endfunction
function [c]=Call_BS_n(s,k,t,r,q,sigma)
  c = (s-k)*(1 -cdfnor("PQ",(k-s)/(15*sqrt(t)),0,1))+15*sqrt(t)*exp(-(s-k)*(s-k)/(2*15*15*t))/(sqrt(2*%pi));
endfunction
function [c]=Put_BS_n(s,k,t,r,q,sigma)
  c = (k-s)*cdfnor("PQ",(k-s)/(15*sqrt(t)),0,1)+15*sqrt(t)*exp(-(s-k)*(s-k)/(2*15*15*t))/(sqrt(2*%pi));
endfunction

function [c] = Brigo(s,k,t,r,q,sigma)
  c = 0.5*Call_BS(s,k,t,r,q,0.2)+0.5*Call_BS(s,k,t,r,q,0.7);
endfunction

function [c] = dens(x,y,t,sigma)
  c = exp(-(log(y/x)+sigma*sigma*t/2)^2/(2*sigma*sigma*t))/(y*sqrt(2*%pi*sigma*sigma*t));
endfunction

function [c]=VOL(s,t)
  c=0.05+0.1*exp(-s/100)+0.01*t;
// $$$   if (s>=90)&(s<110) then
// $$$ 	c=0.3;
// $$$   else
// $$$ 	 c=0.15;
// $$$   end,
//  c = 15*ones(s)./s;
//  c=VO*ones(s);
//  c = sqrt((0.5*0.04*s*s*dens(stock,s,t,0.2)+0.5*0.49*s*s*dens(stock,s,t,0.7))/(0.5*s*s*dens(stock,s,t,0.2)+0.5*s*s*dens(stock,s,t,0.7)+%eps));
endfunction
function [v] = Vega(k,t,x,r,sigma)
  d=(log(x./k)+(r+sigma*sigma/2)*t)/(sigma*sqrt(t));  
  v=sqrt(t)*x*exp(-d*d/2)/sqrt(2*%pi)+%eps;
endfunction
  
function [c]=Implied(price,x,k,t,r);
  c=implied_volatility(price,r,x,t,k,%eps);
endfunction
function [c]=Implied2(price,x,k,t,r);
	deff('[z]=f(sigma)','z=Call_BS(x,k,t,r,0.262,sigma)-price');
	c =fsolve(0.1,f,%eps);
endfunction
function [c]=Implied3(price,x,k,t,r);
	deff('[z]=f(sigma)','z=Put_BS(x,k,t,r,0,sigma)-price');
	c =fsolve(0.1,f,%eps);
endfunction

function [Phi]=Phip2(x,p,vopar)
  pp=p+1;
  b=cdfnor("X",0,1,pas,1-pas);
  for i = 1:nbelement+1
	d = cdfnor("PQ",b(i),0,1)+(1-2*pas)/(nbelement+1);
	b(i+1)=cdfnor("X",0,1,d,1-d);
  end;
  a = vopar*sqrt(1/n)*b;
  if x < a(pp-1) then
	Phi=0;
  end,
  if x >= a(pp+1) then
	Phi=0;
  end,
  if ( x >= a(pp-1)) & (x< a(pp)) then
	Phi=(x-a(pp-1))/(a(pp)-a(pp-1));
  end,
  if ( x >= a(pp)) & (x< a(pp+1)) then
	Phi=(-x+a(pp+1))/(a(pp+1)-a(pp));
  end,  
endfunction

function [Phi]=Phip(x,p,vopar)
  pp=p+1;
  b=cdfnor("X",0,1,pas,1-pas);
  for i = 1:nbelement+1
	d = cdfnor("PQ",b(i),0,1)+(1-2*pas)/(nbelement+1);
	b(i+1)=cdfnor("X",0,1,d,1-d);
  end;
  a =sqrt(1/n)*b*vopar;
  Phi =zeros(x)- ((x-a(pp-1,:)'*ones(1,size(x,'c')))./(a(pp-1,:)'*ones(1,size(x,'c'))-a(pp,:)'*ones(1,size(x,'c')))).*((x>=(a(pp-1,:)'*ones(1,size(x,'c'))))&(x<(a(pp,:)'*ones(1,size(x,'c')))))+(-x+a(pp+1,:)'*ones(1,size(x,'c')))./(a(pp+1,:)'*ones(1,size(x,'c'))-a(pp,:)'*ones(1,size(x,'c'))).*((x>=(a(pp,:)'*ones(1,size(x,'c'))))&(x<(a(pp+1,:)'*ones(1,size(x,'c')))));
endfunction
