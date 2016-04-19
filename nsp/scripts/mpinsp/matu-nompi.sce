t=cputime();
Nt= 60; 
result=[];
Maturities = linspace(1. /12,5,Nt);
for Maturity=Maturities
  exec ('premia.sce');
  result= [result; Maturity,L(1)(3)];
end
t=cputime() - t;
save('matu-nompi.bin',result,CPU=t);








  
  

