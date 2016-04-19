// Etape k -> k+1



//Matrices des Phip
phii = Phip(log((ones(M,1)./Grille)*(Grille')),1,vo);
for p=2:nbelement
  phii = [phii,Phip(log((ones(M,1)./Grille)*(Grille')),p,vo)];
end;



// Matrices liees au Call
tmp =((1+exp(pas))*0.5*ones(NbS,1)*Grille' -Strike*ones(1,M)).*(ones(NbS,1)*Grille' >Strike*ones(1,M));
p=1;
tmp1  = tmp*(phii(:,1+(p-1)*M:p*M)');
q = tmp1*diag(Pik1);
for p=2:nbelement
  tmp1  = tmp*phii(:,1+(p-1)*M:p*M)';
  q = [q, tmp1*diag(Pik1)];
end;
clear tmp;
clear tmp1; 


clear Alpha Beta gam1 gam3 mu1 mu3 Ac Bc;
//matrices liees a l'equation probabiliste
Alpha = [sum(phii(:,1:M),'c'),sum(phii(:,1+M:2*M),'c'),sum(phii(:,2*M+1:3*M),'c')];
//matrices liees a l'equation martingale
Beta  = [(phii(:,1:M)*Grille)./Grille,(phii(:,M+1:2*M)*Grille)./Grille,(phii(:,2*M+1:3*M)*Grille)./Grille];

gam1 =( Alpha(2:M,3) - Beta(2:M,3))./(Beta(2:M,1).*Alpha(2:M,3)-Alpha(2:M,1).*Beta(2:M,3)+%eps);
gam3 =( Alpha(1:M-1,1) - Beta(1:M-1,1))./(Beta(1:M-1,3).*Alpha(1:M-1,1)-Alpha(1:M-1,3).*Beta(1:M-1,1)+%eps);
mu1 =( Alpha(2:M,3).*Beta(2:M,2) -Alpha(2:M,2).* Beta(2:M,3))./(Beta(2:M,1).*Alpha(2:M,3)-Alpha(2:M,1).*Beta(2:M,3)+%eps);
mu3 =( Alpha(1:M-1,1).*Beta(1:M-1,2) -Alpha(1:M-1,2).* Beta(1:M-1,1))./(Beta(1:M-1,3).*Alpha(1:M-1,1)-Alpha(1:M-1,3).*Beta(1:M-1,1)+%eps);

mu1 = [0;mu1];
mu3 = [mu3;0];
gam1 = [1;gam1];
gam3 = [gam3;1];

//p1=p1/1.35;
//p1 =0.0001;
tmp1 =eye(M,M) + diag(mu1);
gam1b = gam1;
clear tmp;
tmp3 = eye(M,M) + diag(mu3);
gam3b = gam3;



lb=zeros(M,1);
ub=%inf*ones(M,1);
x0=diag(PivBS);

// $$$ if k==2 then
// $$$   x0=lb;
// $$$ else
// $$$    x0=xsav;
// $$$ end,

p1=0.00;
//pause;
function [f,g,ind]=costf(x,ind)
  clear Imp_exp grad_imp Call_exp_tmp; 
  Call_exp_tmp = q * [gam1 - x.*mu1;x; gam3 - x.*mu3];
  grad_imp=ones(NbS,M);
  for i=1:NbS
	Imp_exp(1,i) = Implied(Call_exp_tmp(i),stock,Strike(i),k/n,0);

	grad_imp(i,:) = (-q(i,1:M)'.*mu1+q(i,M+1:2*M)'-q(i,2*M+1:3*M)'.*mu3)'/Vega(Strike(i),k/n,stock,0,Imp_exp(1,i));
  end;
  f=p1*(x-x0)'*(x-x0) +(Imp_exp-Imp(k,:))*(Imp_exp-Imp(k,:))'/norm(Imp(k,:),2); 
  g=2*p1*(x-x0)+2*grad_imp'*(Imp_exp-Imp(k,:))'/norm(Imp(k,:),2); 
endfunction

// $$$ function [f,g,ind]=costf(x,ind)
// $$$   clear  grad_imp Call_exp_tmp; 
// $$$   Call_exp_tmp = q * [gam1 - x.*mu1;x; gam3 - x.*mu3];
// $$$   grad_imp=ones(NbS,M);
// $$$   g=0;f=0;
// $$$   for i=1:NbS
// $$$ 	g =g+ (-q(i,1:M)'.*mu1+q(i,M+1:2*M)'-q(i,2*M+1:3*M)'.*mu3)'*(Call_exp_tmp(i)-Call_synt(k,i))/Vega(Strike(i),k/n,stock,0,Imp(k,i));
// $$$ 	f = f+(Call_exp_tmp(i)-Call_synt(k,i))*(Call_exp_tmp(i)-Call_synt(k,i))/q*Vega(Strike(i),k/n,stock,0,Imp(k,i));
// $$$   end;
// $$$   f=f+p1*(tmp1*x-gam1b)'*(tmp1*x-gam1b)/(gam1b'*gam1b)+p1*(tmp3*x-gam3b)'*(tmp3*x-gam3b)/(gam3b'*gam3b);
// $$$   g=g'+2*p1*tmp1'*(tmp1*x-gam1b)+2*p1*tmp3'*(tmp3*x-gam3b);
// $$$ endfunction

[fopt,xopt,gopt] = optim(costf,'b',lb,ub,x0);
////pause;

// $$$ pause;
// $$$ p1 = 1.1*(fopt)/((tmp1*xopt-gam1b)'*(tmp1*xopt-gam1b)/(gam1b'*gam1b));
// $$$ 
// $$$ function [f,g,ind]=costf(x,ind)
// $$$   clear Imp_exp grad_imp Call_exp_tmp; 
// $$$   Call_exp_tmp = q * [gam1 - x.*mu1;x; gam3 - x.*mu3];
// $$$   grad_imp=ones(NbS,M);
// $$$   for i=1:NbS
// $$$ 	Imp_exp(1,i) = Implied(Call_exp_tmp(i),stock,Strike(i),k/n,0);
// $$$ 
// $$$ 	grad_imp(i,:) = (-q(i,1:M)'.*mu1+q(i,M+1:2*M)'-q(i,2*M+1:3*M)'.*mu3)'/Vega(Strike(i),k/n,stock,0,Imp_exp(1,i));
// $$$   end;
// $$$   f=p1*(tmp1*x-gam1b)'*(tmp1*x-gam1b)/(gam1b'*gam1b)+p1*(tmp3*x-gam3b)'*(tmp3*x-gam3b)/(gam3b'*gam3b)+(Imp_exp-Imp(k,:))*(Imp_exp-Imp(k,:))'/norm(Imp(k,:),2); 
// $$$   g=2*p1*tmp1'*(tmp1*x-gam1b)+2*p1*tmp3'*(tmp3*x-gam3b)+2*grad_imp'*(Imp_exp-Imp(k,:))'/norm(Imp(k,:),2); 
// $$$ endfunction
// $$$ 
// $$$ [fopt,xopt,gopt] = optim(costf,'b',lb,ub,xopt);

clear costf;

xsav=xopt;
disp(fopt,"Resolution du systeme");


lambda1  = gam1 - xopt.*mu1;
lambda3  = gam3 - xopt.*mu3;
xopt1 = [lambda1;xopt;lambda3];
disp(norm(q*xopt1-Call_synt(k,:)',1)/norm(Call_synt(k,:),1),"Erreur Call");



// Calcul des poids
tmp = diag(xopt1(1:M))*phii(:,1:M);
for p=2:nbelement
  tmp = tmp+ diag(xopt1(1+(p-1)*M:p*M))*phii(:,1+(p-1)*M:p*M);
end;
Pik = tmp;
clear tmp;
// Calcul de la vol

//Variation quadratique
tmp=ones(Grille);
for i=1:M
  tmp(i) = n*sum((Grille-Grille(i)).*(Grille-Grille(i)).*Pik(i,:)');
  if tmp(i)<0 then
	tmp(i)=vo(i);
  else
	 tmp(i)=sqrt(tmp(i))/Grille(i);
  end,
end;





//Malliavin
// $$$ tmp = (ones(M,1)*Grille'-Grille*ones(1,M)).*(ones(M,1)*Grille'>Grille*ones(1,M));
// $$$ Callinf = Pik*tmp;
// $$$ Callinf2 = Callinf(:,2:M);
// $$$ Callinf = Callinf(:,1:M-1);
// $$$ clear tmp;
// $$$ Alphak =n*(ones(M,1)*Grille'-Grille*ones(1,M));
// $$$ Alphak2 = Alphak(:,2:M);
// $$$ Alphak = Alphak(:,1:M-1);
// $$$ tmp = sum(((Alphak2.*Callinf-Callinf2.*Alphak).^2).*Pik(:,1:M-1),'c');
// $$$ tmp1 =sum((Alphak2.*Callinf-Callinf2.*Alphak).*(Pik(:,1:M-1).^2),'c')+%eps;
// $$$ sxopt=zeros(M,1);
// $$$ sxopt = (sqrt((tmp./tmp1).*(tmp./tmp1>0)))./Grille;
// $$$ //tmp=interpln([Grille';sxopt'],Grille);
// $$$ clear tmp;
// $$$ tmp=sxopt;

//Optimisation sur la Vol
reg = eye(M,M)+diag([-1;-ones(M-2,1)],1);
x=vo';
POIDS=eye(M,M);
p2=10;
p3=0;

if k==2 then
function [f,g,ind]=costf(x,ind)
  f = (POIDS*x - POIDS*tmp)'*(POIDS*x - POIDS*tmp)+p2*(reg*x)'*(reg*x);
  g = 2*POIDS*(POIDS*x - POIDS*tmp)+2*p2*reg'*(reg*x);
endfunction
else
 function [f,g,ind]=costf(x,ind)
  f = (POIDS*x - POIDS*tmp)'*(POIDS*x - POIDS*tmp)+p2*(reg*x)'*(reg*x)+p3*(Vol(k-1,:)-x')*(Vol(k-1,:)'-x);
  g = 2*POIDS*(POIDS*x - POIDS*tmp)+2*p2*reg'*(reg*x)+2*p3*(-Vol(k-1,:)'+x);
 endfunction
end,  
x0 =zeros(vo');


ub=(%eps)^(-1)*ones(M,1);
[sfopt,sxopt,sgopt]=optim(costf,x0,'qn');

clear tmp costf;


Vol(k,:) = sxopt';

// Vraie Volatilite
clear f;
deff('[z]=f(x)','z=VOL(x,(k-1)/n)');
vraievol = feval(Grille,f);
clear f;
// $$$ ii=find((Grille<=0.7*stock)&(Grille*exp(pas)>0.7*stock));
// $$$ jj=find((Grille<=1.3*stock)&(Grille*exp(pas)>1.3*stock));
// $$$ disp(norm(sxopt(ii:jj)-vraievol(ii:jj))/norm(vraievol(ii:jj)),"Erreur Vol 70% et 130%");
// $$$ ii=find((Grille<=0.8*stock)&(Grille*exp(pas)>0.8*stock));
// $$$ jj=find((Grille<=1.2*stock)&(Grille*exp(pas)>1.2*stock));
// $$$ disp(norm(sxopt(ii:jj)-vraievol(ii:jj))/norm(vraievol(ii:jj)),"Erreur Vol 80% et 120%");
// $$$ disp(norm(sxopt-vraievol)/norm(vraievol),"Erreur Vol");

tps = k/n;
tmp=q*xopt1;
// $$$ if k==2
// $$$   tmpr = read("dupiresprices_put.data",n*NbS,3);
// $$$   tmpr1 = tmpr(:,3);
// $$$   Put_synt= matrix(tmpr1,NbS,n)';
// $$$ end,
if k==2
  tmpr = read("put.dat",n*NbS,4);
  tmpr1 = tmpr(:,3);
  Put_synt= matrix(tmpr1,NbS,n)';
end,
Putverif = Put_synt(k,:)';
for i=1:NbS
  //Putverif(i,1) = Put_BS(stock,Strike(i),tps,0,0,VO);
  Imp_exp(1,i) = Implied2(tmp(i),stock,Strike(i),k/n,0);
end;
clear tmp tmpr tmpr1;
Pikp1 = (Pik1' * Pik)';
tmp = (ones(NbS,1)*Grille' <Strike*ones(1,M)).* (-ones(NbS,1)*Grille'+ Strike*ones(1,M));
Put = tmp * Pikp1;
disp(norm(Putverif -Put,1)/norm(Putverif,1),"Erreur sur les Puts");

// $$$ xbasc();
// $$$ xsetech([0,0,0.5,1]);
// $$$ clear tmp tmp1;
// $$$ tmp=interpln([Strike';Imp(k,:)],Grille);
// $$$ tmp1=interpln([Strike';Imp_exp],Grille);
// $$$ xset("use color",0);
// $$$ lege=['Real volatility' 'Approached Volatility' 'Real Implied volatility'  'Approached Implied volatility'];
// $$$ xset("font size",0.1);
// $$$ plot2d(Grille,[vraievol  sxopt tmp' tmp1' ],[1 -1 2 3],rect=[70,0.06,130,0.18]);
// $$$ legends(lege,[1 -1 2 -2],4);
// $$$ 
// $$$ xtitle("Volatility in Dupire model: Pas " + string(k));

// $$$ xsetech([0,0.5,0.5,0.5]);
// $$$ lege = ['real density' 'empirical density'];
// $$$ plot2d(Grille,[Piv(M/2,:)'  Pik(M/2,:)' ], [1 2]);
// $$$ legends(lege,[1 -1 2 -2],4);
// $$$ xtitle("density at the starting point");

// $$$ xsetech([0.5,0,0.5,1]);
// $$$ lege = ['real call' 'approached call' 'real put' 'approached put'];
// $$$ plot2d(Strike,[Call_synt(k,:)' q*xopt1  Put Putverif],[1 -1 2 -2]);
// $$$ legends(lege,[1 -1 2 -2],4);
// $$$ xtitle("Call and Put price in Dupire model");

// $$$ xsetech([0.5,0.5,0.5,0.5]);
// $$$ ii=find((Grille<=0.7*stock)&(Grille*exp(pas)>0.7*stock));
// $$$ jj=find((Grille<=1.3*stock)&(Grille*exp(pas)>1.3*stock));
// $$$ lege = ['real density' 'empirical density' 'real density' 'empirical density'];
// $$$ plot2d(Grille,[Piv(ii,:)'  Pik(ii,:)' Piv(jj,:)'  Pik(jj,:)'], [1 2 3 4]);
// $$$ legends(lege,[1 -1 2 -2],4);
// $$$ xtitle("density at 70 and 130");

// Graph pour Heston
xbasc();
xset("use color",0);
xsetech([0,0,0.5,1]);
lege = ['real call' 'approached call' 'real put' 'approached put'];
plot2d(Strike,[Call_synt(k,:)' q*xopt1  Put Putverif],[1 -1 -2 2]);
legends(lege,[1 -1 -2 2],4);
xtitle("Call and Put price in Heston model");



xsetech([0.5,0,0.5,1]);
lege = ['real call' 'approached call' 'real put' 'approached put'];
clear tmp tmp1 tmp2 tmp3;
callverif=q*xopt1;
for i=1:NbS
  tmp(1,i) = Implied2(Call_synt(k,i),100,Strike(i),k*T/n,0);
  tmp1(1,i) = Implied2(callverif(i),100,Strike(i),k*T/n,0);
  tmp2(1,i) = Implied3(Putverif(i),100,Strike(i),k*T/n,0);
  tmp3(1,i) = Implied3(Put(i),100,Strike(i),k*T/n,0);
end
plot2d(Strike,[tmp' tmp1' tmp2' tmp3'],[1 -1 2 -2]);
legends(lege,[1 -1 2 -2],4);
xtitle("Call and Put price in Heston model");

clear tmp;
//tmp='tmpdir/step'+string(k);
//xsave(tmp,0);
// $$$ if k==2 then
// $$$   pause;
// $$$ end
// $$$ if k==6 then
// $$$   pause;
// $$$ end
if k==12 then 
  ////pause;
end

