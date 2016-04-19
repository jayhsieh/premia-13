mm=20;
nbelement = nbelement+mm;
deff('[z]=f(x,p)','z=Phip(x,p,vo)');
tmp=[1:1:nbelement];
phii=feval(log(Grille/stock),tmp,f);
Alpha=sum(phii,'r')';
Beta = phii'*Grille/stock;

tmp1 = Grille*(exp(pas)+1)/2;
gam = ((ones(NbS,1)*tmp1' - Strike*ones(1,M)).*(ones(NbS,1)*Grille' > Strike*ones(1,M))) *phii;
gam=gam';
reg = eye(nbelement,nbelement)+diag(-ones(1,nbelement-1),-1);

p1 =1e-8;
Q = gam*gam'/(Call_synt(1,:)*Call_synt(1,:)')+p1*reg'*reg;
P = -Call_synt(1,:)*gam'/(Call_synt(1,:)*Call_synt(1,:)');
P=P';
C = [Alpha,Beta]';
b = [1;1];
lb = 0*ones(nbelement,1);
ub=%inf*ones(nbelement,1);
me = 2;

[xopt,lagr,fopt]=quapro(Q,P,C,b,lb,ub,me);
p1 = 1.01*(fopt+0.5)/(0.5*xopt'*reg'*reg*xopt);
Q = gam*gam'/(Call_synt(1,:)*Call_synt(1,:)')+p1*reg'*reg;
P = -Call_synt(1,:)*gam'/(Call_synt(1,:)*Call_synt(1,:)');
P=P';
C = [Alpha,Beta]';
b = [1;1];
lb = 0*ones(nbelement,1);
ub=%inf*ones(nbelement,1);
me = 2;
disp(p1,"Optim 2 avec p1");
[xopt,lagr,fopt]=quapro(Q,P,C,b,lb,ub,me);
Pi1 = phii*xopt;
Volat = sqrt(n*sum((Grille-stock).*(Grille-stock).*Pi1))/stock;
Cexp = (gam'*xopt)';

disp(norm(C*xopt-b)/norm(b),"Les contraintes sont respectees: ");
disp(norm(gam'*xopt-Call_synt(1,:)')/norm(Call_synt(1,:)),"Erreur sur les calls: ");
disp(norm(Volat-VOL(stock,0+%eps))/VOL(stock,0+%eps),"Erreur sur la Vol: ");
disp(norm(Piv(M2,:)-Pi1',%inf)/norm(Pi1,%inf),"Erreur sur les poids");
disp(Volat,"Vol");
//Graphiques
xbasc();
xsetech([0,0,0.5,0.5]);
plot2d(Grille,[Piv(M2,:)' PivBS(M2,:)' Pi1]);
xsetech([0.5,0,0.5,0.5]);
plot2d(Strike,[Call_synt(1,:)' Cexp']);

nbelement = nbelement-mm;
clear mm;

