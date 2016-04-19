//clear;
stacksize(20000000);
// Parametres
NbS=20;
stock=100;
n=12;
T=1;
nbelement=3;
Strike=[80:40/(NbS-1):120]';


//fichiers de declarations
getf ./decl.sci;

//construction de la grille
M=250;
pas=2*log(180/stock)/M;
M2=M/2;
tmp=[1:1:M]';
Grille = stock*exp(pas*(tmp-M2));
Temps=T*[1:1:n]'/n;
disp(pas,"Le pas est de");
disp(Grille(1),"Minimum de la grille");
disp(Grille(M),"Maximum de la grille");

//pause;
VO=0.3;

// Creation de donnees synthetiques Dupire
// $$$ tmp = read("dupiresprices.data",n*NbS,3);
// $$$ tmp1 = tmp(:,3);
// $$$ Call_synt= matrix(tmp1,NbS,n)';

// Creation de donnees synthetiques Heston
tmp = read("Call.dat",n*NbS,4);
tmp1 = tmp(:,3);
Call_synt= matrix(tmp1,NbS,n)';

for tps=1:n
  for i=1:NbS
    //Call_synt(tps,i)=Call_BS(stock,Strike(i),Temps(tps),0,0,VO);
    Imp(tps,i) = Implied(Call_synt(tps,i),stock,Strike(i),Temps(tps),0);
  end;
end;
Impb=interpln([Strike';Imp(1,:)],Grille);


//Pas 1
// Vrai poids

PivBS=cdfnor("PQ",(log((ones(M,1)./Grille)*(Grille*exp(pas))')+(Impb'*ones(1,M)).*(Impb'*ones(1,M))/(2*n))./((Impb*sqrt(1/n))'*ones(1,M)),zeros(M,M),ones(M,M))-cdfnor("PQ",(log((ones(M,1)./Grille)*(Grille)')+(Impb'*ones(1,M)).*(Impb'*ones(1,M))/(2*n))./((Impb*sqrt(1/n))'*ones(1,M)),zeros(M,M),ones(M,M));
Piv = cdfnor("PQ",((ones(M,1)*(Grille*exp(pas))') - Grille*ones(1,M))/(15*sqrt(1/n)),zeros(M,M),ones(M,M))-  cdfnor("PQ",((ones(M,1)*(Grille)') - Grille*ones(1,M))/(15*sqrt(1/n)),zeros(M,M),ones(M,M));
Piv=PivBS;
vo=Impb(1,M2)
j=1;
test = input("Skip step 1?\n yes [y] or no [n]\n","string")
if test == 'n' then
  exec ./BFE_step1.sci;
else
   Volat = vo;
   Pi1 = PivBS(M2,:)';
end,
Vol = Volat*ones(1,M);
Pik1=Pi1;

Pi=hypermat([M M n-1]);
//construction de la grille
erreur=0;
//Boucle en temps
for k=2:n
  Impb=interpln([Strike';Imp(k,:)],Grille);
  PivBS=cdfnor("PQ",(log((ones(M,1)./Grille)*(Grille*exp(pas))')+(Impb'*ones(1,M)).*(Impb'*ones(1,M))/(2*n))./((Impb*sqrt(1/n))'*ones(1,M)),zeros(M,M),ones(M,M))-cdfnor("PQ",(log((ones(M,1)./Grille)*(Grille)')+(Impb'*ones(1,M)).*(Impb'*ones(1,M))/(2*n))./((Impb*sqrt(1/n))'*ones(1,M)),zeros(M,M),ones(M,M));
  disp(k,"Etape:");
  vo=interpln([Strike';Imp(k,:)],Grille);
  exec stepk.sci;
  //Preparation de la boucle suivante
  Pik1 = (Pik1' * Pik)';
  Pi(:,:,k-1)=Pik;
 
end;
//Calcul de l'americaine
//test = input("Calcul de americaine?","string")
test='n';
if test == 'y' then
  K=80:5:120;
  for i=1:9
    Prix = max(K(i)-Grille,0);
    for k=n-1:-1:1
      tmpPrix = Prix;
      clear Prix;
      Prix = max(Pi(:,:,k)*tmpPrix,max(K(i)-Grille,0));
      clear tmpPrix;
    end
    //Pas 0
    tmpPrix = Prix;
    clear Prix;
    Prix =  max(Pi1'*tmpPrix,max(K(i)-stock,0))
    clear tmpPrix;
    PrixSav(i) = Prix;
    clear Prix;
  end;
  clear tmp;
  tmp = read('amer.dat',9,2);
  VraiPrix=tmp(:,2);
  xbasc();
  plot2d(K,[VraiPrix PrixSav]);
  erreur=sqrt(sum((VraiPrix-PrixSav).*(VraiPrix-PrixSav))/sum(VraiPrix.*VraiPrix));
  disp("erreur amer",erreur);
  fic = "amer" + string(NbS);
  write(fic,PrixSav);
end,
clear VraiPrix PrixSav Prix tmpPrix;
//Calcul de la barriere
//test = input("Calcul de la barriere?","string")
test='n';
if test == 'y' then
  K=90:4:110;
  for i=1:6
    Prix = (Grille - K(i));
    for k=n-1:-1:1
      tmpPrix = Prix;
      clear Prix;
      Prix =( Pi(:,:,k)*tmpPrix).*(Grille>=85);
      clear tmpPrix;
    end
    //Pas 0
    tmpPrix = Prix;
    clear Prix;
    Prix =  (Pi1'*tmpPrix)
    clear tmpPrix;
    PrixSav(i) = Prix;
    clear Prix;
  end;
  clear tmp;
  tmp = read('barriere.dat',6,2);
  VraiPrix=tmp(:,2);
  xbasc();
  plot2d(K,[VraiPrix PrixSav]);
  erreur=sqrt(sum((VraiPrix-PrixSav).*(VraiPrix-PrixSav))/sum(VraiPrix.*VraiPrix));
  disp("erreur barriere",erreur);
  fic = "bar" + string(NbS);
  write(fic,PrixSav);
end,
