%
% MATHFI Project, Inria Rocquencourt.
% Jean-Marc Cognet, November 2002.
%
function [hdl] = visu(vecx,vecy,matz,xmin,xmax,ymin,ymax,zmin,zmax,...
    option,coul,langage)
% --------------------------------------------------------------------
% Auteur : Jean-Marc Cognet
% Date de derniere modification : 12/07/02
% Description :
%   A partir de : 
%   - vecx : vecteur d'abscisses X
%   - vecy : vecteur d'ordonnees Y
%   - matz : valeurs de la fonction pour differents couples (X,Y)
%   - xmin, xmax : abscisses min et max telles que : 
%     min(vecx) <= xmin <= xmax <= max(vecx)
%   - ymin, ymax : ordonnees min et max telles que : 
%     min(vecy) <= ymin <= ymax <= max(vecy)
%   - zmin, zmax : valeurs min et max de la fonction telles que : 
%     min(z) <= zmin <= zmax <= max(z)
%   - option : 1 --> projection 2D, 2 --> surface, 3 --> coupe
%   - coul : 'n' --> noir et blanc
%   - langage : 'e' --> english
%   cette fonction permet de visualiser les valeurs de z sous 
%   differentes formes : une projection 2D, une surface ou une coupe 
%   dans une direction.
%   En sortie, on recupere le handle hdl de la figure.
%
% Appel(s) Matlab : figure pcolor caxis colormap colorbar shading
% Appel(s) Matlab : xlabel ylabel surf min max axis error
% Appel(s) : trvind
% --------------------------------------------------------------------
%
% Legendes en francais et en anglais
%
kleg_fr = 'S';
kleg_en = 'S';
tleg_fr = 't (an)';
tleg_en = 't (year)';
%
% Definition de ixmin, ixmax, iymin et iymax
%
valdef = -9999;
%
if xmin==valdef
  ixmin = 1;
else
  ixmin = trvind(vecx,xmin,0);
end
if xmax==valdef
  ixmax = length(vecx);
else
  ixmax = trvind(vecx,xmax,1);
end
if ymin==valdef
  iymin = 1;
else
  iymin = trvind(vecy,ymin,0);
end
if ymax==valdef
  iymax = length(vecy);
else
  iymax = trvind(vecy,ymax,1);
end
%
% Extraction de exvecx, exvecy et exmatz
%
exvecx = vecx(ixmin:ixmax);
exvecy = vecy(iymin:iymax);
exmatz = matz(ixmin:ixmax,iymin:iymax);
%
% Definition de zmin et zmax
%
if zmin==valdef
  fmin = min(min(exmatz));
else
  fmin = zmin;
end
if zmax==valdef
  fmax = max(max(exmatz));
else
  fmax = zmax;
end
%
% Couleur ou N&B ? si coul=='n', on definit la matrice co
%
if coul=='n'
  nc = 16;
  xco = 1.0/(nc-1);
  nt = 16;
  for i = 1:nc
    j = nt*(i-1);
    for k = 1:nt
      co(1,k+j) = 1- xco*(i-1);
      co(2,k+j) = 1- xco*(i-1);
      co(3,k+j) = 1- xco*(i-1);
    end
  end
  co = co';
end
%
% Cas ou fmin=fmax
%
if fmin==fmax
  fmin=fmin*0.8;
  fmax=fmax*1.2;
end
%
% Visualisation en fonction de :
% option : 1 --> projection 2D, 2 --> surface, 3 --> coupe
% coul : 'n' --> noir et blanc
% langage : 'e' --> english
%
figure;
%
exxmin = min(exvecx);exxmax = max(exvecx);
exymin = min(exvecy);exymax = max(exvecy);
if option==1 % projection 2D avec pcolor
  hdl = pcolor(exvecx',exvecy,exmatz');
  caxis([fmin fmax])
  if coul=='n'
    colormap(co);
  end
  colorbar;
  shading flat ;
  if langage=='e'
    xlabel(kleg_en)
    ylabel(tleg_en)
  else
    xlabel(kleg_fr)
    ylabel(tleg_fr)    
  end
elseif option==2 % surface avec surf
  hdl = surf(exvecx',exvecy,exmatz');
  axis([exxmin exxmax exymin exymax fmin fmax])
  caxis([fmin fmax])
  if coul=='n'
    colormap(co);
  end
  colorbar;
  shading flat ;
  if langage=='e'
    xlabel(kleg_en)
    ylabel(tleg_en)
  else
    xlabel(kleg_fr)
    ylabel(tleg_fr)    
  end
elseif option==3 % coupe avec plot
  k_ou_t = input('\nVariable fixee : k ou t ? ','s');
  if k_ou_t == 'k'
    disp(sprintf('Coupe en une valeur proche de (entrer un nombre compris entre %5.2f et %5.2f)',exxmin,exxmax));
    val_k_ou_t = input('  k = ');
    exvec_k_ou_t = exvecy;
    indk = trvind(exvecx,val_k_ou_t,0);
    disp(sprintf('\nVoici une coupe en k = %5.2f',exvecx(indk)));
    exvecz = exmatz(indk,:)';
    if langage=='e'
      legabs = tleg_en;
    else
      legabs = tleg_fr;
    end
    titre = sprintf('k = %5.2f',exvecx(indk));
    absmin = exymin ; absmax = exymax;
    ordmin = fmin ; ordmax = fmax;
  elseif k_ou_t == 't'
    disp(sprintf('Coupe en une valeur proche de (entrer un nombre compris entre %5.2f et %5.2f)',exymin,exymax));
    val_k_ou_t = input('  t = ');    
    exvec_k_ou_t = exvecx;
    indt = trvind(exvecy,val_k_ou_t,0);
    disp(sprintf('\nVoici une coupe en t = %5.2f\n',exvecy(indt)));
    exvecz = exmatz(:,indt);
    if langage=='e'
      legabs = kleg_en;
    else
      legabs = kleg_fr;
    end
    titre = sprintf('t = %5.2f',exvecy(indt));
    absmin = exxmin ; absmax = exxmax;
    ordmin = fmin ; ordmax = fmax;
  else
    error('Taper la lettre k ou la lettre t');
  end
  %
  hdl = plot(exvec_k_ou_t,exvecz);
  %
  xlabel(legabs);
  title(titre);
  axis([absmin absmax fmin fmax])
else
  error('option doit etre egal a 1, 2 ou 3')
end
