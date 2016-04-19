function [hdl]=visu(vecx,vecy,matz,xmin,xmax,ymin,ymax,zmin,zmax,option,coul,langage)
hdl=[];
// 
// MATHFI Project, Inria Rocquencourt.
// Jean-Marc Cognet, November 2002.
// 
// --------------------------------------------------------------------
// Auteur : Jean-Marc Cognet
// Date de derniere modification : 12/07/02
// Description :
//   A partir de : 
//   - vecx : vecteur d'abscisses X
//   - vecy : vecteur d'ordonnees Y
//   - matz : valeurs de la fonction pour differents couples (X,Y)
//   - xmin, xmax : abscisses min et max telles que : 
//     min(vecx) <= xmin <= xmax <= max(vecx)
//   - ymin, ymax : ordonnees min et max telles que : 
//     min(vecy) <= ymin <= ymax <= max(vecy)
//   - zmin, zmax : valeurs min et max de la fonction telles que : 
//     min(z) <= zmin <= zmax <= max(z)
//   - option : 1 --> projection 2D, 2 --> surface, 3 --> coupe
//   - coul : 'n' --> noir et blanc
//   - langage : 'e' --> english
//   cette fonction permet de visualiser les valeurs de z sous 
//   differentes formes : une projection 2D, une surface ou une coupe 
//   dans une direction.
//   En sortie, on recupere le handle hdl de la figure.
// 
// Appel(s) Matlab : figure pcolor caxis colormap colorbar shading
// Appel(s) Matlab : xlabel ylabel surf min max axis error
// Appel(s) : trvind
// --------------------------------------------------------------------
// 
// Legendes en francais et en anglais
// 
kleg_fr = 'S';
kleg_en = 'S';
tleg_fr = 't (an)';
tleg_en = 't (year)';
// 
// Definition de ixmin, ixmax, iymin et iymax
// 
valdef = -9999;
// 
if xmin==valdef then
  ixmin = 1;
else
   
  //!! Unknown function trvind ,the original calling sequence is used
  ixmin = trvind(vecx,xmin,0);
end
if xmax==valdef then
   
  //! unknown arg type, using mtlb_length 
  ixmax = mtlb_length(vecx);
else
   
  //!! Unknown function trvind ,the original calling sequence is used
  ixmax = trvind(vecx,xmax,1);
end
if ymin==valdef then
  iymin = 1;
else
   
  //!! Unknown function trvind ,the original calling sequence is used
  iymin = trvind(vecy,ymin,0);
end
if ymax==valdef then
   
  //! unknown arg type, using mtlb_length 
  iymax = mtlb_length(vecy);
else
   
  //!! Unknown function trvind ,the original calling sequence is used
  iymax = trvind(vecy,ymax,1);
end
// 
// Extraction de exvecx, exvecy et exmatz
// 
exvecx = vecx(ixmin:ixmax);
exvecy = vecy(iymin:iymax);
exmatz = matz(ixmin:ixmax,iymin:iymax);
// 
// Definition de zmin et zmax
// 
if zmin==valdef then
   
  //! mtlb_min(exmatz) may be replaced by 
  //!   min(exmatz) if exmatzis a vector
  //!   min(exmatz,'r') if exmatzis a matrix
   
  //! mtlb_min(mtlb_min(exmatz)) may be replaced by 
  //!   min(mtlb_min(exmatz)) if mtlb_min(exmatz)is a vector
  //!   min(mtlb_min(exmatz),'r') if mtlb_min(exmatz)is a matrix
  fmin = mtlb_min(mtlb_min(exmatz));
else
  fmin = zmin;
end
if zmax==valdef then
   
  //!  mtlb_max(exmatz) may be replaced by 
  //!     max(exmatz) if exmatzis a vector
  //!     max(exmatz,'r') if exmatzis a matrix
   
  //!  mtlb_max(mtlb_max(exmatz)) may be replaced by 
  //!     max(mtlb_max(exmatz)) if mtlb_max(exmatz)is a vector
  //!     max(mtlb_max(exmatz),'r') if mtlb_max(exmatz)is a matrix
  fmax = mtlb_max(mtlb_max(exmatz));
else
  fmax = zmax;
end
// 
// Couleur ou N&B ? si coul=='n', on definit la matrice co
// 
if coul=='n' then
  nc = 16;
  xco = 1/(nc-1);
  nt = 16;
  for i = 1:nc
    j = nt*(i-1);
    for k = 1:nt
      co(1,k+j) = 1-xco*(i-1);
      co(2,k+j) = 1-xco*(i-1);
      co(3,k+j) = 1-xco*(i-1);
    end
  end
  co = co';
end
// 
// Cas ou fmin=fmax
// 
if fmin==fmax then
  fmin = fmin*0.8;
  fmax = fmax*1.2;
end
// 
// Visualisation en fonction de :
// option : 1 --> projection 2D, 2 --> surface, 3 --> coupe
// coul : 'n' --> noir et blanc
// langage : 'e' --> english
// 
xset('window',max(winsid()+1))
max(winsid());
// 
exxmin = min(exvecx);
exxmax = max(exvecx);
exymin = min(exvecy);
exymax = max(exvecy);
if option==1 then
  // projection 2D avec pcolor
   
  //!! Unknown function pcolor ,the original calling sequence is used
  hdl = pcolor(exvecx',exvecy,exmatz');
   
  //!! Unknown function caxis ,the original calling sequence is used
  caxis([fmin,fmax]);
  if coul=='n' then
     
    //! Not enough information using mtlb_colormap
    //! instead of xset('colormap',..
    mtlb_colormap(co);
  end
   
  //!! Unknown function colorbar ,the original calling sequence is used
  colorbar();
   
  //!! Unknown function shading ,the original calling sequence is used
  shading('flat');
  if langage=='e' then
    xtitle(' ',kleg_en,' ');
    xtitle(' ',' ',tleg_en);
  else
    xtitle(' ',kleg_fr,' ');
    xtitle(' ',' ',tleg_fr);
  end
elseif option==2 then
  // surface avec surf
   
  //!! Unknown function surf ,the original calling sequence is used
  hdl = surf(exvecx',exvecy,exmatz');
   
  //!! Unknown function axis ,the original calling sequence is used
  axis([exxmin,exxmax,exymin,exymax,fmin,fmax]);
   
  //!! Unknown function caxis ,the original calling sequence is used
  caxis([fmin,fmax]);
  if coul=='n' then
     
    //! Not enough information using mtlb_colormap
    //! instead of xset('colormap',..
    mtlb_colormap(co);
  end
   
  //!! Unknown function colorbar ,the original calling sequence is used
  colorbar();
   
  //!! Unknown function shading ,the original calling sequence is used
  shading('flat');
  if langage=='e' then
    xtitle(' ',kleg_en,' ');
    xtitle(' ',' ',tleg_en);
  else
    xtitle(' ',kleg_fr,' ');
    xtitle(' ',' ',tleg_fr);
  end
elseif option==3 then
  // coupe avec plot
  k_ou_t = input('\nVariable fixee : k ou t ? ','s');
  if k_ou_t=='k' then
    disp(mtlb_sprintf('Coupe en une valeur proche de (entrer un nombre compris entre %5.2f et %5.2f)',exxmin,exxmax));
    val_k_ou_t = input('  k = ');
    exvec_k_ou_t = exvecy;
     
    //!! Unknown function trvind ,the original calling sequence is used
    indk = trvind(exvecx,val_k_ou_t,0);
    disp(mtlb_sprintf('\nVoici une coupe en k = %5.2f',exvecx(indk)));
    exvecz = exmatz(indk,:)';
    if langage=='e' then
      legabs = tleg_en;
    else
      legabs = tleg_fr;
    end
    titre = mtlb_sprintf('k = %5.2f',exvecx(indk));
    absmin = exymin;
    absmax = exymax;
    ordmin = fmin;
    ordmax = fmax;
  elseif k_ou_t=='t' then
    disp(mtlb_sprintf('Coupe en une valeur proche de (entrer un nombre compris entre %5.2f et %5.2f)',exymin,exymax));
    val_k_ou_t = input('  t = ');
    exvec_k_ou_t = exvecx;
     
    //!! Unknown function trvind ,the original calling sequence is used
    indt = trvind(exvecy,val_k_ou_t,0);
    disp(mtlb_sprintf('\nVoici une coupe en t = %5.2f\n',exvecy(indt)));
    exvecz = exmatz(:,indt);
    if langage=='e' then
      legabs = kleg_en;
    else
      legabs = kleg_fr;
    end
    titre = mtlb_sprintf('t = %5.2f',exvecy(indt));
    absmin = exxmin;
    absmax = exxmax;
    ordmin = fmin;
    ordmax = fmax;
  else
    error('Taper la lettre k ou la lettre t');
  end
  // 
  hdl = mtlb_plot(exvec_k_ou_t,exvecz);
  // 
  xtitle(' ',legabs,' ');
  xtitle(titre,' ',' ');
   
  //!! Unknown function axis ,the original calling sequence is used
  axis([absmin,absmax,fmin,fmax]);
else
  error('option doit etre egal a 1, 2 ou 3');
end
