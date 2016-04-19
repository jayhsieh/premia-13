function []=vprice()
// 
// MATHFI Project, Inria Rocquencourt.
// Jean-Marc Cognet, November 2002.
// 
// VALEURS PAR DEFAUT DES PARAMETRES :
// strfic = 'sigmaest.visu';  nom du fichier contenant les prix
// option = 2;                1 --> projection 2D, 2 --> surface, 3 --> coupe
// coul = 'y';                'n' pour noir et blanc
// langage = 'f';             'e' pour english
// 
// kmin, kmax, tmin, tmax, vmin, vmax peuvent etre definis avant
// l'appel de vprice. Sinon (i.e si ils n'existent pas), ils 
// prennent des valeurs par defaut en fonction de veck, veckt et matv
// 
// EXEMPLES :
// clear ; strfic = 'sigmaest.visu';
// kmin = 75 ; kmax = 125 ; option = 1 ; vprice
// 
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('strfic')==0 then
  // si strfic n'existe pas
  strfic = 'dupiresprices.visu';
end
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('option')==0 then
  option = 2;
end
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('coul')==0 then
  coul = 'y';
end
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('langage')==0 then
  langage = 'f';
end
// 
valdef = -9999;
// 
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('kmin')==0 then
  kmin = valdef;
end
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('kmax')==0 then
  kmax = valdef;
end
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('tmin')==0 then
  tmin = valdef;
end
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('tmax')==0 then
  tmax = valdef;
end
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('vmin')==0 then
  vmin = valdef;
end
 
//! Not enough information using mtlb_exist instead of exists
if mtlb_exist('vmax')==0 then
  vmax = valdef;
end
// 
 
//! Not enough information using mtlb_exist instead of exists
 
//! Not enough information using mtlb_exist instead of exists
 
//! Not enough information using mtlb_exist instead of exists
 
//! Not enough information using mtlb_exist instead of exists
 
//! Not enough information using mtlb_exist instead of exists
if ((((mtlb_exist('n')==0)&(mtlb_exist('m')==0))&(mtlb_exist('veck')==0))&(mtlb_exist('vect')==0))&(mtlb_exist('matv')==0) then
   
  //!! Unknown variable recup ,the original calling sequence is used
  [n,m,veck,vect,matv] = recup(strfic);
end
// 
hdl = visu(veck,vect,matv,kmin,kmax,tmin,tmax,vmin,vmax,option,coul,langage);
