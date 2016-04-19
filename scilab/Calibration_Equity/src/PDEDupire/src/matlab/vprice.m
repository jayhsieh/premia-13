%
% MATHFI Project, Inria Rocquencourt.
% Jean-Marc Cognet, November 2002.
%
% VALEURS PAR DEFAUT DES PARAMETRES :
% strfic = 'sigmaest.visu';  nom du fichier contenant les prix
% option = 2;                1 --> projection 2D, 2 --> surface, 3 --> coupe
% coul = 'y';                'n' pour noir et blanc
% langage = 'f';             'e' pour english
%
% kmin, kmax, tmin, tmax, vmin, vmax peuvent etre definis avant
% l'appel de vprice. Sinon (i.e si ils n'existent pas), ils 
% prennent des valeurs par defaut en fonction de veck, veckt et matv
%
% EXEMPLES :
% clear ; strfic = 'sigmaest.visu';
% kmin = 75 ; kmax = 125 ; option = 1 ; vprice
%
if (exist('strfic')==0) % si strfic n'existe pas
  strfic = 'sigmaest.visu';
end
if (exist('option')==0)
  option = 2; 
end
if (exist('coul')==0)
  coul = 'y';
end
if (exist('langage')==0)
  langage = 'f';
end
%
valdef = -9999;
%
if (exist('kmin')==0)
  kmin = valdef;
end
if (exist('kmax')==0)
  kmax = valdef;
end
if (exist('tmin')==0)
  tmin = valdef;
end
if (exist('tmax')==0)
  tmax = valdef;
end
if (exist('vmin')==0)
  vmin = valdef;
end
if (exist('vmax')==0)
  vmax = valdef;
end
%
if (exist('n')==0 & exist('m')==0 & exist('veck')==0 & ...
      exist('vect')==0 & exist('matv')==0)
  [n,m,veck,vect,matv] = recup(strfic);
end
%
hdl = visu(veck,vect,matv,kmin,kmax,tmin,tmax,vmin,vmax,option,coul,langage);
