%
% MATHFI Project, Inria Rocquencourt.
% Jean-Marc Cognet, November 2002.
%
% VALEURS PAR DEFAUT DES PARAMETRES :
% strfic1 = 'sigmatrue.visu';  nom du premier fichier contenant les prix
% strfic2 = 'sigmaest.visu';   nom du second fichier contenant les prix
% option = 2;                  1 --> projection 2D, 2 --> surface, 3 --> coupe
% coul = 'y';                  'n' pour noir et blanc
% langage = 'f';               'e' pour english
% erreur = 1;                  1 --> absolue, 2 --> relative
%
% kmin, kmax, tmin, tmax, vmin, vmax peuvent etre definis avant
% l'appel de vprice. Sinon (i.e si ils n'existent pas), ils 
% prennent des valeurs par defaut en fonction de veck, veckt et matv
%
% EXEMPLES :
% 
% clear ; strfic1 = 'sigmatrue.visu' ; strfic2 = 'sigmaest0.visu';
% kmin = 75 ; kmax = 125 ; option = 1 ; erreur = 2 ; dprice
% 
if (exist('strfic1')==0) % si strfic1 n'existe pas
  strfic1 = 'sigmatrue.visu';
end
if (exist('strfic2')==0) % si strfic n'existe pas
  strfic2 = 'sigmaest0.visu';
end
if (exist('option')==0) % si strfic n'existe pas
  option = 2;
end
if (exist('coul')==0) % si strfic n'existe pas
  coul = 'y';
end
if (exist('langage')==0) % si strfic n'existe pas
  langage = 'f';
end
if (exist('erreur')==0) % si strfic n'existe pas
  erreur = 1;
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
if (exist('n1')==0 & exist('m1')==0 & exist('veck1')==0 & ...
      exist('vect1')==0 & exist('matv1')==0)
  [n1,m1,veck1,vect1,matv1] = recup(strfic1);
end
if (exist('n2')==0 & exist('m2')==0 & exist('veck2')==0 & ...
      exist('vect2')==0 & exist('matv2')==0)
  [n2,m2,veck2,vect2,matv2] = recup(strfic2);
end
%
if norm(veck1-veck2)~=0
  error('Probleme avec veck');
elseif norm(vect1-vect2)~=0
  error('Probleme avec vect');
else
  veck = veck1;
  vect = vect1;
end
%
if erreur == 1
  matv = abs(matv1 - matv2);
else
  epsilon = 1e-6;
  for i = 1:size(matv1,1)
    for j = 1:size(matv1,2)
      %if matv1(i,j) ~= 0
      if abs(matv1(i,j)) > epsilon
	matv(i,j) = 100 * abs(matv1(i,j) - matv2(i,j)) ./ abs(matv1(i,j));
      else
	matv(i,j) = 100;
      end
    end
  end
end
%
visu(veck,vect,matv,kmin,kmax,tmin,tmax,vmin,vmax,option,coul,langage);
