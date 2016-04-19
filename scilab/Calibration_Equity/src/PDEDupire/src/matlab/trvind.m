%
% MATHFI Project, Inria Rocquencourt.
% Jean-Marc Cognet, November 2002.
%
function ind = trvind(vec,val,mode)
% --------------------------------------------------------------------
% Auteur : Jean-Marc Cognet
% Date de derniere modification : 23/05/02
% Description :
%   A partir du vecteur vec, du nombre val et du mode (0 ou 1), 
%   cette fonction permet de determiner un indice ind tel que :
%   si mode = 0 (min) : ind correspond au plus grand indice 
%                       tel que vec(ind) <= val
%   si mode = 1 (max) : ind correspond au plus petit indice 
%                       tel que vec(ind) >= val
%   On suppose au prealable que le vecteur vec contient des nombres 
%   ranges par ordre croissant et que : min(vec) <= val <= max(vec), 
%   i.e que vec(1) <= val <= vec(length(vec)).
%
% Appel(s) Matlab : length min max error break
% Appel(s) : 
% --------------------------------------------------------------------
%
l = length(vec);
%
minvec = vec(1);
maxvec = vec(l);
%
if val<minvec | val>maxvec
  error('Probleme avec la fonction trvind : vecmin <= val <= vecmax')
end
%
if mode==0 % min
  ind = 1;
  while ind < l
    if vec(ind)>val
      ind = ind - 1;
      break
    else
      ind = ind+1;
    end
  end
elseif mode==1 % max
  ind = l;
  while ind > 1
    if vec(ind)<val
      ind = ind + 1;
      break
    else
      ind = ind-1;
    end
  end
else
  error('Probleme avec la fonction trvind : mode = 0 ou 1')
end
