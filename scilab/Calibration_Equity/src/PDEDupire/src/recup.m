%
% MATHFI Project, Inria Rocquencourt.
% Jean-Marc Cognet, November 2002.
%
function [nx,ny,vecx,vecy,matz] = recup(sfic)
% --------------------------------------------------------------------
% Auteur : Jean-Marc Cognet
% Date de derniere modification : 23/05/02
% Description :
%   A partir de la chaine de caracteres sfic contenant le nom d'un 
%   fichier ASCII, cette fonction permet de recuperer les nombres 
%   entiers nx et ny, les vecteurs vecx(1:nx+1) et vecy(1:ny+1), 
%   ainsi que la matrice matz(1:nx+1,1:ny+1) avec : 
%     matz(i,j) = z(vecx(i),vecy(j)) pour i=1:nx+1 et j=1:ny+1.
%   On suppose que les (nx+2)(ny+2)+1 nombres contenus dans le 
%   fichier sfic sont ranges de la facon suivante :
%          nx
%          ny
%          vecx(1)
%          ...
%          vecx(nx+1)
%          vecy(1)
%          ...
%          vecy(ny+1)
%          matz(1,1)
%          ...
%          matz(1,ny+1)
%          matz(2,1)
%          ...
%          matz(2,ny+1)
%          ...
%          ...
%          matz(nx+1,1)
%          ...
%          matz(nx+1,ny+1)
%
% Appel(s) Matlab : eval load break
% Appel(s) : 
% --------------------------------------------------------------------
%
% Obtention du vecteur vecdata a partir du fichier ASCII sfic
%
eval(['load ',sfic,' -ascii ;']);
for ind = 1:length(sfic)
  if sfic(ind)=='.'
    break
  else
    temporaire(ind) = sfic(ind);
  end
end
eval(['vecdata = ',temporaire,';']);
%
% Recuperation de nx, ny, vecx, vecy et matz a partir du vecteur 
% vecdata
%
nx = vecdata(1);
ny = vecdata(2);
%
vecx = vecdata(3:nx+3);
vecy = vecdata(nx+4:nx+ny+4);
%
nbr = nx+ny+4;
for i=1:nx+1
  for j=1:ny+1
    matz(i,j) = vecdata(nbr+(i-1)*(ny+1)+j);
  end
end
