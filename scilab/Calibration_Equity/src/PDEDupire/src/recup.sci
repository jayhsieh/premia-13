function [nx,ny,vecx,vecy,matz]=recup(sfic)
nx=[];ny=[];vecx=[];vecy=[];matz=[];
// 
// MATHFI Project, Inria Rocquencourt.
// Jean-Marc Cognet, November 2002.
// 
// --------------------------------------------------------------------
// Auteur : Jean-Marc Cognet
// Date de derniere modification : 23/05/02
// Description :
//   A partir de la chaine de caracteres sfic contenant le nom d'un 
//   fichier ASCII, cette fonction permet de recuperer les nombres 
//   entiers nx et ny, les vecteurs vecx(1:nx+1) et vecy(1:ny+1), 
//   ainsi que la matrice matz(1:nx+1,1:ny+1) avec : 
//     matz(i,j) = z(vecx(i),vecy(j)) pour i=1:nx+1 et j=1:ny+1.
//   On suppose que les (nx+2)(ny+2)+1 nombres contenus dans le 
//   fichier sfic sont ranges de la facon suivante :
//          nx
//          ny
//          vecx(1)
//          ...
//          vecx(nx+1)
//          vecy(1)
//          ...
//          vecy(ny+1)
//          matz(1,1)
//          ...
//          matz(1,ny+1)
//          matz(2,1)
//          ...
//          matz(2,ny+1)
//          ...
//          ...
//          matz(nx+1,1)
//          ...
//          matz(nx+1,ny+1)
// 
// Appel(s) Matlab : eval load break
// Appel(s) : 
// --------------------------------------------------------------------
// 
// Obtention du vecteur vecdata a partir du fichier ASCII sfic
// 
 
//! mtlb_eval can be replaced by evstr if 'load '+sfic+' -ascii ;'
//! is a valid scilab instruction
mtlb_eval('load '+sfic+' -ascii ;');
 
//! unknown arg type, using mtlb_length 
for ind = 1:mtlb_length(sfic)
  if sfic(ind)=='.' then
    break
     ;
     
  else
    temporaire(1,ind) = sfic(ind);
  end
end
 
//! mtlb_eval can be replaced by evstr if 'vecdata = '+temporaire+';'
//! is a valid scilab instruction
mtlb_eval('vecdata = '+temporaire+';');
// 
// Recuperation de nx, ny, vecx, vecy et matz a partir du vecteur 
// vecdata
// 
 
//!! Unknown function vecdata ,the original calling sequence is used
nx = vecdata(1);
 
//!! Unknown function vecdata ,the original calling sequence is used
ny = vecdata(2);
// 
 
//!! Unknown function vecdata ,the original calling sequence is used
vecx = vecdata(3:nx+3);
 
//!! Unknown function vecdata ,the original calling sequence is used
vecy = vecdata(nx+4:nx+ny+4);
// 
nbr = nx+ny+4;
for i = 1:nx+1
  for j = 1:ny+1
     
    //!! Unknown function vecdata ,the original calling sequence is used
    matz(i,j) = vecdata(nbr+(i-1)*(ny+1)+j);
  end
end
