%
% MATHFI Project, Inria Rocquencourt.
% Jean-Marc Cognet, November 2002.
%
function [xmin,xmax,ymin,ymax,zmin,zmax] = valdefaut(vecx,vecy,matz)
% --------------------------------------------------------------------
% Auteur : Jean-Marc Cognet
% Date de derniere modification : 23/05/02
% Description :
%   A partir des vecteurs vecx, vecy et de la matrice matz, cette 
%   fonction donne des valeurs par defaut au nombres xmin, xmax, 
%   ymin, ymax, zmin et zmax tes que :
%     xmin et xmax correspondent a la premiere et a la derniere 
%     valeur de vecx,
%     ymin et ymax correspondent a la premiere et a la derniere 
%     valeur de vecy,
%     zmin et zmax correspondent a la valeur min et a la valeur max 
%     de la matrice matz.
%
% Appel(s) Matlab : min max
% Appel(s) : 
% --------------------------------------------------------------------
%
xmin = vecx(1);
xmax = vecx(length(vecx));
%
ymin = vecy(1);
ymax = vecy(length(vecy));
%
zmin = min(min(matz));
zmax = max(max(matz));
