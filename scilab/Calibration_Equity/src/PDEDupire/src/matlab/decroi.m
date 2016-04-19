%
% MATHFI Project, Inria Rocquencourt.
% Jean-Marc Cognet, November 2002.
%
% Visualisation de la decroissance du critere a partir du fichier ficj
%
% EXEMPLE :
% ficj = 'sigmaest_0.out.j' ; decroi

legx_fr = 'Nombre d''itérations';
legy_fr = 'Fonction coût';
legx_en = 'Number of iterations';
legy_en = 'Cost function';
format  = '%6.4f'; % 4.2f

if (exist('ficj')==0) % si ficj n'existe pas
  ficj = input('\nEntrer le nom du fichier : ','s');
  disp(' ');
end

% Exemple :
% ficj = 'sigmaest_nm1_0.out.j';
% tabj = load('sigmaest_nm1_0.out.j');
eval(['tabj = load(''',ficj,''');']);

number_iter = tabj(:,1);
val_j       = tabj(:,2);

figure

plot(number_iter,val_j)

xlabel(legx_fr);
ylabel(legy_fr);
jinit = val_j(1);
jest  = val_j(length(val_j));
strtitle_j = ['F : ',sprintf(format,jinit),' --> ',sprintf(format,jest),'.'];
  title(strtitle_j);
