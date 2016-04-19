//prod_txt = ['No companies';'Nominal';'Maturite';'Frequence';'Tranche(s)'];
Method = ['Copula';'Implied_method']
Method = x_choose(Method, 'Which method ?')
prod_text=[];
prod_def=[];

select Method, 
case 1 then 
prod_txt = ['No companies';'Nominal';'Maturite';'Frequence';'Tranche(s)'];
if ((exists('product'))&(size(product,'*')==5)) then
    prod_def = ['100';'0.01';'5';'0.25';'[0;0.03;0.06;0.1;1]'];
else
    prod_def = ['100';'0.01';'5';'0.25';'[0;0.03;0.06;0.1;1]'];
end;
product = x_mdialog('Produit', prod_txt, prod_def);
if (product == []) then
    product = prod_def;
    abort;
else 
    n_comp = evstr(product(1));
    test_nom = execstr('evstr(product(2))','errcatch');
    if test_nom > 0 then
	nominal = read('../datas/'+product(2), n_comp, 1);
    else
	nominal = evstr(product(2));
    end;
    if (size(nominal,'*') == 1) then nominal=nominal*ones(n_comp,1); end;
    dates = [evstr(product(4)):evstr(product(4)):evstr(product(3))]';
    tranches = evstr(product(5));
end;

function print_product()
  printf("Product:\n");
  printf("   Number of companies: %s\n", prod_def(1));
  printf("   Nominal of each company: %s\n", prod_def(2));
  printf("   Maturity: %s\n", prod_def(3));
  printf("   Time step: %s\n", prod_def(4));
endfunction

case 2 then 
prod_txt = ['No companies';'Nominal';'Maturite';'Frequence';'Tranche(s)'];
if ((exists('product'))&(size(product,'*')==5)) then
    prod_def = ['125';'0.008';'5';'0.25';'[0;0.03]'];//prod_def = product;
else
    prod_def = ['125';'0.008';'5';'0.25';'[0;0.03]'];
end;
product = x_mdialog('Produit', prod_txt, prod_def);
if (product == []) then
    product = prod_def;
    abort;
else 
    n_comp = evstr(product(1));
    test_nom = execstr('evstr(product(2))','errcatch');
    if test_nom > 0 then
	nominal = read('../datas/'+product(2), n_comp, 1);
    else
	nominal = evstr(product(2));
    end;
    if (size(nominal,'*') == 1) then nominal=nominal*ones(n_comp,1); end;
    dates = [evstr(product(4)):evstr(product(4)):evstr(product(3))]';
    tranches = evstr(product(5));
end;

function print_product()
  printf("Product:\n");
  printf("   Number of companies: %s\n", prod_def(1));
  printf("   Nominal of each company: %s\n", prod_def(2));
  printf("   Maturity: %s\n", prod_def(3));
  printf("   Time step: %s\n", prod_def(4));
endfunction


end 