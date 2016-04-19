prod_txt = ['No companies';'Nominal';'Maturité';'Fréquence';'Tranche(s)'];
if exists('product') && size(product,'*')==5 then
    prod_def = product;
else
    prod_def = ['100';'0.01';'5';'0.25';'[0;0.03;0.06;0.1;1]'];
end;
prod_def = ['100';'0.01';'5';'0.25';'[0;0.03;0.06;0.1;1]'];
product = x_mdialog('Produit', prod_txt, prod_def);
if (product == []) then
    product = prod_def;
    abort;
else 
    n_comp = evstr(product(1));
    test_nom = execstr('evstr(product(2))',errcatch=%t);
    if ~test_nom then
	nominal = read(product(2), n_comp, 1);
    else
	nominal = evstr(product(2));
    end;
    if (size(nominal,'*') == 1) then nominal=nominal*ones(n_comp,1); end;
    dates = [evstr(product(4)):evstr(product(4)):evstr(product(3))]';
    tranches = evstr(product(5));
end;
