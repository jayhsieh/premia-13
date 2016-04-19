printf("Computing...\n");
[price, pl, dl] = price_cdo(n_comp, nominal, dates, evstr( product(5)), intensity, [0.,5.], [0.03,0.15], model1(2), evstr(model2(2)), model1(4), evstr(model2(4)), method, [10, 100, 10000]);

if (method <= 2) 
    setmenu('CDO',5);
else 
    if (method <= 5)
	setmenu('CDO',6);
    end;
end;

nt=size(tranches,'*');
print([tranches(1:nt-1), tranches(2:nt), price, dl, pl]);

