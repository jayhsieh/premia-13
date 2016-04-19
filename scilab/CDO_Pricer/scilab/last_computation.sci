//winId=progressionbar('Calcul en cours...');


printf("***************************************************\n");
print_product();
print_model();
print_method();
printf("***************************************************\n");
printf("Computing...\n");

select Method,

case 1 then 

[price, dl, pl] = price_cdo(n_comp, nominal, dates, evstr(product(5)), intensity, taux(:,1), taux(:,2), evstr(model1(2)), evstr(model2(2)), evstr(model1(4)), evstr(model2(4)), method, evstr(params_method));
//winclose(winId);

if (method <= 2) 
    setmenu('CDO',5);
else 
    if (method <= 5)
	setmenu('CDO',6);
    end;
end;

nt=size(tranches,'*');

if (method <= 5) 
// Affichage
for (i=1:nt-1),
    printf("Tranche %i: %.2f\% - %.2f\% \n", i, tranches(i), tranches(i+1));
    printf("         Price: %.5f \n", price(i));
    printf("   Default Leg: %.5f \n", dl(i));
    printf("   Payment Leg: %.5f \n\n", pl(i));
end;

printf("Results are in variables: price, dl, and pl.\n");
else 
   if(method==6)
// Affichage pour Monte-Carlo avec IC
for (i=1:nt-1),
    printf("Tranche %i: %.2f\% - %.2f\% \n", i, tranches(i), tranches(i+1));
    icpl = 1.96 * sqrt(pl(i+nt-1) / evstr(params_method(1)));
    icdl = 1.96 * sqrt(dl(i+nt-1) / evstr(params_method(1)));
    printf("         Price: %.5f \t[%.5f, %.5f]\n", price(i), ((dl(i)-icdl)/(pl(i)+icpl))*10000, ((dl(i)+icdl)/(pl(i)-icpl))*10000);
    printf("   Default Leg: %.5f \t[%.5f, %.5f]\n", dl(i), dl(i)-icdl, dl(i)+icdl);
    printf("   Payment Leg: %.5f \t[%.5f, %.5f]\n\n", pl(i), pl(i)-icpl, pl(i)+icpl);
end
 
else 
   if(method==7)
// Affichage pour Monte-Carlo avec IC
for (i=1:nt-1),
    printf("Tranche %i: %.2f\% - %.2f\% \n", i, tranches(i), tranches(i+1));
    icpl = 1.96 * sqrt(pl(i+nt-1) / evstr(params_method(1)));
    icdl = 1.96 * sqrt(dl(i+nt-1) / evstr(params_method(1)));
    printf("         Price: %.5f \t[%.5f, %.5f]\n", price(i), ((dl(i)-icdl)/(pl(i)+icpl))*10000, ((dl(i)+icdl)/(pl(i)-icpl))*10000);
    printf("   Default Leg: %.5f \t[%.5f, %.5f]\n", dl(i), dl(i)-icdl, dl(i)+icdl);
    printf("   Payment Leg: %.5f \t[%.5f, %.5f]\n\n", pl(i), pl(i)-icpl, pl(i)+icpl);

end;

else 
   if(method==8)

for (i=1:nt-1),
    printf("Tranche %i: %.2f\% - %.2f\% \n", i, tranches(i), tranches(i+1));
    printf("         Price: %.5f \n", price(i));
    printf("   Default Leg: %.5f \n", dl(i));
    printf("   Payment Leg: %.5f \n\n", pl(i));
end;

printf("Results are in variables: price, pl, and dl.\n");
  end;
  end;
 end;
end;

mode(0);


case 2 then   

select model1(4),

case 1 then 

method=9;

[price, dl, pl] = price_cdo(n_comp, nominal, dates, evstr(product(5)), intensity, taux(:,1), taux(:,2), evstr(model1(2)), evstr(model2(2)), evstr(model1(4)), evstr(model2(4)), method, evstr(params_method));
//winclose(winId);


nt=size(tranches,'*');
x=ones(5,1);
y=ones(5,1);

 for (i=1:1),
  if(tranches(i+1)<=0.03) then  
    printf("Tranche %i: %.2f\% - %.2f\% \n", i, tranches(i), tranches(i+1));
    printf("         Price: %.5f \n",100*(dl(i)-0.05*pl(i))/((tranches(i+1)-tranches(i))*n_comp*nominal));
    printf("   Default Leg: %.5f \n", dl(i)   );
    printf("   Payment Leg: %.5f \n\n",pl(i)  );
    printf("spread is given in percent if detachment point is less than 0.03.\n");
    
 else
    printf("Tranche %i: %.2f\% - %.2f\% \n", i, tranches(i), tranches(i+1));
    printf("         Price: %.5f \n", price(i));
    printf("   Default Leg: %.5f \n", dl(i));
    printf("   Payment Leg: %.5f \n\n", pl(i));

  end; 
 end;

printf("***************************************************\n");

 for(i=nt:nt+4),
   printf("   Base_correlation: %.5f \n", pl(i));
 end;

for(i=1:5),
 x(i,1)=dl(nt+i-1);
 y(i,1)=pl(nt+i-1); 
end;

plot(x',y);
xtitle('Base correlation','att_det','correlation');

mode(0);

case 2 then 

method=10;

[price, dl, pl] = price_cdo(n_comp, nominal, dates, evstr(product(5)), intensity, taux(:,1), taux(:,2), evstr(model1(2)), evstr(model2(2)), evstr(model1(4)), evstr(model2(4)), method, evstr(params_method));
//winclose(winId);
x=ones(50,1);
y=ones(50,1);


nt=size(tranches,'*');

 for (i=1:nt-1),
  if (tranches(i+1)<=0.03) then   
    printf("Tranche %i: %.2f\% - %.2f\% \n", i, tranches(i), tranches(i+1));
    printf("         Price: %.5f \n", 100*(dl(i)-0.05*pl(i))/((tranches(i+1)-tranches(i))*n_comp*nominal));
    printf("   Default Leg: %.5f \n", dl(i));
    printf("   Payment Leg: %.5f \n\n",pl(i)  );
    printf("spread is given in percent if detachment point is less than 0.03.\n");
    
 else  
    printf("Tranche %i: %.2f\% - %.2f\% \n", i, tranches(i), tranches(i+1));
    printf("         Price: %.5f \n", price(i));
    printf("   Default Leg: %.5f \n", dl(i));
    printf("   Payment Leg: %.5f \n\n", pl(i));

  end; 
 end;


for(i=1:50),
 x(i,1)=dl(nt+i-1);
 y(i,1)=pl(nt+i-1);
  
end;

plot(y',x);
xtitle('Implied intensity density curve ','X_axis','Y_axis');

mode(0);
 end;
end;

 


