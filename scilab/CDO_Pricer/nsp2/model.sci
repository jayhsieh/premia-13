intensity=list('combo','Intensity',1,['Constant','From a file']);
recovery=list('combo', 'Recovery',1,['Constant','Uniforme','Gaussian']);
interest=list('combo', 'Interest rate',1,['Constant','From a file']);
copule=list('combo', 'One factor copula',1,['Gaussian','Clayton','NIG']);
if exists('model1') && length(model1)==4 then
     intensity(3) = model1(1);
     recovery(3) = model1(2);
     interest(3) = model1(3);
     copule(3) = model1(4);
end;

model1_def = list(intensity, recovery, interest, copule);
[Lres, L, model1_new] = x_choices('Model -1-', model1_def);
if model1_new == [] then 
  abort
end
modif_model = [%t, %t, %t, %t];
if ((exists('model1'))&&(length(model1)==4)) then 
    if (model1_new(1) == model1(1)) then modif_model(1) = %f; end;
    if (model1_new(2) == model1(2)) then modif_model(2) = %f; end;
    if (model1_new(3) == model1(3)) then modif_model(3) = %f; end;
    if (model1_new(4) == model1(4)) then modif_model(4) = %f; end;
end;
model1 = model1_new;

// initiliaze empty smat
if (~ exists('model2_txt')) then
  model2_txt=smat_create(0,0);
end

if (~ exists('model2_def')) then
  model2_def=smat_create(0,0);
end


if (exists('model2')) then bool_model2 = ~(model2 == []); 
    else bool_model2 = %f; end;
modif_model = modif_model||~bool_model2
if modif_model(1) then
    select model1(1),
	case 1 then model2_txt(1) = 'Intensity value';
		    model2_def(1) = '0.01';
	case 2 then model2_txt(1) = 'Name of the file for the intensity';
		    model2_def(1) = '../datas/intensity.dat';
    end
end;
if modif_model(2) then
    select model1(2),
	case 1 then model2_txt(2) = 'Recovery value';
		    model2_def(2) = '0.4';
	case 2 then model2_txt(2) = 'Bounds [a, b] of the uniform recovery';
		    model2_def(2) = '[0.3 , 0.5]';
	case 3 then model2_txt(2) = 'Mean and variance [m, s] of the Gaussian recovery';
		    model2_def(2) = '[0.4, 0.2]';
    end
end;
if modif_model(3) then
    select model1(3),
	case 1 then model2_txt(3) = 'Enter the interest rate';
		    model2_def(3) = '0.04';
	case 2 then model2_txt(3) = 'Enter the name of the file for the rate';
		    model2_def(3) = 'taux.dat';
    end
end;
if modif_model(4) then
    select model1(4),
	case 1 then model2_txt(4) = 'Gaussian Copula - Correlation parameter';
		    model2_def(4) = '0.03';
	case 2 then model2_txt(4) = 'Clayton Copula - Theta parameter';
		    model2_def(4) = '0.2';
	case 3 then model2_txt(4) = 'NIG Copula - [Correlation, alpha, beta]';
		    model2_def(4) = '[0.06, 1.2558, -0.2231]';
	case 4 then model2_txt(4) = 'DoubleT Copula - [Correlation, T1, T2]';
		    model2_def(4) = '[0.2, 5, 5]';
    end
end;
model2 = x_mdialog('Model -2-', model2_txt, model2_def);
if (model2 == []) then 
    abort;
else
    model2_def = model2; 
end;
test_int = execstr('evstr(model2(1))',errcatch=%t);
if  test_int > 0 then
    intensity = read(model2(1), n_comp, 1);
else
    intensity = evstr(model2(1));
end;
if (size(intensity,'*')==1) then intensity=intensity*ones(n_comp,1); end

