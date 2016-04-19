select Method,
//model2_txt=[];
//model2_def=[];

case 1 then 
intensity_l=list('Intensity',1,['Constant','From a file']);
recovery=list('Recovery',1,['Constant','Uniforme','Gaussian']);
interest=list('Interest rate',1,['Constant','Linear step function from a file']);
copule=list('copula',1,['Gaussian','Clayton','NIG','Student','Double_T']);


if ((exists('model1'))&(length(model1)==4)) then
     intensity_l(2) = model1(1);
     recovery(2) = model1(2);
     interest(2) = model1(3);
     copule(2) = model1(4);
     //Implied_Methods(2)=model1(5);
end;

model1_def = list(intensity_l, recovery, interest, copule);
model1_new = x_choices('Model -1-', model1_def);
if (model1_new == []) then 
    abort;
//    model1_new = [1, 1, 1, 1];
end;
modif_model = [ %T, %T, %T,%T];
if ((exists('model1'))&(length(model1)==5)) then 
    if (model1_new(1) == model1(1)) modif_model(1) = %F; end;
    if (model1_new(2) == model1(2)) modif_model(2) = %F; end;
    if (model1_new(3) == model1(3)) modif_model(3) = %F; end;
    if (model1_new(4) == model1(4)) modif_model(4) = %F; end;
    //if (model1_new(5) == model1(5)) modif_model(5) = %F; end;
end;
model1 = model1_new;

//model2_txt = [];
if (exists('model2')) then bool_model2 = ~(model2 == []); 
    else bool_model2 = %F; end;
modif_model = modif_model|~bool_model2
if modif_model(1) then
    select model1(1),
	case 1 then model2_txt(1) = 'Intensity value';
		    model2_def(1) = '0.01';
	case 2 then model2_txt(1) = 'Name of the file for the intensity';
		    model2_def(1) = 'intensity.dat';
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
	case 4 then model2_txt(4) = 'Student Copula - [Correlation, dlib T1]';
		    model2_def(4) = '[0.02, 5]';
       case 5 then model2_txt(4) = 'Double T Copula - [Correlation, dlib T1,dlib T2]';
		    model2_def(4) = '[0.03 ,5 ,7]';
     end 
end;
 
model2 = x_mdialog('Model -2-', model2_txt, model2_def);

if (model2 == []) then 
    abort;
else
    model2_def = model2; 
end;

// intensity
test_int = execstr('evstr(model2(1))','errcatch');
if  test_int > 0 then
    intensity = read('../datas/'+model2(1), n_comp, 1);
else
    intensity = evstr(model2(1));
end;
if (size(intensity,'*')==1) then intensity=intensity*ones(n_comp,1); end

// taux 
test_taux = execstr('evstr(model2(3))','errcatch');
if  test_taux > 0 then
    taux = read('../datas/'+model2(3), -1, 2); 
else
    taux = evstr(model2(3));
end;
if (size(taux,'*')==1) then taux=[0,0;evstr(prod_def(3)),evstr(model2(3))*evstr(prod_def(3))]; end

function print_model()
   printf("Model:\n");
   printf("   Intensity: %s = %s\n", intensity_l(3)(evstr(model1(1))), model2_def(1));
   printf("   Recovery: %s =  %s\n", recovery(3)(evstr(model1(2))), model2_def(2));
   printf("   Interest rate: %s =  %s\n", interest(3)(evstr(model1(3))), model2_def(3));
   printf("   Copula: %s =  %s\n", copule(3)(evstr(model1(4))), model2_def(4));
   
  endfunction

case 2 then 
intensity_l=list('Intensity',1,['Constant','From a file']);
recovery=list('Recovery',1,['Constant','Uniforme','Gaussian']);
interest=list('Interest rate',1,['Constant','Linear step function from a file']);
Implied_Methods=list('Methodes implicites',1,['Base_correlation','Implied_copula']);

if ((exists('model1'))&(length(model1)==4)) then
     intensity_l(2) = model1(1);
     recovery(2) = model1(2);
     interest(2) = model1(3);
     Implied_Methods(2) = model1(4);
     
end;
model1_def = list(intensity_l, recovery, interest, Implied_Methods);
model1_new = x_choices('Model -1-', model1_def);
if (model1_new == []) then 
    abort;
//    model1_new = [1, 1, 1, 1];
end;
modif_model = [ %T, %T, %T,%T];
if ((exists('model1'))&(length(model1)==4)) then 
    if (model1_new(1) == model1(1)) modif_model(1) = %F; end;
    if (model1_new(2) == model1(2)) modif_model(2) = %F; end;
    if (model1_new(3) == model1(3)) modif_model(3) = %F; end;
    if (model1_new(4) == model1(4)) modif_model(4) = %F; end;
    
end;
model1 = model1_new;

//model2_txt = [];
if (exists('model2')) then bool_model2 = %F; // bool_model2 = ~(model2 == []); 
    else bool_model2 = %F; end;
modif_model = modif_model|~bool_model2
if modif_model(1) then
    select model1(1),
	case 1 then model2_txt(1) = 'Intensity value';
		    model2_def(1) = '0.01';
	case 2 then model2_txt(1) = 'Name of the file for the intensity';
		    model2_def(1) = 'intensity.dat';
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
		    model2_def(3) = '0.08';
	case 2 then model2_txt(3) = 'Enter the name of the file for the rate';
		    model2_def(3) = 'taux.dat';
    end
end;
if modif_model(4) then
    select model1(4),
	case 1 then model2_txt(4) = 'Base Correlation - itraxx_CDX';
		    model2_def(4) = '1';
	case 2 then model2_txt(4) = 'Implied Copula - itraxx_CDX';
		    model2_def(4) = '1';
    end 	
end;
 
model2 = x_mdialog('Model -2-', model2_txt, model2_def);

if (model2 == []) then 
    abort;
else
    model2_def = model2; 
end;

// intensity
test_int = execstr('evstr(model2(1))','errcatch');
if  test_int > 0 then
    intensity = read('../datas/'+model2(1), n_comp, 1);
else
    intensity = evstr(model2(1));
end;
if (size(intensity,'*')==1) then intensity=intensity*ones(n_comp,1); end

// taux 
test_taux = execstr('evstr(model2(3))','errcatch');
if  test_taux > 0 then
    taux = read('../datas/'+model2(3), -1, 2); 
else
    taux = evstr(model2(3));
end;
if (size(taux,'*')==1) then taux=[0,0;evstr(prod_def(3)),evstr(model2(3))*evstr(prod_def(3))]; end

   function print_model()
   printf("Model:\n");
   //printf("   Intensity: %s = %s\n", intensity_l(3)(evstr(model1(1))), model2_def(1));
   printf("   Recovery: %s =  %s\n", recovery(3)(evstr(model1(2))), model2_def(2));
   printf("   Interest rate: %s =  %s\n", interest(3)(evstr(model1(3))), model2_def(3));
   printf("   Implied Method: %s =  %s\n", Implied_Methods(3)(evstr(model1(4))), model2_def(4));
   
endfunction

end 

//function print_model()
 //   printf("Model:\n");
 //   printf("   Intensity: %s = %s\n", intensity_l(3)(evstr(model1(1))), model2_def(1));
 //   printf("   Recovery: %s =  %s\n", recovery(3)(evstr(model1(2))), model2_def(2));
 //   printf("   Interest rate: %s =  %s\n", interest(3)(evstr(model1(3))), model2_def(3));
 //   printf("   Copula: %s =  %s\n", copule(3)(evstr(model1(4))), model2_def(4));
   
//endfunction
