select Method,
case 1 then 
if (~(isdef('product'))|(product == [])) then x_message('Product does not exist !'); abort; end;
if (~(isdef('model2'))|(model2 == [])) then x_message('Model does not exist !'); abort; end;

methods = ['Number of defaults - Hull & White';'Number of defaults - Laurent & Gregory';'Recurrence (homogenous case only) - Hull & White';'Recurrence - Hull & White';'FFT - Laurent & Gregory';'Monte Carlo';'Monte Carlo + Control Variable';'Saddlepoint']

method = x_choose(methods, 'Which method ?')
params_method_txt=[];
params_method_def=[];

select method, 
    case 1 then params_method_txt(1) = 'Subdivisions of the time';
		params_method_def(1) = '4';
    case 2 then params_method_txt(1) = 'Subdivisions of the time';
		params_method_def(1) = '4';
    case 3 then params_method_txt(1) = 'Subdivisions of the time';
		params_method_def(1) = '4';
		params_method_txt(2) = 'Subdivisions of the losses';
		params_method_def(2) = '100';
    case 4 then params_method_txt(1) = 'Subdivisions of the time';
		params_method_def(1) = '4';
		params_method_txt(2) = 'Subdivisions of the losses';
		params_method_def(2) = '100';
    case 5 then params_method_txt(1) = 'Subdivisions of the time';
		params_method_def(1) = '4';
		params_method_txt(2) = 'Subdivisions of the losses';
		params_method_def(2) = '100';
    case 6 then params_method_txt(1) = 'Number of Monte-Carlo iterations';
		params_method_def(1) = '10000';
    case 7 then params_method_txt(1) = 'Number of Monte-Carlo iterations';
		params_method_def(1) = '10000';
		params_method_txt(2) = 'Subdivisions of the time for the control variate';
		params_method_def(2) = '4';
   
   case 8 then params_method_txt(1) = 'Subdivisions of the time';
		params_method_def(1) = '4';
		
	

end
params_method = x_mdialog('Parameters', params_method_txt, params_method_def);
if (params_method == []) then 
    abort;
else
    params_method_def = params_method; 
end;


case 2 then 
params_method_txt=[];
params_method_def=[];
if((model2_def(4)=='1')) then 
  params_method_txt(1)='spread [0.00 0.03]';
  params_method_def(1)='24';
  params_method_txt(2)='spread [0.03 0.06]';
  params_method_def(2)='81';
  params_method_txt(3)='spread [0.06 0.09]';
  params_method_def(3)='39';
  params_method_txt(4)='spread [0.09 0.12]';
  params_method_def(4)='12';
  params_method_txt(5)='spread [0.12 0.22]';
  params_method_def(5)='9';
  params_method_txt(6)='spread - index CDS';
  params_method_def(6)='36.375';

 else 
  params_method_txt(1)='spread [0.00 0.03]';
  params_method_def(1)='35';
  params_method_txt(2)='spread [0.03 0.07]';
  params_method_def(2)='292';
  params_method_txt(3)='spread [0.07 0.1]';
  params_method_def(3)='85';
  params_method_txt(4)='spread [0.1 0.15]';
  params_method_def(4)='39';
  params_method_txt(5)='spread [0.15 0.3]';
  params_method_def(5)='12';
  params_method_txt(6)='spread -CDX index CDS';
  params_method_def(6)='41';
 
end;
params_method = x_mdialog('Parameters', params_method_txt, params_method_def);
if (params_method == []) then 
    abort;
else
    params_method_def = params_method; 
end;

//function print_method()
   
    //for (i = 1:size(params_method_def,'*')), 
	//printf("   %s: %s\n", params_method_txt(i), params_method_def(i));
  // end;
//endfunction
 
//exec('last_computation.sci',-1);

end; 
function print_method()
   
    for (i = 1:size(params_method_def,'*')), 
	printf("   %s: %s\n", params_method_txt(i), params_method_def(i));
    end;
endfunction

exec('last_computation.sci',-1);