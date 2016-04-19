if (~(exists('product'))|(product == [])) then x_message('Product does not exist !'); abort; end;
if (~(exists('model2'))|(model2 == [])) then x_message('Model does not exist !'); abort; end;

method = x_choose(['Number of defaults - Hull & White';'Number of defaults - Laurent & Gregory';'Recurrence (homogenous case only) - Hull & White';'Recurrence - Hull & White';'FFT - Laurent & Gregory';'Monte Carlo';'Monte Carlo + Control Variable'],'Which method ?')


params_method_txt=smat_create(0,0);
params_method_def=smat_create(0,0);
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
                params_method_def(2) = '10';
end
params_method = x_mdialog('Parameters', params_method_txt, params_method_def);
if (params_method == []) then
    abort;
else
    params_method_def = params_method;
end;
exec('last_computation.sci');
