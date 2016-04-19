//-*- fill-column: 90 -*-

// This script must be launched from the directory scripts.
// It creates the test files for every Premia pb.
//

exec('../libpremia/loader.sce');

function Lnew = relative_path_data_aux (L, outdir)
  Lnew = list();
  if is(L, %types.List) && length(L)>0 && ~is(L(1), %types.List); then
    Lnew = L;
    if Lnew(1) == "FILENAME"; then
      file ("copy", [Lnew(3), outdir]);
      Lnew(3) = file("join", [outdir, file ("tail", Lnew(3))]);
    end
  else
    for i=(1:size(L))
      Lnew.add_last[relative_path_data_aux (L(i), outdir)];
    end
  end
endfunction


function relative_path_data (P, outdir)
  L = P.get_model_values[];
  Lnew = relative_path_data_aux (L, outdir);
  P.set_model_values[Lnew];
  L = P.get_option_values[];
  Lnew = relative_path_data_aux (L, outdir);
  P.set_option_values[Lnew];
  L = P.get_method_values[];
  Lnew = relative_path_data_aux (L, outdir);
  P.set_method_values[Lnew];
endfunction

function gener_test_files (outdir)
  P=premia_create();
  assets = premia_get_assets();
  
  for asset = assets'
    models = premia_get_models ( asset=asset );
    nmodels = size ( models, '*' );
    families = premia_get_families ( asset=asset );
    nfamilies = size ( families, '*' );
    P.set_asset[str=asset];
    for imodel=(1:nmodels)
       P.set_model[str=models(imodel)];
       for ifamily=(1:nfamilies)
         options = premia_get_family ( ifamily, imodel, asset=asset );
         noptions = size(options, '*');
         for ioption=(1:noptions);
           P.set_option[str=options(ioption)];
           methods=P.get_methods[];
           nmethods=size ( methods, '*' );
           for imethod=(1:nmethods)
             names = [ asset; models(imodel); options(ioption); methods(imethod) ];
             P.set_method[imethod];
             filename = catenate ( names, sep="@" ) + ".bin";
             fullname = file ( "join", [ outdir, filename ] )
             file ( "mkdir", outdir ); 
             relative_path_data (P, outdir);
             printf(filename + "\n");
             save(fullname,P);
           end
        end
      end
    end
  end
endfunction

gener_test_files ("../portfolio-1");
