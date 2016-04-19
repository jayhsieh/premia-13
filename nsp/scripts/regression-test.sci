//-*- fill-column: 90 -*-

// This script must be launched from the directory scripts.
// It creates the test files for every Premia pb.

output_dir = '../tests';
exec('../libpremia/loader.sce');

function res = detect_nan (L)
  res = 0;
  if is(L, %types.List)
  then
    for i=(1:size(L))
      res = res + detect_nan (L(i));
    end
  else
    if is(L, %types.Mat)
      res = res + sum (isnan (L));
    end
  end
endfunction

function gener_test_files ()
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
             filename = catenate ( names, sep="@" ) + ".tst";
             dirname = catenate ( names, sep="/" );
             outdir = file ( "join", [ output_dir, dirname ] );
             fullname = file ( "join", [ outdir, filename ] )
             file ( "mkdir", outdir );
             printf(filename + "\n");
             if ~file("exists", fullname) 
             then
               P.compute[];
               if (detect_nan(P.get_method_results[]) > 0) 
               then printf("\tNaN values detected !\n");
               end
               save(fullname,P);
              end
           end
        end
      end
    end
  end
endfunction

gener_test_files ();
