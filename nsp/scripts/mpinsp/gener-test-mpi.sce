// Jérôme Lelong, January 2011
//
// This script is intended to generate missing test files and compute their
// solutions on a parallel architecture. Because this script will mainly be run
// on an 8-core computer, we do not bother sending Premia objects but we rather
// access them directly on the disk.
//
// This script can be run through 
// mpirun -np 8 nsp -nw -f gener-test-mpi.sce -e "quit"

// Some global variables
output_dir = '../../tests';
MPINSP="/usr/local/src/nsp2_dev/contribs/mpinsp/"
exec(MPINSP + 'src/loader.sce')
VERBOSE = %t;

// Initialize Premia stuff
exec('../../libpremia/loader.sce');
exec ('../compare_eps.sci');
premia_init()

// Initialize MPI
if ~MPI_Initialized() then   MPI_Init();end 
MPI_COMM_WORLD=mpicomm_create('WORLD');
[mpi_rank] = MPI_Comm_rank (MPI_COMM_WORLD);    
[mpi_size] = MPI_Comm_size (MPI_COMM_WORLD);

// Global variables for MPI
SLV = (mpi_rank <> 0)
MST = ~ SLV;
TAGWORK = 4;
TAGFINISH = 5;
TAGEXIT = 6;

printf('rank=%d,mpi_size=%d\n",mpi_rank,mpi_size);


// Walks through all the possible Premia problems and creates one instance for
// each of them it it does not already exist.
// The solutions of the created problems are NOT computed. 
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
            if ~file("exists", fullname) then save(fullname,P); end
          end
        end
      end
    end
  end
endfunction


function solve_problems (list_pb)
// All slaves send result back
  if SLV then 
    while %t then 
      name = MPI_Recv_Obj(0,TAGWORK,MPI_COMM_WORLD); // receives the name
      if name == '' then 
        if VERBOSE; then printf("I stop working (%d)\n",mpi_rank); end
        break;
      end
      load (name);
      if (length(P.get_method_results[]) == 0)
      then
        P.compute[];
        save(name, P);
      end
      MPI_Send(1,0,TAGFINISH,MPI_COMM_WORLD); // tells the master it has finished his work
      if VERBOSE; then printf("Sending answer to master from %d\n",mpi_rank); end
    end
  else // Here at master
    Lpb = list_pb;
    slv = 1;
    for pb=Lpb(1:mpi_size-1)'   // send
      if VERBOSE then printf("Sending job to %d %s\n",slv,pb); end
      MPI_Send_Obj (pb,slv,TAGWORK,MPI_COMM_WORLD); // send name
      slv = slv + 1;
    end
    Lpb(1:mpi_size-1)=[];
    for pb=Lpb'
      [stat] = MPI_Probe(-1,-1,MPI_COMM_WORLD);
      fake = zeros(1,1);
      MPI_Recv(fake,stat.src,TAGFINISH,MPI_COMM_WORLD);
      if VERBOSE; then printf("Received answer from %d\n",stat.src); end
      if VERBOSE; then printf("Sending job to %d %s\n",stat.src,pb); end
      MPI_Send_Obj (pb,stat.src,TAGWORK,MPI_COMM_WORLD); // send name
    end
    // we still have mpi_size-1 Revc to perform
    for slv=1:mpi_size-1      // collect in any order
      [stat] = MPI_Probe(-1,-1,MPI_COMM_WORLD);
      MPI_Recv(fake,stat.src,TAGFINISH,MPI_COMM_WORLD);
      if VERBOSE; then printf("Received answer from %d\n",stat.src); end
    end
    // tell slaves to stop working 
    for slv=1:mpi_size-1
      MPI_Send_Obj([''],slv,TAGWORK,MPI_COMM_WORLD);
    end
  end
endfunction


if MST
then
    gener_test_files ();
end
PBLIST=glob(output_dir + '/*/*/*/*/*.tst');
solve_problems (PBLIST);

// finalize slaves and master 
MPI_Finalize();
