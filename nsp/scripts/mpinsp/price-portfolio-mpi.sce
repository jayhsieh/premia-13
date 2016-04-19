// This file can be run through 
// /usr/local/openmpi-1.2/bin/mpirun -hostfile ~/bhost  -c 10 nsp -nw \
// -f price-portfolio.sce -e "quit"
// The list of pricing problems must be set to the variable PBLIST

MPINSP="/usr/local/src/nsp2_dev/contribs/mpinsp/"
exec(MPINSP + 'src/loader.sce')
VERBOSE = %t;

if ~MPI_Initialized() then   MPI_Init();end 
MPI_COMM_WORLD=mpicomm_create('WORLD');
[mpi_rank] = MPI_Comm_rank (MPI_COMM_WORLD);	
[mpi_size] = MPI_Comm_size (MPI_COMM_WORLD);
    
SLV = (mpi_rank <> 0) //
MST = ~ SLV;	  // 
TAG = 4;

printf('rank=%d,mpi_size=%d\n",mpi_rank,mpi_size);
// //
// // au lieu de serialize + pack on peut utiliser MPI_Send_Obj
// // regarder l'envoi de Smat.
// //
// function send_premia_pb( name, slv )
//   load(name);
//   MPI_Send_Obj (name,slv,TAG,MPI_COMM_WORLD);
//   MPI_Send_Obj (P,slv,TAG,MPI_COMM_WORLD);  
// endfunction

// Charge un object Premia sous forme serialisée et l'envoie précédé de son nom au
// processus slv.
function send_premia_pb( name, slv )
  ser_obj = sload(name);
  MPI_Send_Obj (name,slv,TAG,MPI_COMM_WORLD); // send name
  pack_obj = MPI_Pack (ser_obj, MPI_COMM_WORLD); // pack
  MPI_Send (pack_obj,slv,TAG,MPI_COMM_WORLD);  // send the packed object
endfunction


function [sl, result] = receive_res ()
  [stat] = MPI_Probe(-1,-1,MPI_COMM_WORLD);
  sl = stat.src;
  result = MPI_Recv_Obj(sl,TAG,MPI_COMM_WORLD);
  if VERBOSE; then printf("Received answer from %d\n",sl); end
endfunction

exec('../../libpremia/loader.sce');
premia_init()

// PARALLEL // computation (depends on mpi_rank/mpi_size)

if SLV				// All slaves send result back
  while %t then 
    name = MPI_Recv_Obj(0,TAG,MPI_COMM_WORLD); // receives the name
    if name == '' then 
      if vERBOSE; then printf("I stop working (%d)\n",mpi_rank); end
      break;
    end 
    [stat]=MPI_Probe(-1,-1,MPI_COMM_WORLD)
    [elems]=MPI_Get_count(stat);
    pack_obj=mpibuf_create(elems);                 // enough size 
    stat=MPI_Recv (pack_obj, 0, TAG, MPI_COMM_WORLD); // receives the packed object
    ser_obj = MPI_Unpack (pack_obj, MPI_COMM_WORLD); // unpack
    P = unserialize(ser_obj); // unserialize
    P.compute[];
    L = P.get_method_results[];
    result= L(1)(3);
    MPI_Send_Obj(result,0,TAG,MPI_COMM_WORLD); // sends the results back
    if VERBOSE; then printf("Sending answer to master from %d\n",mpi_rank); end
  end
else				// Here at master

  pb_list = PBLIST;

  t=cputime();
  Nt= size(pb_list, '*'); 
  nb_per_node = floor (Nt / (mpi_size-1));  
  prm=grand(1,'perm',Nt);
  Lpb=pb_list(prm);
  for pb=Lpb(1:mpi_size-1)'			// send 
    if VERBOSE then printf("Sending job to %d %s\n",slv,pb); end
    send_premia_pb (pb, slv);
    slv = slv + 1;
  end
  res=list();
  Lpb(1:mpi_size-1)=[];
  for pb=Lpb'
    [sl, result] = receive_res ();
    res.add_last[list(sl, result)];
    if VERBOSE then printf("Sending job to %d %s\n",sl,pb); end
    send_premia_pb (pb, sl);
  end
  for slv=1:mpi_size-1			// send 
    if VERBOSE; then printf("Sending job to %d %s\n",slv,Lpb(slv)); end
    send_premia_pb (Lpb(slv), slv);
  end
  // we still have mpi_size-1 Revc to perform
  for slv=1:mpi_size-1		// collect in any order
    [sl, result] = receive_res ();
    res.add_last[list(sl, result)];
  end
  // tell slaves to stop working 
  for slv=1:mpi_size-1
    MPI_Send_Obj([''],slv,TAG,MPI_COMM_WORLD);
  end
  t=cputime()-t;
  save('pb-res.bin',res,CPU=t);
  printf ("CPU time %f \n", t);
end

// finalize slaves and master 

MPI_Finalize();
