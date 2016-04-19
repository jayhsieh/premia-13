// This file can be run through 
// /usr/local/openmpi-1.2/bin/mpirun -c 2 /usr/bin/nsp -nw -f premia-pb.sce -e "quit"
// /usr/local/openmpi-1.2/bin/mpirun -hostfile ~/bhost  -c 10 /usr/bin/nsp -nw -f matu-pb.sce -e "quit"

MPINSP="/usr/local/src/nsp2_dev/contribs/mpinsp/"
exec(MPINSP + 'src/loader.sce')
PBDIR='../portfolio-1/'; // do not forget the final slash

if ~MPI_Initialized() then   MPI_Init();end 
MPI_COMM_WORLD=mpicomm_create('WORLD');
[mpi_rank] =	MPI_Comm_rank (MPI_COMM_WORLD);	
[mpi_size] =	MPI_Comm_size (MPI_COMM_WORLD);
    
SLV = (mpi_rank <> 0) //
MST = ~ SLV;	  // 
TAG=4;

printf('rank=%d,mpi_size=%d\n",mpi_rank,mpi_size);


//
// au lieu de serialize + pack on peut utiliser MPI_Send_Obj
// regarder l'envoi de Smat.
//
function send_premia_pb( name, slv )
  load(name);
  MPI_Send_Obj (name,slv,TAG,MPI_COMM_WORLD);
  MPI_Send_Obj (P,slv,TAG,MPI_COMM_WORLD);  
endfunction

exec('../../libpremia/loader.sce');
premia_init()

// PARALLEL // computation (depends on mpi_rank/mpi_size)
////////////////////////	// vectorized code, equivalent to

if SLV				// All slaves send result back
  while %t then 
    name = MPI_Recv_Obj(0,TAG,MPI_COMM_WORLD); // receives the name
    if name == '' then 
      printf("I stop working (%d)\n",mpi_rank);
      break;
    end 
    P=MPI_Recv_Obj (0,TAG,MPI_COMM_WORLD); // receives the packed object
    P.compute[];
    L = P.get_method_results[];
    result= list(mpi_rank,name,L(1)(3));
    MPI_Send_Obj(result,0,TAG,MPI_COMM_WORLD); // sends the results back
    printf("Sending answer to master from %d\n",mpi_rank);
  end
else				// Here at master

  pb_list = glob(PBDIR + '*.bin');

  t=cputime();
  Nt= size(pb_list, '*'); 
  nb_per_node = floor (Nt / (mpi_size-1));  
  prm=grand(1,'perm',Nt);
  Lpb=pb_list(prm);
  for slv=1:mpi_size-1			// send 
    printf("Sending job to %d %s\n",slv,Lpb(slv))
    send_premia_pb (Lpb(slv), slv);
  end
  Lpb(1:mpi_size-1)=[];
  res=list();
  while %t 
    printf("Receiving job stack %d\n",size(Lpb, '*'));
    result = MPI_Recv_Obj(-1,TAG,MPI_COMM_WORLD);
    sl=result(1);
    printf("Received answer from %d\n",sl);
    res.add_last[result];
    if ~isempty(Lpb) then 
      printf("Sending job to %d %s\n",sl,Lpb(1))
      send_premia_pb (Lpb(1), sl);
      Lpb(1)=[];
    else
      break;
    end
  end
  // we still have mpi_size-2 Revc to perform
  for slv=1:mpi_size-2		// collect in any order
    printf("Receiving job stack %d\n",size(Lpb, '*'));
    result = MPI_Recv_Obj(-1,TAG,MPI_COMM_WORLD);
    printf("Received answer from %d\n",result(1));
    res.add_last[result];
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

//MPI_Finalize();
