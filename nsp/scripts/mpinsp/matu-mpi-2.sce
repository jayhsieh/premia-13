// This file can be run through 
// /usr/local/openmpi-1.2/bin/mpirun -c 2 /usr/bin/nsp -nw -f premia-mpi-2.sce -e "quit"
// /usr/local/openmpi-1.2/bin/mpirun -hostfile ~/bhost  -c 10 /usr/bin/nsp -nw -f matu-mpi-2.sce -e "quit"

MPINSP="/usr/local/src/nsp2_dev/contribs/mpinsp/"
exec(MPINSP + 'src/loader.sce')

if ~MPI_Initialized() then   MPI_Init();end 
MPI_COMM_WORLD=mpicomm_create('WORLD');
[rank] =	MPI_Comm_rank (MPI_COMM_WORLD);	
[size] =	MPI_Comm_size (MPI_COMM_WORLD);
    
SLV = (rank <> 0) //
MST = ~ SLV;	  // 
TAG=4;

printf('rank=%d,size=%d\n",rank,size);

// PARALLEL // computation (depends on rank/size)
////////////////////////	// vectorized code, equivalent to

if SLV				// All slaves send result back
  while %t then 
    Maturity=0;
    MPI_Recv (Maturity,-1,TAG,MPI_COMM_WORLD); // receives the vector
    if Maturity < 0 then 
      printf("I stop working (%d)\n",rank);
      break;
    end 
    result=[];
    exec ('premia.sce');
    result= [rank,Maturity,L(1)(3)];
    MPI_Send(result,0,TAG,MPI_COMM_WORLD); // sends the results back
    printf("Sending answer to master from %d\n",rank);
  end
else				// Here at master
  t=cputime();
  Nt= 60; 
  nb_per_node = floor (Nt / (size-1));  
  Maturities = linspace(1. /12,5,Nt);
  prm=grand(1,'perm',length(Maturities));
  M=Maturities(prm);
  for slv=1:size-1			// send 
    printf("Sending job to %d %f\n",slv,M(slv))
    MPI_Send (M(slv),slv,TAG,MPI_COMM_WORLD);
  end
  M(1:size-1)=[];
  res=[];
  result=ones(1,3);
  while %t 
    printf("Receiving job stack %d\n",length(M));
    MPI_Recv(result,-1,TAG,MPI_COMM_WORLD);
    printf("Received answer from %d\n",result(1));
    sl=result(1);
    res=[res;result];
    if ~isempty(M) then 
      printf("Sending job to %d %f\n",sl,M(1))
      MPI_Send (M(1),sl,TAG,MPI_COMM_WORLD);
      M(1)=[];
    else
      break;
    end
  end
  // we still have size-2 Revc to perform
  for slv=1:size-2		// collect in any order
    printf("Receiving job stack %d\n",length(M));
    MPI_Recv(result,-1,TAG,MPI_COMM_WORLD);
    printf("Received answer from %d\n",result(1));
    res=[res;result];
  end
  // tell slaves to stop working 
  for slv=1:size-1
    MPI_Send([-1],slv,TAG,MPI_COMM_WORLD);
  end
  t=cputime()-t;
  save('matu.bin',res,CPU=t);
  printf ("CPU time %f \n", t);
end

// finalize slaves and master 

MPI_Finalize();
