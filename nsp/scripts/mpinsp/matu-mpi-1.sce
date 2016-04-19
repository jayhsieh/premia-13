// This file can be run through 
// /usr/local/openmpi-1.2/bin/mpirun -c 2 /usr/bin/nsp -nw -f matu-mpi-1.sce -e "quit"
// /usr/local/openmpi-1.2/bin/mpirun -hostfile ~/bhost  -c 10 /usr/bin/nsp -nw -f matu-mpi-1.sce -e "quit"

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
  nb = 0; // receives the size of the vector to come
  MPI_Recv (nb,-1,TAG,MPI_COMM_WORLD);
  Maturities=zeros(1,nb);
  MPI_Recv (Maturities,-1,TAG,MPI_COMM_WORLD); // receives the vector
  result=[];
  for Maturity=Maturities // iterate over the vector
    exec ('premia.sce');
    result= [L(1)(3), result];
  end
  MPI_Send (result,0,TAG,MPI_COMM_WORLD); // sends the results back

else				// Here at master
  t=cputime();
  Nt= 24*3; 
  nb_per_node = floor (Nt / (size-1));  
  Maturities = linspace(0.25,5,Nt);
  for slv=1:size-1			// send 
    MPI_Send (nb_per_node,slv,TAG,MPI_COMM_WORLD);
    MPI_Send (Maturities((slv-1)*nb_per_node+1:slv*nb_per_node),slv,TAG,MPI_COMM_WORLD);
  end

// iterate over the part of the vector that has not been sent
  res=[];
  for Maturity=Maturities((size-1) * nb_per_node + 1:$)     
    exec ('premia.sce');
    res= [L(1)(3), res];
  end
  
// preparing to collect results
  result=zeros(1,nb_per_node);
  for slv=1:size-1			// collect in any order
    MPI_Recv(result,-1,TAG,MPI_COMM_WORLD);
    res=[res,result];
  end
  t=cputime()-t;
  save('matu.bin',res,CPU=t);
end

// finalize slaves and master 

MPI_Finalize();








  
  

