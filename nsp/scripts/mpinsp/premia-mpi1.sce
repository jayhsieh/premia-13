//// Copyright (C) 2007 Javier Fernández Baldomero
//// This program is free software
// myPi:	Classic PI computation by numeric integration of arctan'(x) in [0..1]
//
//  N	[1E7]	//subdivisions of the [0, 1] interval
//  mod	['s']	communication modality:  (s)end (r)educe
//
//  printed results struct contains
//	pi	estimated pi value
//	err	error
//	time	from argument xmit to pi computed
//
// Examples:
//
// This file can be run through 
// /usr/local/openmpi-1.2/bin/mpirun -c 2 /usr/bin/nsp -nw -f premia-mpi1.sce -e "quit"

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
  Maturity=0;
  Strike=0;
  MPI_Recv(Maturity,-1,TAG,MPI_COMM_WORLD);
  MPI_Recv(Strike,-1,TAG,MPI_COMM_WORLD);
  exec('premia.sce');
  result= L(1)(3);
  MPI_Send(result,0,TAG,MPI_COMM_WORLD);
else				// Here at master
  for slv=1:size-1			// collect in any order
    Maturity=slv;
    Strike=100;
    MPI_Send(Maturity,slv,TAG,MPI_COMM_WORLD);
    MPI_Send(Strike,slv,TAG,MPI_COMM_WORLD);
  end
  result=0;
  res=[];
  for slv=1:size-1			// collect in any order
    MPI_Recv(result,-1,TAG,MPI_COMM_WORLD);
    res=[res,result];
  end
  save('poo.bin',res);
end

// finalize slaves and master 

MPI_Finalize();








  
  

