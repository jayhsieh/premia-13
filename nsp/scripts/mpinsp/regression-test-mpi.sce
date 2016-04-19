// This file can be run through 
// /usr/local/openmpi-1.2/bin/mpirun -hostfile ~/bhost  -c 10 nsp -nw \
// -f regression-test-mpi.sce -e "quit"
// The list of pricing problems must be set to the variable PBLIST

MPINSP="/usr/local/src/nsp2_dev/contribs/mpinsp/"
exec(MPINSP + 'src/loader.sce')
VERBOSE = %t;

if ~MPI_Initialized() then   MPI_Init();end 
MPI_COMM_WORLD=mpicomm_create('WORLD');
[mpi_rank] = MPI_Comm_rank (MPI_COMM_WORLD);    
[mpi_size] = MPI_Comm_size (MPI_COMM_WORLD);
    
SLV = (mpi_rank <> 0)
MST = ~ SLV;
TAG = 4;

if MST
  PBLIST=glob('../../tests/*/*/*/*/*.tst');
end
printf('rank=%d,mpi_size=%d\n",mpi_rank,mpi_size);


// Charge un object Premia sous forme serialisée et 
// l'envoie précédé de son nom au processus slv.
function send_premia_pb( name, slv )
  ser_obj = sload(name);
  MPI_Send_Obj (name,slv,TAG,MPI_COMM_WORLD); // send name
  pack_obj = MPI_Pack (ser_obj, MPI_COMM_WORLD); // pack
  MPI_Send (pack_obj,slv,TAG,MPI_COMM_WORLD);  // send the packed object
endfunction

exec('../../libpremia/loader.sce');
exec ('../compare_eps.sci');
premia_init()

if SLV // All slaves send result back
  while %t then 
    name = MPI_Recv_Obj(0,TAG,MPI_COMM_WORLD); // receives the name
    if name == '' then 
      if VERBOSE; then printf("I stop working (%d)\n",mpi_rank); end
      break;
    end
    [stat]=MPI_Probe(-1,-1,MPI_COMM_WORLD)
    [elems]=MPI_Get_count(stat)       // 
    pack_obj=mpibuf_create(elems);                 // enough size 
    stat=MPI_Recv (pack_obj, 0, TAG, MPI_COMM_WORLD); // receives the packed object
    ser_obj = MPI_Unpack (pack_obj, MPI_COMM_WORLD); // unpack
    clear P;
    if test_unserialize(ser_obj) == %t; then
        P = unserialize(ser_obj); // unserialize
        Lold = P.get_method_results[];
        P.compute[];
        Lnew = P.get_method_results[];
        [status, res] =compare_eps_diff (Lold,Lnew);
        if status ~= 0 then
            msg = [name + " FAIL"; m2s(res, "%.12f") ];
        else
            msg = name + " OK";
        end
    else
        msg =  name + " UNLOADABLE";
    end
    MPI_Send_Obj(msg,0,TAG,MPI_COMM_WORLD); // sends the results back
    if VERBOSE; then printf("Sending answer to master from %d\n",mpi_rank); end
  end
else // Here at master
  Lpb = PBLIST;
  slv = 1;
  for pb=Lpb(1:mpi_size-1)'			// send 
    if VERBOSE then printf("Sending job to %d %s\n",slv,pb); end
    send_premia_pb (pb, slv);
    slv = slv + 1;
  end
  Lpb(1:mpi_size-1)=[];
  res=smat_create(0,0);
  for pb=Lpb'
    [stat] = MPI_Probe(-1,-1,MPI_COMM_WORLD);
    result = MPI_Recv_Obj(stat.src,TAG,MPI_COMM_WORLD);
    if VERBOSE; then printf("Received answer from %d\n",stat.src); end
    res = [res ; result];
    if VERBOSE; then printf("Sending job to %d %s\n",stat.src,pb); end
    send_premia_pb (pb, stat.src);
  end
  // we still have mpi_size-1 Revc to perform
  for slv=1:mpi_size-1      // collect in any order
    [stat] = MPI_Probe(-1,-1,MPI_COMM_WORLD);
    result = MPI_Recv_Obj(stat.src,TAG,MPI_COMM_WORLD);
    if VERBOSE; then printf("Received answer from %d\n",stat.src); end
    res = [res ; result];
  end
  // tell slaves to stop working 
  for slv=1:mpi_size-1
    MPI_Send_Obj([''],slv,TAG,MPI_COMM_WORLD);
  end
  Htime = localtime();
  out = "regression-test-" + m2s(Htime.year + 1900) + m2s(Htime.mon + 1, format="%02.0f") ...
    + m2s(Htime.mday, format="%02.0f")+ m2s(Htime.hour, format="%02.0f") ...
    + m2s(Htime.min, format="%02.0f") + ".txt";
  F = fopen (out, mode ="w");
  F.put_smatrix[res];
  F.putstr["\n"];
  F.close[]

  save('regression-test.bin',res);
end

// finalize slaves and master 
MPI_Finalize();
