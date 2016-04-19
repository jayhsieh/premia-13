//-------------------------------------------------------
//  spawning a nsp child then send/recv 
//-------------------------------------------------------

MPI_Init(); 
COMM_WORLD=mpicomm_create('WORLD');
COMM_SELF=mpicomm_create('SELF');
INFO_NULL=mpiinfo_create('NULL');

MPINSP="/usr/local/src/nsp2_dev/contribs/mpinsp/"
nsp_exe = getenv('SCI')+'/bin/nsp';
args=["-e", "exec(MPINSP + ''src/loader.sce'');MPI_Init();","-name","nsp-child"];
[children errs] = MPI_Comm_spawn (nsp_exe,args,1,INFO_NULL,0, COMM_SELF);
                       // a new nsp is started and MPI_Init called 
		       [parent]=MPI_Comm_get_parent();
		       [NEWORLD] = MPI_Intercomm_merge (parent, 1)

// %%% parent nsp %%%%%%%%%%%%%%
[NEWORLD] = MPI_Intercomm_merge (children, 0)	        // Intercomm unblocks
A=[1,5,89];
MPI_Send_Obj(A,1,7,NEWORLD); 		        // dst=1, tag=7
                        // The proble unblocks and we get from parent 
			A=MPI_Recv_Obj(0,7,NEWORLD);
			A == [1,5,89]
                        // send to parent 
                        A=[6,7,89];
			MPI_Send_Obj(A,0,7,NEWORLD)	// dst=0, tag=7
// %%%% recieve from child 
B=MPI_Recv_Obj(1,7,NEWORLD)		        // src=1, tag=7, stat.len=178!!!
                        info=MPI_Finalize()	// 
			quit			// quit child 

						
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Point-to-point send modes
// MPI_Bsend, _Ssend, _Rsend, _Buffer_attach, _Buffer_detach
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exec(MPINSP + 'src/loader.sce')
MPI_Init(); 
COMM_SELF=mpicomm_create('SELF');
INFO_NULL=mpiinfo_create('NULL');
cmd = "exec(MPINSP + ''src/loader.sce'');MPI_Init();";
cmd = cmd + "parent=MPI_Comm_get_parent();[NEWORLD]=MPI_Intercomm_merge(parent,1);";
nsp_exe = getenv('SCI')+'/bin/nsp';
args=["-name","nsp-child","-e", cmd];
[children errs] = MPI_Comm_spawn (nsp_exe,args,1,INFO_NULL,0, COMM_SELF);
// child will execute cmd 
[NEWORLD] = MPI_Intercomm_merge (children, 0);

// -------------------------------------------------------
// We'll first send a "big" message (greater than tcp_short=64KB)
// so that the different modes are better understood later
// -------------------------------------------------------

A=ones(8,1024);                           	// 64KB
MPI_Send(A,1,7,NEWORLD)			// returns immediately (rpi tcp)
			[stat]=MPI_Probe(-1,-1,NEWORLD)	// (lamd blocks)
			B=zeros(8,1024);// 
			[stat]=MPI_Recv(B,0,7,NEWORLD);

A=ones(9,1024);                                //greater than 64KB
MPI_Send(A,1,7,NEWORLD)		       // returns when recvd (rpi tcp)
			[stat]=MPI_Probe(-1,-1,NEWORLD)	// not yet
			B=zeros(9,1024);
			[stat]=MPI_Recv(B,0,7,NEWORLD); // 

//-------------------------------------------------------
// User can provide buffer, instead of system buffer
//-------------------------------------------------------

B='example of string'; B=strcat(repmat(B,4,1024));	// long string
sb=length(B);					        // 69632 bytes> 64KB
A=B;sa=sb;
snd=sa+sb+2*40;                                         // MPI_BSEND_OVERHEAD is 40 

// take care that buf must be kept alive until detach is called 

buf = mpibuf_create(snd);
MPI_Buffer_attach(buf)				// user buffer
info=MPI_Bsend (A,1,7,NEWORLD)				// won't block
info=MPI_Bsend (B,1,7,NEWORLD)				// won't block
							// won't block,will fail
// MPI_Errhandler_set(NEWORLD,'RETURN')			// avoid default abort
MPI_Bsend (A,1,7,NEWORLD)				// this one
                                                        // should fail.

		[stat]=MPI_Probe(-1,-1,NEWORLD)	        // 
		[elems]=MPI_Get_elements(stat,'a')      // elems= 69633
							// room for A
		A='1234567 of string'; A=strcat(repmat(A,4,1024));
		[stat]=MPI_Recv(A,0,7,NEWORLD)	// blocks !?!? :-(

[flag]=MPI_Iprobe (-1,-1,NEWORLD)			// must be awaken

		A(1)					// now unblocks
		[stat]=MPI_Probe(-1,-1,NEWORLD)	// stat.len=69632
		B='plenty of room for B'
		B=strcat(repmat(B,4,1024));		// looong string
		[info stat]=MPI_Recv(B,0,7,NEWORLD)	// blocked again !!!

[flag]=MPI_Iprobe (-1,-1,NEWORLD)			// again same trick



MPI_Buffer_detach()				// back to sysbuf

//-------------------------------------------------------
// User can send synchronous, waiting for starting Recv (Send may finish sooner)
// LAM may buffer Sends (specially in tcp mode) and make MPI_Send finish sooner
// Ssend is a way of not allowing that
// you need msgsize<=64KB in tcp rpi to see any difference
//-------------------------------------------------------
help modes

help MPI_Ssend
help MPI_Send

A=ones(8,1024); s=whos('A'), sa=s.bytes		// 64KB
info=MPI_Send (A,1,7,NEWORLD)			// returns immediately (rpi tcp)
info=MPI_Send (A,1,7,NEWORLD)			// returns immediately
info=MPI_Ssend(A,1,7,NEWORLD)			// blocks until Recv

			[info stat]=MPI_Probe(-1,-1,NEWORLD)
			B=zeros(8,1024); whos B, B(1)
			[info stat]=MPI_Recv(B,0,7,NEWORLD)	// matches Send
			B(1), B(1)=0; B(1)
			[info stat]=MPI_Probe(-1,-1,NEWORLD)
			[info stat]=MPI_Recv(B,0,7,NEWORLD)	// matches 2nd
			B(1), B(1)=0; B(1)
			[info stat]=MPI_Probe(-1,-1,NEWORLD)	// matches Ssend
			[info stat]=MPI_Recv(B,0,7,NEWORLD)	// & unblocks it
			B(1)

//-------------------------------------------------------
// User can send to an already waiting recv ("ready" send)
// DOUBT: no error? defaults to MPI_Send when no recv posted?
//-------------------------------------------------------
help MPI_Rsend
			[info stat]=MPI_Probe(-1,-1,NEWORLD)	// blocks

info=MPI_Rsend(A,1,7,NEWORLD)		// returns immediately, unblocks probe
					// defaults to MPI_Send?
info=MPI_Rsend(A,1,7,NEWORLD)		// returns immediately ?!?

			[info stat]=MPI_Recv(B,0,7,NEWORLD)	// matches 1st
			B(1), B(1)=0; B(1)
			[info stat]=MPI_Probe(-1,-1,NEWORLD)	// doesn't blk
			[info stat]=MPI_Recv(B,0,7,NEWORLD)	// matches 2nd
			B(1), B(1)=0; B(1)
			[info stat]=MPI_Recv(B,0,7,NEWORLD)	// blocks
info=MPI_Rsend(A,1,7,NEWORLD)					// unblocks
			B(1)

			info=MPI_Finalize
			quit
info=MPI_Finalize
quit


//-------------------------------------------------------
// Starting all over within Nsp
//-------------------------------------------------------

MPI_Finalize
[info flag]=MPI_Initialized
[M,MEX,J]=inmem
clear(MEX{:})					// clear all MPITB from mem
[info flag]=MPI_Initialized			// love that
getenv('LAM_MPI_SSI_rpi')			// should be 'lamd'
putenv('LAM_MPI_SSI_rpi=tcp')			// if you followed instructions
getenv('LAM_MPI_SSI_rpi')			// should be 'tcp' now
MPI_Init

			////////// ouch!!! forgot about the other child Nsp
			////////// and I only have 1 Nsp license for ox1
			MPI_Finalize
			quit

[info children errs] = MPI_Comm_spawn ('/usr/X11R6/bin/xterm',args,...
					1,'NULL',0,'SELF')
[info NEWORLD] = MPI_Intercomm_merge (children, 0)	// child unblocks

			////////// type this on new child //////////////////////////////////
			getenv('LAM_MPI_SSI_rpi')	// should be 'tcp'

//-------------------------------------------------------
// Checking for tcp behavior
//-------------------------------------------------------
A=ones(8,1024); s=whos('A'), sa=s.bytes		// 64KB
info=MPI_Send(A,1,7,NEWORLD)			// won't block
			B=zeros(8,1024); B(1)
			[info stat]=MPI_Recv(B,0,7,NEWORLD)
			B(1), B(1)=0; B(1)

A=ones(9,1024); s=whos('A'), sa=s.bytes		// >64KB
info=MPI_Send(A,1,7,NEWORLD)			// blocks
			B=zeros(9,1024); B(1)
			[info stat]=MPI_Recv(B,0,7,NEWORLD)	// unblocks
			B(1), B(1)=0; B(1)
			[info stat]=MPI_Recv(B,0,7,NEWORLD)	// ready recv

info=MPI_Send(A,1,7,NEWORLD)			// won't block, unblocks recv
			B(1)
			MPI_Finalize
			quit
MPI_Finalize
quit

// ====================================================
// Non-blocking (immediate) point-to-point
// MPI_Isend, _Ibsend, _Issend, _Irsend, _Irecv, _Iprobe
// MPI_Test, _Testall, _Testany, _Testsome
// MPI_Wait, _Waitall, _Waitany, _Waitsome
// ====================================================

exec src7.x/loader.sce
MPI_Init(); 
COMM_SELF=mpicomm_create('SELF');
INFO_NULL=mpiinfo_create('NULL');
cmd = "exec(''src7.x/loader.sce'');MPI_Init();";
cmd = cmd + "parent=MPI_Comm_get_parent();[NEWORLD]=MPI_Intercomm_merge(parent,1);";
nsp_exe = getenv('SCI')+'/bin/nsp';
args=["-name","nsp-child","-e", cmd];
[children errs] = MPI_Comm_spawn (nsp_exe,args,1,INFO_NULL,0, COMM_SELF);
// child will execute cmd 
[NEWORLD] = MPI_Intercomm_merge (children, 0);

// The same but whe use a nsp macro 

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(1);

//-------------------------------------------------------
// _Isend: Immediate (non-blocking) send
// must _Test/_Wait for it to progress (at least under tcp rpi)
//-------------------------------------------------------
			// child 
			stat=MPI_Probe(-1,-1,NEWORLD)	// blocks

A=7;B="poo";C=ones(70,1024);
[rqa]=MPI_Isend(A,1,7,NEWORLD);		// returns immediately/unblocks
[rqb]=MPI_Isend(B,1,8,NEWORLD);		// that doesn't mean
[rqc]=MPI_Isend(C,1,9,NEWORLD);		// xmit has progressed

			[stata]=MPI_Probe(0,7,NEWORLD)	// tag 7 queued
			[statb]=MPI_Probe(0,8,NEWORLD)	// tag 8 too
			[statc]=MPI_Probe(0,9,NEWORLD)	// all them

stats = MPI_Testall ({rqa,rqb,rqc});
[idx,stat] = MPI_Testany ({rqa,rqb,rqc})// ANY of them? yes, progress
help MPI_PROC_NULL				// flag==1, idx=0 (1st msg)
help MPI_UNDEFINED				// notice len==4 ?!?

[idxs,stats] = MPI_Testsome ({rqa,rqb,rqc})	// rqa 
[idxs,stats] = MPI_Testsome ({rqa,rqb,rqc})	// rqb 

[idx, stats] = MPI_Waitany({rqa,rqb,rqc})	// blocks due to rqc

			// nsp child 
			[stat]=MPI_Probe(-1,-1,NEWORLD)	// tag 7 
			A=10;[stat]=MPI_Recv(A,0,7,NEWORLD)
			A==7 
			
			[stat]=MPI_Probe(-1,-1,NEWORLD)	// tag 8 
			B="abc";[stat]=MPI_Recv(B,0,8,NEWORLD)
			B=="poo"
			
			[stat]=MPI_Probe(-1,-1,NEWORLD)	// tag 9 
			C=0*ones(70,1024);
			[stat]=MPI_Recv(C,0,9,NEWORLD)	// unblocks the MPI_Waitany
			A=0;B="abc";C=0*ones(70,1024);			
			
			// XXX Note that A,B,C should be protected against access 
			// when a Irecv is initialized to prevent nsp crash. 
			
			[rqa]=MPI_Irecv (A,0,7,NEWORLD)	// 3 non-blcking
			[rqb]=MPI_Irecv (B,0,8,NEWORLD)	// concurrent
			[rqc]=MPI_Irecv (C,0,9,NEWORLD)	// requests
			
			[flag]=MPI_Iprobe(-1,-1,NEWORLD)	// empty queue
			[idx,stat]=MPI_Waitany({rqa,rqb,rqc}) // blocks

MPI_Send(A,1,7,NEWORLD)			// won't block, unblocks Waitany

			idx==0, stat.src==0, stat.tag==7,
			A-7 
			[idxs stats]=MPI_Waitsome({rqa,rqb,rqc}); // blocks

MPI_Send(C,1,9,NEWORLD)			// won't block, unblks Waitsome

			idxs==2, stats.src==0, stats.tag==9, //
			C(1)==1
			
			[stats]=MPI_Waitall({rqa,rqb,rqc}); // blocks for b 

MPI_Send(A,1,7,NEWORLD)	// won't block(tcp rpi), neither unblock(tag 7)
MPI_Send(B,1,8,NEWORLD)	// won't block, unblocks Waitall

			B == "poo";
			
			// only stats(2) should be relevant 
			// 

			[stats]=MPI_Waitall({rqa,rqb,rqc}) // done
			stats 
			
			[stats]=MPI_Testall({rqa,rqb,rqc}) // flag==1
			stats 
						
			[stat]=MPI_Iprobe(-1,-1,NEWORLD) // recall A
			[stat]=MPI_Recv(A,0,7,NEWORLD)
			[stat]=MPI_Iprobe(-1,-1,NEWORLD) // emptyqueue

			MPI_Finalize()		// too large for just 1 session
			quit			// jump over the next...
MPI_Finalize()					// ...Init cycle if you want
quit						// (not tired yet? great! :-)

//-------------------------------------------------------
// _Ibsend: Immediate (non-blocking) buffered (user provided) send
// User provides buffer, instead of system buffer
// Init cycle again... skip it if you didn't Finalize just above
//-------------------------------------------------------

putenv(['LAM_MPI_SSI_rpi=tcp']), MPI_Init	// better understood with tcp

args={}					// choose one of these
args={'-display',getenv('DISPLAY')}	// depending on your rsh/ssh environ

args={args{:},...				// xterm switches for DISPLAY
	'-e','matlab','-nosplash','-nojvm',...	// matlab switches for NoJava
	'-r','startup_mergeParent'}		// matlab switches to run script

[info children errs] = MPI_Comm_spawn ('/usr/X11R6/bin/xterm',args,...
					1,'NULL',0,'SELF')
[info NEWORLD] = MPI_Intercomm_merge (children, 0)	// child unblocks
MPI_Errhandler_set(NEWORLD,'RETURN')			// safer this way

						// some assorted vars
A='Repeat this string until _8KB_! '		//    8KB
A=repmat(A,1,128);
		s=whos('A'), sa=s.bytes
B=ones(8,1024); s=whos('B'), sb=s.bytes		//   64KB
C=ones(9,1024); s=whos('C'), sc=s.bytes		// > 64KB (73728 bytes)

			// //////////// ON CHILD ////////////////////////////////////////////////////////////
			A=repmat(' ',1,4096);	// recv counterparts
			B=zeros(8,1024);
			C=zeros(9,1024); whos

help nonblocking

//-------------------------------------------------------
// Init cycle done... jump here if you didn't Finalize just above
//-------------------------------------------------------
help MPI_Ibsend

whos, sa, sb, sc				// A,B,C, -- 8,64,72 KB
sz=sa+sb+2*MPI_BSEND_OVERHEAD			// room for A & B, later C
sc+MPI_BSEND_OVERHEAD < sz
snd=repmat(' ',1,sz/2); s=whos('snd')		// around 72KB
info = MPI_Buffer_attach(snd)			// for Ibsend

			// //////////////// under child Nsp ////////////////////////////////////////
			whos
			[info flag stat]=MPI_Iprobe(-1,-1,NEWORLD) // emptyqueue

[info rqa]=MPI_Ibsend(A,1,7,NEWORLD)		// doesn't block or progress
[info rqb]=MPI_Ibsend(B,1,8,NEWORLD)		// neither
[info rqc]=MPI_Ibsend(C,1,9,NEWORLD)		// neither
[info      stats]=MPI_Waitall(rqa,rqb,rqc)	// nothing to wait
stats(1), stats(2), stats(3)			// len==0
stats(1).src==MPI_PROC_NULL
stats(2).tag==MPI_UNDEFINED
stats(3).err==MPI_SUCCESS

			[info      stat]=MPI_Probe(0,7,NEWORLD)	// tag 7 queued
			[info flag stat]=MPI_Iprobe(0,8,NEWORLD)// tag 8 too
			[info flag stat]=MPI_Iprobe(0,9,NEWORLD)// tag 9 as well

			A(1:30)=deal(' '); A(1:60)	// cleared to
			B(:,1:2)=deal(0); B(:,1:3)	// detect when
			C(:,1:2)=deal(0); C(:,1:3)	// it's recvd

			[info rqa]=MPI_Irecv (A,0,7,NEWORLD)	// 3 non-blcking
			[info rqb]=MPI_Irecv (B,0,8,NEWORLD)	// concurrent
			[info rqc]=MPI_Irecv (C,0,9,NEWORLD)	// requests

			A(1:60), B(:,1:3), C(:,1:3)	// C not yet... why?
							// yup, snd buffer
			[info idxs stats]=MPI_Testsome(rqa,rqb,rqc) // idxs 0/1
			stats(1), stats(2)			    // tags 7/8

			[info flag stats]=MPI_Testall(rqa,rqb,rqc) // flag==0
			[info      stats]=MPI_Waitall(rqa,rqb,rqc) // blocks!!!

[info flag stat]=MPI_Test(rqc)				// must push it :-)

			stats(1), stats(2), stats(3)	// ok, unblocked rqc
			[info flag stat]=MPI_Iprobe(-1,-1,NEWORLD)	// done
			clear flag idxs info rqa rqb rqc stat stats
			whos

info = MPI_Buffer_detach
clear errs flag info rqa rqb rqc s snd stat stats sz
whos

MPI_Finalize					// too large for just 1 session
quit						// either Init cycle or jump
			MPI_Finalize		// over to the next session
			quit

//-------------------------------------------------------
// _Issend: Immediate (non-blocking) synchronous (waits for matching Recv) send
// Rather self-contradictory :-) basically sending without buffers
// paradoxically, far easier to understand ;-)
// Init cycle again... skip it if you didn't Finalize just above
//-------------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(1);

						// some assorted vars
A='Repeat this string until _8KB_! '		//    8KB
A=repmat(A,1,128);
		s=whos('A'), sa=s.bytes
B=ones(8,1024); s=whos('B'), sb=s.bytes		//   64KB
C=ones(9,1024); s=whos('C'), sc=s.bytes		// > 64KB (73728 bytes)

			// //////////// ON CHILD ////////////////////////////////////////////////////////////
			A=repmat(' ',1,4096);	// recv counterparts
			B=zeros(8,1024);
			C=zeros(9,1024); whos

help nonblocking

//-------------------------------------------------------
// Init cycle done... jump here if you didn't Finalize just above
//-------------------------------------------------------
help MPI_Issend
			// //////////////// under child Nsp ////////////////////////////////////////
			whos
			[info flag stat]=MPI_Iprobe(-1,-1,NEWORLD) // emptyqueue

[info rqa]=MPI_Issend(A,1,7,NEWORLD)		// doesn't block or progress
[info rqb]=MPI_Issend(B,1,8,NEWORLD)		// neither
[info rqc]=MPI_Issend(C,1,9,NEWORLD)		// neither
[info idx  stat ]=MPI_Waitany(rqa,rqb,rqc)	// blocks, of course

			[info flag stat]=MPI_Iprobe(0,7,NEWORLD)// tag 7 queued
			[info flag stat]=MPI_Iprobe(0,8,NEWORLD)// tag 8 too
			[info flag stat]=MPI_Iprobe(0,9,NEWORLD)// tag 9 as well

			A(1:30)=deal(' '); A(1:60)	// cleared to
			B(:,1:2)=deal(0); B(:,1:3)	// detect when
			C(:,1:2)=deal(0); C(:,1:3)	// it's recvd

			[info stat]=MPI_Recv(C,0,9,NEWORLD)	// unblks idx=2
			C(:,1:3)			// any order would do
			B(:,1:3), A(1:60)

[info idxs stats]=MPI_Waitsome(rqa,rqb,rqc)		// blocks rqa/b

			[info rqa]=MPI_Irecv (A,0,7,NEWORLD)	// unblks idxs=0
			B(:,1:3), A(1:60)

[info      stats]=MPI_Waitall(rqa,rqb,rqc)		// blocks due to rqb

			[info stat]=MPI_Recv (B,0,8,NEWORLD)	// unblocks
			B(:,1:3)

stats(1), stats(2), stats(3)

			[info flag stat]=MPI_Iprobe(-1,-1,NEWORLD)	// done
			clear flag info rqa stat
			whos

clear errs idx idxs info rqa rqb rqc s snd stat stats
whos

MPI_Finalize					// too large for just 1 session
quit						// either Init cycle or jump
			MPI_Finalize		// over to the next session
			quit

//-------------------------------------------------------
// _Irsend: Immediate (non-blocking) ready ("requires" matching Recv) send
// Init cycle again... skip it if you didn't Finalize just above
//-------------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(1);
						// some assorted vars
A='Repeat this string until _8KB_! '		//    8KB
A=repmat(A,1,128);
		s=whos('A'), sa=s.bytes
B=ones(8,1024); s=whos('B'), sb=s.bytes		//   64KB
C=ones(9,1024); s=whos('C'), sc=s.bytes		// > 64KB (73728 bytes)

			// //////////// ON CHILD ////////////////////////////////////////////////////////////
			A=repmat(' ',1,4096);	// recv counterparts
			B=zeros(8,1024);
			C=zeros(9,1024); whos

//-------------------------------------------------------
// Init cycle done... jump here if you didn't Finalize just above
//-------------------------------------------------------
help nonblocking

help MPI_Irsend
			// //////////////// under child Nsp ////////////////////////////////////////
			whos
			[info flag stat]=MPI_Iprobe(-1,-1,NEWORLD) // emptyqueue

[info rqa]=MPI_Irsend(A,1,7,NEWORLD)		// doesn't block or progress
[info rqb]=MPI_Irsend(B,1,8,NEWORLD)		// neither
[info rqc]=MPI_Irsend(C,1,9,NEWORLD)		// neither
[info idx  stat ]=MPI_Waitany (rqa,rqb,rqc)	// idx ==0 done
[info idxs stats]=MPI_Waitsome(rqa,rqb,rqc)	// idxs==1 done
[info      stats]=MPI_Waitall (rqa,rqb,rqc)	// blocks due to rqc

			[info flag stat]=MPI_Iprobe(0,7,NEWORLD)// tag 7 queued
			[info flag stat]=MPI_Iprobe(0,8,NEWORLD)// tag 8 too
			[info flag stat]=MPI_Iprobe(0,9,NEWORLD)// tag 9 as well

			A(1:30)=deal(' '); A(1:60)	// cleared to
			B(:,1:2)=deal(0); B(:,1:3)	// detect when
			C(:,1:2)=deal(0); C(:,1:3)	// it's recvd

			[info stat]=MPI_Recv(C,0,9,NEWORLD)	// unblks Wtall
			C(:,1:3)
			B(:,1:3), A(1:60)

stats(1), stats(2), stats(3)
[info      stats]=MPI_Waitall (rqa,rqb,rqc)		// can't wait anymore

			[info stat]=MPI_Recv(A,0,7,NEWORLD)
			B(:,1:3), A(1:60)
			[info stat]=MPI_Recv(B,0,8,NEWORLD)
			B(:,1:3)

			[info flag stat]=MPI_Iprobe(-1,-1,NEWORLD)	// done
			MPI_Finalize
			quit
MPI_Finalize
quit


// ---------------------------------------------------
// Concurrent Send/Recv
// MPI_Sendrecv, _Sendrecv_replace
// 	*** WE ARE USING 3 COMPUTERS THIS TIME :-) ***
// ---------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(2);

					// string hostname same size everywhere
[info rank]=MPI_Comm_rank(NEWORLD)
[stat, hnam]=system('hostname')

hnam(1,60)=' '
A(1:60)=deal(' ')
s=whos('hnam'); if s.bytes~=120 disp(['take it easy, dude!!!' char(7)]), end

		// //////// ON FIRST CHILD //////////////
		[info rank]=MPI_Comm_rank(NEWORLD)
		[stat, hnam]=system('hostname')

		hnam(1,60)=' '
		A(1:60)=deal(' ')
		s=whos('hnam'); if s.bytes~=120 disp(['slower!' char(7)]), end

				// //////// ON SECOND CHILD //////////////
				[info rank]=MPI_Comm_rank(NEWORLD)
				[stat, hnam]=system('hostname')

				hnam(1,60)=' '
				A(1:60)=deal(' ')
				s=whos('hnam')
				if s.bytes~=120 disp(['linewise' char(7)]), end

help modes

//-------------------------------------------------------
// Can Sendrecv from same child
//-------------------------------------------------------
help MPI_Sendrecv
[info stat]=MPI_Sendrecv(hnam,1,7, A,1,8, NEWORLD), A	// blocks (recv part)

		// //////// ON FIRST CHILD //////////////
		[info stat]=MPI_Probe(-1,-1,NEWORLD)
		[info stat]=MPI_Recv   (A,0,7,NEWORLD), A
		 info      =MPI_Send(hnam,0,8,NEWORLD)	// now it unblocks

//-------------------------------------------------------
// It is ...recv what really blocks Sendrecv, send... gets buffered
//-------------------------------------------------------
			////////// ON SECOND CHILD //////////////
			 info      =MPI_Send(hnam,0,9,NEWORLD)
			[info stat]=MPI_Recv(   A,0,7,NEWORLD), A	// blcks
A
[info stat]=MPI_Sendrecv(hnam,2,7, A,2,9, NEWORLD)	// doesn't block, unblks
A

//-------------------------------------------------------
// Can Sendrecv from/to different children (doesnot cascade :-)
//-------------------------------------------------------
			////////// ON SECOND CHILD //////////////
			[info stat]=MPI_Recv(A,0,7,NEWORLD), A	// blocks

[info stat]=MPI_Sendrecv(hnam,2,7, A,1,8, NEWORLD), A	// blocks(recv)
							// _and_ unblocks 2nd
		////////// ON FIRST CHILD //////////////
		info=MPI_Send(hnam,0,8,NEWORLD)		// unblocks, no cascade

//-------------------------------------------------------
// Can Sendrecv against Sendrecv
//-------------------------------------------------------
A(1:10)=deal(' ')
[info stat]=MPI_Sendrecv(hnam,1,7, A,1,8, NEWORLD), A	// blocks (recv)

		////////// ON FIRST CHILD //////////////
		A(1:15)=deal(' ')			// unblocks
		[info stat]=MPI_Sendrecv(hnam,0,8, A,0,7, NEWORLD), A

//-------------------------------------------------------
// Circular Sendrecv
//-------------------------------------------------------
A(1:10)=deal(' ')
[info stat]=MPI_Sendrecv(hnam,1,7, A,2,9, NEWORLD), A	// blocks (recv)

		////////// ON FIRST CHILD //////////////
		A(1:15)=deal(' ')			// doesn't block
		[info stat]=MPI_Sendrecv(hnam,2,8, A,0,7, NEWORLD), A

			////////// ON SECOND CHILD //////////////
			A(1:15)=deal(' ')			// unblocks
			[info stat]=MPI_Sendrecv(hnam,0,9, A,1,8, NEWORLD), A

//-------------------------------------------------------
// Can Sendrecv in-place (same buffer), as long as there is room enough for recv
//-------------------------------------------------------
A(1:10)=deal(' ')						// blocks (recv)
[info stat]=MPI_Sendrecv_replace(hnam,1,7, 2,9, NEWORLD), hnam

		////////// ON FIRST CHILD //////////////
		A(1:15)=deal(' ')				// doesn't block
		[info stat]=MPI_Sendrecv_replace(hnam,2,8, 0,7, NEWORLD), hnam

			////////// ON SECOND CHILD //////////////
			A(1:15)=deal(' ')
			[info stat]=MPI_Recv(   A,1,8,NEWORLD), A
			 info      =MPI_Send(hnam,0,9,NEWORLD)	// unblocks
			MPI_Finalize
			quit
		MPI_Finalize
		quit
MPI_Finalize
quit

// ============================================
// Recv cancellation (no-op for Send)
// MPI_Cancel, _Test_cancelled
// ============================================

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(1);

						// some assorted vars
A='Repeat this string until _8KB_! '		//    8KB
A=repmat(A,1,128);
		s=whos('A'), sa=s.bytes
B=ones(8,1024); s=whos('B'), sb=s.bytes		//   64KB
C=ones(9,1024); s=whos('C'), sc=s.bytes		// > 64KB (73728 bytes)
whos

			////////////// ON CHILD ////////////////////////////////////////////////////////////
			A=repmat(' ',1,4096);	// recv counterparts
			B=zeros(8,1024);
			C=zeros(9,1024); whos

//-------------------------------------------------------
// Can cancel only receives, not sends
//-------------------------------------------------------
help cancel

help MPI_Cancel
help MPI_Test_cancelled
			[info rqa]=MPI_Irecv (A,0,7,NEWORLD)	// 3 non-blcking
			[info rqb]=MPI_Irecv (B,0,8,NEWORLD)	// concurrent
			[info rqc]=MPI_Irecv (C,0,9,NEWORLD)	// requests

			[info flag]=MPI_Iprobe(-1,-1,NEWORLD)	// empty queue
			 info      =MPI_Cancel(rqb)		// cancel B
			[info idx  stat ]=MPI_Waitany(rqa,rqb,rqc) // idx==1
			stat.src==MPI_CANCEL_SOURCE
			stat.tag==MPI_UNDEFINED
			stat.err==MPI_SUCCESS
			[info flag statb]=MPI_Test_cancelled(stat)

			//////////////////////// notice how cancel info "washes out"
			stat, statb
			[info      stat ]=MPI_Wait(rqb)		// won't wait
			stat.src==MPI_ANY_SOURCE
			stat.tag==MPI_ANY_TAG
			stat.err==MPI_ERR_PENDING
			[info flag stat2]=MPI_Test_cancelled(stat ) // flag==0
			[info flag stat2]=MPI_Test_cancelled(statb) // flag==1

			[info idx  stat ]=MPI_Waitany(rqa,rqb,rqc) // blocks

 info     =MPI_Send (A,1,7,NEWORLD)		// won't block, unblks Waitany
 info     =MPI_Send (B,1,8,NEWORLD)		// not expected at receiver
[info rqc]=MPI_Isend(C,1,9,NEWORLD)		// can't be cancelled? sure?
 info     =MPI_Cancel(rqc)			// no-op for send ?!? info==0
[info flag statc] = MPI_Test(rqc)		// certainly, not completed
[info flag stat ] = MPI_Test(rqc)		// info wanishing out
[info flag stat ] = MPI_Test_cancelled(statc)	// nope, not cancelled
//[info     stat ] = MPI_Wait(rqc)		// don't do that, would block

			[info flag stat]=MPI_Iprobe(0,-1,NEWORLD) // tag 8 1st
			[info flag stat]=MPI_Iprobe(0, 9,NEWORLD) // no progress
			[info flag stat]=MPI_Iprobe(0, 7,NEWORLD) // noprogress

			 info      =MPI_Cancel(rqa)		// info==13
			 info     ==MPI_ERR_ARG
			 info      =MPI_Cancel(rqc)		// cancelled!!!
			[info flag statc] = MPI_Test(rqc)	// not completed
			[info flag stat ] = MPI_Test_cancelled(statc) // nope
		//////////	[info      stat ] = MPI_Wait(rqc)	// would block

 info     =MPI_Cancel(rqc)			// double-cancelled ?!? info==0

			[info flag stats] = MPI_Testall (rqa,rqb,rqc) // not ALL
			stats(1), stats(2), stats(3)
			[info flag stat ] = MPI_Test(rqa)	// survivor
			A(1:60)					// there it is
			B(:,1:3), C(:,1:3)			// untouched

			 clear stat stats
			[info flag stat ] = MPI_Test(rqa)	// completed
			[info flag stat ] = MPI_Test_cancelled(stat) // nope
			[info      stat ] = MPI_Wait(rqa)	// won't wait

			[info flag stat ] = MPI_Test(rqb)	// cancel info
			[info flag stat ] = MPI_Test_cancelled(stat) // vanished
			[info flag stat ] = MPI_Test_cancelled(statb)
			[info      stat ] = MPI_Wait(rqb)	// won't wait

			[info flag stat ] = MPI_Test(rqc)	// not completed
			[info flag stat ] = MPI_Test_cancelled(stat)  // not
			[info flag stat ] = MPI_Test_cancelled(statc) // really
		////////////	[info      statc] = MPI_Wait(rqc)	// would block

			[info flag stat ] = MPI_Iprobe(-1,-1,NEWORLD)
			[info      stat ] = MPI_Wait(rqb)	// didn't block
			B(:,1:3)			// recall unexpected
			[info      stat ] = MPI_Recv(B,0,8,NEWORLD)
			B(:,1:3)

			MPI_Finalize
			quit
MPI_Finalize
quit

// ===========================================================
// Persistent requests
// MPI_[-BSR]send_init, _Recv_init, _Start[all], _Request_free
// MPI_Address useful at last! :-)
// ===========================================================

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(1);

						// some assorted vars
A='Repeat this string until _8KB_! '		//    8KB
A=repmat(A,1,128);
		s=whos('A'), sa=s.bytes
B=ones(8,1024); s=whos('B'), sb=s.bytes		//   64KB
C=ones(9,1024); s=whos('C'), sc=s.bytes		// > 64KB (73728 bytes)
whos

			////////////// ON CHILD ////////////////////////////////////////////////////////////
			A=repmat(' ',1,4096);	// recv counterparts
			B=zeros(8,1024);
			C=zeros(9,1024); whos

//-------------------------------------------------------
// Persistent requests: concept (Recv)
//-------------------------------------------------------
help persistent

help MPI_Recv_init
help MPI_Start
help MPI_Request_free

help eval					// get 3 copies of A
for n = 1:3					// called A1, A2 and A3
  eval(['A' num2str(n) ' = A;'])		// we will receive them
end						// repeatedly over CHILD's A
whos('A*')

[info addr] = MPI_Address (A)			// MPI_Address useful at last!
for n = 1:3
  m = num2str(n);
  eval( [ '[info addr' m ']=MPI_Address(A' m ')' ] )
end
isequal(addr,addr1,addr2,addr3)			// Nsp lazy copy mechanism
dump(A)						// Peter Boettcher's dump prog.
						// notice link->A3->A2->A1->A
A1(1:10)='1st copy! '; A1(1:30)			// new unshared var on writes
addr==addr1
[info addr1] = MPI_Address (A1)			// A1 address changed on write
addr==addr1
dump(A1)					// link -> (nil)
dump(A2)					// link -> A -> A3 -> A2 itself

A2(1:10)='2nd copy! '; A2(1:30)			// Ok, 3 unshared vars now
[info addr2] = MPI_Address (A2)			// A3 still shared with A
addr==addr2

			////////////////// ON CHILD //////////////////////////////////////////////////////
			[info req]=MPI_Recv_init(A,0,7,NEWORLD)
			[info flag stat]=MPI_Test(req)	// flag==1, done
			stat.src==MPI_ANY_SOURCE
			stat.tag==MPI_ANY_TAG
			stat.err==MPI_ERR_PENDING
			[info flag stat]=MPI_Test(req)	// nothing happens
			[info      stat]=MPI_Wait(req)	// doesn't block

			for n = 1:3
			  info           =MPI_Start(req)	// active
			 [info flag stat]=MPI_Test (req)	// flag==0
			 [info      stat]=MPI_Wait (req)	// blocks
			  A(1:30)				// see results
			end

info = MPI_Send(A1,1,7,NEWORLD)		// won't block, unblks Wait, A='1st...'
info = MPI_Send(A2,1,7,NEWORLD)		// won't block, unblks Wait, A='2nd...'
info = MPI_Send(A3,1,7,NEWORLD)		// won't block, unblks Wait, A='Repeat'

			info = MPI_Request_free (req)		// cool, isn't?
			info = MPI_Request_free (req)		// info==7
			info== MPI_ERR_REQUEST
			help   MPI_ERR_REQUEST

//-------------------------------------------------------
// Persistent requests: concept (Send)
//-------------------------------------------------------
help persistent

help MPI_Send_init

[info req]=MPI_Send_init(C,1,9,NEWORLD)	// C > 64KB (tcp_short)
[info flag stat]=MPI_Test(req)		// flag==1, done
stat.src==MPI_ANY_SOURCE
stat.tag==MPI_ANY_TAG
stat.err==MPI_ERR_PENDING
[info flag stat]=MPI_Test(req)		// nothing happens
[info      stat]=MPI_Wait(req)		// doesn't block

			[info flag stat ] = MPI_Iprobe(-1,-1,NEWORLD) // clean
for n = 1:3
  C(:,1:3)=deal(n);			// C > 64KB (tcp_short)
  info           =MPI_Start(req)	// activate, will block (tcp rpi)
  info           =MPI_Start(req)	// can't start next, wait for this
  info          ==MPI_ERR_REQUEST
 [info msg]=MPI_Error_string(info)
 [info flag stat]=MPI_Test (req)	// flag==0, not completed
 [info      stat]=MPI_Wait (req)	// blocks
end
			[info stat]=MPI_Probe(-1,-1,NEWORLD)	// envelope here
					      C(:,1:3)		// C==0
			[info stat]=MPI_Recv (C,0,9,NEWORLD)	// unblocks
					      C(:,1:3)		// C==1
			[info stat]=MPI_Recv (C,0,9,NEWORLD)
					      C(:,1:3)		// C==2
			[info stat]=MPI_Recv (C,0,9,NEWORLD)
					      C(:,1:3)		// C==3
			[info flag stat ] = MPI_Iprobe(-1,-1,NEWORLD) // clean

info = MPI_Request_free (req)		// how about B==64KB ???
[info req]=MPI_Send_init(B,1,8,NEWORLD)	// B = 64KB (tcp_short)

for n = 1:3
  B(:,1:3)=deal(n);			// C > 64KB (tcp_short)
  info           =MPI_Start(req)	// activate, will block (tcp rpi)
  info           =MPI_Start(req)	// can't start next, wait for this
  info          ==MPI_ERR_REQUEST	// 1st turn / 2nd turn / 3rd turn (blks)
 [info flag stat]=MPI_Test (req)	// flag==1  / flag==0  / flag=0 (blocks)
  stat.src==MPI_PROC_NULL
  stat.tag==MPI_UNDEFINED
  stat.err==MPI_SUCCESS
 [info      stat]=MPI_Wait (req)	// won't blk / won't   /)
  stat.src==[MPI_ANY_SOURCE,	MPI_PROC_NULL]
  stat.tag==[MPI_ANY_TAG,	MPI_UNDEFINED]
  stat.err==[MPI_ERR_PENDING,	MPI_SUCCESS]
end
			[info stat]=MPI_Probe(-1,-1,NEWORLD)	// unblocks!!!
					      B(:,1:3)		// B==0
			[info stat]=MPI_Recv (B,0,8,NEWORLD)	// unblocks
					      B(:,1:3)		// B==1
			[info stat]=MPI_Recv (B,0,8,NEWORLD)
					      B(:,1:3)		// B==2
			[info stat]=MPI_Recv (B,0,8,NEWORLD)
					      B(:,1:3)		// B==3
			[info flag stat ] = MPI_Iprobe(-1,-1,NEWORLD) // clean

info = MPI_Request_free (req)

MPI_Finalize					// too large for just 1 session
quit						// either Init cycle or jump
			MPI_Finalize		// over to the next session
			quit

//-------------------------------------------------------
// Persistent requests: modes ([BSR]send)
// Several requests at once: _Startall
// Init cycle again... skip it if you didn't Finalize just above
//-------------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(1);
						// some assorted vars
A='Repeat this string until _8KB_! '		//    8KB
A=repmat(A,1,128);
		s=whos('A'), sa=s.bytes
B=ones(8,1024); s=whos('B'), sb=s.bytes		//   64KB
C=ones(9,1024); s=whos('C'), sc=s.bytes		// > 64KB (73728 bytes)
whos

			////////////// ON CHILD ////////////////////////////////////////////////////////////
			A=repmat(' ',1,4096);	// recv counterparts
			B=zeros(8,1024);
			C=zeros(9,1024); whos

//-------------------------------------------------------
// Init cycle done... jump here if you didn't Finalize just above
//-------------------------------------------------------
help persistent

help MPI_Bsend_init
help MPI_Rsend_init
help MPI_Ssend_init
help MPI_Startall

sz=sa+MPI_BSEND_OVERHEAD				// 8232 bytes
snd=' '; snd=repmat(snd,1,sz/2); s=whos('snd')		// room for Bsend(A)
info=MPI_Buffer_attach(snd)				// user buffer

[info rqa]=MPI_Bsend_init(A,1,7,NEWORLD)
[info rqb]=MPI_Rsend_init(B,1,8,NEWORLD)
[info rqc]=MPI_Ssend_init(C,1,9,NEWORLD)

			[info flag stat ] = MPI_Iprobe(-1,-1,NEWORLD) // clean

info=MPI_Startall(rqa,rqb,rqc)				// envelopes snt eagerly

			[info flag stat ] = MPI_Iprobe(-1,-1,NEWORLD) // A,tag 7
			[info flag stat ] = MPI_Iprobe(-1, 8,NEWORLD) // B too
			[info flag stat ] = MPI_Iprobe(-1, 9,NEWORLD) // and C

[info idxs stats]=MPI_Testsome(rqa,rqb,rqc)		// A-B can progress
idxs==[0 1]						// A - it's buffered
[stats.src]==MPI_PROC_NULL				// B - it's <= 64KB
[stats.tag]==MPI_UNDEFINED
[stats.err]==MPI_SUCCESS
[info idxs stats]=MPI_Testsome(rqa,rqb,rqc)		// idxs=[], C can't
[info flag stats]=MPI_Testall (rqa,rqb,rqc)		// not ALL

			A(1:60), B(:,1:3), C(:,1:3)	// none recv yet

			[info rqa]=MPI_Irecv (A,0,7,NEWORLD)	// 3 non-blcking
			[info rqb]=MPI_Irecv (B,0,8,NEWORLD)	// concurrent
			[info rqc]=MPI_Irecv (C,0,9,NEWORLD)	// requests
			
			A(1:60), B(:,1:3), C(:,1:3)	// A, B, just recvd

			[info flag stats]=MPI_Testall(rqa,rqb,rqc) // "pushing"
			stats(1), stats(2), stats(3)
			[info idxs stats]=MPI_Testsome(rqa,rqb,rqc) // idxs[0 1]
			stats(1), stats(2)			    // A & B

			A(1:60), B(:,1:3), C(:,1:3)	// nothing changed (A/B)

			[info idxs stats]=MPI_Waitsome(rqa,rqb,rqc)
							// blocks!!! :-)

[info stats]=MPI_Waitall(rqa,rqb,rqc)		// gosh, hard to push ;-)
stats(1), stats(2), stats(3)

			C(:,1:3)
			[info flag stat ] = MPI_Iprobe(-1,-1,NEWORLD) // clean

			info=MPI_Request_free(rqa)	// not usually done
			info=MPI_Request_free(rqb)	// for Irecv/I[brs]send
			info=MPI_Request_free(rqc)	// see help page
			help MPI_Request_free
			 info    ==MPI_ERR_REQUEST
			[info msg]=MPI_Error_string(info)
			MPI_Finalize
			quit

info = MPI_Buffer_detach

info=MPI_Request_free(rqa)			// usually done for MPI_*_init
info=MPI_Request_free(rqb)
info=MPI_Request_free(rqc)

if rqc==MPI_REQUEST_NULL, disp('NULL'), end
[info flag stat]=MPI_Test (rqc)			// graceful testing on NULL
[info      stat]=MPI_Wait (rqc)			// at least doesn't block
[info flag stat]=MPI_Test_cancelled(stat)	// no, cancel ~= free

MPI_Finalize
quit

// ==============================================================
// Collective operation
// MPI_Barrier, _Bcast, _Scatter, _Scatterv, _Gather, _Gatherv
// MPI_Allgather, _Allgatherv, _Alltoall, _Alltoallv  
// MPI_Reduce, _Allreduce, _Scan, _Reduce_scatter 
// 	*** WE ARE USING 3 COMPUTERS FROM NOW ON :-) ***
// ==============================================================

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(2);

[rank]=MPI_Comm_rank(NEWORLD)

		// ////////// ON FIRST CHILD //////////////
		[rank]=MPI_Comm_rank(NEWORLD)

				// ////////// ON SECOND CHILD //////////////
				[rank]=MPI_Comm_rank(NEWORLD)

//-------------------------------------------------------
// Can synchronize all ranks in communicator group
//-------------------------------------------------------
				//////////// ON SECOND CHILD //////////////
				MPI_Barrier(NEWORLD)	// 2nd child blocks
		// //////// ON FIRST CHILD //////////////
		MPI_Barrier(NEWORLD)		// 1st child blocks
// ////////// ON PARENT Nsp //////////////
MPI_Barrier(NEWORLD)				// unblocks both

//-------------------------------------------------------
// Can broadcast a msg to all ranks in communicator group
//-------------------------------------------------------
				//////// ON SECOND CHILD //////////////
				C=[1 3 5 7]
				MPI_Bcast(C,0,NEWORLD) // blocks
				
				
////////// ON PARENT Nsp //////////////
A=[2 4;8 16]
MPI_Bcast(A,0,NEWORLD)				// unblocks

		////////// ON FIRST CHILD //////////////
		B=[2;4;6;8]
		MPI_Bcast(B,0,NEWORLD) //	// doesn't block

//-------------------------------------------------------
// The three nsp get the same data [2 4;8 16] stored in 
// matrices whose shape depends on the given argument 

//-------------------------------------------------------

//-------------------------------------------------------
// Can scatter a msg among all ranks in communicator group
// variant with displacements, too
// send blocks of size R to each group member from an element of root 
// giving a proper A is important only for the sender. 
// 
//-------------------------------------------------------

help MPI_Scatter
				//////////// ON SECOND CHILD //////////////
				R=[7 8]
				MPI_Scatter([],R,0,NEWORLD);	// blocks

//////////// ON PARENT Nsp //////////////
A=[1 2;3 4;5 6], R=[7;8]
MPI_Scatter(A,R,0,NEWORLD);	// unblocks

		//////////// ON FIRST CHILD //////////////
		R=[7;8]
		MPI_Scatter([],R,0,NEWORLD);		// doesn't blk

//-------------------------------------------------------
// Parent had	[7;8]	gets	[1;3]	from	[1 2;
// child1 had	[7;8]	gets	[5;2]		 3 4;
// child2 had	[7 8]	gets	[4 6]		 5 6]	notice row-major order
//-------------------------------------------------------
// variant with displacements, too
//-------------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(2);

help MPI_Scatterv

		// ON FIRST CHILD
		R=ones(1,4);
		MPI_Scatterv([],[],[],R,0,NEWORLD), R	// blocks

// On parent nsp 
A=[0 1 2 3 4 5 6 7]
R=ones(1,4);
MPI_Scatterv(A,[1 2 3],[0 3 5],R,0,NEWORLD)		// unblocks
R							

				// ON SECOND CHILD 
				R=ones(1:4);
				MPI_Scatterv([],[],[],R,0,NEWORLD), R

//-------------------------------------------------------
// Parent gets 1 value from [0 1 2 3 4 5 6 7] at position 0 
// child1 gets 2 values from [0 1 2 3 4 5 6 7] at position 3 
// child2 gets 3 values from [0 1 2 3 4 5 6 7] at position 5 
//-------------------------------------------------------

//-------------------------------------------------------
// Can gather a msg from among all ranks in communicator group
// variant with displacements, too
//-------------------------------------------------------

if %t then 
  exec src7.x/loader.sce 
  exec utils/NSP_spawn.sci 
  NEWORLD=NSP_spawn(2);
end


help MPI_Gather

		S=[4 5 6]					// doesn't
		MPI_Gather(S,[],0,NEWORLD)

A=ones(3,3), S=[1;2;3]
MPI_Gather(S,A,0,NEWORLD), A				// blocks

				S=[7 8 9]
				MPI_Gather(S,[],0,NEWORLD)	// unbl

//-------------------------------------------------------
// Parent sends	[1;2;3]	gets	[1 4 7
// child1 sends	[4 5 6]		 2 5 8
// child2 sends	[7 8 9]		 3 6 9]			notice row-major order
//-------------------------------------------------------
// variant with displacements, too
//-------------------------------------------------------

if %t then 
  exec src7.x/loader.sce 
  exec utils/NSP_spawn.sci 
  NEWORLD=NSP_spawn(2);
end

help MPI_Gatherv

A=ones(3,3), S=[1;2;3];
MPI_Gatherv(S,A,[3 3 6],[0 2 3],0,NEWORLD);		// blocks

		S=[4 5 6];					// doesn't
		MPI_Gatherv(S,[],[],[],0,NEWORLD);

				S=[7 9 7;8 8 6]
				MPI_Gatherv(S,[],[],[],0,NEWORLD);	// unbl

//-------------------------------------------------------
// Parent sends	[1;2;3]	gets	[1 7 8		notice overlapping parnt/chld1
// child1 sends	[4 5 6]		 2 8 7		(elm#2) and chld1/2 (elms#3,4)
// child2 sends	[7 9 7		 4 9 6]		higher rank overwrites lower's
//		 8 8 6]					notice row-major order
//-------------------------------------------------------

MPI_Finalize()						// ok, too large
quit							// cycle-Init or
		MPI_Finalize()				// skip to next session
		quit
				MPI_Finalize()
				quit

//-------------------------------------------------------
// Init cycle ... skip it if you didn't Finalize just above
//-------------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(2);

//-------------------------------------------------------
// Can gather in all ranks, not only at root
// variant with displacements, too
//-------------------------------------------------------

XX OK 
help MPI_Allgather

                rank=MPI_Comm_rank(NEWORLD);
		R=ones(2,3); S=[1,2]*(rank+1);
		MPI_Allgather(S,R,NEWORLD);		// blocks

rank=MPI_Comm_rank(NEWORLD);
R=ones(2,3); S=[1,2]*(rank+1);
MPI_Allgather(S,R,NEWORLD);		// blocks

                        rank=MPI_Comm_rank(NEWORLD);
		        R=ones(2,3); S=[1,2]*(rank+1);
		        MPI_Allgather(S,R,NEWORLD); // unblock all
			
//-------------------------------------------------------
// variant with displacements, too
//-------------------------------------------------------
help MPI_Allgatherv

A=ones(size(3)), S=[1;2;3]
info=MPI_Allgatherv(S,A,[3 3 6],[0 2 3],NEWORLD), A		// blocks

		B=linspace(0,0,9), S=[4 5 6]			// blocks
		info=MPI_Allgatherv(S,B,[3 3 6],[0 2 3],NEWORLD), B

			C=repmat(0,2,5), S=[7 9 7;8 8 6]	// unblocks
			info=MPI_Allgatherv(S,C,[3 3 6],[0 2 3],NEWORLD), C

//-------------------------------------------------------
// Parent had	[1 1 1	sends	[1;2;3]	gets	[1 7 8		// overlapping
//		 1 1 1				 2 8 7		// as seen in
//		 1 1 1]				 4 9 6]		// MPI_Gatherv
// child1 had	[0 0...	sends	[4 5 6]	gets	[1 2 4 7 8 9 8 7 6]
// child2 had	[0 0...	sends	[7 9 7	gets	[1 4 8 8 6
//		 0 0...		 8 8 6]		 2 7 9 7 0]	// greater room
//-------------------------------------------------------

A=zeros_deprectaed(3), S=[1 2 3]
info=MPI_Allgatherv(S,A,[3 2 3],[0 3 6],NEWORLD), A		// blocks

		B=zeros_deprectaed(3), S=[4 5]				// blocks
		info=MPI_Allgatherv(S,B,[3 2 3],[0 6 3],NEWORLD), B

			C=zeros_deprectaed(3), S=[7 8 9]			// unblocks
			info=MPI_Allgatherv(S,C,[3 2 3],[3 0 6],NEWORLD), C

//-------------------------------------------------------
// sizes (cnts in help file) must be equal in all ranks (or somebody is lying)
// in particular, varsize must correspond to declared size in cnts
// positions (disps in help file) need not be equal, can be altered at will
//-------------------------------------------------------

//-------------------------------------------------------
// Can scatter-gather from all to all
// variant with displacements, too
//-------------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(2);

// OK 

help MPI_Alltoall

		R=ones(1,3), S=[1 2 3]
		MPI_Alltoall(S,R,2,NEWORLD),R		// blocks
				
R=ones(1,3), S=[1 2 3]*10
MPI_Alltoall(S,R,2,NEWORLD), R				// blocks

                                R=ones(1,3), S=[1 2 3]*100
				MPI_Alltoall(S,R,2,NEWORLD),R	// unblocks

//-------------------------------------------------------
// variant with displacements, too
//-------------------------------------------------------
help MPI_Alltoallv

R=[0 0 0], S=[1 2 2 3 3 3]'
MPI_Alltoallv ( S,[1 2 3],[0 1 3],...
		R,[1 1 1],[2 1 0],NEWORLD), R			// blocks

		R=zeros_deprectaed(2), S=[1 4 3;4 3 3]			// blocks
		MPI_Alltoallv ( S,[1 2 3],[0 1 3],...
				R,[2 2 2],[0 2 1],NEWORLD), R

				R=zeros(4,3), S=[5 3;2 3;2 3]	// unblocks
				MPI_Alltoallv ( S,[1 2 3],[0 1 3],...
						R,[3 3 3],[1 5 9],NEWORLD),R

//-------------------------------------------------------
// Parent had	[0 0 0]	sends	[1 2 2 3 3 3]	gets	[5 1 1]
// child1 had	[0 0	sends	[1 2 3		gets	[2 2
//		 0 0]		 2 3 3]			 2 4]
// child2 had	[0 0 0	sends	[4 3		gets	[0 0 0
//		 0 0 0		 5 3			 3 3 3
//		 0 0 0		 5 3]			 3 3 3
//		 0 0 0]					 3 3 3]
//-------------------------------------------------------
// as with Allgatherv,
// sizes sent (scnts in help file) must be coherent with sizes recvd (rcnts)
// (or somebody is lying)
// positions (disps in help file) need not be equal, can be altered at will
// again, there can be overlapping, and higher ranks overwrite lower ones
//-------------------------------------------------------

MPI_Finalize						// ok, too large
quit							// cycle-Init or
		MPI_Finalize				// skip to next session
		quit
				MPI_Finalize
				quit

//-------------------------------------------------------
// Init cycle ... skip it if you didn't Finalize just above
//-------------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(2);

[info rank]=MPI_Comm_rank(NEWORLD)

		////////// ON FIRST CHILD //////////////
		[info rank]=MPI_Comm_rank(NEWORLD)

				////////// ON SECOND CHILD //////////////
				[info rank]=MPI_Comm_rank(NEWORLD)

//-------------------------------------------------------
// Init cycle done... jump here if you didn't Finalize just above
//-------------------------------------------------------

help collective

//-------------------------------------------------------
// Can reduce distributed values to its single (sum,prod,max,min...)
//-------------------------------------------------------

////////// ON PARENT Nsp //////////////

OK 

A=[-6 ,11 ,22]; R=ones(1,3);
MPI_Reduce(A,R,'MIN',0,NEWORLD);				// blocks

		// ON FIRST CHILD
		B=[20 ,-8 ,12]
		MPI_Reduce(B,[],'MIN',0,NEWORLD);		// doesn't

				// ON SECOND CHILD 
				C=[10 ,21 ,-10]
				MPI_Reduce(C,[],'MIN',0,NEWORLD);// unblocks

//-------------------------------------------------------
// Can perform the reduction in all, not just root
//-------------------------------------------------------

help MPI_Allreduce
OK 

R=ones(1,3);
MPI_Allreduce(A,R,'MAX',NEWORLD);			// blocks

		R=ones(1,3);
		MPI_Allreduce(B,R,'MAX',NEWORLD);	// DOES block

				R=ones(1,3);
				MPI_Allreduce(C,R,'MAX',NEWORLD); // unblocks

//-------------------------------------------------------
// NOTICE that rank 1 DOES get blocked... all them like root blocked before
//-------------------------------------------------------
// Can perform partial reductions in each rank (with its lesser ranks)
//-------------------------------------------------------

help MPI_Scan
OK 

R=ones(1,3);A=ones(1,3);
MPI_Scan(A,R,'SUM',NEWORLD);				// doesn't hang

				R=ones(1,3);C=ones(1,3);
				MPI_Scan(C,R,'SUM',NEWORLD);	// depends on...

		R=ones(1,3);B=ones(1,3);
		MPI_Scan(B,R,'SUM',NEWORLD);		// unblocks

//-------------------------------------------------------
// NOTICE that each rank depends on its lower ranks, blocking if not available
//-------------------------------------------------------
// And scatter the reduction after computing it
//-------------------------------------------------------
help MPI_Reduce_scatter

XXXXX 

A=[0 1 3;1 8 1], R=[1 1 1]
info=MPI_Reduce_scatter(A,R,[1 1 4],'PROD',NEWORLD), R		// blocks

		B=[1 4 1;0 3 2], R=[2;2;2]
		MPI_Reduce_scatter(B,R,[1 1 4],'PROD',NEWORLD)	// blocks too
		R

				C=[1 6 6;1 1 9], R=[3 3;3 3]	// unblocks
				MPI_Reduce_scatter(C,R,[1 1 4],'PROD',NEWORLD)
				R

//-------------------------------------------------------
// Parent had	[1 1 1]	sends	[0  1  3	gets	[0 1 1]
//				 1  8  1]
// child1 had	[2 2 2]'sends	[1  4  1	gets	[0 2 2]'
//				 0  3  2]
// child2 had	[3 3	sends	[1  6  6	gets	[24 18
//		 3 3]		 1  1  9]		 24 18]
//			prod is	[0 24 18
//				 0 24 18]
// NOTICE that counts must be coherent, as with Allgatherv / Alltoallv
//-------------------------------------------------------

R=[1 1 1]
MPI_Reduce_scatter(A,R,[2 2 2],'PROD',NEWORLD), R		// blocks

		R=[2;2;2]					// blocks!!!
		info=MPI_Reduce_scatter(B,R,[1 2 3],'PROD',NEWORLD)
		R

				R=[3 3;3 3]			// unblocks
				MPI_Reduce_scatter(C,R,[1 1 4],'PROD',NEWORLD)
				R

		info						// info==0 !!!
	
//-------------------------------------------------------
// Parent had	[0 1 1]	sends	[0  1  3	gets	[0 0 1]	the two reqstd
//				 1  8  1]			from self
// child1 had	[2 2 2]'sends	[1  4  1	gets	[24 24]'
//				 0  3  2]
// child2 had	[3 3	sends	[1  6  6	gets	[18 3
//		 3 3]		 1  1  9]		 18 3]
//			prod is	[0 24 18
//				 0 24 18]
// try [1 1 3] as counts in 1st child
//-------------------------------------------------------
// Done
//-------------------------------------------------------

MPI_Finalize						// cycle-Init or
quit							// skip to next session
		MPI_Finalize
		quit
				MPI_Finalize
				quit

//-------------------------------------------------------
// Init cycle ... skip it if you didn't Finalize just above
//-------------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(2);

[info rank]=MPI_Comm_rank(NEWORLD)
!hostname

		////////// ON FIRST CHILD //////////////
		[info rank]=MPI_Comm_rank(NEWORLD)
		!hostname

				////////// ON SECOND CHILD //////////////
				[info rank]=MPI_Comm_rank(NEWORLD)
				!hostname

//-------------------------------------------------------
// Init cycle done... jump here if you didn't Finalize just above
//-------------------------------------------------------

help MPI_Reduce

//-------------------------------------------------------
// suitable Types for each Operator
//-------------------------------------------------------
help MPI_REPLACE				// replaces all types

////////// ON PARENT Nsp //////////////
A=[10 11 12], R=[00 00 00]
info=MPI_Reduce(A,R,'REPLACE',0,NEWORLD), R

		////////// ON FIRST CHILD //////////////
		B=[20 21 22]
		info=MPI_Reduce(B,[],'REPLACE',0,NEWORLD)

				////////// ON SECOND CHILD //////////////
				C=[30 31 32]
				MPI_Reduce(C,[],'REPLACE',0,NEWORLD)

//-------------------------------------------------------
// Parent had	[10 11 12]	gets	[10 11 12]
// child1	irrelevant
// child2	irrelevant
//-------------------------------------------------------
// Same for doubles/singles, ints (includes chars), bytes (includes logicals)
//-------------------------------------------------------
A=single([10 11 12]), R=single([00 00 00])
info=MPI_Reduce(A,R,'REPLACE',0,NEWORLD), R
		B=single([20 21 22])
		info=MPI_Reduce(B,[],'REPLACE',0,NEWORLD)
				C=single([30 31 32])
				MPI_Reduce(C,[],'REPLACE',0,NEWORLD)

A=int64([10 11 12]), R=int64([00 00 00])
info=MPI_Reduce(A,R,'REPLACE',0,NEWORLD), R
		B=int64([20 21 22])
		info=MPI_Reduce(B,[],'REPLACE',0,NEWORLD)
				C=int64([30 31 32])
				MPI_Reduce(C,[],'REPLACE',0,NEWORLD)

A=uint32([10 11 12]), R=uint32([00 00 00])
info=MPI_Reduce(A,R,'REPLACE',0,NEWORLD), R
		B=uint32([20 21 22])
		info=MPI_Reduce(B,[],'REPLACE',0,NEWORLD)
				C=uint32([30 31 32])
				MPI_Reduce(C,[],'REPLACE',0,NEWORLD)

A=int16([10 11 12]), R=int16([00 00 00])
info=MPI_Reduce(A,R,'REPLACE',0,NEWORLD), R
		B=int16([20 21 22])
		info=MPI_Reduce(B,[],'REPLACE',0,NEWORLD)
				C=int16([30 31 32])
				MPI_Reduce(C,[],'REPLACE',0,NEWORLD)

A=uint8([10 11 12]), R=uint8([00 00 00])
info=MPI_Reduce(A,R,'REPLACE',0,NEWORLD), R
		B=uint8([20 21 22])
		info=MPI_Reduce(B,[],'REPLACE',0,NEWORLD)
				C=uint8([30 31 32])
				MPI_Reduce(C,[],'REPLACE',0,NEWORLD)

A=logical([00 11 12]), R=logical([00 00 00])
info=MPI_Reduce(A,R,'REPLACE',0,NEWORLD), R
		B=logical([20 00 22])
		info=MPI_Reduce(B,[],'REPLACE',0,NEWORLD)
				C=logical([30 31 00])
				MPI_Reduce(C,[],'REPLACE',0,NEWORLD)

//-------------------------------------------------------
// arithmetic operators MAX/MIN/SUM/PROD, not for bytes
//-------------------------------------------------------
help MPI_SUM

A=[10 11 12], R=[00 00 00]			// SUM is [60 63 66]
info=MPI_Reduce(A,R,'SUM',0,NEWORLD), R
		B=[20 21 22]
		info=MPI_Reduce(B,[],'SUM',0,NEWORLD)
				C=[30 31 32]
				MPI_Reduce(C,[],'SUM',0,NEWORLD)

A=single([10 11 12]), R=single([00 00 00])
info=MPI_Reduce(A,R,'SUM',0,NEWORLD), R
		B=single([20 21 22])
		info=MPI_Reduce(B,[],'SUM',0,NEWORLD)
				C=single([30 31 32])
				MPI_Reduce(C,[],'SUM',0,NEWORLD)

A=int64([10 11 12]), R=int64([00 00 00])
info=MPI_Reduce(A,R,'SUM',0,NEWORLD), R
		B=int64([20 21 22])
		info=MPI_Reduce(B,[],'SUM',0,NEWORLD)
				C=int64([30 31 32])
				MPI_Reduce(C,[],'SUM',0,NEWORLD)

A=uint32([10 11 12]), R=uint32([00 00 00])
info=MPI_Reduce(A,R,'SUM',0,NEWORLD), R
		B=uint32([20 21 22])
		info=MPI_Reduce(B,[],'SUM',0,NEWORLD)
				C=uint32([30 31 32])
				MPI_Reduce(C,[],'SUM',0,NEWORLD)

A=int16([10 11 12]), R=int16([00 00 00])
info=MPI_Reduce(A,R,'SUM',0,NEWORLD), R
		B=int16([20 21 22])
		info=MPI_Reduce(B,[],'SUM',0,NEWORLD)
				C=int16([30 31 32])
				MPI_Reduce(C,[],'SUM',0,NEWORLD)

// int8, uint8 and logical are all bytes, we'll abort on them later
// let' try now chars

A='zbc', R='   '				// MAX will be zyx
info=MPI_Reduce(A,R,'MAX',0,NEWORLD), R		// SUM would be non-didactical
		B='dyf'
		info=MPI_Reduce(B,[],'MAX',0,NEWORLD)
				C='ghx'
				MPI_Reduce(C,[],'MAX',0,NEWORLD)

// Panic! Panic! not for Bytes!!!
help MPI_Reduce				// read carefully the last 2 paragraphs
MPI_Errhandler_set('WORLD','RETURN')	// on errors. See the result:

A=uint8([10 11 12]), R=uint8([00 00 00])	// wrong SUM [30 31 32]
info=MPI_Reduce(A,R,'SUM',0,NEWORLD), R
		B=uint8([20 21 22])
		info=MPI_Reduce(B,[],'SUM',0,NEWORLD)
				C=uint8([30 31 32])
				MPI_Reduce(C,[],'SUM',0,NEWORLD)

MPI_Errhandler_set('WORLD','FATAL')	// if default errhandler, abort:

A=uint8([10 11 12]), R=uint8([00 00 00])
info=MPI_Reduce(A,R,'SUM',0,NEWORLD), R
		B=uint8([20 21 22])
		info=MPI_Reduce(B,[],'SUM',0,NEWORLD)
				C=uint8([30 31 32])
				MPI_Reduce(C,[],'SUM',0,NEWORLD)

//-------------------------------------------------------
// We get:
// Rank (0, MPI_COMM_WORLD): Call stack within LAM:
// Rank (0, MPI_COMM_WORLD):  - MPI_Reduce()
// Rank (0, MPI_COMM_WORLD):  - main()
//-------------------------------------------------------
reset
lamclean
matlab
//-------------------------------------------------------
// Ok, restart again, go on with byte-wise operators
//-------------------------------------------------------

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(2);


[info rank]=MPI_Comm_rank(NEWORLD)

		////////// ON FIRST CHILD //////////////
		[info rank]=MPI_Comm_rank(NEWORLD)

				////////// ON SECOND CHILD //////////////
				[info rank]=MPI_Comm_rank(NEWORLD)

//-------------------------------------------------------
// bitewise operators BAND/BOR/BXOR, not for floats
//-------------------------------------------------------
help MPI_BXOR

A=int64([01 02 04]), R=int64([00 00 00])
info=MPI_Reduce(A,R,'XOR',0,NEWORLD), R
		B=int64([08 16 32])
		info=MPI_Reduce(B,[],'XOR',0,NEWORLD)
				C=int64([64 128 32])
				MPI_Reduce(C,[],'XOR',0,NEWORLD)

//-------------------------------------------------------
// Parent had	[0000 0001	0000 0010	0000 0100]
// child1 had	[0000 1000	0001 0000	0010 0000]
// child2 had	[0100 0000	1000 0000	0010 0000]
// XOR    is	[0100 1001	1001 0010	0000 0100]
//		 64+9=73	128+18=146		4
//-------------------------------------------------------

A=uint32([01 02 04]), R=uint32([00 00 00])
info=MPI_Reduce(A,R,'XOR',0,NEWORLD), R
		B=uint32([08 16 32])
		info=MPI_Reduce(B,[],'XOR',0,NEWORLD)
				C=uint32([64 128 32])
				MPI_Reduce(C,[],'XOR',0,NEWORLD)

A=int16([01 02 04]), R=int16([00 00 00])
info=MPI_Reduce(A,R,'XOR',0,NEWORLD), R
		B=int16([08 16 32])
		info=MPI_Reduce(B,[],'XOR',0,NEWORLD)
				C=int16([64 128 32])
				MPI_Reduce(C,[],'XOR',0,NEWORLD)

A=uint8([01 02 04]), R=uint8([00 00 00])
info=MPI_Reduce(A,R,'XOR',0,NEWORLD), R
		B=uint8([08 16 32])
		info=MPI_Reduce(B,[],'XOR',0,NEWORLD)
				C=uint8([64 128 32])
				MPI_Reduce(C,[],'XOR',0,NEWORLD)

A='abc', R='   '				// XOR accordig to ASCII codes
info=MPI_Reduce(A,R,'XOR',0,NEWORLD), R		// ok, rather unverifiable :-)
		B='def'				// I get 'bol', hope you too
		info=MPI_Reduce(B,[],'XOR',0,NEWORLD)
				C='ghi'
				MPI_Reduce(C,[],'XOR',0,NEWORLD)

//-------------------------------------------------------
// Done
//-------------------------------------------------------
MPI_Finalize
quit
		MPI_Finalize
		quit
				MPI_Finalize
				quit


////////////////////////////////////////////////////////////
// Groups, Communicators, etc //
////////////////////////////////////////////////////////////

help comms
help groups
help process

// =====================================================
// Communicators:
// MPI_Comm_spawn, _Comm_get_parent, _Intercomm_merge
// =====================================================
// 		*** WE KEEP ON USING 3 COMPUTERS :-) ***
//-------------------------------------------------------
// Communicators:
// MPI_Comm_rank, _size, _remote_size, _test_inter
// MPI_Comm_dup, _compare, _free
// MPI_COMM_WORLD, _SELF, _NULL
// MPI_IDENT, _CONGRUENT, _SIMILAR, _UNEQUAL
// =======================================================

exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
[NEWORLD,children]=NSP_spawn(2);

MPI_COMM_WORLD=mpicomm_create('WORLD');
MPI_COMM_NULL=mpicomm_create('NULL');
MPI_COMM_SELF=mpicomm_create('SELF');

		////////// ON FIRST CHILD //////////////
		size=MPI_Comm_size(NEWORLD)	// size==3
		rank=MPI_Comm_rank(NEWORLD)	// rank==1

				////////// ON SECOND CHILD //////////////
				[size]=MPI_Comm_size(NEWORLD)
				[rank]=MPI_Comm_rank(NEWORLD)

[size]=MPI_Comm_size(MPI_COMM_SELF)
[size]=MPI_Comm_size(MPI_COMM_WORLD)
[size]=MPI_Comm_size(NEWORLD)			// size==3
[rank]=MPI_Comm_rank(NEWORLD)			// rank==0

//-------------------------------------------------------
// Can check for intercommunicators and their remote size
//-------------------------------------------------------

[flag] = MPI_Comm_test_inter (MPI_COMM_WORLD)		// flag==0
[flag] = MPI_Comm_test_inter (MPI_COMM_SELF)		// flag==0
[flag] = MPI_Comm_test_inter (NEWORLD)		// flag==0
[flag] = MPI_Comm_test_inter (children)		// flag==1
[rsiz] = MPI_Comm_remote_size(children)		// rem.size==2 children
[lsiz] = MPI_Comm_size       (children)		// loc.size==1 parent

[size] = MPI_Comm_remote_size (NEWORLD)		// Error, info==5



		////////// ON FIRST CHILD //////////////
		MPI_COMM_WORLD=mpicomm_create('WORLD');
		[flag] = MPI_Comm_test_inter (NEWORLD)	// 0
		[flag] = MPI_Comm_test_inter (MPI_COMM_WORLD)	// 0
		[size] = MPI_Comm_size       (MPI_COMM_WORLD)	// size==2
		[flag] = MPI_Comm_test_inter (parent)	// flag==1
		[rsiz] = MPI_Comm_remote_size(parent)	// rsiz==1 parnt
		[lsiz] = MPI_Comm_size       (parent)	// lsiz==2 chldr

//-------------------------------------------------------
// Can duplicate, compare and free communicators
//-------------------------------------------------------

[WORLDNEW] = MPI_Comm_dup (NEWORLD)			// blocks

		////////// ON FIRST CHILD //////////////
		[WORLDNEW] = MPI_Comm_dup (NEWORLD)	// blocks

			////////// ON SECOND CHILD //////////////
			[WORLDNEW] = MPI_Comm_dup (NEWORLD)// unblocks
			[flag]=MPI_Comm_test_inter(WORLDNEW) // 0
			MPI_Comm_free (WORLDNEW)

		[flag]=MPI_Comm_test_inter(WORLDNEW)	// flag==0
		MPI_Comm_free (WORLDNEW)

[flag]=MPI_Comm_test_inter(WORLDNEW)			// flag==0
								// keep it

[child2] = MPI_Comm_dup (children)				// blocks
		[par2] = MPI_Comm_dup (parent)		// blocks
			[par2] = MPI_Comm_dup (parent)	// unblocks
			[flag]=MPI_Comm_test_inter(par2)	// flag==1
			MPI_Comm_free (par2)
		[flag]=MPI_Comm_test_inter(par2)		// flag==1
		MPI_Comm_free (par2)
[flag]=MPI_Comm_test_inter(child2)				// flag==1
MPI_Comm_free (child2)

WORLD3 = WORLDNEW					// copy by ref.
// WORLD3 and WORLDNEW points to same object 

ZZZZZ 

[result] = MPI_Comm_compare (WORLD3, WORLDNEW)	// ident
result==MPI_IDENT
[result] = MPI_Comm_compare (WORLD3, NEWORLD)	// congruent
result==MPI_CONGRUENT					// obtained with _dup
[result] = MPI_Comm_compare (WORLD3, MPI_COMM_WORLD)	// unequal
result==MPI_UNEQUAL

help MPI_Comm_compare
help MPI_SIMILAR

//-------------------------------------------------------
// We will check for similar when we build similar groups (changing order)
//-------------------------------------------------------
// Can't Init-cycle right now, sorry, we must keep both communicators
// for later use, remember they were obtained using
// [ WORLDNEW]=MPI_Comm_dup(NEWORLD)
// WORLD3=WORLDNEW
//-------------------------------------------------------

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Groups:
// MPI_Comm_group, MPI_Group_size, _rank,
// MPI_Group_compare, _free
// MPI_GROUP_NULL, _EMPTY
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
help groups

//-------------------------------------------------------
// Groups: can extract them from communicators, check for rank & size,
//-------------------------------------------------------
help MPI_Comm_group
help MPI_Group_size
help MPI_Group_rank

[size]= MPI_Comm_size (NEWORLD)				// 3
[rank]= MPI_Comm_rank (NEWORLD)				// 0
[grp ]= MPI_Comm_group(NEWORLD)
[size]= MPI_Group_size(grp)				// 3
[rank]= MPI_Group_rank(grp)				// 0

		[size]= MPI_Comm_size (NEWORLD)		// 3
		[rank]= MPI_Comm_rank (NEWORLD)		// 1
		[grp ]= MPI_Comm_group(NEWORLD)
		[size]= MPI_Group_size(grp)		// 3
		[rank]= MPI_Group_rank(grp)		// 1

				[size]= MPI_Comm_size (NEWORLD)
				[rank]= MPI_Comm_rank (NEWORLD)
				[grp ]= MPI_Comm_group(NEWORLD)
				[size]= MPI_Group_size(grp)	// 3
				[rank]= MPI_Group_rank(grp)	// 2

//-------------------------------------------------------
// Groups: can compare and free them
//-------------------------------------------------------
help MPI_Group_compare
help MPI_Group_free

[gr2 ]=MPI_Comm_group(WORLDNEW)	// congruent to NEWORLD
[size]=MPI_Group_size(gr2)		// size==3
[rslt]=MPI_Group_compare(grp,gr2)	// result==1
rslt==MPI_IDENT				// ident. groups, same order & members

gr3 = gr2				// new identical group
dump(gr3)				// linked

[rsl1]=MPI_Group_compare(gr2, gr3)	// result1==ident
[rsl2]=MPI_Group_compare(grp, gr3)	// result2==ident
isequal(MPI_IDENT,rsl1,rsl2)

 info= MPI_Group_free (gr2)		// info==MPI_SUCCESS
//info= MPI_Group_free (gr3)		// eek!!! same address!!!
 help  MPI_Group_free			// on output, MPI_GROUP_NULL
strcmp(gr2,MPI_GROUP_NULL)		// indeed it is
[adr1]=MPI_Address(gr3)
[adr2]=MPI_Address(gr2)
adr1==adr2				// same address
strcmp(gr3,MPI_GROUP_NULL)		// for certain
help MPI_GROUP_NULL
help MPI_GROUP_EMPTY

MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')	// Errhandler on NEWORLD not enough
info=MPI_Group_free(gr3)		// it was a copy, already freed
[nfo msg]=MPI_Error_string(info)	// invalid group
info==MPI_ERR_GROUP

= MPI_Comm_free (WORLD3)
= MPI_Comm_free (WORLDNEW)		// copy already freed
[nfo msg]=MPI_Error_string(info)	// invalid communicator
info==MPI_ERR_COMM

//-------------------------------------------------------
// Too long for a single session? wanted to do this for a long :-)
//-------------------------------------------------------

MPI_Errhandler_set(MPI_COMM_WORLD,'FATAL')	// original Errhandler on 'WORLD'
info=MPI_Group_free(MPI_GROUP_EMPTY)	// any of these will do :-)
info=MPI_Group_free(gr3)		// copy already freed!!!
info=MPI_Group_free(MPI_GROUP_NULL)
= MPI_Comm_free (WORLDNEW)		// copy already freed

//-------------------------------------------------------
// Rank (0, MPI_COMM_WORLD): Call stack within LAM:
// Rank (0, MPI_COMM_WORLD):  - MPI_Group/Comm_free()
// Rank (0, MPI_COMM_WORLD):  - main()
// now you must reset & lamclean
//-------------------------------------------------------
reset					# now keyboard echoes again
lamclean				# children die
matlab					# start over

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Groups:
// MPI_Group_incl, _excl, _range_incl, _range_excl
// MPI_Comm_remote_group, _Comm_create
// MPI_Group_translate_ranks
// MPI_Group_intersection, _union, _difference
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////


exec src7.x/loader.sce 
exec utils/NSP_spawn.sci 
NEWORLD=NSP_spawn(2);


		////////// ON FIRST CHILD //////////////
		!hostname
		[rank]=MPI_Comm_rank(NEWORLD)	// rank==1

				////////// ON SECOND CHILD //////////////
				!hostname
				[rank]=MPI_Comm_rank(NEWORLD)

help groups

//-------------------------------------------------------
// Groups: can operate with them, including & excluding ranks
//-------------------------------------------------------
help MPI_Group_incl
help MPI_Group_excl

[grp ]=MPI_Comm_group (NEWORLD)		// grp ==[par ch1 ch2]
[size]=MPI_Group_size (grp )		// size==3
[grp1]=MPI_Group_excl (grp, [0])		// grp1==    [ch1 ch2]
[size]=MPI_Group_size (grp1)		// size==2
[grp2]=MPI_Group_incl (grp, [2 1])		// grp2==    [ch2 ch1]
[size]=MPI_Group_size (grp2)
[grp3]=MPI_Group_incl (grp, [1 2])		// grp3==    [ch1 ch2]
[size]=MPI_Group_size (grp3)
[grp0]=MPI_Comm_group (MPI_COMM_WORLD)		// grp0==[par]
[size]=MPI_Group_size (grp0)		// size==1

//-------------------------------------------------------
// Groups: can extract groups from remote part of intercomm
//-------------------------------------------------------
help MPI_Comm_remote_group

[grpc]=MPI_Comm_group   (children)		// children: local  [parent]
[size]=MPI_Group_size   (grpc)		//         : remote [ch1 ch2]
[size]=MPI_Comm_size    (children)		// returns local group size 1
[rslt]=MPI_Group_compare(grp0, grpc)	// grp0 we got it from WORLD
      rslt==MPI_IDENT				// children's local side is par
      =MPI_Group_free   (grpc)		// done
strcmp(grpc,MPI_GROUP_NULL)			// grpc NULL

[grpc]=MPI_Comm_remote_group(children)	// remote side [ch1 ch2]
[size]=MPI_Group_size       (grpc)		// size==2
[size]=MPI_Comm_remote_size (children)	// ok, understood
[rslt]=MPI_Group_compare  (grp3, grpc)	// grp3 we built it from [1 2]
      rslt==MPI_IDENT
[rslt]=MPI_Group_compare  (grp2, grpc)	// grp3 we built it from [2 1]
      rslt==MPI_SIMILAR
[rslt]=MPI_Group_compare  (grp1, grpc)	// grp1 from [0 1 2] removing 0
      rslt==MPI_IDENT
[rslt]=MPI_Group_compare  (grp0, grpc)	// grp1 from WORLD, it's [0]
      rslt==MPI_UNEQUAL
      =MPI_Group_free           (grpc)	// done
     strcmp(MPI_GROUP_NULL,	      grpc)	// grpc NULL

      =MPI_Comm_free	 (children)	// done with intercomm as well
     strcmp(MPI_GROUP_NULL,	  children)	// children NULL

//-------------------------------------------------------
// Groups: can translate ranks in different groups
//-------------------------------------------------------
help MPI_Group_translate_ranks

MPI_IDENT, MPI_SIMILAR, MPI_UNEQUAL, MPI_UNDEFINED
[result] = MPI_Group_compare	  (grp,		 grp1)	// unequal
[ranks ] = MPI_Group_translate_ranks (grp, [0 1 2], grp1)	// missing 0
      ranks == [MPI_UNDEFINED, 0, 1]
[result] = MPI_Group_compare	  (grp1,        grp2)	// similar
[ranks ] = MPI_Group_translate_ranks (grp1, [0 1], grp2)	// reversed
      ranks == [1, 0]
[result] = MPI_Group_compare	  (grp1,        grp3)	// ident
[ranks ] = MPI_Group_translate_ranks (grp1, [0 1], grp3)	// same order
      ranks == [0, 1]
[result] = MPI_Group_compare	  (grp2,        grp3)	// similar
[ranks ] = MPI_Group_translate_ranks (grp2, [0 1], grp3)	// reversed
      ranks == [1, 0]

//-------------------------------------------------------
// Groups: can compute intersection/union/difference among groups
//-------------------------------------------------------
help MPI_Group_intersection
help MPI_Group_union
help MPI_Group_difference

[grpi]=MPI_Group_intersection(grp, grp0)		// grp =[par ch1 ch2]
[size]=MPI_Group_size(grpi)			// grp0=[par] = grpi
[rnks]=MPI_Group_translate_ranks(grp,[0 1 2],grpi)	// only rank 0
      rnks==[0 MPI_UNDEFINED MPI_UNDEFINED]
[rslt]=MPI_Group_compare(grp0, grpi)		// ident
      rslt==MPI_IDENT
      =MPI_Group_free(grpi)			// done
     strcmp(MPI_GROUP_NULL,grpi)

[grpd]=MPI_Group_difference(grp, grp0)		// grpd=[ch1 ch2]
[size]=MPI_Group_size(grpd)			// size==2
[rnks]=MPI_Group_translate_ranks(grp,[0 1 2],grpd)	// missing 0
      rnks==[MPI_UNDEFINED 0 1]
[rslt]=MPI_Group_compare(grp1,grpd)		// ident
      =MPI_Group_free        (grpd)

 = MPI_Group_free (grp2)
 = MPI_Group_free (grp3)

[grps]=MPI_Group_union    (grp1, grp0)		// grp0=[par]
[size]=MPI_Group_size           (grps)		// grp1=[ch1 ch2]
[rslt]=MPI_Group_compare   (grp, grps)		// similar
      rslt==MPI_SIMILAR
[rnks]=MPI_Group_translate_ranks(grps,[0 1 2],grp)	// 0 is last
      rnks==[1 2 0]

//-------------------------------------------------------
// Groups: can operate with them, including & excluding rank ranges
//-------------------------------------------------------

help MPI_Group_range_incl
help MPI_Group_range_excl

      =MPI_Group_free(grp0)		// we're creating a new grp0/1
      =MPI_Group_free(grp1)		// free former ones/no mem leak
MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')
      =MPI_Group_free(grp1)		// can't do it twice
[info2 msg]=MPI_Error_string(info)		// invalid group
     ==MPI_ERR_GROUP

[grp0]=MPI_Group_range_excl     (grp,[1 2 1]'    )	// grp=[par ch1 ch2]
[size]=MPI_Group_size			 (grp0)	// exclude 1 to 2 step 1
[rnks]=MPI_Group_translate_ranks(grp,[0 1 2],grp0)	// only parent left
      rnks==[0 MPI_UNDEFINED MPI_UNDEFINED]

[grp1]=MPI_Group_range_incl     (grp,[1;2;1]     )	// include 1 to 2 step 1
[size]=MPI_Group_size (grp1)			// size==2
[rnks]=MPI_Group_translate_ranks(grp,[0 1 2],grp1)	// [ch1 ch2] left
      rnks==[MPI_UNDEFINED 0 1]
[rnks]=MPI_Group_translate_ranks(grp1,[0 1 2],grp)	// [ch1 ch2] left
      rnks==[1 2 0]					// there is no rank 2
[info2 msg]=MPI_Error_string(info)			// invalid rank in grp1
     ==MPI_ERR_RANK

[grps]=MPI_Group_union(grp1, grp0)		// notice order: [ch1 ch2 par]
[size]=MPI_Group_size (grps)			// size==3
[rnks]=MPI_Group_translate_ranks(grps,[0 1 2],grp)	// ranks 0 1 2 in grps
      rnks==[1 2 0]				// are the ranks 1 2 0 in grp

= MPI_Group_free (grp0)			// done
= MPI_Group_free (grp1)			// keep grps to create comm

//-------------------------------------------------------
// Communicators: can create them. Recall we said
// "We will check for similar when we build similar groups (changing order)"
//-------------------------------------------------------
help MPI_Comm_create

[comms] = MPI_Comm_create (NEWORLD, grps)	// blocks
						// ouch! no grps in children!!!

		////////// ON FIRST CHILD //////////////
		[grp0]=MPI_Comm_group (parent)	// grp0=[par]
		[size]=MPI_Group_size (grp0)	// local size 2
		      =MPI_Group_free (grp0)	// not interested
		      =MPI_Comm_free (parent)	// anymore on that side
		 strcmp    (MPI_COMM_NULL, parent)
		      =MPI_Comm_free (parent)	// can't do it twice
		[info2 msg]=MPI_Error_string(info)	// invalid communicator
		     ==MPI_ERR_COMM

		[grp ]=MPI_Comm_group (NEWORLD)	// grp=[par ch1 ch2]
		[grps]=MPI_Group_range_incl(grp,...
				[1 2 1 ;...		// include 1 to 2 step 1
				 0 0 1]')		// then    0 to 0 step 1
		[rnks]=MPI_Group_translate_ranks(grps,[0 1 2],grp)
		      rnks==    [1 2 0]			// grps==[ch1 ch2 par]

		[comms] = MPI_Comm_create (NEWORLD, grps)

				////////// ON SECOND CHILD //////////////
				[grp ]=MPI_Comm_group(NEWORLD)
				[grps]=MPI_Group_range_incl(grp,...
					[1 0;2 0;1 1])
				[comms]=MPI_Comm_create(NEWORLD, grps)
							////////////////////////////
							// unblocks all
							////////////////////////////
[rslt]=MPI_Comm_compare(NEWORLD, comms)		// similar
      rslt==MPI_SIMILAR

//-------------------------------------------------------
// cleaning
//-------------------------------------------------------
= MPI_Group_free (grps)
= MPI_Comm_free (comms)
MPI_Finalize
quit

		////////// ON FIRST CHILD //////////////
		= MPI_Group_free (grp )
		= MPI_Group_free (grps)
		= MPI_Comm_free (comms)
		MPI_Finalize
		quit

				////////// ON SECOND CHILD //////////////
				= MPI_Group_free (grp )
				= MPI_Group_free (grps)
				= MPI_Comm_free (comms)
				MPI_Finalize
				quit


// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Communicators:
// MPI_Intercomm_create, MPI_Comm_split
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////

putenv(['LAM_MPI_SSI_rpi=tcp']), MPI_Init

args={}					// choose one of these
args={'-display',getenv('DISPLAY')}	// depending on your rsh/ssh environ

args={args{:},...
	'-e','matlab','-nosplash','-nojvm',...
	'-r','startup_mergeParent'}

[children errs] = MPI_Comm_spawn ('/usr/X11R6/bin/xterm',args,...
				2,'NULL',0,'SELF')
[NEWORLD] = MPI_Intercomm_merge (children, 0)
MPI_Errhandler_set(NEWORLD,'RETURN')


!hostname

		////////// ON FIRST CHILD //////////////
		!hostname
		[size]=MPI_Comm_size(NEWORLD)	// size==3
		[rank]=MPI_Comm_rank(NEWORLD)	// rank==1

				////////// ON SECOND CHILD //////////////
				!hostname
				[size]=MPI_Comm_size(NEWORLD)
				[rank]=MPI_Comm_rank(NEWORLD)

help comms			// we're almost done
help MPI_Intercomm_create
help MPI_Comm_split

//-------------------------------------------------------
// Groups: one more experiment on creating Comms
//-------------------------------------------------------
help MPI_Comm_create
help MPI_Comm_group
help MPI_Group_incl

[grp ]=MPI_Comm_group(NEWORLD)			// [par ch1 ch2]
[grp1]=MPI_Group_incl(grp,[1])			// [ch1]
[com1]=MPI_Comm_create(NEWORLD, grp1)		// blocks

		[grp ]=MPI_Comm_group(NEWORLD)
		[grp1]=MPI_Group_incl(grp,[1])		// [ch1]
		[com1]=MPI_Comm_create(NEWORLD, grp1)	// blcks

				[grp ]=MPI_Comm_group(NEWORLD)
				[grp1]=MPI_Group_incl(grp,[1]) // unbl
				[com1]=MPI_Comm_create(NEWORLD, grp1)

// naturally, blocks in NEWORLD, but only meaningful in rank 1

MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')			// avoid abort
[size]=MPI_Comm_size(com1)				// we're not in com1
[inf2 msg ]=MPI_Error_string(info)			// invalid comm
     ==MPI_ERR_COMM

		[size]=MPI_Comm_size(com1)		// no problem, size==1

				[size]=MPI_Comm_size(com1)
				 info==MPI_ERR_COMM	// didn't need Errhndlr

//nfo = MPI_Group_free(grp1)				// keep grp1/later use
= MPI_Comm_free (com1)				// invalid
==MPI_ERR_COMM

//		= MPI_Group_free(grp1)		// keep it
		= MPI_Comm_free (com1)
		==MPI_SUCCESS

				= MPI_Group_free(grp1)
				==MPI_SUCCESS
				= MPI_Comm_free (com1)
			 	==MPI_ERR_COMM	// not in it
				= MPI_Group_free(grp1)
			 	==MPI_ERR_GROUP	// already freed
				= MPI_Group_size(grp1)
			 	==MPI_ERR_GROUP

		[size] = MPI_Group_size(grp1)	// still valid here
		       = MPI_Group_free(grp1)	// done

[size] = MPI_Group_size(grp1)			// funny, size==1
       = MPI_Group_free(grp1)			// done as well
[size] = MPI_Comm_size (com1)			// not in it
       ==MPI_ERR_COMM

//-------------------------------------------------------
// Can create intercomms with local/remote parts
//-------------------------------------------------------
help MPI_Intercomm_create

[icom01]=MPI_Intercomm_create('SELF',0, ...	// local part: parent
				 children,0, 666)	// remote: children
							// remote peer: ch1
		[icom01]=...
		MPI_Intercomm_create(MPI_COMM_WORLD,0, ...	// local: children
				      parent,0, 666)	// remote: parent

				////////// ON SECOND CHILD //////////////
				[icom01]=...
				MPI_Intercomm_create('WORLD',0, ...
						      parent,0, 666)

[flag] = MPI_Comm_test_inter  (icom01)		// flag==1
[size] = MPI_Comm_remote_size (icom01)		// remote size==2 chld
[grpr] = MPI_Comm_remote_group(icom01)		// remote group=children
[size] = MPI_Group_size(grpr)			// size==2

[size] = MPI_Group_size(grp )			// size==3 from NEWORLD
[grpc] = MPI_Group_excl(grp,0)			// grpc=[ch1 ch2]
[size] = MPI_Group_size(grpc)			// size==2 now
[rslt] = MPI_Group_compare(grpr, grpc)		// identical
      rslt == MPI_IDENT

 = MPI_Group_free (grpc)				// done
 = MPI_Group_free (grpr)

[size] = MPI_Comm_size (icom01)			// local size==1 parent
[grpl] = MPI_Comm_group(icom01)			// local group=myself
[size] = MPI_Group_size(grpl)			// size==2

[grps] = MPI_Comm_group(MPI_COMM_SELF)		// myself
[size] = MPI_Group_size(grps)			// size==1
[rslt] = MPI_Group_compare(grpl, grps)		// identical
      rslt == MPI_IDENT

 = MPI_Group_free (grps)				// done
 = MPI_Group_free (grpl)

		[flag]=MPI_Comm_test_inter(icom01)	// flag==1
		[size]=MPI_Comm_size (icom01)	// local group 2 chldrn
		[grpl]=MPI_Comm_group(icom01)
		[size]=MPI_Group_size(grpl)

		[grpw]=MPI_Comm_group(MPI_COMM_WORLD)
		[rslt]=MPI_Group_compare(grpl,grpw)// identical
		      rslt == MPI_IDENT

		= MPI_Group_free(grpw)	// done
		= MPI_Group_free(grpl)

		= MPI_Comm_free(icom01)	// done

				////////// ON SECOND CHILD //////////////
				= MPI_Comm_free(icom01)

= MPI_Comm_free(icom01)
strcmp(MPI_COMM_NULL,icom01)

//-------------------------------------------------------
// Can split communicators in parts
//-------------------------------------------------------

help MPI_Comm_split

[black]=MPI_Comm_split(NEWORLD, 0, 69)	// color 0 ("black") order 69
								// blocks
		[white]=MPI_Comm_split(NEWORLD, 1, 69)	// blocks

				////////// ON SECOND CHILD //////////////	// unblocks
				[black]=MPI_Comm_split(NEWORLD, 0, 66)
						// color 0 order 66 before 69

[size]=MPI_Comm_size (black)		// size==2
[grpb]=MPI_Comm_group(black)		// [child2 parent]
[rnks]=MPI_Group_translate_ranks(grp,[0 1 2],grpb) // new ranks are:
 rnks == [1 MPI_UNDEFINED 0]
  =    MPI_Group_free(grpb)
  =    MPI_Comm_free (black)

		[size]=MPI_Comm_size (white)	// size==1
		[grpw]=MPI_Comm_group(white)	// group=[ch1]
		[rnks]=MPI_Group_translate_ranks(grp,[0 1 2],grpw)
		     rnks==[MPI_UNDEFINED 0 MPI_UNDEFINED]
		     = MPI_Group_free(grpw)
		     = MPI_Comm_free (white)

				////////// ON SECOND CHILD //////////////
				     = MPI_Comm_free (black)

//-------------------------------------------------------
// cleaning
//-------------------------------------------------------
= MPI_Group_free(grp)
= MPI_Comm_free (NEWORLD)
MPI_Finalize
quit

		= MPI_Group_free (grp)
		= MPI_Comm_free (NEWORLD)
		MPI_Finalize
		quit

				= MPI_Group_free (grp)
				= MPI_Comm_free (NEWORLD)
				MPI_Finalize
				quit


// =============================================================
// The MPI_Info object (MPI-2.0)
// MPI_Info_create, _dup, _free
// MPI_Info_set, _get, _get_valuelen, _get_nkeys, _get_nthkey
// MPI_INFO_NULL
// no children this time
// ==============================================================

exec src7.x/loader.sce
MPI_Init();

help info

//-------------------------------------------------------
// Info: create an object, add entries, lookup entries, query length
//-------------------------------------------------------

help MPI_Info_create
help MPI_Info_set
help MPI_Info_get
help MPI_Info_get_valuelen

INFO_NULL = mpiinfo_create('NULL');
inf = mpiinfo_create(); // nsp object creation 

MPI_Info_set (inf, 'Key1', 'first element''s value')
MPI_Info_set (inf, 'Key2', 'second element''s value')

[val] = MPI_Info_get (inf, 'Key3')		// flag==0
[val] = MPI_Info_get (inf, 'Key2')		// flag==1

[len] = MPI_Info_get_valuelen (inf, 'Key2')	// len==22

[val] = MPI_Info_get (inf, 'Key1')		// flag==1

[len] = MPI_Info_get_valuelen (inf, 'Key1')	// len==21
[len] = MPI_Info_get_valuelen (inf, 'Key0')	// flag==0

//-------------------------------------------------------
// Info: query no. entries, get indexed entry
//-------------------------------------------------------
help MPI_Info_get_nkeys
help MPI_Info_get_nthkey

[nkeys] = MPI_Info_get_nkeys  (inf)		// nkeys==2
[key] = MPI_Info_get_nthkey (inf, 1)		// key=='Key2'
[len] = MPI_Info_get_valuelen (inf, key)	// len==22
[val] = MPI_Info_get          (inf, key)	// val=='second...'
[key ] = MPI_Info_get_nthkey (inf, 2)		// no such key

//-------------------------------------------------------
// Info: duplicate / free objects
//-------------------------------------------------------

help MPI_Info_dup
help MPI_Info_free
help MPI_INFO_NULL

[inf2] = MPI_Info_dup       (inf)
[nkeys] = MPI_Info_get_nkeys (inf2)		// same 2 keys

[key ]   =MPI_Info_get_nthkey(inf2, 1)		// key=='Key2'
[val]=MPI_Info_get       (inf2, key)		// val=='second...'

[key ]   =MPI_Info_get_nthkey(inf2, 0)		// Key1
[val]=MPI_Info_get       (inf2, key)		// first...


//-------------------------------------------------------
// Info: delete entries
//-------------------------------------------------------

help info						// we're done
help MPI_Info_delete					// this is the last one

MPI_Info_delete (inf, 'Key3')
[nkeys] = MPI_Info_get_nkeys (inf)			// 2 keys
[key ] = MPI_Info_get_nthkey(inf, 0)		// 1st one keyed 'Key1'
[val]=MPI_Info_get      (inf, key)		// BTW, value 'first...'
MPI_Info_delete    (inf, key)		// delete it
[nkeys] = MPI_Info_get_nkeys (inf)			// 1 key left

MPI_Finalize();

quit

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Process management
// MPI_Comm_spawn[_multiple], MPI_Comm[_accept,_connect,_disconnect],
// MPI_[Open_,Close_]port, MPI_[ [Un] publish_, Lookup_] name
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// funny mixing of xterms, matlabs, etc
// we require 3 hosts in $LAMBHOST, your host computer last in list
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////

help process

//-------------------------------------------------------
// can spawn program copies, as long as they MPI_Init
//-------------------------------------------------------
help MPI_Comm_spawn

exec src7.x/loader.sce;
MPI_Init()

COMM_SELF=mpicomm_create('SELF');
INFO_NULL=mpiinfo_create('NULL');
cmd = "exec(''src7.x/loader.sce'');MPI_Init();";
cmd = cmd + "MPI_Finalize();quit;";
nsp_exe = getenv('SCI')+'/bin/nsp';
args=["-name","nsp-child","-e", cmd];
[children errs] = MPI_Comm_spawn (nsp_exe,args,1,INFO_NULL,0, COMM_SELF);
// child will execute cmd 

//ZZZ   MPI_Comm_free(children)
// we can accept MPI_Comm_free(children)
// which could free children and set children to MPI_COMM_NULL 

children =[];

// cannot start more children 

[children errs] = MPI_Comm_spawn (nsp_exe,args,1,INFO_NULL,0, COMM_SELF); 
MPI_Comm_free(children)

//-------------------------------------------------------
// can spawn arbitrary processes, as long as they MPI_Init
//-------------------------------------------------------
help MPI_Comm_spawn_multiple

[children errs] = MPI_Comm_spawn_multiple (...
	{'/usr/X11R6/bin/xterm', '/usr/X11R6/bin/xterm'},...
	{{args{:},'-e','matlab','-nojvm','-nosplash'}...// args 1st xterm
	  args       },...				// args 2nd xterm
	{2,2},...					// 2 copies of each
	{MPI_INFO_NULL,MPI_INFO_NULL},...		// no additional info
	 0,'SELF')					// and i'm the parent

		////////// ON xterms  //////////////
		$MPITB_ROOT/stdlones/MPI_Init		# unblock spawn
		exit

				////////// ON matlabs  //////////////
				!hostname
				MPI_Init		// unblock spawn
				quit			// dirty quitting
MPI_Comm_free(children)

//-------------------------------------------------------
// How could we communicate 2 MPI processes not related by spawn?
// ie, not mpirun together, not spawned one from the other
//-------------------------------------------------------

		#---- open a new xterm and log into cluster node  -------

		ssh ox1		# ox1 is the name of my slave cluster node
		hostname	# should be ox1 by now
		matlab -nojvm

		//---- Nsp mode from now on ---------------------------

		MPI_Init
		[size]=MPI_Comm_size(MPI_COMM_WORLD)	// size==1
		[rank]=MPI_Comm_rank(MPI_COMM_WORLD)	// rank==0
		
//-------------------------------------------------------
// we have no inter-communicator to merge... see the point?
//-------------------------------------------------------
// can open ports, accept connections on them, connect to them
//-------------------------------------------------------

exec src7.x/loader.sce;
MPI_Init()

COMM_WORLD=mpicomm_create(MPI_COMM_WORLD);
COMM_SELF=mpicomm_create('SELF');
INFO_NULL=mpiinfo_create('NULL');

help MPI_Open_port
help MPI_Comm_accept
help MPI_Comm_connect

[port] = MPI_Open_port(INFO_NULL)	// write down port name

		////////// ON cluster node Nsp  //////////////
		
		port = 'n2:i11:323'		// or whatever you got above
		whos port			// requires knowing port name
		COMM_WORLD=mpicomm_create(MPI_COMM_WORLD);
		INFO_NULL=mpiinfo_create('NULL');
		[svcomm]=MPI_Comm_connect(port,INFO_NULL,0,COMM_WORLD)
						// blocks, of course

[clcomm]=MPI_Comm_accept (port, INFO_NULL, 0, COMM_SELF)	// unblocks
[size] = MPI_Comm_size        (clcomm)			// 1 local
[size] = MPI_Comm_remote_size (clcomm)			// 1 remote

		[size] = MPI_Comm_size        (svcomm)	// 1 local
		[size] = MPI_Comm_remote_size (svcomm)	// 1 remote

//-------------------------------------------------------
// having inter-communicator, it's the same old stuff now
//-------------------------------------------------------
[NEWORLD] = MPI_Intercomm_merge(clcomm,0)			// blocks

		[NEWORLD] = MPI_Intercomm_merge(svcomm,1)	// unblocks
		[size]    = MPI_Comm_size      (NEWORLD)	// 2 ranks

[rank] = MPI_Comm_rank(NEWORLD)				// server rnk 0
= MPI_Send (port, 1, 0, NEWORLD)				// send sthing

		whos port
		port(1:10)=deal(' ')				// room for recv
		[stat] = MPI_Recv (port, 0, 0, NEWORLD)	// recv sthing
		port						// if only bfore

		MPI_Comm_free(NEWORLD)				// done with it

MPI_Comm_free(NEWORLD)

//-------------------------------------------------------
// done, disconnecting. Can disconnect from ports, and clos'em
//-------------------------------------------------------
help MPI_Comm_disconnect
help MPI_Close_port

		////////// ON cluster node Nsp  //////////////
		= MPI_Comm_disconnect (svcomm)
		strcmp(MPI_COMM_NULL,       svcomm)

= MPI_Comm_disconnect (clcomm)
strcmp(MPI_COMM_NULL,       clcomm)

// = MPI_Close_port (port)		// not yet, keep it for next try
// port

//-------------------------------------------------------
// nice, but required knowing port name at client side
// can publish a port name and look it up from client side
// this only requires according previously on a "well-known" service name
//-------------------------------------------------------

help process					// we're almost done

help MPI_Publish_name
help MPI_Lookup_name
help MPI_Unpublish_name

INFO_NULL=mpiinfo_create('NULL');

MPI_Publish_name('ServName', INFO_NULL, port)	// same port
port

		////////// ON cluster node Nsp  //////////////
		port
		port(1:10)=deal(' ')				// clean it
		[port] = MPI_Lookup_name('ServName','NULL')// there it is

		//---------------------
		// same old stuff now
		//---------------------
		[svcomm]=MPI_Comm_connect(port,'NULL',0,MPI_COMM_WORLD)
								// blocks

[clcomm] =MPI_Comm_accept (port,'NULL',0,'SELF')		// unblocks
[NEWORLD]=MPI_Intercomm_merge(clcomm,0)			// blocks

		[NEWORLD] = MPI_Intercomm_merge(svcomm,1)	// unblocks
		[rank]    = MPI_Comm_rank      (NEWORLD)	// rank 1
		port(1:10)=deal(' ')				// room for recv
		[stat] = MPI_Recv (port, 0, 0, NEWORLD)	// from rank 0

= MPI_Send (port,1,0,NEWORLD)				// unblocks
= MPI_Comm_free     (NEWORLD)				// done

		MPI_Comm_free(NEWORLD)				// done with it
		MPI_Comm_disconnect (svcomm)
		strcmp(MPI_COMM_NULL,svcomm)
		strcmp(MPI_COMM_NULL,NEWORLD)
		port						// it worked

= MPI_Comm_disconnect (clcomm)
strcmp(MPI_COMM_NULL,       clcomm)
strcmp(MPI_COMM_NULL,      NEWORLD)

= MPI_Unpublish_name ('ServName', 'NULL', port)		// done with it
= MPI_Close_port                         (port)

MPI_Errhandler_set('SELF','RETURN')				// trying hard
MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')				// to survive

=    MPI_Unpublish_name ('ServName', 'NULL', port)		// info==25
[nfo msg]=MPI_Error_string(info)				// publish svc
help      MPI_ERRCODES
 ==  MPI_ERR_SERVICE

=    MPI_Close_port  (port)				// info==16
[nfo msg]=MPI_Error_string(info) 				// unclassified
help      MPI_ERR_OTHER

//-------------------------------------------------------
// it's funny to be naughty sometimes
//-------------------------------------------------------
= MPI_Comm_disconnect (clcomm)		// can't disconnect twice
						// no matter how hard we tried
//-------------------------------------------------------
// MPI process rank 0 (n2, p11315) caught a SIGSEGV in MPI_Comm_disconnect.
// Rank (0, MPI_COMM_WORLD): Call stack within LAM:
// Rank (0, MPI_COMM_WORLD):  - MPI_Comm_disconnect()
// Rank (0, MPI_COMM_WORLD):  - main()
//-------------------------------------------------------
		////////// ON cluster node xterm  //////////////
		[flag]=MPI_Initialized
		
reset						# in xterm
lamclean					# the other dies

		////////// ON cluster node xterm  //////////////
		reset
		lamclean

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Errhandlers & Miscellaneous
// MPI_Errhandler_ [create,free,set,get]
// MPI_Error_ [class|string]
// MPI_ERRORS_ARE_FATAL / RETURN
// MPI_ERRHANDLER_NULL / Nsp_ERRORS_RETURN
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Just 1 computer this time
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////

MPI_Init

help errors

//-------------------------------------------------------
// can get/set error handler
//-------------------------------------------------------
help MPI_Errhandler_get
help MPI_Errhandler_set

[ eh]=MPI_Errhandler_get (MPI_COMM_WORLD)
strcmp(eh, MPI_ERRORS_ARE_FATAL)

    = MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')
[eh]= MPI_Errhandler_get (MPI_COMM_WORLD)
strcmp(eh, MPI_ERRORS_RETURN)

help MPI_ERRORS_ARE_FATAL
help MPI_ERRORS_RETURN

//-------------------------------------------------------
// can obtain informational error string
//-------------------------------------------------------
help MPI_Error_string

help MPI_ERRCODES
help MPI_ERR_BASE
help MPI_ERR_LASTCODE
[msg] = MPI_Error_string(MPI_ERR_BASE)	// notice info==0
[errc msg] = MPI_Error_string(MPI_ERR_LASTCODE)	// notice errc==13
[msg] = MPI_Error_string(errc)		// LASTCODE was invalid
 errc     == MPI_ERR_ARG			// invalid arg for Error_string
[cls] = MPI_Error_class (errc)		// cls==errc
[msg] = MPI_Error_string(cls)
 cls      == MPI_ERR_ARG

//-------------------------------------------------------
// can play with Nsp error handlers
//-------------------------------------------------------
help MPI_Errhandler_create

which Nsp_ERRORS_RETURN			// suitable Nsp cmd
type  Nsp_ERRORS_RETURN			// for Errhandler_create

[EHND]=MPI_Errhandler_create('Nsp_ERRORS_RETURN')
      =MPI_Errhandler_set(MPI_COMM_WORLD,EHND)
[ eh ]=MPI_Errhandler_get(MPI_COMM_WORLD)
strcmp(eh,  EHND)
 
[errc msg] = MPI_Error_string(MPI_ERR_LASTCODE)	// causing a 'dummy' error
//-----------------------------------------------// Whoa! lots of output
// returned error:
// errclass==13, error==22, errmsg==
// in MPI_COMM_WORLD
// won't take any action on it, except for writing its MPI_Error_string:
// MPI_Error_string: invalid argument
//-------------------------------------------------------
[cls] = MPI_Error_class (errc)		// same previous error
 errc     == MPI_ERR_ARG			// invalid arg for Error_string
 cls      == MPI_ERR_ARG

//-------------------------------------------------------
// undoing Nsp errhandler
//-------------------------------------------------------
help MPI_Errhandler_free

 info= MPI_Errhandler_free(EHND)		// wrong!!! still set in WORLD!
strcmp(MPI_ERRHANDLER_NULL,EHND)

[ eh2] = MPI_Errhandler_get (MPI_COMM_WORLD)	// if we get an error right now
strcmp(eh2,eh)					// LAM will call freed EHND

help errors					// well, we have seen it all

[errc msg] = MPI_Error_string(MPI_ERR_LASTCODE)	// stubborn as only we can
//-------------------------------------------------------
// MPI process rank 0 (n2, p12215) caught a SIGSEGV in MPI_Error_string.
// Rank (0, MPI_COMM_WORLD): Call stack within LAM:
// Rank (0, MPI_COMM_WORLD):  - MPI_Error_string()
// Rank (0, MPI_COMM_WORLD):  - main()
//-------------------------------------------------------
// we should have done this to avoid that
//      = MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')
// [ eh]= MPI_Errhandler_get(MPI_COMM_WORLD)
// strcmp(eh,  MPI_ERRORS_RETURN)
// MPI_Finalize
// quit
//-------------------------------------------------------
// now we are forced to this instead
//-------------------------------------------------------
reset
lamclean


// =========================================================
// Attribute Caching
// MPI_Keyval_ [create|free|], MPI_Attr_ [put|get|delete]
// MPE_Attr_put, Nsp_COPY_FN, Nsp_DEL_FN
// =========================================================
// Again, just 1 computer, although the interesting part
// of attributes is that they're stored in the communicator
// and thus immediately accessible to all ranks in comm
// =========================================================

help environ
help caching

MPI_Init()
A={1 2;3 4};
COMM_WORLD=mpicomm_create(MPI_COMM_WORLD);

//-------------------------------------------------------
// Creating, setting and accessing attributes
//-------------------------------------------------------
help MPI_Keyval_create
help MPI_Attr_put
help MPI_Attr_get

[kv]     = MPI_Keyval_create ('NULL', 'NULL', {})	// kv==19 or so
         = MPI_Attr_put (COMM_WORLD, kv, A)		// store attr. in WORLD
[B flag] = MPI_Attr_get (COMM_WORLD, kv)		// flag 1, B==A

[C flag] = MPI_Attr_get (COMM_WORLD, 666)		// info==28
[info2   msg] = MPI_Error_string(info)			// invalid key value
         ==MPI_ERR_KEYVAL

//-------------------------------------------------------
// Overwriting, deleting attributes
//-------------------------------------------------------
help MPI_Attr_delete

	    B = [1 2;3 4]

         = MPI_Attr_put (MPI_COMM_WORLD, kv, B)		// overwrite
[C flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// flag==1, C==B

	    B = 'hello'
         = MPI_Attr_put (MPI_COMM_WORLD, kv, B)
[C flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// flag==1, C==B

         = MPI_Attr_delete(MPI_COMM_WORLD, kv)
[B flag] = MPI_Attr_get   (MPI_COMM_WORLD, kv)		// flag==0, B=={}

//-------------------------------------------------------
// The copy-function stuff
//-------------------------------------------------------
	    C = {A B;A B}
         = MPI_Attr_put (MPI_COMM_WORLD, kv, C)	// notice we put C
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)	// flag==1, B==C

[WORLD2] = MPI_Comm_dup (MPI_COMM_WORLD)		// copy_fn was NULL
[B flag] = MPI_Attr_get ( WORLD2, kv)	// attr not dup, flag==0
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)	// exists in WORLD, flag==1

//-------------------------------------------------------
// The persistent copy problem
//-------------------------------------------------------
	    C = 'hello'
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)	// it is just C, just!
       dump(C)

      clear C
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)	// accessing freed memory?
//-------------------------------------------------------
// ??? Error using ==> display			// Nsp tried to display B
// Unknown command option.
//-------------------------------------------------------
      B						// see the problem?   
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv);	// using ; Nsp won't display
	    C = []
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)	// ouch!! we shouldn't have
//-------------------------------------------------------
// MPI process rank 0 (n2, p13287) caught a SIGSEGV.
//-------------------------------------------------------
reset
lamclean
matlab

//-------------------------------------------------------
// Ok, again from the beginning, copy-function revisited
//-------------------------------------------------------

MPI_Init
MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')			// try to survive
[kv]     = MPI_Keyval_create ('NULL', 'NULL', {})	// kv==19 or so
	    A = [1 2;3 4]
	    B = 'hello'
	    C = {A B;A B}
         = MPI_Attr_put (MPI_COMM_WORLD, kv, C)	// notice we put C
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)	// flag==1, B==C
	    C = 'hello'
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)	// it's just C
       dump(C)

[WORLD2] = MPI_Comm_dup (MPI_COMM_WORLD)		// copy_fn was NULL
[B flag] = MPI_Attr_get ( WORLD2, kv)	// attr not dup, flag==0
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)	// exists in WORLD, flag==1

//-------------------------------------------------------
// starting over, this time with some copy-function
//-------------------------------------------------------
help MPI_Keyval_free
help MPI_KEYVAL_INVALID

         = MPI_Comm_free (WORLD2)		// first Comm, then Key
         = MPI_Keyval_free (kv), kv	// kv==-1, info==0, ok
 kv           ==MPI_KEYVAL_INVALID
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)	// invalid key
[info2   msg] = MPI_Error_string(info)		// invalid key value
         ==MPI_ERR_KEYVAL

//-------------------------------------------------------
// standard LAM/MPI copy function
//-------------------------------------------------------
help MPI_Keyval_create
help MPI_NULL_COPY_FN
help MPI_DUP_FN

      A       = [1 2;3 4]
[kv]     = MPI_Keyval_create ('DUP','NULL',{})	// now we dup
         = MPI_Attr_put (MPI_COMM_WORLD, kv, A)		// using MPI callback
[WORLD2] = MPI_Comm_dup (MPI_COMM_WORLD)			// kv==20 or so
[B flag] = MPI_Attr_get ( WORLD2, kv)		// now it dups, B==A
[C flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// have both

      A       = {1 2;3 4}
[B flag] = MPI_Attr_get  ( WORLD2, kv)		// it is just A
[C flag] = MPI_Attr_get  (MPI_COMM_WORLD, kv)		// in both
         = MPI_Attr_delete(WORLD2, kv)		// deleted from WORLD2
[B flag] = MPI_Attr_get  ( WORLD2, kv)		// here
[C flag] = MPI_Attr_get  (MPI_COMM_WORLD, kv)		// not here

         = MPI_Keyval_free (kv), kv		// first Key, then Com.
 kv           ==MPI_KEYVAL_INVALID			// same result
[B flag] = MPI_Attr_get ( WORLD2, kv)		// finished
         ==MPI_ERR_KEYVAL				// dirty quitting
         = MPI_Comm_free ( WORLD2)	// Keyval_free while attr. put in WORLD
         = MPI_Comm_free (MPI_COMM_WORLD)			// Ok, starting all over
//-------------------------------------------------------
// MPI process rank 0 (n2, p13392) caught a SIGSEGV in MPI_Comm_free.
// Rank (0, MPI_COMM_WORLD): Call stack within LAM:
// Rank (0, MPI_COMM_WORLD):  - MPI_Comm_free()
// Rank (0, MPI_COMM_WORLD):  - main()
//-------------------------------------------------------
reset
lamclean
matlab
//-------------------------------------------------------
// Ok, again from the beginning, persistent copy revisited
//-------------------------------------------------------

MPI_Init()
MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')			// try to survive
[kv]     = MPI_Keyval_create('Nsp_COPY_FN',...
				  'Nsp_DEL_FN',{})	// using Nsp skeletns
	A     = [1 2;3 4]
         = MPI_Attr_put (MPI_COMM_WORLD, kv, A)		// attr NOT PERSISTENT!!
[WORLD2] = MPI_Comm_dup (MPI_COMM_WORLD)			// Nspcallback called
//-------------------------------------------------------
// Copying attribute keyed 19				// we get this message
//      1     2						// from Nsp_COPY_FN
//      3     4
//  
// Nsp copy wrapper0x0: 1 pers. copies
//-------------------------------------------------------
type Nsp_COPY_FN

[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// this not persistent
[C flag] = MPI_Attr_get ( WORLD2, kv)		// this is  persistent
      A       = 'hello'
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// changes, it's A
[C flag] = MPI_Attr_get ( WORLD2, kv)		// notice persistence

clear A
[C flag] = MPI_Attr_get   ( WORLD2, kv)		// notice persistence
         = MPI_Attr_delete( WORLD2, kv)		// Nspcallback called
//-------------------------------------------------------
// Currently, callback in MPI_KEYVAL.mexglx will delete this:
//      1     2						// we get this message
//      3     4						// from Nsp_DEL_FN
//  
// Nsp dlt wrapper0x0: 0 pers. copies
//-------------------------------------------------------
         = MPI_Comm_free (WORLD2)			// exit as clean as pos

//B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// no, won't try this
         = MPI_Attr_delete(MPI_COMM_WORLD, kv)		// ouch! this also hurts
//        = MPI_Keyval_free (kv), kv		// had no time to free
//-------------------------------------------------------
// MPI process rank 0 (n2, p13604) caught a SIGSEGV in MPI_Comm_delete_attr.
// Rank (0, MPI_COMM_WORLD): Call stack within LAM:
// Rank (0, MPI_COMM_WORLD):  - MPI_Comm_delete_attr()
// Rank (0, MPI_COMM_WORLD):  - MPI_Attr_delete()
// Rank (0, MPI_COMM_WORLD):  - main()
//-------------------------------------------------------
reset
lamclean
matlab
//-------------------------------------------------------
// To avoid that, we created MPE_Attr_put, which puts attr as persistent
//-------------------------------------------------------
MPI_Init
MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')			// try to survive
[kv]     = MPI_Keyval_create('Nsp_COPY_FN',...
				  'Nsp_DEL_FN',{})	// using Nsp skeletns
	A     = [1 2;3 4]
         = MPE_Attr_put (MPI_COMM_WORLD, kv, A)		// PERSISTENT this time
//-------------------------------------------------------// we get this msg
// MPE_Attr_put: 1 persist. copies
//-------------------------------------------------------
[WORLD2] = MPI_Comm_dup (MPI_COMM_WORLD)			// Nspcallback called
//-------------------------------------------------------// we get this msg
// Copying attribute keyed 19
//      1     2
//      3     4
//  
// Nsp copy wrapper0x0: 2 pers. copies
//-------------------------------------------------------
[C flag] = MPI_Attr_get ( WORLD2, kv)		// pers. by Nsp_COPY
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// pers. by MPE_Attr

      A       = 'hello'
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// both persistent
[C flag] = MPI_Attr_get ( WORLD2, kv)

clear A							// no problem
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// wanted to do this
[C flag] = MPI_Attr_get ( WORLD2, kv)		// for a long time :-)

         = MPI_Attr_delete(MPI_COMM_WORLD, kv)		// dissapears from WORLD
//-------------------------------------------------------// we get this msg
// Currently, callback in MPI_KEYVAL.mexglx will delete this:
//      1     2
//      3     4
//  
// Nsp dlt wrapper0x0: 1 pers. copies
//-------------------------------------------------------
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// flag==0, B={}
[C flag] = MPI_Attr_get ( WORLD2, kv)		// flag==1, C==[1 2;3 4]

         = MPI_Attr_delete(WORLD2,kv)		// real clean exit now
//-------------------------------------------------------// we get this msg
// Currently, callback in MPI_KEYVAL.mexglx will delete this:
//      1     2
//      3     4
//  
// Nsp dlt wrapper0x0: 0 pers. copies
//-------------------------------------------------------
         = MPI_Keyval_free (kv), kv		// operations commute
 kv           ==MPI_KEYVAL_INVALID			// same result
         = MPI_Comm_free (WORLD2)			// can make it both ways
 strcmp(WORLD2, MPI_COMM_NULL)
[B flag] = MPI_Attr_get (MPI_COMM_WORLD, kv)		// Key finished
         ==MPI_ERR_KEYVAL
[C flag] = MPI_Attr_get ( WORLD2, kv)		// Comm finished
         ==MPI_ERR_COMM

MPI_Finalize						// really wanted a
quit							// clean exit :-)


// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// More attribute Caching
// MPI_ TAG_UB, _HOST, _IO, _WTIME_IS_GLOBAL, _UNIVERSE_SIZE
// MPI_APPNUM, _WIN_BASE, _WIN_SIZE, _WIN_DISP_UNIT
// LAM_UNIVERSE_NCPUS, _NNODES
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Try with a lamboot of say 9 computers
// if some of them are multiprocessor, better (difference nnodes/ncpus)
// just LAM/MPI required on remaining 8 computers
// and Nsp license on 1st computer (in addition to local license)
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////

putenv(['LAM_MPI_SSI_rpi=tcp']), MPI_Init

args={}					// choose one of these
args={'-display',getenv('DISPLAY')}	// depending on your rsh/ssh environ

args={args{:},'-e','matlab','-nosplash','-nojvm','-r','startup_mergeParent'}

[children errs] = MPI_Comm_spawn ('/usr/X11R6/bin/xterm',args,...
					1,'NULL',0,'SELF')
[NEWORLD] = MPI_Intercomm_merge (children, 0)


		////////// ON CHILD ////////////////////
		MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')
		[rank]=MPI_Comm_rank(NEWORLD)	// rank==1

MPI_Errhandler_set( MPI_COMM_WORLD,'RETURN')			// avoid future aborts
MPI_Errhandler_set(NEWORLD ,'RETURN')
[rank]=MPI_Comm_rank(NEWORLD)			// rank==0

help environ
help MPI_ATTRKVALS

//-------------------------------------------------------
// Those 18 predefined attributes are the reason that
// newly created attributes are keyed from 19 on
//-------------------------------------------------------

//-------------------------------------------------------
// We already used processor_name
//-------------------------------------------------------

help MPI_Get_processor_name
[name] = MPI_Get_processor_name

//-------------------------------------------------------
// Accessing predefined attributes
//-------------------------------------------------------

help MPI_TAG_UB
[tub flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_TAG_UB)	// tub==2147483647
format bank						// to see more digits
tub
format							// back to short mode
tub
fprintf('//d\n',tub)					// this is better

		////////// ON CHILD ////////////////////
		[tub flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_TAG_UB)
		fprintf('//d\n',tub)			// same 2147483647
		MPI_Send(tub,0,tub,NEWORLD)

[stat]=MPI_Recv(tub ,1,tub ,NEWORLD)
tub1       =tub+1
      =MPI_Send(tub1,1,tub1,NEWORLD)		// info==4
[nfo msg]  =MPI_Error_string(info)			// invalid tag
     ==MPI_ERR_TAG

[tub flag]=MPI_Attr_get(NEWORLD,MPI_TAG_UB)	// tub=={}, flag==0

//-------------------------------------------------------

help MPI_HOST
[hst flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_HOST)		// hst==0, rank0 here

		////////// ON CHILD ////////////////////
		[hst flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_HOST)	// hst==-2
		hst==MPI_PROC_NULL			// not in local procssor

[hst flag]=MPI_Attr_get(NEWORLD,MPI_HOST)		// hst=={} not found

//-------------------------------------------------------

help MPI_IO
[mio flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_IO)		// mio==-1
mio==MPI_ANY_SOURCE

		////////// ON CHILD ////////////////////
		[mio flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_IO)	// mio==-2
		mio==MPI_PROC_NULL

[mio flag]=MPI_Attr_get(NEWORLD,MPI_IO)		// mio not found

//-------------------------------------------------------

help MPI_WTIME_IS_GLOBAL
[wtg flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_WTIME_IS_GLOBAL) // wtg==0, implies

		////////// ON CHILD ////////////////////
		[wtg flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_WTIME_IS_GLOBAL) //==0

[wtg flag]=MPI_Attr_get(NEWORLD,MPI_WTIME_IS_GLOBAL) // wtg not found

//-------------------------------------------------------

help MPI_UNIVERSE_SIZE
[mus flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_UNIVERSE_SIZE)	// mus==9 lambooted

	[mus flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_UNIVERSE_SIZE)	// mus==9
	[mus flag]=MPI_Attr_get(NEWORLD,MPI_UNIVERSE_SIZE)	// not found

[mus flag]=MPI_Attr_get(NEWORLD,MPI_UNIVERSE_SIZE)	// mus not found

//-------------------------------------------------------

help MPI_APPNUM
[apn flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_APPNUM)	// apn not found
[apn flag]=MPI_Attr_get(NEWORLD,MPI_APPNUM)	// anywhere

	[apn flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_APPNUM) // apn==0!!! found!
	[apn flag]=MPI_Attr_get(NEWORLD,MPI_APPNUM) // not here

//-------------------------------------------------------

help MPI_WIN_BASE
[mwb flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_WIN_BASE)	// info==28
[mwb flag]=MPI_Attr_get(NEWORLD,MPI_WIN_BASE)	// info==28

	[mwb flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_WIN_BASE) // info==28
	[mwb flag]=MPI_Attr_get(NEWORLD,MPI_WIN_BASE) // info==28

[nfo msg]  =MPI_Error_string(info)			// invalid key value
     ==MPI_ERR_KEYVAL				// !?! it's documented

//-------------------------------------------------------

help MPI_WIN_SIZE
[mws flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_WIN_SIZE)	// info==28
[mws flag]=MPI_Attr_get(NEWORLD,MPI_WIN_SIZE)	// info==28

	[mws flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_WIN_SIZE) // info==28
	[mws flag]=MPI_Attr_get(NEWORLD,MPI_WIN_SIZE) // info==28

[nfo msg]  =MPI_Error_string(info)			// invalid key value
     ==MPI_ERR_KEYVAL				// !?! it's documented

//-------------------------------------------------------

help MPI_WIN_DISP_UNIT
[mwu flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_WIN_DISP_UNIT)	// info==28
[mwu flag]=MPI_Attr_get(NEWORLD,MPI_WIN_DISP_UNIT)	// info==28

	[mwu flag]=MPI_Attr_get(MPI_COMM_WORLD,MPI_WIN_DISP_UNIT)	// info==28
	[mwu flag]=MPI_Attr_get(NEWORLD,MPI_WIN_DISP_UNIT)	// info==28

[nfo msg]  =MPI_Error_string(info)			// invalid key value
     ==MPI_ERR_KEYVAL				// !?! it's documented

//-------------------------------------------------------

help LAM_UNIVERSE_NCPUS
[lnc flag]=MPI_Attr_get(MPI_COMM_WORLD,LAM_UNIVERSE_NCPUS) // gosh! crashes!!!

//-------------------------------------------------------
// on xterm prompt
//-------------------------------------------------------
reset
lamclean
matlab
//-------------------------------------------------------
// trying again
//-------------------------------------------------------

putenv(['LAM_MPI_SSI_rpi=tcp']), MPI_Init

args={}					// choose one of these
args={'-display',getenv('DISPLAY')}	// depending on your rsh/ssh environ

args={args{:},'-e','matlab','-nosplash','-nojvm','-r','startup_mergeParent'}

[children errs] = MPI_Comm_spawn ('/usr/X11R6/bin/xterm',args,...
					1,'NULL',0,'SELF')
[NEWORLD] = MPI_Intercomm_merge (children, 0)


		////////// ON CHILD ////////////////////
		MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')
		[rank]=MPI_Comm_rank(NEWORLD)	// rank==1

MPI_Errhandler_set( MPI_COMM_WORLD,'RETURN')			// trying to survive
MPI_Errhandler_set(NEWORLD ,'RETURN')
[rank]=MPI_Comm_rank(NEWORLD)			// rank==0

//-------------------------------------------------------

help LAM_UNIVERSE_NCPUS

	////////// ON CHILD ////////////////////
	[lnc flag]=MPI_Attr_get(MPI_COMM_WORLD,LAM_UNIVERSE_NCPUS) // ouch!!

[lnc flag]=MPI_Attr_get(NEWORLD,LAM_UNIVERSE_NCPUS) // not found, good

//-------------------------------------------------------

help LAM_UNIVERSE_NNODES
[lnn flag]=MPI_Attr_get(MPI_COMM_WORLD,LAM_UNIVERSE_NNODES) // ouch!!

//-------------------------------------------------------
// on xterm prompt
//-------------------------------------------------------
reset
lamclean
matlab
//-------------------------------------------------------
// again trying again
//-------------------------------------------------------

putenv(['LAM_MPI_SSI_rpi=tcp']), MPI_Init

args={}					// choose one of these
args={'-display',getenv('DISPLAY')}	// depending on your rsh/ssh environ

args={args{:},'-e','matlab','-nosplash','-nojvm','-r','startup_mergeParent'}

[children errs] = MPI_Comm_spawn ('/usr/X11R6/bin/xterm',args,...
					1,'NULL',0,'SELF')
[NEWORLD] = MPI_Intercomm_merge (children, 0)


		////////// ON CHILD ////////////////////
		MPI_Errhandler_set(MPI_COMM_WORLD,'RETURN')
		[rank]=MPI_Comm_rank(NEWORLD)	// rank==1

MPI_Errhandler_set( MPI_COMM_WORLD,'RETURN')			// trying to survive
MPI_Errhandler_set(NEWORLD ,'RETURN')
[rank]=MPI_Comm_rank(NEWORLD)			// rank==0

//-------------------------------------------------------

help LAM_UNIVERSE_NCPUS

	////////// ON CHILD ////////////////////
	[lnc flag]=MPI_Attr_get(NEWORLD,LAM_UNIVERSE_NCPUS) // not found

[lnc flag]=MPI_Attr_get(NEWORLD,LAM_UNIVERSE_NCPUS) // not found, good

//-------------------------------------------------------

help LAM_UNIVERSE_NNODES
[lnn flag]=MPI_Attr_get(NEWORLD,LAM_UNIVERSE_NNODES)	// not found

	[lnn flag]=MPI_Attr_get(NEWORLD,LAM_UNIVERSE_NNODES) // flag==0

//-------------------------------------------------------
// crashing it all
//-------------------------------------------------------

	[lnc flag]=MPI_Attr_get(MPI_COMM_WORLD,LAM_UNIVERSE_NCPUS) // gosh!

[lnn flag]=MPI_Attr_get(MPI_COMM_WORLD,LAM_UNIVERSE_NNODES)	// crashed!
//-------------------------------------------------------
// MPI process rank 0 (n8, p17814) caught a SIGSEGV.
//-------------------------------------------------------
reset
lamclean


// ====================================================================
// Topologies
// MPI_[Dims/Cart/Graph]_create, Topo_test, [Graphdims/Cartdim]_get
// MPI_[Graph/Cart]_get, Cart_[rank/coords], Graph_neighbors[_count]
// MPI_Cart_[shift/sub], MPI_[Cart/Graph]_map
// ====================================================================
// requires 9-machines-wide LAMBOOT :-) your computer last in $LAMBHOST
// ok, it may work with less hosts, as long as you don't have license problems
// anyways, for real work, it's not a big deal 2 Nsps on same host
// (even on a biprocessor ;-)
// ===================================================================

putenv(['LAM_MPI_SSI_rpi=tcp']), MPI_Init

args={}					// choose one of these
args={'-display',getenv('DISPLAY')}	// depending on your rsh/ssh environ

args={args{:},...
	'-e','matlab','-nosplash','-nojvm',...
	'-r','startup_mergeParent'}

[children errs] = MPI_Comm_spawn ('/usr/X11R6/bin/xterm',args,...
					8,'NULL',0,'SELF')
[NEWORLD] = MPI_Intercomm_merge (children, 0)


MPI_Errhandler_set( MPI_COMM_WORLD,'RETURN')			// avoid future aborts
MPI_Errhandler_set(NEWORLD ,'RETURN')
[Nsiz]=MPI_Comm_size(NEWORLD)			// Nsiz==9
[Nrnk]=MPI_Comm_rank(NEWORLD)			// Nrnk==0
!hostname

		////////// ON CHILDREN //////////////
		MPI_Errhandler_set(MPI_COMM_WORLD ,'RETURN')
		[Nrnk]=MPI_Comm_rank(NEWORLD)	// Nrnk==i

//-------------------------------------------------------
// CAVEAT: take your time to order your xterm/Nsp windows based on rank
// (that's why we asked for it right now)
// or you'll end up in a mess and you won't understand anything at all
// Yet better if you "cascade" them such that you can copy-paste easily
// In this session, you will paste the same text to every child Nsp
// You have been warned :-)
//-------------------------------------------------------

help topo

//-------------------------------------------------------
// Topologies: computing dims (not so easy if you have 1,354,297 computers)
//-------------------------------------------------------
help MPI_Dims_create

ndims = 2					// for grid?
[dims] = MPI_Dims_create (Nsiz, ndims)	// (3x3)
ndims = 3					// for hypercube?
[dims] = MPI_Dims_create (Nsiz, ndims)	// (3x3x1), because 9=3x3
[dims] = MPI_Dims_create (8   , ndims)	// (2x2x2), for 8 in 3-D grid
[dims] = MPI_Dims_create (8   ,     2)	// (4x2),   for 8 in 2-D


//-------------------------------------------------------
// Topologies: creating cartesian dims (low level call, DO NOT USE)
//-------------------------------------------------------
//  <->	0  <->	1  <->	2  <->		(3x3) grid
//	^	^	^		periodic (toroidal)
//	|	|	|		on 2nd dim (less significant)
//	v	v	v
//  <->	3  <->	4  <->	5  <->
//	^	^	^
//	|	|	|
//	v	v	v
//  <->	6  <->	7  <->	8  <->
//-------------------------------------------------------
help MPI_Cart_map

dims    = [3 3]					// (3x3) grid
periods = [0 1]					// periodic in 2nd dim
[rank]=MPI_Cart_map(NEWORLD,dims,periods)	// this would be rank 0
Nrnk==rank					// not too difficult :-)

		////////// ON CHILDREN //////////////
		[rank]=MPI_Cart_map(NEWORLD,[3 3],[0 1])
		Nrnk==rank			// not too difficult :-)


//-------------------------------------------------------
// Topologies: creating graph dims (low level call, DO NOT USE)
//-------------------------------------------------------
//			   ----	3  <->	4 -	This time we want a ring
//			 /		    \	  like this one with 8 ranks
//	0  <->	1  <->	2 		     5     and rank 8 out again
//			 \		    /	  We indicate neighbors and
//			   ----	7  <->	6 -	 cummulative neighbor count
//-------------------------------------------------------
// nodes	   0	1	2	3	4	5	6	7	8
// neighbs  1	0-2	1-3-7	2-4	3-5	4-6	5-7	6-2	-
// num-ngb  1	2	3	2	 2	 2	 2	 2	 0
// cummul   1	3	6	8	10	12	14	16	16
// index = [1	3	6	8	10	12	14	16  ]
// edges = [1	0 2	1 3 7	2 4	3 5	4 6	5 7	6 2 ]
//-------------------------------------------------------
help MPI_Graph_map

index =  [1 3   6     8  10  12  14  16  ]	// see above
edges =  [1 0 2 1 3 7 2 4 3 5 4 6 5 7 6 2]	// rnk8 makes graph unconnected

[grnk]=MPI_Graph_map(NEWORLD,index,edges)	// grnk==0
Nrnk==grnk					// again, not too complex

		////////// ON CHILDREN //////////////
		index =  [1 3   6     8  10  12  14  16 0]
		edges =  [1 0 2 1 3 7 2 4 3 5 4 6 5 7 6 2]
		[grnk]=MPI_Graph_map(NEWORLD,index,edges)
		Nrnk==grnk			// not too complex either :-)

		// Subtle difference between undefined and unconnected

		////////// ON CHILD #8 //////////////
		index =  [1 3   6     8  10  12  14  16  ]	// missing 0 end
		[grnk]=MPI_Graph_map(NEWORLD,index,edges)
		      grnk==MPI_UNDEFINED


//-------------------------------------------------------
// Topologies: creating cartesian communicator
//-------------------------------------------------------
//  <->	0  <->	1  <->	2  <->		(3x3) grid
//	^	^	^		periodic (toroidal)
//	|	|	|		on 2nd dim (less significant)
//	v	v	v
//  <->	3  <->	4  <->	5  <->
//	^	^	^
//	|	|	|
//	v	v	v
//  <->	6  <->	7  <->	8  <->
//-------------------------------------------------------
help MPI_Cart_create
help MPI_Topo_test

[GRID]=MPI_Cart_create(NEWORLD,[3 3],[0 1],1)	// reorder flag ignored
[topo]=MPI_Topo_test  (GRID)			// collective call,
      topo==MPI_CART					// got blocked

		////////// ON CHILDREN //////////////
		[GRID]=MPI_Cart_create(NEWORLD,[3 3],[0 1],1)
		[Gsiz]=MPI_Comm_size  (GRID)	// Gsiz==9
		[Grnk]=MPI_Comm_rank  (GRID)
		Nrnk==Grnk				// easy to understand

[topo] = MPI_Topo_test(NEWORLD)
      topo == MPI_UNDEFINED
[Gsiz] = MPI_Comm_size(GRID)
[Grnk] = MPI_Comm_rank(GRID)
Nrnk==Grnk

		////////// ON ANY CHILD ////////////// Just 1 of them, just to check
		[eh]=MPI_Errhandler_get(GRID)	// inherited from
		strcmp(eh,MPI_ERRORS_RETURN)		// NEWORLD?

[eh]=MPI_Errhandler_get(GRID)
strcmp(eh,MPI_ERRORS_RETURN)
//   =MPI_Errhandler_set(GRID,'RETURN')		// not required then


//-------------------------------------------------------
//
//	  |/	  |/
//	--4 ----- 5 --	Yup, I know you can't see anything
//	|/|	|/|	Well, imagine it's a 3-D (not-so-hyper-)cube
//    --	0 +----	1 +-		  4 -- 5
//      /|-6 ---/+ 7 --		0 -- 1 |
//	|/|	|/|		| 6 -|-7
//    --	2 -----	3 --		2 -- 3
//      /       /	and it's periodic (cyclic, toroidal) on each dim
//
//-------------------------------------------------------

[CUBE]=MPI_Cart_create(NEWORLD,[2 2 2],[1 0 1],1)	// rank 8 is out!
[topo]=MPI_Topo_test  (CUBE)
      topo==MPI_CART

		////////// ON CHILDREN //////////////
		[CUBE]=MPI_Cart_create(NEWORLD,[2 2 2],[1 0 1],1)
		[Csiz]=MPI_Comm_size  (CUBE)
		[Crnk]=MPI_Comm_rank  (CUBE)
		Nrnk==Crnk				// boring to understand

		// If you unadvertently pasted the previous text to child 8,
		// simply skip next instruction in child 8 (it's a collective)

		////////// ON CHILD #8 //////////////
		[CUBE]=MPI_Cart_create(NEWORLD,[2 2 2],[1 0 1],1)
		strcmp(CUBE,MPI_COMM_NULL)		// rank 8 out, CUBE=NULL
		[topo]=MPI_Topo_test  (CUBE)	// info==5
		     ==MPI_ERR_COMM		// can't ask topo(NULL)

[Csiz] = MPI_Comm_size(CUBE)
[Crnk] = MPI_Comm_rank(CUBE)
Nrnk==Crnk						// boringer and boringer

		////////// ON ANY CHILD ////////////// Just 1 of them, just to check
		[eh]=MPI_Errhandler_get(CUBE)	// inherited from
		strcmp(eh,MPI_ERRORS_RETURN)		// NEWORLD?

[eh]=MPI_Errhandler_get(CUBE)
strcmp(eh,MPI_ERRORS_RETURN)
//   =MPI_Errhandler_set(CUBE,'RETURN')		// not required then


//-------------------------------------------------------
// Topologies: creating graph communicator
//-------------------------------------------------------
//			   ----	3  <->	4 -	The same previous ring
//			 /		    \
//	0  <->	1  <->	2 		     5
//			 \		    /
//			   ----	7  <->	6 -
//-------------------------------------------------------
// nodes	   0	1	2	3	4	5	6	7	8
// neighbs  1	0-2	1-3-7	2-4	3-5	4-6	5-7	6-2	-
// num-ngb  1	2	3	2	 2	 2	 2	 2	 0
// cummul   1	3	6	8	10	12	14	16	16
//-------------------------------------------------------
help MPI_Graph_create

index =  [1 3   6     8  10  12  14  16  ]	// see above
edges =  [1 0 2 1 3 7 2 4 3 5 4 6 5 7 6 2]	// rnk8 makes graph unconnected

[RING] = MPI_Graph_create (NEWORLD,index,edges,1)	// reorder flag ignored
[topo] = MPI_Topo_test    (RING)
      topo == MPI_GRAPH

		////////// ON CHILD #8 //////////////
		index =  [1 3   6     8  10  12  14  16  ]
		edges =  [1 0 2 1 3 7 2 4 3 5 4 6 5 7 6 2]
		[RING]=MPI_Graph_create(NEWORLD,index,edges,1)
		strcmp(RING,MPI_COMM_NULL)
		[topo]=MPI_Topo_test   (RING)	// rank 8 out

		////////// ON CHILDREN //////////////
		index =  [1 3   6     8  10  12  14  16  ]
		edges =  [1 0 2 1 3 7 2 4 3 5 4 6 5 7 6 2]
		[RING]=MPI_Graph_create (NEWORLD,index,edges,1)
		[Rsiz]=MPI_Comm_size(RING)		// it says 8
		[Rrnk]=MPI_Comm_rank(RING)
		Nrnk==Rrnk				// speechless

[Rsiz] = MPI_Comm_size(RING)			// 8 didn't enter graph
[Rrnk] = MPI_Comm_rank(RING)
Nrnk==Rrnk						// no comment

		////////// ON ANY CHILD ////////////// Just 1 of them, just to check
		[eh]=MPI_Errhandler_get(RING)
		strcmp(eh,MPI_ERRORS_RETURN)

[eh]=MPI_Errhandler_get(RING)			// protected to aborts
strcmp(eh,MPI_ERRORS_RETURN)


//-------------------------------------------------------
// Topologies: retrieving graph connectivity (topology) and neighbors
//-------------------------------------------------------
help MPI_Graphdims_get
help MPI_Graph_get

[nnodes nedges] = MPI_Graphdims_get (CUBE)		// nope! not a graph!
[info2       string] = MPI_Error_string  (info)		// invalid topology
               == MPI_ERR_TOPOLOGY
[nnodes nedges] = MPI_Graphdims_get (RING)		// 8 nodes, 16 edges

		////////// ON ANY CHILD ////////////// Just 1 of them, just to check
		[nnodes nedges] = MPI_Graphdims_get (RING)

clear index edges
[index edges] = MPI_Graph_get (RING)
nnodes==length(index)
nedges==length(edges)

		////////// ON ANY CHILD ////////////// Just 1 of them, just to check
		index, edges
		[index2 edges2] = MPI_Graph_get (RING)
		isequal(index,index2)
		isequal(edges,edges2)

[nneighs]= MPI_Graph_neighbors_count (RING, Rrnk)	// rank in RING (same)
[neighs] = MPI_Graph_neighbors       (RING, Rrnk)	// can ask for other rnk

//-------------------------------------------------------
//			   ----	3  <->	4 -	Recall rank 0 is neighbor to 1
//			 /		    \
//	0  <->	1  <->	2 		     5
//			 \		    /
//			   ----	7  <->	6 -	rank 2 neighbor to 1,3,7
//-------------------------------------------------------

		////////// ON ANY CHILD ////////////// Just 1 of them, just to check
		[nneighs]= MPI_Graph_neighbors_count (RING, 2)
		[neighs] = MPI_Graph_neighbors       (RING, 2)


//-------------------------------------------------------
// Topologies: retrieving cartesian connectivity
//-------------------------------------------------------
help MPI_Cartdim_get
help MPI_Cart_get

[ndims] = MPI_Cartdim_get (RING)		// nope! it's not cartesian
[inf2 strng] = MPI_Error_string(info)		// invalid topology
       == MPI_ERR_TOPOLOGY
[ndims] = MPI_Cartdim_get (CUBE)

clear dims periods
[dims periods coords]=MPI_Cart_get(CUBE)	// each one asks for theirs

		////////// ON EVERY CHILD ////////////// In order to understand dim signif
		[dims periods coords]=MPI_Cart_get(CUBE)

//-------------------------------------------------------
// Particularly, for this simple case (3-D cube)
// the 3 dims are of length 2 and coords==rank in binary
// ie: rank 0 is at (0 0 0), rank 1 at (0 0 1), 2 @ (0 1 0)
//-------------------------------------------------------
//-------------------------------------------------------
//
//	  |/	  |/
//	--4 ----- 5 --	In these drawings, less significant dimension is
//	|/|	|/|	horizontal (X), next is vertical (Y), Z is depth
//    --	0 +----	1 +-		  4 -- 5
//      /|-6 ---/+ 7 --		0 -- 1 |
//	|/|	|/|		| 6 -|-7
//    --	2 -----	3 --		2 -- 3
//      /       /
//
//-------------------------------------------------------
// Topologies: rank <-> coords translations
//-------------------------------------------------------
help MPI_Cart_rank
help MPI_Cart_coords

[rank  ]=MPI_Cart_rank  (CUBE, coords)
isequal(rank,Crnk,Nrnk,Grnk)
[coords]=MPI_Cart_coords(CUBE, 0)		// can ask for other's

		////////// ON ANY CHILD ////////////// Just 1 of them, just to check
		[coords]=MPI_Cart_coords(CUBE, 0     )
		[rank] = MPI_Cart_rank  (CUBE, coords)

		[coords]=MPI_Cart_coords(CUBE, 1     )
		[rank] = MPI_Cart_rank  (CUBE, coords)

		[coords]=MPI_Cart_coords(CUBE, 2     )
		[rank] = MPI_Cart_rank  (CUBE, coords)

		[coords]=MPI_Cart_coords(CUBE, 4     )
		[rank] = MPI_Cart_rank  (CUBE, coords)


//-------------------------------------------------------
// Topologies: retrieving cartesian routing
//-------------------------------------------------------
help MPI_Cart_shift

disp=1, horz=2, vert=1, dpth=0		// 1 hop along each dim Horz/Vert/Depth

[Hsrc Hdst] = MPI_Cart_shift (CUBE, horz, disp)	// LSB routing
[Vsrc Vdst] = MPI_Cart_shift (CUBE, vert, disp)	// middle-bit
[Dsrc Ddst] = MPI_Cart_shift (CUBE, dpth, disp)	// MSB routing

//-------------------------------------------------------
//	  4 -- 5	Who is routed from/to node 0?
//	0 -- 1 |	Horizontal: from left (none) / to right (1)
//	| 6 -|-7	Vertical:   from up   (none) / to down  (2)
//	2 -- 3		Depth:      from fore (none) / to back  (4)
//			but recall all dims are periodic, so 1,2,4 again
//-------------------------------------------------------
		////////// ON ANY CHILD ////////////// 7 is a nice dual of 0, all bits set
		disp=1, horz=2, vert=1, dpth=0
		[Hsrc Hdst]=MPI_Cart_shift(CUBE,horz,disp)	// LSB routing
		[Vsrc Vdst]=MPI_Cart_shift(CUBE,vert,disp)	// middle-bit
		[Dsrc Ddst]=MPI_Cart_shift(CUBE,dpth,disp)	// MSB routing

MPI_PROC_NULL				// when no neighbor, PROC_NULL
help MPI_PROC_NULL			// pity I set all periodic
help MPI_Cart_shift			// ok, I'm setting vert not periodic


//-------------------------------------------------------
// Topologies: sub-dividing cartesian grid in cartesian sub-grids
//-------------------------------------------------------
help MPI_Cart_sub

[NEW] = MPI_Cart_sub (RING, [0 1])		// nope! not cartesian
     == MPI_ERR_TOPOLOGY			// we'll drop depth (most signf)
[LMS] = MPI_Cart_sub (CUBE, [0 1 1])	// two halves: MSHalf, LSHalf

		////////// ON EVERY CHILD but 8 //////////////
		// Blocks only in CUBE, rank 8 not needed
		//-------------------------------------------------------
		[LMS] = MPI_Cart_sub (CUBE, [0 1 1])

//-------------------------------------------------------
//	  4 -- 5	From depth point of view, there is a
//	0 -- 1 |	less significant half 0,1,2,3, and a
//	| 6 -|-7	most significant half 4,5,6,7
//	2 -- 3		(the MSHalf has the least depth coordinate possible)
// In more general cases, several communicators ("slices") are produced
// with ranks sharing the same coordinate in the dropped dimension(s)
//-------------------------------------------------------
[Lsiz]= MPI_Comm_size(LMS)			// 2 halves sized (2x2)
[Lrnk]= MPI_Comm_rank(LMS)			// of course we're on LSHalf

		////////// ON EVERY CHILD but 8 //////////////
		[Lsiz]= MPI_Comm_size(LMS)	// Lsiz==4
		[Lrnk]= MPI_Comm_rank(LMS)	// first 4 are 0,1,2,3 and again

[gLMS] = MPI_Comm_group (LMS)
[gsiz] = MPI_Group_size (gLMS)
[gCUB] = MPI_Comm_group (CUBE)
[info gsiz] = MPI_Group_size (gCUB)
[info ranks]= MPI_Group_translate_ranks(gCUB, 0:gsiz-1, gLMS)
isequal(ranks,...
 [0 1 2 3 MPI_UNDEFINED MPI_UNDEFINED MPI_UNDEFINED MPI_UNDEFINED])

		////////// ON SOME CHILD of the other half  ////////////// say, 7
		[info gLMS] = MPI_Comm_group (LMS)
		[info gsiz] = MPI_Group_size (gLMS)
		[info gCUB] = MPI_Comm_group (CUBE)
		[info gsiz] = MPI_Group_size (gCUB)
		[info ranks]= MPI_Group_translate_ranks(gCUB,0:gsiz-1,gLMS)
		isequal(ranks,[MPI_UNDEFINED*ones(1,4) 0 1 2 3])

//-------------------------------------------------------
// Recall the 2-D GRID we built at the beginning
//
//	0 - 1 - 2	We can cut it in 3 slices
//	|   |   |	along the... say... most significant direction
//	3 - 4 - 5	(vertical), so that it goes 0,1,2, then 3,4,5,
//	|   |   |	and finally 6,7,8 (in three different comms)
//	6 - 7 - 8
//
//-------------------------------------------------------

[info ROW] = MPI_Cart_sub (GRID, [0 1])		// row-sliced, throwing row dim
[info size]= MPI_Comm_size(ROW)			// 3 rows sized (1x3)

		////////// ON EVERY CHILD //////////////
		[info ROW] = MPI_Cart_sub (GRID, [0 1])
		[info size]= MPI_Comm_size(ROW)

[info gROW] = MPI_Comm_group (ROW)
[info gsiz] = MPI_Group_size (gROW)
[info gRID] = MPI_Comm_group (GRID)
[info gsiz] = MPI_Group_size (gRID)
[info ranks]= MPI_Group_translate_ranks(gRID, 0:gsiz-1, gROW)
isequal(ranks, [0 1 2 MPI_UNDEFINED*ones(1,6)])

		////////// ON EVERY CHILD //////////////
		[info gROW] = MPI_Comm_group (ROW)
		[info gsiz] = MPI_Group_size (gROW)
		[info gRID] = MPI_Comm_group (GRID)
		[info gsiz] = MPI_Group_size (gRID)
		[info ranks]= MPI_Group_translate_ranks(gRID,0:gsiz-1,gROW)
		isequal(ranks, [0 1 2 MPI_UNDEFINED*ones(1,6)])
		isequal(ranks, [MPI_UNDEFINED*ones(1,3) 0 1 2,...
				MPI_UNDEFINED*ones(1,3)])
		isequal(ranks, [MPI_UNDEFINED*ones(1,6) 0 1 2])


//-------------------------------------------------------
// way too tired to try to clean
//-------------------------------------------------------
////// MPI_Finalize
////// quit
//////		////////// ON EVERY CHILD //////////////
//////		MPI_Finalize
//////		quit

quit
lamhalt		# in the shell. This was faster ;-)

-------------------------------------------------------------------------------
*******************************************************************************
-------------------------------------------------------------------------------

help mpi
help missing

//////////////////////////////////////
// WON'T IMPLEMENT //	//////////// datatypes //////////////
//////////////////////////////////////

	REASON: Nsp users have neither 1) info on how does Nsp
store variables, nor 2) way of controlling data layout

1) have no C definition of mxArray type
2) no constructions in Nsp languaje to impose data layout

For MPI layout-compatible types (mxDOUBLE but not mxCOMPLEX, and mxCHAR,
basically) some MPI type can be used (MPI_DOUBLE, _UNSIGNED_SHORT).

For the remaining (mxSTRUCT and mxCELL, basically), MPI_Pack/Unpack
must be used. These routines use MPI_PACKED datatype and MEX calls
(mxGetNumberOfElements, mxGetClassID, mxGetName...)
to analyze (at source) and reconstruct (at destination) the Nsp var.
Need or usefulness of datatypes is thus obviated

