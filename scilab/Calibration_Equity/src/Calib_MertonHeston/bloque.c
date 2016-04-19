gradFunction(dimx,x00,grad0);
for(i=0;i<dimx;i++)  if(grad0[i]==0.)  grad0[i] = 1.e-8;
//
//
//
tmp = fabs(grad[0]/grad0[0]);
ibidon = 0;
for(i=0;i<dimx;i++) if(fabs(grad[i]/grad0[i]) > tmp)
{
  tmp = fabs(grad[i]/grad0[i]);
  ibidon = i;
}
if(tmp == 0.) tmp = 1;
bloqueheston = 1;
for(un=0;un<3;un++){
  for(i=0;i<dimx;i++) {
    xmin[i] = xmin0[i];
    xmax[i] = xmax0[i];
    nbd[i]  = nbd0[i];
  }
  
  if(bloqueheston==1 && TypeModel!=2 && un == 0)
    {
      //	  un = -1*un;
      //  un = optim%3;	  
      printf(" On bloque les variables V0, theta et rho\n");
      for(i=0;i<dimx;i++) {
	if(i!=1 && i!=3 && i!=4) {
	  printf(" On bloque la variable %d \n",i);
	  xmin[i] = x0[i];
	  xmax[i] = xmin[i];
	  nbd[i] = 2;
	}
      }
    }
  else if(bloqueheston==1 && TypeModel!=2 && un == 1)
    {
      //	  un = -1*un;
      //	  un = optim%3;	  
      printf(" On bloque Heston \n");
      for(i=0;i<dimx;i++) {
	if(i!=0 && i!=5 && i!=6 && i!=7) {
	  printf(" On bloque la variable %d \n",i);
	  xmin[i] = x0[i];
	  xmax[i] = xmin[i];
	  nbd[i] = 2;
	}
      }
    }
  else if(bloqueheston==1 && TypeModel!=2 && un == 2)
    {
      //	  un = -1*un;
      // un = optim%3;	  
      printf(" On bloque tou sauf ================================== KAPPA \n");
      for(i=0;i<dimx;i++) {
	if(i!=1) {
	  printf(" On bloque la variable %d \n",i);
	  xmin[i] = x0[i];
	  xmax[i] = xmin[i];
	  nbd[i] = 2;
	}
      }
      
    }
  else if(TypeModel==2)
    {
      printf(" On bloque Heston \n");
      for(i=0;i<dimx;i++) {
	if(i!=1 && i!=2 && i!=3 & i!=4) {
	  printf(" On bloque la variable %d \n",i);
	  xmin[i] = x0[i];
	  xmax[i] = xmin[i];
	  nbd[i] = 2;
	}
      }
      
    }  else 
      {
	
	for(i=0;i<dimx;i++)  if(fabs(grad[i]/grad0[i])/tmp < 0.1)  //if(i!=1)  //if(fabs(grad[i]/grad0[i])/tmp < 0.1)
	  {
	    printf(" On bloque la variable %d \n",i);
	    //	  x0[i] = sol[i];
	    
	    xmin[i] = x0[i];
	    xmax[i] = xmin[i];
	    nbd[i] = 2;
	  }
      }
}



      // basculement de Merton a Merton + Heston 
  
  if( (TypeModelPrix==3) && (TypeModel==2) && ( (optim>2) || (fx < 1.e-2) ) )
	{
	  printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ On bascule de Merton a Merton + Heston \n");
	  LogNorme = 0;
	  TypeModel=3;
	  init_sol(TypeModel,&dimx,sol,dimx1,dimx2,dimx3,V0,kappa,theta,sigmav,rho,lambda,m0,v);
	  init_bornes(TypeModel,dimx,nbd,xmin,xmax);
	  init_bornes(TypeModel,dimx,nbd0,xmin0,xmax0);
  // for(i=0;i<dimx;i++) printf("nbd[%d] = %d\n",i,nbd[i]);
  // printf("press_a_key_to_continue...\n");
  // getchar();
	  for(i=1;i<4;i++) x0[i+4] = x0[i];
	  for(i=1;i<5;i++)
		{
		  rand = uniforme();
		  // voir ici pour modifier les xmin et xmax
		  x0[i] = xmin[i] + rand*(xmax[i]-xmin[i]);
		}
	  fprintf(ftest,"x0 bis  = ");
	  for(i=0;i<dimx;i++) fprintf(ftest,"%f, ",x0[i]); fprintf(ftest,"\n");

	}
