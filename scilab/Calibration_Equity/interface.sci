
flag=1;
color=2;
while flag==1
 
  l1=list('Calibration Method',1,[' PDE Dupire simulator', 'Monte Carlo Weighted', 'Stochastic Control', 'Markovian Calibration']);
  rep1=x_choices('Methods Choice',list(l1));
  
  if rep1==1 then
    l2=list('OPTIONS',2,['PUT','CALL']);
    callorput=x_choices('Options Choiche',list(l2))-1; 
    x0=[590.0     //S_0
    0.06;        // r
    0.0262;     //q
    0.;        //time origin
       ];

    y0=[2.0         //maturity
       ];
    
    z0 = [400 //Space Step
          100; // Time Step
          ];
    ystr_models=['s'; 'r';'divid';'Time origin'];
    ystr_options=['Maturity'];
    ystr_methods=['Space Step';'Time Step'];

    x0=evstr(x_mdialog(['Model Parameters'],ystr_models,string(x0)));
    y0=evstr(x_mdialog(['Option Parameters'],ystr_options,string(y0)));
    z0=evstr(x_mdialog(['Method Parameters'],ystr_methods,string(z0)));
    chdir("./src/PDEDupire/src");
    exec builder.sce;
    exec loader.sce;
    
    if MSDOS then unix('del calib.in');
    else unix('rm -f calib.in');
    end,
   
    u = file('open','calib.in','unknown')
    fprintf(u,'%f\n',x0(1));
    fprintf(u,"%f\n",x0(2));
    fprintf(u,"%f\n",x0(3));
    fprintf(u,"%d\n",callorput);
    fprintf(u,"%f\n",x0(4));
    fprintf(u,"%f\n",y0(1));
    fprintf(u,"%f\n",-7);
    fprintf(u,"%f\n",7);
    fprintf(u,"%d\n",z0(1));
    fprintf(u,"%d\n",z0(2));
    fprintf(u,"%d\n",0);
    fprintf(u,"%f\n",.5);
    fprintf(u,"%d\n",2);

    fprintf(u,"optim2.in");

    fprintf(u,"sp500prices.data");

    fprintf(u,"sigmainit_n1m1_0_013.ddl");

    fprintf(u,"v2_008_1_sigmaest_n1m1_0_013.ddl");

    fprintf(u,"visusigma");

    fprintf(u,"sigmainit_n1m1_0_013.visu");

    fprintf(u,"v2_008_1_sigmaest_n1m1_0_013.visu");
    file('close',u);

    calib()

    flag=x_choose(['Stop';'Continue'],'Another Test ?');
    chdir("../..");
    end,


  if rep1==2 then
    x0=[120.0     //S_0
    0.06;        // r
    0.0262;     //q
       ];

    y0=[1 // Option Volatility
    1 // Option Constraint
    1.0;        //maturity
       ];
    
    z0 = [0.3 //Volatility of generated path
          0.001; // Time Step
          20; // Number of Calibration Put
          1; //Perturbation coefficient of volatility
          5.;// Return to mean
          0.7; // Correlation
          1;//Option constraint
          ];
    ystr_models=['s'; 'r';'divid'];
    ystr_options=['Option Volatility'; 'Option Constraint';'Maturity'];
    ystr_methods=['Volatility of generated path';'Time Step';'Number of Calibration Put';'Perturbation coefficient of volatility';
            'Return to mean';' Correlation';'Option constraint'];
    x_dialog("Message","The datas have to be in file ./src/CalibMc/datas.in") ;
    x0=evstr(x_mdialog(['Model Parameters'],ystr_models,string(x0)));
    y0=evstr(x_mdialog(['Option Parameters'],ystr_options,string(y0)));
    z0=evstr(x_mdialog(['Method Parameters'],ystr_methods,string(z0)));
    chdir("./src/Calib_Mc");
    exec builder.sce;
    exec loader.sce;
    if MSDOS then unix('del parameters.in');
    else unix('rm -f parameters.in');
    end,
   
    u = file('open','parameters.in','unknown')
    fprintf(u,'%f\n',x0(1));
    fprintf(u,"%f\n",x0(2));
    fprintf(u,"%f\n",x0(3));
    fprintf(u,"%f\n",y0(1));
    fprintf(u,"%f\n",y0(2));
    fprintf(u,"%f\n",y0(3));
    fprintf(u,"%d\n",z0(1));
    fprintf(u,"%d\n",z0(2));
    fprintf(u,"%d\n",z0(3));
    fprintf(u,"%d\n",z0(4));
    fprintf(u,"%d\n",z0(5));
    fprintf(u,"%d\n",z0(6));
    fprintf(u,"%d\n",z0(7));
    calib()

    flag=x_choose(['Stop';'Continue'],'Another Test ?');
    chdir("../..");
    end,

  if rep1==3 then
    x0=[100.0     //S_0
    0.05;        // r
    0.0;     //q
       ];

    y0=[
        100; //    N : number of space steps of the fine grid
        0.21; //  sigma_0 : --> prior
        0.1;//  sigma_min
        0.46; // sigma_max 
        0.48; // sigma_bar : -->
     ];
    z0 = [0.00001 // gradtol : tolerance on the relative gradient
    0.00001 // steptol : tolerance on the relative change of x
    1; // verbosity : level of printed information (0 --> 3)
    0; // saveSuccessiveXinFile : save successive x0 in the file data.out (0 or 1)
    100; // maxCounter : maximum number of iterations
    0; //  lambda : Initial default value of lagrange parametre
    1; // alpha : Tune this parametre when the program do not converge
          ];
    ystr_models=['s'; 'r';'divid'];
    ystr_options=['number of space steps';'sigma_0 : --> prior';'sigma_min';'sigma_max';'sigma_bar'];
    ystr_methods=['tolerance on gradient';'tolerance on change';'Verbosity';'saveSuccessiveXinFile';'maximum number of iterations'
        'Initial default value of lagrange parameters';'Tune this parameter when no convergence'];

    x0=evstr(x_mdialog(['Model Parameters'],ystr_models,string(x0)));
    y0=evstr(x_mdialog(['Option Parameters'],ystr_options,string(y0)));
    z0=evstr(x_mdialog(['Method Parameters'],ystr_methods,string(z0)));
    chdir("./src/Calib_Ave/src");
    exec builder.sce;
    exec loader.sce;
    
    if MSDOS then unix('del calib_Avellaneda.in');
    else unix('rm -f calib_Avellaneda.in');
    end,
    x_dialog("Message","The datas have to be in file ./src/CalibMc/VolStoPut.data") ;
  
    u = file('open','calib_Avellaneda.in','unknown')
    fprintf(u,'%f\n',x0(1));
    fprintf(u,"%f\n",x0(2));
    fprintf(u,"%f\n",x0(3));
    fprintf(u,"%f\n",y0(1));
    fprintf(u,"%f\n",y0(2));
    fprintf(u,"%f\n",y0(3));
    fprintf(u,"%f\n",y0(4));
    fprintf(u,"%f\n",y0(5));
    fprintf(u,"%d\n",z0(1));
    fprintf(u,"%d\n",z0(2));
    fprintf(u,"%d\n",z0(3));
    fprintf(u,"%d\n",z0(4));
    fprintf(u,"%d\n",z0(5));
    fprintf(u,"%d\n",z0(6));
    fprintf(u,"%d\n",z0(7));
    fprintf(u,"VolStoPut.data");
    fprintf(u,"VolStoLoc.out");
    file('close',u);

    calib_ave_scilab()

    flag=x_choose(['Stop';'Continue'],'Another Test ?');
    chdir("../../..");
    end,

  if rep1==4 then
    x0=[100.0     //S_0
       ];

    y0=[1.0         //maturity
    ];
    
    z0 = [ 12; //    N : number of time steps 
        20; //  M Number of space steps (large grid)
        3;//  Number of finite elements
        250; // Number of space steps (fine grid)
     ];

    ystr_models=['s'];
    ystr_options=['Maturity'];
    ystr_methods=['number of time steps ';'Number of space steps';'Number of finite elements';'Number of space steps'];

    x0=evstr(x_mdialog(['Model Parameters'],ystr_models,string(x0)));
    y0=evstr(x_mdialog(['Option Parameters'],ystr_options,string(y0)));
    z0=evstr(x_mdialog(['Method Parameters'],ystr_methods,string(z0)));
    chdir("./src/Markov/src");
    
    exec("builder.sce");
    exec("loader.sce");
    
    if MSDOS then unix('del param.in');
    else unix('rm -f param.in');
    end,
    x_dialog("Message","The datas have to be in file ./src/Markov/Call.dat") ;
  
    u = file('open','param.in','unknown')
    fprintf(u,'%f\n',x0(1));
    fprintf(u,"%f\n",y0(1));
    fprintf(u,"%d\n",z0(1));
    fprintf(u,"%d\n",z0(2));
    fprintf(u,"%d\n",z0(3));
    fprintf(u,"%d\n",z0(4));
    file('close',u);
    exec calib_ef.sci;



    flag=x_choose(['Stop';'Continue'],'Another Test ?');
    chdir("../../..");
    end,




end;
