exec loader.sce;


flag=1;
color=2;
while flag==1
  
  l1=list('Calibration Method',1,[' PDE Dupire simulator']);
  rep1=x_choices('Methods Choice',list(l1));
  
  if rep1==1 then

    l2=list('OPTIONS',2,['PUT','CALL']);
    callorput=x_choices('Options Choiche',list(l2))-1; 


     

    x0=[590.0     //S_0
	0.06;        // r
	0.0262;    	//q
    0.;        //time origin
	   ];

    y0=[2.0     	//maturity
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


     
  end;


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

end,
