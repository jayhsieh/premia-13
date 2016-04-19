


flag=2;
color=2;
while flag==2
  
  l1=list('Calibration Method',1,[' PDE Dupire simulator']);
  rep1=x_choices('Methods Choice',list(l1));
  
  if rep1==1 then

    l2=list('OPTIONS',1,['PUT','CALL']);
    callorput=x_choices('Options Choiche',list(l2)); 


     

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

 
