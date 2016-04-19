M=premia_create()
// Choose a model 
models=premia_get_models();
// n=x_choose(models,['Choose a model']);
n=14;
M.set_model[n];
M.get_model[]
B=M.get_model_values[];

b=list( list( "A double" ,   4, %t ),...
	list( "A1 unsetable" ,   7.890, %f ),...
	list( "B double array" , 0:9, %t),...
	list( "C a numfunc" ,...
	      list(list( "F1" ,10, %t),...
		   list( "F2 unsetable" ,11, %f),...
		   list( "F3" ,124.500, %t))));

if or(B<>b) then pause;end 

// change the model variables 

b=list( list( "A double" ,   6, %t ),...
	list( "A1 unsetable" ,   7.890,%f ),...
	list( "B double array" , 0:9, %t ),...
	list( "C a numfunc" ,...
	      list(list( "F1" ,10*3, %t ),...
		   list( "F2 unsetable" ,11, %f ),...
		   list( "F3" ,-124.500*5, %t ))));

M.set_model_values[b];
 
B=M.get_model_values[];

if or(B<>b) then pause;end 

S=M.model_check[];

exec(getenv('NSP')+'/src/libpremia/macros.sci');

premia_model_values(M) 


  



