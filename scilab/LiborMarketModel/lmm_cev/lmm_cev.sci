function [result]=lmm_cev_price(model, product,pricing,discr_scheme,alpha,Fspot,strike,expiry,swp_lnght,int_pts)
  disp("linking library.....");
  link(["lmm_cev_sci.o","lmm_cev_pricer.o"],"lmm_cev_sci","C");
  disp(".....done");
  _model=[model];
  _product=[product];
  _pricing=[pricing];
  _discr_scheme=[discr_scheme];
  _alpha=[alpha];
  _Fspot=[Fspot];
  _strike=[strike];
  _expiry=[expiry];
  _swp_lnght=[swp_lnght];
  _int_pts=[int_pts];
  c=fort("lmm_cev_sci",_model,2,"i",_product,3,"i",_pricing,4,"i",_discr_scheme,5,"i",_alpha,6,"d",_Fspot,7,"d",_strike,8,"d",_expiry,9,"d",_swp_lnght,10,"d",_int_pts,11,"i","out",[2,1],1,"d");
  result=c;
endfunction


// $$$ cev_price(size_t type_model,size_t type_product,size_t type_pricing, 
// $$$ 		 size_t type_scheme, double alpha, double fzero, double H,
// $$$ 		 double expiry,double swp_lenght,int int_pts,double out[2])
		 
disp("  //////////////////////////////////////////////");
disp("  //CEV/LCEV Pricer for Caplets and Swaptions //");
disp("  //////////////////////////////////////////////");
model=input("Choose model:CEV=0, LCEV=1  ");
product=input("Choose product:Caplet=0,  Swaption=1 ");
method=input("Choose pricing method: Closed formula=0, Monte Carlo=1 ");
if((model==1)&(product==0)&(method==0)) then,
  disp("Sorry, NO closed formula for LCEV caplets, I will use MC instead");
  method=1;
else, end  
if (method==1) then,
  if (product==1) then,
    disp("Sorry, Swaption MC pricing unavalaible, using closed formula instead");
    method=0;discret=3;num_MC=10000; 
  else,
    discret=input("Choose discretization scheme: Log-Euler=0, Euler=1 ");
    num_MC=input("Choose the number of MC samples ");
  end
else, discret=3;num_MC=10000; 
end

alpha=input("Choose CEV exponent ");
disp("Default time structure: MAX 20 years, fix and  floating ");
disp("legs  payments each 0.5 years");
fzero=input("Choose spot  FLAT Libor curve value ");
H=input("Choose the strike ");

check=%f;
while (check==%f) do 
 expiry=input("Choose the expiry of the option (years,<20)");
 if(expiry<20) then check=%t; else, end
end   
 
if(product==1) then, 
  check=%f;
  while (check==%f) do, 
    swp_lenght=input("Choose the duration of the swap (years)");
    if (expiry+swp_lenght<21) then check=%t; else, end
  end    
else swp_lenght=1;, end  
tmp=lmm_cev_price(model,product,method,discret,alpha,fzero,H,expiry,swp_lenght,num_MC);
printf("Price %f\n Error %f",tmp(1),tmp(2));
  
