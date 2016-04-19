

function resu=lmm_swpt_stovol_sci(period , nb_fac , swpt_mat , swp_mat , perct)

list_file=[ "lmm_swpt_stovol_sci.o", "lmm_volatility.o" , "lmm_zero_bond.o" , "lmm_stochastic_volatility.o" , "lmm_random_generator.o" , "lmm_products.o" , "lmm_numerical.o" , "lmm_mathtools.o" , "lmm_libor.o" , "premia_files/complex.o " ]

link( list_file ,"lmm_swpt_stovol_sci","C")

tenor=[ period ];
numFac=[ nb_fac ];
swaptionMat=[ swpt_mat ];
swapMat=[ swp_mat ];
percent=[ perct ];
c=fort("lmm_swpt_stovol_sci",tenor,2,"d",numFac , 3, "i" , swaptionMat , 4 , "d", swapMat , 5 , "d" , percent , 6 , "d" , "out",[1,1],1,"d");
resu=c*100*100;

endfunction

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////         martingale X
////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function resu=lmm_cap_martX_sci(period , nb_fac , maturity , strike )

list_file=[ "lmm_martingaleX_sci.o", "lmm_martingaleX.o" , "lmm_volatility.o" , "lmm_zero_bond.o" ,  "lmm_random_generator.o" , "lmm_products.o" , "lmm_numerical.o" , "lmm_mathtools.o" , "lmm_libor.o" , "cumulfunc.o" ,  "dcdflib.o" , "ipmpar.o" ]

link( list_file ,"lmm_cap_martX_sci","C")

tenor=[ period ];
numFac=[ nb_fac ];
Mat=[ maturity ];
K=[strike];
numcap=round(maturity/period)

c=fort("lmm_cap_martX_sci",tenor,3,"d",numFac , 2, "i" , Mat , 5 , "d", strike , 4 , "d" , "out",[numcap , 1],1,"d");
m=zeros(numcap,2);
for i=1:numcap,
  m(i,1)=c(i),
  m(i,2)=(i-1)*period,
end;

resu=m;

endfunction

function resu=lmm_swpt_martX_sci(period , nb_fac , swpt_maturity , swp_maturity , strike )

list_file=[ "lmm_martingaleX_sci.o", "lmm_martingaleX.o" , "lmm_volatility.o" , "lmm_zero_bond.o" ,  "lmm_random_generator.o" , "lmm_products.o" , "lmm_numerical.o" , "lmm_mathtools.o" , "lmm_libor.o" , "cumulfunc.o" ,  "dcdflib.o" , "ipmpar.o" ]

link( list_file ,"lmm_swpt_martX_sci","C")

tenor=[ period ];
numFac=[ nb_fac ];
swpt_mat = [ swpt_maturity ];
swp_mat = [ swp_maturity ];
K=[ strike];

c=fort("lmm_swpt_martX_sci",tenor,6 ,"d",numFac , 4, "i" ,swpt_mat , 2 , "d", swp_mat , 3 , "d", K , 5 , "d" , "out",[1 , 1],1,"d");

resu=c;

endfunction


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////         martingale V
////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function resu=lmm_cap_spotV_sci(period , nb_fac , maturity , strike )

list_file=[ "lmm_martingaleV_sci.o", "lmm_martingaleV.o" , "lmm_volatility.o" , "lmm_zero_bond.o" ,  "lmm_random_generator.o" , "lmm_products.o" , "lmm_numerical.o" , "lmm_mathtools.o" , "lmm_libor.o" , "cumulfunc.o" ,  "dcdflib.o" , "ipmpar.o" ]

link( list_file ,"lmm_cap_spotV_sci","C")

tenor=[ period ];
numFac=[ nb_fac ];
Mat=[ maturity ];
K=[strike];
numcap=round(maturity/period)

c=fort("lmm_cap_spotV_sci",tenor,3,"d",numFac , 2, "i" , Mat , 5 , "d", strike , 4 , "d" , "out",[numcap , 1],1,"d");
m=zeros(numcap,2);
for i=1:numcap,
  m(i,1)=c(i),
  m(i,2)=(i-1)*period,
end;
resu=m;

endfunction

function resu=lmm_swpt_spotV_sci(period , nb_fac , swpt_maturity , swp_maturity , strike )

list_file=[ "lmm_martingaleV_sci.o", "lmm_martingaleV.o" , "lmm_volatility.o" , "lmm_zero_bond.o" ,  "lmm_random_generator.o" , "lmm_products.o" , "lmm_numerical.o" , "lmm_mathtools.o" , "lmm_libor.o" , "cumulfunc.o" ,  "dcdflib.o" , "ipmpar.o" ]

link( list_file ,"lmm_swpt_spotV_sci","C")

tenor=[ period ];
numFac=[ nb_fac ];
swpt_mat = [ swpt_maturity ];
swp_mat = [ swp_maturity ];
K=[strike];

c=fort("lmm_swpt_spotV_sci",tenor,6 ,"d",numFac , 4, "i" ,swpt_mat , 2 , "d", swp_mat , 3 , "d", K , 5 , "d" , "out",[1 , 1],1,"d");

resu=c;

endfunction

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////         Jump
////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
function resu=lmm_swpt_jump_sci(period , nb_fac , swpt_maturity , swp_maturity , strike )

list_file=[ "lmm_jump_sci.o", "lmm_jump.o" , "lmm_volatility.o" , "lmm_zero_bond.o" ,  "lmm_random_generator.o" , "lmm_products.o" , "lmm_numerical.o" , "lmm_mathtools.o" , "lmm_libor.o" , "cumulfunc.o" ,  "dcdflib.o" , "ipmpar.o" ]

link( list_file ,"lmm_swpt_jump_sci","C")

tenor=[ period ];
numFac=[ nb_fac ];
swpt_mat = [ swpt_maturity ];
swp_mat = [ swp_maturity ];
K=[ strike];

c=fort("lmm_swpt_jump_sci",tenor,6 ,"d",numFac , 4, "i" ,swpt_mat , 2 , "d", swp_mat , 3 , "d", K , 5 , "d" , "out",[1 , 1],1,"d");

resu=c;

endfunction

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///  Bermuda swaption Pedersen interface
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function resu=lmm_bermuda_LS_sci(period , nb_fac , swpt_maturity , swp_maturity , strike , payoff_Reg ,Regr_basis_dim)

list_file=[ "lmm_bermuda_LS_sci.o", "lmm_basis.o" , "lmm_bermudaprice.o" , "lmm_volatility.o" , "lmm_zero_bond.o" ,  "lmm_random_generator.o" , "lmm_products.o" , "lmm_numerical.o" , "lmm_mathtools.o" , "lmm_libor.o"  ]

link( list_file ,"lmm_bermuda_LS_sci","C")

tenor=[ period ];
numFac=[ nb_fac ];
swpt_mat = [ swpt_maturity ];
swp_mat = [ swp_maturity ];
K=[ strike];

payoff_as_regressor=[ payoff_Reg ];
Regr_Basis_Dimension=[ Regr_basis_dim ];

c=fort("lmm_bermuda_LS_sci",tenor ,2 ,"d",numFac , 3, "i" ,swpt_mat , 4 , "d", swp_mat , 5 , "d", K , 6 , "d" ,  payoff_as_regressor , 7 , "d"  ,  Regr_Basis_Dimension , 8 , "i" , "out",[1 , 1],1,"d");

resu=c*100*100;

endfunction

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///  Bermuda swaption Andersen interface
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function resu=lmm_bermuda_andersen_sci(period , nb_fac , swpt_maturity , swp_maturity , strike )

list_file=[ "lmm_bermuda_andersen_sci.o", "lmm_basis.o" , "lmm_bermudaprice_andersen.o" , "lmm_volatility.o" , "lmm_zero_bond.o" ,  "lmm_random_generator.o" , "lmm_products.o" , "lmm_numerical.o" , "lmm_mathtools.o" , "lmm_libor.o"  ]

link( list_file ,"lmm_bermuda_andersen_sci","C")

tenor=[ period ];
numFac=[ nb_fac ];
swpt_mat = [ swpt_maturity ];
swp_mat = [ swp_maturity ];
K=[ strike];

c=fort("lmm_bermuda_andersen_sci",tenor ,2 ,"d",numFac , 3, "i" ,swpt_mat , 4 , "d", swp_mat , 5 , "d", K , 6 , "d", "out",[1 , 1],1,"d");

resu=c*100*100;

endfunction

