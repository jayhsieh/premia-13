exec ("../libpremia/loader.sce");
premia_init ();

// Compares the numerical entries of the two lists P1 and P2.  Returns
// a matrix containing the difference between the numerical entries of
// P1 and P2
//
// type(L)=='Mat' does not actually do the comparison you think,
// instead use is(L, %types.Mat)

function res=compare_eps (P1,P2)
  res = [];  
  if is(P1, %types.List)
  then       
    for i=(1:size(P1))
      res = [res ; compare_eps (P1(i),P2(i)) ];
    end
  elseif is(P1, %types.Mat)
  then 
    res = [res ; max (abs(P1-P2)) ];
  end
endfunction



function [res] = premia_diff (P)
  Lold = P.get_method_results[]
  P.compute[]; 
  Lnew = P.get_method_results[];
  res = compare_eps(Lold, Lnew);
endfunction
