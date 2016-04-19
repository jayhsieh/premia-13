epsilon_comp = 1e-9;

// Compares the numerical entries of the two lists P1 and P2. Returns
// status=0 if the numerical entries are equal up to epsilon, returns
// a strictly positive value otherwise. res is a matrix containing the
// diffrences between P1 and P2
//
// type(L)=='Mat' does not actually do the comparison you think,
// instead use is(L, %types.Mat)

function [status, res] =compare_eps_diff (P1,P2)
  res = [];
  status = 0;
  if is(P1, %types.List)
  then		
	for i=(1:size(P1))
      [x, y] = compare_eps_diff (P1(i),P2(i));
      status = status + x;
      res = [ res ; y];
	end
  elseif is(P1, %types.Mat)
  then
    if norm(P1-P2,'inf')>epsilon_comp
    then
      status = 1;
      res = [res ; max (abs(P1-P2)) ];
    end
  end
endfunction

function [status] =compare_eps (P1,P2, out, tst)
  [status, res] =compare_eps_diff (P1,P2);
  if status ~= 0
  then
    F = fopen (out, mode ="a");
    F.putstr[tst + "\n"]
    F.put_matrix[res, format='%.12f'];
    F.putstr["\n"];
    F.close[]
  end
endfunction
  
