// Computes the implied volatility for a given surface of prices
// call the function compute_vol (P, spot, interest)
// P is a Premia object (the computation must have already been done).
// spot is the spot price
// 1+interest/100 is the annual interest rate to be coherent with Premia

function p = cdf_nor(x)
  [p, q] = cdfnor("PQ", x, 0, 1);
endfunction

function y = bs (x, r, sigma, k, t)
  d1 = (log(x/k) + (r+sigma^2/2)*t)/(sigma*sqrt(t));
  d2 = d1 - sigma*sqrt(t);
  y = x * cdf_nor(d1) - k*exp(-r*t)*cdf_nor(d2);
endfunction

function y = vega(x, r, sigma, k, t)
  d1 = (log(x/k) + (r+sigma^2/2)*t)/(sigma*sqrt(t));
  d2 = d1 - sigma*sqrt(t);
  y = x * exp(-d1*d1/2)/sqrt(2*%pi)* (-d1/sigma + sqrt(t)) - ...
      k * exp(-r*t -d2*d2/2)/sqrt(2*%pi)* (-d2/sigma - sqrt(t));
endfunction

function [sig] = vol_implicit(x, r, p, k, t)
  sigma(1) = 0.1;
  i =1;
  while %t && i<30
    f = p - bs(x,r,sigma(i), k, t);
    veg = 2 * vega(x,r,sigma(i), k, t);
    if (abs(f * veg) / p < 10^(-12)); then break; end
    sigma(i+1) = sigma(i) + f / veg;
    i = i+1;
  end
  sig = sigma(i);
endfunction

function [vol] = vol_impli(x, r, p_mat, k_mat, t_mat)
  i = 1; 
  for k = k_mat'
    j = 1;
    for t = t_mat'
      vol(i,j)=vol_implicit(x,r,p_mat(i,j), k, t);
      j = j+1;
    end
    i = i+1;
  end
endfunction


function vol = compute_vol (P, spot, interest)
  L = P.get_method_results_iter[];
  vol = vol_impli(spot, log(1+interest/100), L(3), L(1),L(2));
  xset('window', xget('window')+1);
  plot3d( L(1), L(2), vol, theta=5, alpha=89.95);
endfunction
