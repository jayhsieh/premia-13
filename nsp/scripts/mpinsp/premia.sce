[test,ilib] = c_link( "libpremiamodel_Interf");
if ~test  then
  exec('../../loader.sce');
  premia_init();
end
load('test-premia-1.bin');
L=P.get_option_values[];
// Maturity
if exists('Maturity') then L(2)(3)=Maturity;end 
if exists('Strike') then L(3)(3)=Strike;end 
P.set_option_values[L];
P.compute[];
L=P.get_method_results[];
save('/tmp/test-res.bin',L);
