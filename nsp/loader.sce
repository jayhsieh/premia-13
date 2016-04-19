
// test if lib already loaded, unload first
[test,ilib] = c_link( "libpremiamodel_Interf");
if test  then
  printf("   library alreday loaded : nothing done\n");
else
  exec('libpremia/loader.sce')
end
