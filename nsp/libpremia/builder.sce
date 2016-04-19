// -*- mode: scilab -*-

ilib_name='libpremiamodel';
files=['premiamod.o'; 'premia_wrap.o'; "premia_vars.o"];
libs = [];
table = [ "premia_create", "int_premiamodel_create";
          "premia_init", "int_premiamodel_init";
          "premia_get_models", "int_premia_get_models";
          "premia_get_family", "int_premia_get_family";
          "premia_get_families", "int_premia_get_families";
          "premia_get_methods", "int_premia_get_methods";
          "premia_get_assets", "int_premia_get_assets",
          "test_unserialize", "int_serial_test_unserialize"];

cflags=["-I../../include/"];
ldflags = [ "-L../../lib/ -lpremia" ];
ilib_build(ilib_name,table,files,libs,cflags=cflags,ldflags=ldflags);
