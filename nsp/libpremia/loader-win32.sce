// generated by builder.sce: Please do not edit this file
// ------------------------------------------------------

libpremia_path=file('join',['.','../../lib/libpremia.dll']);
libpnl_path=file('join',['.','../../lib/libpnl.dll']);
link(libpnl_path);
link(libpremia_path);
libpremiamodel_path=file('join',['.','libpremiamodel.dll']);
addinter(libpremiamodel_path,'libpremiamodel');
