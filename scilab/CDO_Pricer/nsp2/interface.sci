try 
    type(cdo_factor);
    delmenu('CDO');
catch
    exec builder.sce;
    exec loader.sce; 
end;
addmenu('CDO',['Product';'Model';'Computation';'Last computation';'Numdef';'Losses']);
CDO=['exec(''product.sci'')';'exec(''model.sci'')';'exec(''computation.sci'')';'exec(''last_computation.sci'')';'exec(''animate_numdef'')';'exec(''animate_losses.sci'')'];
unsetmenu('CDO',5);
unsetmenu('CDO',6);

exec product.sci
exec model.sci
exec computation.sci
