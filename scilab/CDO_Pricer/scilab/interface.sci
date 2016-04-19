try 
    type(cdo_factor);
    delmenu('CDO');
catch
//    exec builder.sce;
    exec loader.sce; 
end;
addmenu('CDO',['Product';'Model';'Computation';'Last computation';'Numdef';'Losses']);
CDO=['exec(''product.sci'',-1)';'exec(''model.sci'',-1)';'exec(''computation.sci'',-1)';'exec(''last_computation.sci'',-1)';'exec(''animate_numdef'',-1)';'exec(''animate_losses.sci'',-1)'];
unsetmenu('CDO',5);
unsetmenu('CDO',6);

exec('product.sci',-1);
exec('model.sci',-1);
exec('computation.sci',-1);
