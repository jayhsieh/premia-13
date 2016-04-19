// vim: set tw=90 sw=2 ts=2 sts=2: 

function clear_table (win, typeobj)
  table = win.get_data[typeobj + '_table'];
  table.destroy[];
  table = gtktable_new(rows=1, columns=1, homogeneous=%f);
  select typeobj 
  case "model" then win.set_data[model_table=table];
  case "option" then win.set_data[option_table=table];
  case "method" then win.set_data[method_table=table];
  end
endfunction


//
// data = list (win, combo_model)
// cb = asset_combo
//
function asset_combo_changed(cb,data)
  win=data(1);
  model_cb=data(2);
  option_cb=data(3);
  method_cb=data(4);
  premia_obj = win.get_data['premia_model'];
  ts_model = cb.get_model[]
  iter=cb.get_active_iter[];
  asset_name = ts_model.get_value[iter,0];
  if ~isempty(asset_name)
  then
    premia_obj.set_asset[str=asset_name];
    ts_model = gtktreestore_new(premia_get_models(asset=asset_name));
    model_cb.set_model[model=ts_model];
    model_cb.set_active[-1];
    option_cb.set_active[-1];
    method_cb.set_active[-1];
    clear_table(win, "model");
    clear_table(win, "option");
    clear_table(win, "method");
  end
  win.show_all[];
endfunction

//
// data = list (win, combo_option)
// cb = model_combo
//
function model_combo_changed(cb,data)
  win=data(1);
  opt_cb=data(2);
  method_cb=data(3);
  premia_obj = win.get_data['premia_model'];

  if (cb.get_active[]~=-1  && isempty(premia_obj.get_asset[]))
  then
    throw_message(GTK.MESSAGE_ERROR, "select the asset type first");
    cb.set_active[-1];
    return;
  end

  if (cb.get_active[] == -1)
    I = [];
  else
    iter=cb.get_active_iter[];
    I = cb.get_active[] +1 ;  
  end

  clear_table(win, "model");
  clear_table(win, "option");
  clear_table(win, "method");
  if  ~isempty(I) 
    premia_obj.set_model[I];
    [ts_model,opts,family_index] = create_premia_option_tree(I, premia_obj.get_asset[])
    opt_cb.set_model[model=ts_model];
    [lab, entr, table] = premia_model_values(win,premia_obj);
  end
  opt_cb.set_active[-1];
  method_cb.set_active[-1];
  win.show_all[];
endfunction

//
// data = list (win, combo_method)
// cb = option_combo
//
function opt_combo_changed(cb, data)
  win = data(1);
  meth_cb = data(2);
  premia_obj = win.get_data['premia_model'];
  it = cb.get_active[];
  
  if (it ~= -1 && (isempty(premia_obj.get_asset[]) || isempty(premia_obj.get_model[]))) 
  then
    throw_message(GTK.MESSAGE_ERROR, ...
                  "select the asset type and the model first");
    cb.set_active[-1];
    return;
  end
  
  clear_table(win, "option");
  clear_table(win, "method");
  if it ~= -1 then
    asset_name = premia_obj.get_asset[];
    [option, family] = get_value_from_comb_opt(cb, premia_obj);
    
    if ~isempty(option) then 
      premia_obj.set_option[family, option];
      meth=premia_obj.get_methods[];
      set_combo_method(win, meth_cb, meth);
      [lab, entr, table ] = premia_option_values(win,premia_obj);
    end
  end
  meth_cb.set_active[-1];
  win.show_all[];
endfunction

//
// data = list (win)
// cb = method_combo
//
function meth_combo_changed(cb, data)
  win = data(1);
  
  premia_obj = win.get_data['premia_model'];
  asset_name = premia_obj.get_asset[];
  
  it = cb.get_active[]
  
  if (it~=-1 && (isempty(premia_obj.get_asset[]) || isempty(premia_obj.get_model[]) ...
                 || isempty(premia_obj.get_option[])))
  then
    throw_message(GTK.MESSAGE_ERROR, ...
                  "select the asset type, the model and the option first");
    cb.set_active[-1];
    return;
  end

  if it<>-1 then
    iter=cb.get_active_iter[];
    ts_model = cb.get_model[];
    meth_name = ts_model.get_value[iter,0];
    meth = premia_obj.get_methods[];
    i_meth = find(meth==meth_name);
  else
    i_meth = [];
  end
  clear_table(win, "method");
  if ~isempty(i_meth) then
    premia_obj.set_method[i_meth];
    [lab, entr, table ] = premia_method_values(win, premia_obj);
  end
  win.show_all[];
endfunction


//
// Creates the list of entries and labels appearing in the GUI
//
function premia_get_values (win, M, typeobj)
  select typeobj
  case "model" then premia_model_values (win, M);
  case "option" then premia_option_values (win, M);
  case "method" then premia_method_values (win, M);
  end
endfunction

function [L] = get_object_values(M,typeobj)
  select typeobj
  case "model" then L=M.get_model_values[];
  case "option" then L=M.get_option_values[];
  case "method" then L=M.get_method_values[];
  end
endfunction

function set_object_values(M,typeobj,data)
  select typeobj
  case "model" then M.set_model_values[data];
  case "option" then M.set_option_values[data];
  case "method" then M.set_method_values[data];
  end
endfunction

//
// Goes through the list data to search for the key label and
// set the corresponding value to val
//
//
function [L] = tweak_object_values(data,label,val)
  L = list();
  // Beware that label may have been wrapped, so we need to ensure that it is a single line
  label1 = strsubst(label, "\n", "");
  for it = data
    Li = it;
    select it(1)
    case "combo" then
      if it(2) == label1; then
        Li(3) = it(6)(val+1)
      else
        subList = list ();
        for l = it(7)
          subList.add_last[tweak_object_values(l,label1,val)];
        end
       Li(7) = subList; 
      end
    end
    L.add_last[Li];
  end
endfunction

//
// data = list(win,M,typeobj,label)
// cb = an enumeration combo
//
function enum_combo_changed (cb, data)
  win = data(1);
  M = data(2);
  typeobj = data(3);
  label = data(4)
  it = cb.get_active[]
  if it == -1; then return; end

  // get the new value
  iter=cb.get_active_iter[];
  ts_model = cb.get_model[];
  path = ts_model.get_path[iter];
  choice = evstr(path.to_string[]);

  // copy the old value of cb back to store all entries
  L = get_object_values(M,typeobj);
  Lc = premia_var_to_xchoices(L);
  [oldlabels,oldentries] = recursive_table(win,M,Lc,typeobj);
  entries = win.get_data[typeobj + "_entries"];

  // sets the old value of the combo back
  i = 1;
  for l=oldlabels; if l.get_text[] == label.get_text[]; then break; else i=i+1; end; end
  entries(i) = oldentries(i);
  data = get_all_data(Lc, entries);
  set_object_values(M,typeobj,data);
  L = tweak_object_values (data, label.get_text[], choice);
  set_object_values(M,typeobj,L);
  premia_get_values(win,M,typeobj);
  win.show_all[];

endfunction

