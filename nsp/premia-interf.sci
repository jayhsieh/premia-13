// vim: set tw=90 sw=2 ts=2 sts=2: 
// -*- fill-column: 90 -*-

//
// Make sure v is a column vector
//
function r = vect_to_col(v)
  [m,n] = size(v);
  r = v
  if m==1 && n>1 then r = v'; end
endfunction


//
// Make sure v is a row vector
//
function r = vect_to_row(v)
  [m,n] = size(v);
  r = v
  if n==1 && m>1 then r = v'; end
endfunction

//
// wrap too long labels
//
function w = wrap(str)
  if length(str) < 20 then w = str; return; end

  B=strindex (str, ' ');
  if ~isempty(B) then
    i = max(B(B < 20));
    if isempty(i) then
      i = min(B(B >= 20));
    end
    str1 = part(str, (1:i))
    str2 =  part(str, (i+1: length(str)))
    w = [ str1; "\n"; wrap(str2)];
  else
    w = str;
  end
endfunction

//
// Creates a list suitable for building a table of Gtk widgets
//
//  L: is a list as returned by get_{model,option,method}_values[]
//  
// Returns a list whose entries are of the type
//      list (type, name, value, [choice_labels, ids])
//          - type: is a string 'combo', 'selectfile', 'entry' or 'ignore'
//          - name: a label to draw next to the Gtk widget
//          - value: value of the parameter (my be of different type)
//          - choice_labels (only when type=='combo', ie. enumaration case): a vector of the different
//            labels for the enumeration
//          - ids (only when type=='combo', ie. enumaration case): a vector of the different
//            keys for the enumeration
//          - a list of lists of the current type (only when type='combo').
//            Each entry of the list is created by calling premia_var_to_xchoices on the
//            list of parameters for the corresponding choice in the combo
//          Note that the key for choice_labels(i) is ids(i) 
//
function Lc=premia_var_to_xchoices(L)
  // convert a premia var list to a format suitable for x_choices 
  Lc=list();
  for i=1:length(L)
    select L(i)(1)
    case "ENUM" then
      if L(i)(4)== %t then  tag = 'combo'; else  tag='ignore';  end
      str=sprint(L(i)(3),as_read=%t);
      str=str(2); // first line is a header 
      subLc = list ();
      for it=L(i)(7)
        subLc.add_last[premia_var_to_xchoices(it)]
      end
      Lc.add_last[list(tag,L(i)(2),str, L(i)(5), L(i)(6),subLc)];
    case "FILENAME" then 
      if L(i)(4)== %t then  tag = 'selectfile'; else  tag='ignore';  end
      Lc.add_last[list(tag,L(i)(2),L(i)(3))];
    else
      if L(i)(4)== %t then tag = 'entry' ;else tag='ignore'; end
      str=sprint(vect_to_row(L(i)(3)),as_read=%t);
      str=str(2:size(str,'*'));  // first line is a header 
      str=strsubst(str, '...', ''); // only one line matrices are accepted
      str=catenate(str); // no ... at end of line.
      Lc.add_last[list(tag,L(i)(2),str)];
    end
  end
endfunction

function [ts_model,options,family_index]= create_premia_option_tree(pmodel,asset_name)
  // select options which are compatible with given
  // pmodel. Options are grouped by families. 
  n=1;
  while ~isempty(premia_get_family(n, asset = asset_name)) then n=n+1;end 
  // n-1 gives the total number of families 
  // extract the sub list for pmodel 
  options = list();
  Top = [];
  Top_names = [];
  for i=1:n-1 
    S=premia_get_family(i,pmodel, asset=asset_name);
    if ~isempty(S) then 
      options($+1)=S;
      Top = [Top;i];
      Top_names= [Top_names; premia_get_families(i,asset=asset_name )];
    end 
  end
  family_index=Top;
  ts_model = gtktreestore_new(Top_names); 
  // fill the ts_model at next level 
  it = ts_model.get_iter_first[] 
  i=1;

  while %t 
    ts_model.insert[it,0,list(options(i))];
    i = i+1;
    if ts_model.iter_next[it] == %f 
      break;
    end
  end
endfunction

function [option,  family] = get_value_from_comb_opt(cb, premia_obj)
  ts_model = cb.get_model[];

  if (cb.get_active[] == -1) then
    option = []; family = [];
    return
  end
  iter=cb.get_active_iter[];
  path = ts_model.get_path[iter];
  i = path.get_indices[];
  family=i(1)+1; option=i(2)+1;
endfunction

function set_combo_method(win, cb, meth)
  ts =  gtktreestore_new(meth);
  cb.set_model[model=ts];
  iter=ts.get_iter[[0]];
  cb.set_active_iter[iter];
  win.show_all[];
endfunction

function [table] = premia_object_values(table,labels,entries,vbox)
  n = length(labels);
  table.destroy[];
  table = gtktable_new(rows=n, columns=2, homogeneous=%f);
  for i=1:n
    table.attach_defaults[ labels(i),  0, 1, i-1, i]
    table.attach_defaults[ entries(i), 1, 2, i-1, i];
  end
  vbox.pack_start[table, expand=%f, fill=%t, padding=0];
endfunction

function  [labels, entries, table] = premia_model_values(win,M)
  B=M.get_model_values[];
  Bx=premia_var_to_xchoices(B);
  table = win.get_data['model_table'];
  vbox_m=win.get_data['vbox_model'];
  [labels, entries] = recursive_table(win,M,Bx,"model");
  [table] = premia_object_values(table,labels,entries,vbox_m);
  win.set_data[model_labels=labels];
  win.set_data[model_entries=entries];
  win.set_data[model_table=table];
endfunction

function  [labels, entries, table ] = premia_option_values(win,M)
  B=M.get_option_values[];
  Bx=premia_var_to_xchoices(B);
  table= win.get_data['option_table'];
  vbox_opt=win.get_data['vbox_option'];
  [labels, entries] = recursive_table(win,M,Bx,"option");
  [table] = premia_object_values(table,labels,entries,vbox_opt);
  win.set_data[option_labels=labels];
  win.set_data[option_entries=entries];
  win.set_data[option_table=table];
endfunction

function  [labels, entries, table] = premia_method_values(win,M)
  table= win.get_data['method_table'];
  vbox_meth=win.get_data['vbox_method'];
  B=M.get_method_values[];
  Bx=premia_var_to_xchoices(B);
  [labels, entries] = recursive_table(win,M,Bx,"method");
  [table] = premia_object_values(table,labels,entries,vbox_meth);
  win.set_data[method_labels=labels];
  win.set_data[method_entries=entries];
  win.set_data[method_table=table];
endfunction

// 
// win is the main window object
// L is a list as returned by premia_var_to_xchoices
// M is the Premia object
// typeobj is a string "model", "option", "method"
//
// From this list, a new gtktable is built to reflect the list of setable parameters
//
//  - labels is a list of labels for the different parameters. Each parameter is associated to a
//    gtkwidget
//  - entries is a list of gtkwidgets. Each widget is associated to the corresponding value of
//    labels
//
function [labels, entries] = recursive_table(win,M,L,typeobj)

  labels=list(); entries=list();

  for it = L
    select it(1) 
    case "entry" then
      labels_i = gtklabel_new(str=it(2));
      labels_i.set_line_wrap[%t];
      labels_i.set_text[catenate(wrap(it(2)))];
      entries_i = gtkentry_new ();
      entries_i.set_text[it(3)];
      labels.add_last[labels_i];
      entries.add_last[entries_i];
    case "combo" then
      cell_renderer = gtkcellrenderertext_new ();
      labels_i = gtklabel_new(str=it(2));
      entries_i = gtkcombobox_new();
      labels_i.set_text[catenate(wrap(it(2)))];
      ts_gener = gtktreestore_new(it(4));
      entries_i.set_model[model=ts_gener];
      entries_i.pack_start[ cell_renderer,expand= %t];
      entries_i.add_attribute[cell_renderer,"text",0];
      index = find(it(5)==evstr(it(3)));
      iter=ts_gener.get_iter[index-1];
      entries_i.set_active_iter[iter];
      entries_i.connect["changed", enum_combo_changed, list(win,M,typeobj,labels_i)]
      labels.add_last[labels_i];
      entries.add_last[entries_i];
      if length(it(6)(index)) > 0; then
        [sublabels, subentries] = recursive_table(win,M,it(6)(index),typeobj);
        labels.concat[sublabels];
        entries.concat[subentries];
      end
    case "selectfile" then
      labels_i = gtklabel_new(str=it(2));
      entries_i = gtkfilechooserbutton_new['Select File', 'open'];
      entries_i.set_filename[it(3)];
      labels_i.set_text[catenate(wrap(it(2)))];
      labels.add_last[labels_i];
      entries.add_last[entries_i];
    else
    end
  end
endfunction


//
// L is a list as returned by premia_var_to_xchoices
//
// data is a list suitable for being given to .set_{model,option,method}
//
function [data] = copy_all_data (L)
  data=list();

  for it = L
    data_i = list();
    data_i = it;
    select it(1)
    case "ignore" then
      data_i(4) = %f;
    case "combo" then
      data_i(4) = %t;
      data_i(5) = it(4);
      data_i(6) = it(5);
      subList = list ();
      for l = it(6)
        subList.add_last[copy_all_data(l)];
      end
      data_i(7) = subList;
    case "entry" then
      data_i(3) = evstr(it(3));
    else
      data_i(4) = %t
    end
    data.add_last[data_i];
  end
endfunction

function [data] = get_all_data(L, entries)
  data=list();
  i=1;
  for it = L
    data_i = list();
    select it(1)
    case "entry" then
      data_i(1) = it(1);
      data_i(2) = it(2)
      data_i(3) = evstr(entries(i).get_text[]);
      data_i(4) = %t;
      i = i+1;
    case "ignore" then
      data_i(1) = it(1);
      data_i(2) = it(2);
      data_i(3) = it(3);
      data_i(4) = %f;
    case "combo" then
      data_i(1) = it(1);
      data_i(2) = it(2);
      Labels = it(4);
      ts = entries(i).get_model[]
      iter=entries(i).get_active_iter[];
      name = ts.get_value[iter,0];
      choice = (find(Labels == name));
      data_i(3) = it(5)(choice);
      data_i(4) = %t;
      data_i(5) = it(4);
      data_i(6) = it(5);
      subList = list();
      for l = it(6)
        subList.add_last[copy_all_data(l)];
      end
      data_i(7) = subList;
      data_i(7)(choice) = get_all_data(it(6)(choice), list(entries(i+1:$)));
      i = i + length(data_i(7)(choice)) + 1;
    case "selectfile" then
      data_i(1) = it(1);
      data_i(2) = it(2);
      data_i(3) = entries(i).get_filename[];
      data_i(4) = %t;
      i = i+1;
    end
    data.add_last[data_i];
  end
endfunction

function ok=premia_check_basic(L,Lold)
  // check that values in L matches types in Lold 
  // we assume here that L and Lold have the same shapes 
  ok=%t;
  for i=1:length(L)
    if L(i)(4)== %t then 
      // we have to check here that values match
      if L(i)(1) == "entry" &  type(L(i)(3),'short') <> 'm' then 
        x_message(L(i)(2)+' should be of Matrix type');
        ok=%f;
        return;
      end
    end
  end
endfunction

// gets the data from win and sets the corresponding parameters of the PremiaObject.
function set_premia_obj_values (win, M)
  // model
  B=M.get_model_values[];
  L=premia_var_to_xchoices(B);
  entries = win.get_data['model_entries'];
  data = get_all_data(L, entries);

  if premia_check_basic(data,B) then 
    M.set_model_values[data];
    // Now we have to check if the model accepts the parameters 
    S=M.model_check[];
    if ~isempty(S) then 
      M.set_model_values[B];
      x_message(S(1)+': '+S(2));
      return
    end
  end

  // option
  B=M.get_option_values[];
  L=premia_var_to_xchoices(B);
  entries = win.get_data['option_entries'];
  data = get_all_data(L, entries)
  if premia_check_basic(data,B) then 
    M.set_option_values[data];
    // Now we have to check if the option accepts the parameters 
    S=M.option_check[];
    if ~isempty(S) then 
      M.set_option_values[B];
      x_message(S(1)+': '+S(2));
      return
    end
  end

  // method
  B=M.get_method_values[];
  L=premia_var_to_xchoices(B);
  entries = win.get_data['method_entries'];
  data = get_all_data(L, entries);
  if premia_check_basic(data,B) then 
    // Now we have to check if the method  accepts the parameters 
    M.set_method_values[data];
    S=M.method_check[];
    if ~isempty(S) then 
      M.set_method_values[B];
      x_message(S(1)+': '+S(2));
      return
    end
  end
endfunction

function compute_results(win)
  M = win.get_data['premia_model'];
  if isempty(M.get_asset[]) then
    throw_message(GTK.MESSAGE_ERROR, "Asset is not set");
    return;
  end

  if isempty(M.get_model[]) then
    throw_message(GTK.MESSAGE_ERROR, "Model is not set");
    return;
  end

  if isempty(M.get_option[]) then
    throw_message(GTK.MESSAGE_ERROR, "Option is not set");
    return;
  end

  if isempty(M.get_method[]) then
    throw_message(GTK.MESSAGE_ERROR, "Method is not set");
    return;
  end

  set_premia_obj_values (win, M);
  // compute price and delta and show them
  M.compute[];
  compute_err = M.get_compute_err[];
  if (~ isempty(compute_err)) then
    throw_message(GTK.MESSAGE_ERROR, compute_err);
    return;
  end
endfunction


