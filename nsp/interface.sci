// vim: set tw=90 sw=2 ts=2 sts=2: 
// -*- fill-column: 90 -*-

exec('premia-interf.sci');
exec('callbacks.sci');

function y=macro_set_capital_sensitive(cell_view,cell,model,iter,data)
// unset senstitivity to the X -- Y fields 
  sensitive = model.iter_has_child[iter];
  cell.set_property['sensitive',~sensitive]
  y=%t;
endfunction

//
// Create a vbox either for asset, model, option or method.
// The vbox holds a combo to choose in the list and the table
// to be used for setting parameters
// 
//  name is the label of the combo
//  ts is a gtk treestore used to build the combo
//
function [vbox,combo,table] = create_vbox (name, ts)
  vbox = gtkvbox_new(homogeneous=%f, spacing=8);
  label = gtklabel_new(str=name);
  vbox.pack_start[label, expand=%f, fill=%f, padding=0];
  combo = gtkcombobox_new();
  combo.set_model[model=ts];
  combo.set_add_tearoffs[%t];
  cell_renderer = gtkcellrenderertext_new ();
  combo.pack_start[ cell_renderer,expand= %t];
  combo.add_attribute[cell_renderer,"text",0];
  combo.set_active[-1];
  vbox.pack_start[combo, expand=%f, fill=%t, padding=0];
  table=gtktable_new(rows=1, columns=1, homogeneous=%f);
  vbox.pack_start[table, expand=%t, fill=%t, padding=0];
endfunction


function [M] = premia()

  M=premia_create();

  title='Premia';
  description='An option pricer';
  assets = premia_get_assets();

  models=premia_get_models();

  // we start with bs1d CallEuro and BS formula
  [ts_model,opts,family_index] = create_premia_option_tree(1, assets(1))

  // main window
  win = gtkwindow_new();
  win.set_title[title];
  vbox_main = gtkvbox_new(homogeneous=%f,spacing=0);
  win.add[vbox_main];
  
  // add properties to the main window
  win.set_data[premia_model=M];
  win.set_data[model_entries=list()];
  win.set_data[option_entries=list()];
  win.set_data[method_entries=list()];
  win.set_data[model_labels=list()];
  win.set_data[option_labels=list()];
  win.set_data[method_labels=list()];
  // properties for the graphs 
  win.set_data[color_style=1];
  win.set_data[graph_x=[]];
  win.set_data[graph_y=[]];
  win.set_data[leg=''];
  
  // logo
  stock = gtkimage_new("file","premia.gif");
  //help menu
  menubar = gtkmenubar_new();
  menu_file = gtkmenuitem_new("File");
  menu = gtkmenuitem_new("Help");
  menu_tools = gtkmenuitem_new("Tools");
  about = gtkmenuitem_new("About");
  about.set_right_justified[%t];
  // sub menu
  
  sub_menu_file = gtkmenu_new ();
  load_file = gtkmenuitem_new(label='Load');
  save_file = gtkmenuitem_new(label='Save');
  quit_file = gtkmenuitem_new(label='Quit');
  sub_menu_file.append[load_file];
  sub_menu_file.append[save_file];
  sub_menu_file.append[quit_file];
  menu_file.set_submenu[sub_menu_file];
  
  sub_menu_help = gtkmenu_new();
  model_help = gtkmenuitem_new(label='Model Help');
  option_help = gtkmenuitem_new(label='Option Help');
  method_help = gtkmenuitem_new(label='Method Help');
  sub_menu_help.append[model_help];
  sub_menu_help.append[option_help];
  sub_menu_help.append[method_help];
  menu.set_submenu[sub_menu_help];
  
  sub_menu_tools = gtkmenu_new();
  draw_graph = gtkmenuitem_new(label='Draw');
  redraw_graph = gtkmenuitem_new(label='Redraw');
  legend_graph = gtkmenuitem_new(label='Legend');
  clear_graph = gtkmenuitem_new(label='Clear Graph');
  clear_data = gtkmenuitem_new(label='Clear Data');
  sub_menu_tools.append[draw_graph];
  sub_menu_tools.append[redraw_graph];
  sub_menu_tools.append[legend_graph];  
  sub_menu_tools.append[clear_graph];
  sub_menu_tools.append[clear_data];
  menu_tools.set_submenu[sub_menu_tools];
  
  sub_menu_about = gtkmenu_new();
  premia_about = gtkradiomenuitem_new(label='Premia');
  sub_menu_about.append[premia_about];
  about.set_submenu[sub_menu_about];
  menubar.append[menu_file];
  menubar.append[menu];
  menubar.append[menu_tools];
  menubar.append[about];
  vbox_main.pack_start[ menubar,expand=%f,fill=%t,padding=0];

  // 2 tabs
  notebook = gtknotebook_new ();
  
  // first we put the interface in a scrollable window
  scrolled_window = gtkscrolledwindow_new();
  scrolled_window.set_policy[ GTK.POLICY_AUTOMATIC, GTK.POLICY_AUTOMATIC]
  scrolled_window.set_shadow_type[ GTK.SHADOW_IN]
  notebook.append_page[scrolled_window, gtklabel_new(mnemonic="_Interface")]; 
  hb_inter = gtkvbox_new();
  scrolled_window.add_with_viewport[hb_inter];
  

  vbox_main.pack_start[ notebook,expand=%t,fill=%t,padding=0]
  hb_graph = gtkhbox_new();
  notebook.append_page[hb_graph, gtklabel_new(mnemonic="_Graphics")]; 

  // be sure to pack hb_graph before creating the graph window !
  // don't use opengl=%t otherwise closing the main window makes Nsp
  // crash
  winid=nsp_graphic_new(win,hb_graph,opengl=%f);
  // we must keep track of the winid somewhere 
  // just in case other graphics windows are also activated.
  win.set_data[graphid=winid];
  
  // hbox for combos to define the pb 
  hbox1 = gtkhbox_new(homogeneous=%f,spacing=8);
  hbox1.set_border_width[8]
  // hbox for the textview
  hbox_textview = gtkhbox_new(homogeneous=%f,spacing=8);
  hbox_textview.set_border_width[8]
  // hbox for the buttons "compute" and "close"
  hbox_but = gtkhbox_new(homogeneous=%f,spacing=8);
  hbox_but.set_border_width[8]
  
  hb_inter.pack_start[ hbox1,expand=%f,fill=%f,padding=0];
  hb_inter.pack_end[ hbox_textview,expand=%t,fill=%t,padding=0];
  hb_inter.pack_end[ hbox_but,expand=%f,fill=%t,padding=0];

  // asset type
  ts_assets = gtktreestore_new(assets);
  [vbox_asset,asset_combo,asset_table] = create_vbox ('Asset type', ts_assets);
  hbox1.pack_start[vbox_asset, expand=%f, fill=%t, padding=0];
  
  // model choice
  ts_models = gtktreestore_new(models); 
  [vbox_model,model_combo,model_table] = create_vbox ('Model', ts_assets);
  hbox1.pack_start[vbox_model, expand=%f, fill=%t, padding=0];
  win.set_data[model_table=model_table];

  // family & option choice
  [vbox_option,option_combo,option_table] = create_vbox ('Option', ts_assets);
  hbox1.pack_start[vbox_option, expand=%f, fill=%t, padding=0];
  win.set_data[option_table=option_table];

  // method choice
  methods=premia_get_methods(1,1,1);
  ts_meth=gtktreestore_new(methods); 
  [vbox_method,method_combo,method_table] = create_vbox ('Method', ts_assets);
  hbox1.pack_start[vbox_method, expand=%f, fill=%t, padding=0];
  win.set_data[method_table=method_table];

  // changed methods
  asset_combo.connect["changed", asset_combo_changed, list(win, model_combo, option_combo, method_combo)]
  model_combo.connect["changed", model_combo_changed, list(win, option_combo, method_combo)]
  option_combo.connect["changed", opt_combo_changed, list(win, method_combo)]
  method_combo.connect["changed", meth_combo_changed, list(win)]


  // vboxes need to be updated on interactive choices.  
  // need to link them to the main window to have an asynchronous interface.
  win.set_data[vbox_asset= vbox_asset];
  win.set_data[vbox_model=vbox_model];
  win.set_data[vbox_option=vbox_option];
  win.set_data[vbox_method=vbox_method];
  

  // hanbdler for the menu
  model_help.connect["activate", help_activated, list(win, M, 'model')]
  option_help.connect["activate", help_activated, list(win, M, 'option')]
  method_help.connect["activate", help_activated, list(win, M, 'method')]
  premia_about.connect["activate", raise_logo,list(win)];
  draw_graph.connect["activate", draw_graph_callback,list(win, M, notebook)];
  redraw_graph.connect["activate", redraw_graph_callback,list(win, M)];
  clear_graph.connect["activate", clear_graph_callback,list(win, M)];
  clear_data.connect["activate", clear_data_callback,list(win, M)];
  legend_graph.connect["activate", legend_graph_callback,list(win, M)];
  load_file.connect["activate", load_premia_obj_from_file,list(win, asset_combo, model_combo, option_combo, method_combo, M)];
  save_file.connect["activate", save_premia_obj_to_file, list(win, M)];
  quit_file.connect["activate", close_win, list(win)];

  // display output
  view = gtktextview_new[]
  sw = gtkscrolledwindow_new();
  sw.set_policy[ GTK.POLICY_AUTOMATIC,GTK.POLICY_AUTOMATIC];
  sw.add[view];

  hbox_textview.pack_start[sw, expand=%t, fill=%t]
  buffer = view.get_buffer[];
  win.set_data[buffer=buffer];
  
  
  // buttons
  compute_but=gtkbutton_new[label='Compute'];
  close_but=gtkbutton_new[label='Close'];
  hbox_but.pack_end[close_but, expand=%f, fill=%f, padding=8];
  hbox_but.pack_end[compute_but, expand=%f, fill=%f, padding=8];
  close_but.connect["clicked", close_win, list(win)];
  compute_but.connect["clicked", compute_but_callback, list(win)];
  
  win.set_default_size[  600, 400]
  win.show_all[];

endfunction 

// handler for the "close" button
function close_win(but, data)
  win = data(1);
  win.destroy[];
endfunction

// handler for the "compute" button
function compute_but_callback(but, data)
  win=data(1);
  M=win.get_data['premia_model'];
  buffer = win.get_data['buffer'];
  compute_results(win);
  iter = buffer.get_iter_at_offset [0];
  buffer.insert[iter, M.get_model[] + " / " + M.get_option[] + " / " + ...
                M.get_method[] + "\n"];
  if (M.is_with_iter[] == %t) then
    [L, Vnames] = M.get_method_results_iter[];
    if length(L) == 2 then
      x = L(1);
      y = L(2);
      x_old = win.get_data['graph_x'];
      y_old = win.get_data['graph_y'];
      x_more = vect_to_col(x);
      y_more = vect_to_col(y);
      if ((~isempty(x) && size(x_old)(1) == size(x_more)(1)) ... 
          && (~isempty(y) && size(y_old)(1) == size (y_more)(1))) then
        win.set_data[graph_x=[x_old x_more]];
        win.set_data[graph_y=[y_old y_more]];
      else
        win.set_data[graph_x=x_more];
        win.set_data[graph_y=y_more];
        // clear graph
        xbasc(win.get_data['graphid']);
        win.set_data[color_style = 1];
      end
      res = [ vect_to_col (L(1))'; vect_to_col(L(2))' ];
      res = m2s([[0.; 0.] res]);
      res(1,1) = Vnames(1); res(2,1) = '';
      indent = length(Vnames(1));
    else
      res = [ vect_to_col(L(1)) L(3) ]; res = [ [0.0, vect_to_col(L(2))'] ;res]
      res = m2s(res); res(1,1) = Vnames(1) + "/" + Vnames(2);
      indent = max(length(res(1,1)) - length(res(2,1)), 0.)
    end
    space = ' ';
    tab = '\n' + strcat(space(ones(1,indent+2)))
    buffer.insert[iter, catenate(res, col=" ", row=tab) + "\n"]
  else
    L=M.get_method_results[]
    for i=1:length(L)
      it = L(i);
      len = length(it(3));
      buffer.insert[iter, string(it(2)) + ' ']
      buffer.insert[iter, catenate(string(it(3)), sep=" ") + "\n"]
    end
  end
  buffer.insert[iter, "\n\n"];
endfunction

// handler to draw the graph. Compute must have been run before
function draw_graph_callback(but, data)
  win=data(1);
  M = data(2);
  notebook = data(3);
  color_style = win.get_data['color_style'];
  L  = M.get_method_results_iter[];
  wincur=xget('window');
  xset('window',win.get_data['graphid']);
  if  length(L) == 2 then
    plot2d(L(1),L(2), style=color_style);
    win.set_data[color_style = color_style+1];
  else
    plot3d(L(1), L(2), L(3));
  end
  xset('window',wincur);
  notebook.next_page[];
endfunction

// handler to redraw a graph
function redraw_graph_callback(but, data)
  win=data(1);
  M = data(2);
  x = win.get_data['graph_x'];
  y = win.get_data['graph_y'];
  wincur=xget('window');
  xset('window',win.get_data['graphid']);
  xbasc();
  plot2d(x,y);
  xset('window',wincur);
endfunction

// handler to clear the graph
function clear_graph_callback(but, data)
  win=data(1);
  M = data(2);
  xbasc(win.get_data['graphid']);
  win.set_data[color_style = 1];
endfunction

// handler to clear the data
function clear_data_callback(but, data)
  win=data(1);
  M = data(2);
  win.set_data[color_style = 1];
  win.set_data[graph_x=[]];
  win.set_data[graph_y=[]];
endfunction

// handler to add a legend. Actually redraws everything
function legend_graph_callback(but, data)
  win=data(1);
  M = data(2);
  x = win.get_data['graph_x'];
  y = win.get_data['graph_y'];
  pos_list = ['ul', 'ur', 'dl', 'dr'];
  leg = win.get_data['leg'];
  l1 = list('entry', 'legend', 1, leg);
  l2 = list('combo', 'legend position', 1, pos_list);
  l3 = list('entry', 'style (optional)', 1, '');
  l4 = list('entry', 'rect (optional)', 1, '');
  [Lres, L] = x_choices('Legend options', list(l1,l2,l3,l4), %t);
  if isempty (Lres) then return; end
  leg = Lres(1);
  pos = pos_list(Lres(2));
  style = evstr(Lres(3));
  rect = evstr(Lres(4));  
  wincur=xget('window');
  win.set_data[leg=leg];
  xset('window',win.get_data['graphid']);
  if isempty(style) && isempty(rect)
  then
    plot2d(x,y, leg=leg, leg_pos = pos);
  end
  if ~isempty(style) && isempty(rect)
  then
    plot2d(x,y, style=style, leg=leg, leg_pos=pos);
  end
  if isempty(style) && ~isempty(rect)
  then
    plot2d(x,y, rect=rect, leg=leg, leg_pos=pos);
  end
  if ~isempty(style) && ~isempty(rect)
  then
    plot2d(x,y, style=style, rect=rect, leg=leg, leg_pos=pos);
  end
  xset('window',wincur);
endfunction

// loads a PremiaObject into the interface from a file.
function load_premia_obj_from_file (but, data)
  win = data(1);
  asset_cb = data(2);
  model_cb = data(3);
  opt_cb = data(4);
  meth_cb = data (5);
  M = data(6);
  
  fic = xgetfile();
  if fic == "" then return; end
  load(fic);
  // after the load, the newly loaded object is called P
  asset = P.get_asset[];
  model = P.get_model[];
  opt = P.get_option[];
  method = P.get_method[];
  
  // set the asset
  assets = premia_get_assets();
  i_asset = find( assets == asset);
  if isempty (i_asset) then return; end
  asset_cb.set_active[i_asset-1];
   
  // set the model
  models = premia_get_models(asset=asset);
  i_model = find (models == model);
  if isempty (i_model) then return; end
  model_cb.set_active[i_model-1];  
  
  [lab, entr, table] = premia_model_values(win,P);
  win.set_data[model_labels=lab];
  win.set_data[model_entries=entr];
  win.set_data[model_table=table];
  
  //set the option. Need to construct the correct path through the
  //opt_combo
  n=1; 
  while ~isempty(premia_get_family(n, asset = asset)) then n=n+1; end 
  // n-1 gives the total number of families 
  // extract the sub list for model 

  // i_family is the index of the current family used to build the
  // path. Empty families are not counted. This index refers to the
  // accepted families given the chosen model.
  found = %f; i_local_family=1; i_family=1;
  while ~found & i_family<n
  then
    S=premia_get_family(i_family,i_model, asset=asset);
    if ~isempty(S) then 
      i_opt = find (S == opt);
      if ~isempty (i_opt); 
        ts = opt_cb.get_model[];
        path = string(i_local_family-1) + ":" + string(i_opt-1);
        it=ts.get_iter_from_string[path];
        opt_cb.set_active_iter[it];
        found = %t;
        break
      end
        i_local_family = i_local_family+1;
    end
    i_family=i_family+1;
  end
  if ~found then throw_message(GTK.MESSAGE_ERROR,"Cannot find option type"); end

  [lab, entr, table] = premia_option_values(win,P);
  win.set_data[option_labels=lab];
  win.set_data[option_entries=entr];
  win.set_data[option_table=table];
  
  // set the method
  methods = P.get_methods[];
  i_meth = find (methods == method);
  if isempty (i_meth) then throw_message(GTK.MESSAGE_ERROR,"Cannot find method"); end
  meth_cb.set_active[i_meth-1];
  [lab, entr, table] = premia_method_values(win,P);
  win.set_data[method_labels=lab];
  win.set_data[method_entries=entr];
  win.set_data[method_table=table];
  win.show_all[]
endfunction

// save a premia object defined from the interface to a file using the
// standard save method for PremiaObject
function save_premia_obj_to_file (but, data)
  win = data(1);
  M = data(2);
  if (isempty(M.get_asset[]) || isempty(M.get_model[]) ...
      || isempty(M.get_option[]) || isempty(M.get_method[]))
    throw_message(GTK.MESSAGE_ERROR, ...
                  "Problem not completely defined\nbe sure to select the asset type,"+...
                  " the model, the option and the method");
    return;
  end
  
  set_premia_obj_values (win, M);
  fic = xgetfile();
  if fic == "" then return; end
  save(fic, P=M);
endfunction

// handler for Help button
function help_activated(menuitem, L)
  win=L(1); 
  M = L(2); tag=L(3);
  M.get_help[tag];
  printf("display help on " + tag + "\n");  
endfunction


// mess_type = [GTK.MESSAGE_INFO;
//                GTK.MESSAGE_WARNING;
//                GTK.MESSAGE_QUESTION;
//                GTK.MESSAGE_ERROR]
function throw_message(t,mess)
  dialog = gtkmessagedialog_new (flags= GTK.DIALOG_MODAL, type= t, ...
                                 buttons= GTK.BUTTONS_OK,...
                                 message = mess );
  dialog.run[];
  dialog.destroy[];
endfunction

function raise_logo(menuitem, data)
  dialog = gtkdialog_new(title= "Interactive Dialog",parent=data(1),
                         flags = ior(GTK.DIALOG_MODAL,GTK.DIALOG_DESTROY_WITH_PARENT),
                         buttons = "gtk-ok");
  
  stock = gtkimage_new("file","premia.gif");
  dialog.vbox.pack_start[stock,expand=%f,fill=%f,padding=0];
  dialog.show_all[];
  reponse=dialog.run[];
  
  if reponse==1 then dialog.destroy[]; end
  
endfunction


