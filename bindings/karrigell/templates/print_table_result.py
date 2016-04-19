"""
def script(s):
    print "<script type=\"text/javascript\">"
    print s
    print "</script>"

def todata(v):
   s = "["
   for i in range(len(v)):
      s += "[" + str(res_t[0][1][i]) + "," + str(v[i]) + "],"
   s += "]"
   return s

#print REQUEST

result_labels = []
def printVectorResult(table_result_i, res_t):
    script("var datas = []")
    idx = 1
    for k,v in res_t:
        idx = idx + 1
        is_header = k == res_t[0][0]
        is_selected = ('showGraph_' + k in REQUEST)

        if type(v[0]) == list:

            table_result_i <= TR(TD(k, align='right',rowspan=len(v[0])) + Sum([TD(L[0]) for L in v]),bgcolor=clr(res_colors,idx))

            if not is_header:
               result_labels.extend(map(lambda i: (REQUEST['showGraphLabel_' + k]+"["+str(i)+"]", ('showGraph_' + k in REQUEST)), range(len(v[0]))))          

            if not is_header and is_selected:               
                script("datas.push({ data : "+todata(map(lambda L: L[0], v))+", label : \'"+REQUEST['showGraphLabel_' + k]+"[0]\' })")
            for i in range(1, len(v[0])):
                table_result_i <= TR(Sum([TD(L[i])  for L in v]),bgcolor=clr(res_colors,idx))
                if not is_header and is_selected:               
                   script("datas.push({ data : "+todata(map(lambda L: L[i], v))+", label : \'"+REQUEST['showGraphLabel_' + k]+"["+str(i)+"]\' })")          
        else:
            table_result_i <= TR(TD(k,align='right') + Sum([TD(L) for L in v]),bgcolor=clr(res_colors,idx))
            if not is_header:
               result_labels.append((REQUEST['showGraphLabel_' + k], ('showGraph_' + k in REQUEST)))
               
            if not is_header and is_selected:               
               script("datas.push({ data : "+todata(v)+", label : \'"+REQUEST['showGraphLabel_' + k]+"\' })")          

"""
