def printScalarResult(table_method, res):
    idx = 1
    table_method <= TR(TD(B("Result:"),align='right') + TD(B("")))
    for k,v in res:
        idx = idx + 1
        if type(v) == list:
            table_method <= TR(TD(k, align='right',rowspan=len(v)) + TD(v[0]),bgcolor=clr(res_colors,idx))
            for i in range(1, len(v)):
                table_method <= TR(TD(v[i]),bgcolor=clr(res_colors,idx))
        else:
            table_method <= TR(TD(k,align='right') + TD(v),bgcolor=clr(res_colors,idx))

