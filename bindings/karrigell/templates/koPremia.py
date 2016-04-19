from HTMLTags import *
Include('../import.py')

from premia import interop

def js_include(filename):
   return SCRIPT(type="text/javascript", src=filename)

print """
<head>
   <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
   <link rel="stylesheet" type="text/css" href="style.css" />
   <link rel="stylesheet" type="text/css" href="hsd-flotr2.css" />
   <script type="text/javascript" src="jquery-1.8.2.min.js"></script>
   <script type="text/javascript" src="jquery.json-2.3.min.js"></script>
   <script type="text/javascript" src="knockout-2.1.0.js"></script>
   <script type="text/javascript" src="flotr2.min.js"></script>
</head>
"""
Include("shaders.py")
Include("param_table.html")
Include("current.html")

print SCRIPT("var data_dir = '" + interop.data_dir + "';")

if len(REQUEST):
   request = REQUEST.iterkeys().__iter__().next()
   print SCRIPT("var s = '"+ request + "';", type=r"text/javascript")
   print SCRIPT("var request = $.parseJSON(s);", type=r"text/javascript") 
else:
   print SCRIPT("var request = undefined;", type=r"text/javascript") 


print js_include("script.js")
