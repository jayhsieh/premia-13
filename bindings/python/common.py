def getParametersAsList(x): return [x.__dict__['_' + y] for y in x.parameters()] 

def getParametersAsDict(x):  
	d = dict()
	for y in x.parameters():
	   d[y] = x.__dict__['_' + y]
	return d

def getRepr(x, label):
	return label + ' ' + x.__class__.__name__ + ' ' + getParametersAsDict(x).__repr__() 

def extend(lst, size): lst.extend((size - len(lst)) * [lst[0]]) 

def writearray(x): 
   import interop 
   interop.write_array(len(x)) 
   for y in x: interop.write_array_double(y)
   interop.write_array_end() 

def ignorearray(x): 
   import interop 
   interop.ignore_array(len(x)) 
   for y in x: interop.ignore_array_double(y)
   interop.ignore_array_end() 

def readarray(): 
   import interop 
   L = []
   sz = interop.read_array_size() 
   for i in range(0,sz): L.append(interop.read_array_double())
   interop.read_array_end() 
   return L

def getResultArray(idx):
   import interop 
   N = interop.getResultArraySize(idx) 
   R = []
   for i in range(0, N): 
      R.append(interop.getResultArrayItem(idx, i)) 
   return R
