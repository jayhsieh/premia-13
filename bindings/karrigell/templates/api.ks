import json
import sys
import traceback
from Queue import Empty
from multiprocessing import Process,Queue
from time import time

Include('../import.py')
import premia.assets
import premia.enum

def _return(x):
   print json.dumps(x)
   
def _lookupAsset(a):
   return eval('premia.assets.'+a)
   
def _lookupModel(m):
   return eval('premia.models.'+m)
   
def _lookupFamily(f):
   return eval('premia.options.'+f)
   
def _lookupOption(f, o):
   return eval('premia.options.'+f+'.'+o)
   
def _lookupMethod(m, f, meth):
   m_id = _lookupModel(m).ID()
   return eval('premia.pricings.'+ m_id + '.' + m_id + '_' + f + '.' + meth)

def assets():
   _return(map(lambda a: a.__name__, premia.assets.all()))
   
def models(a):
   _return(map(lambda m: m.__name__, _lookupAsset(a).models()))
   
def families(m):
   _return(_lookupModel(m).families())

def options(m, f):
   m_id = _lookupModel(m).ID()
   opts = _lookupFamily(f).all()
   pricing = eval('pricings.' + m_id + '.' + m_id + '_' + f + '.methods_for_options()')
   _return([o.__name__ for o in opts if o in pricing])
   
def methods(m, f, o):
   m_id = _lookupModel(m).ID()
   opt = _lookupOption(f, o)
   pricing = eval('pricings.' + m_id + '.' + m_id + '_' + f + '.methods_for_options()')
   methods = [x.__name__ for x in pricing[opt]]
   _return(methods)

def _params(obj):
   return obj.meta()
   
def model_params(m):
   _return(_params(_lookupModel(m)()))

def option_params(f, o):
   _return(_params(_lookupOption(f, o)()))

def method_params(m, f, meth):
   _return(_params(_lookupMethod(m, f, meth)()))

def method_results(m, f, meth):
   _return(_lookupMethod(m, f, meth).results() + [('Computation time', -1, False)])

def myImport(m):
   exec 'import ' + m
   return eval(m)

def _ksModelModule(m_id): 
   return myImport('kspremia.mod.' + m_id + '.model')

def _ksOptionModule(f, o): 
      return myImport('kspremia.opt.'+f+'.'+o)    

def _ksMethodModule(m_id,f,meth): 
   return myImport('kspremia.mod.'+m_id+'.'+m_id+'_'+f+'.'+meth)

def model_html(m):
   m_id = _lookupModel(m).ID()
   _return(_ksModelModule(m_id).html())

def option_html(f,o):
   _return(_ksOptionModule(f,o).html())

def family_html(f,o):
   _return(_ksOptionModule(f,o).familyHtml())

def method_html(m,f,meth):
   m_id = _lookupModel(m).ID()
   _return(_ksMethodModule(m_id, f, meth).html())

def _parse(kwargs):
   request = kwargs.iterkeys().__iter__().next()
   return json.loads(request)

def roundIfNeeded(dstType):
   if dstType == type(1):
      return lambda x: int(round(x))
   if dstType == type(1L):
      return lambda x: long(round(x))
   return lambda x: x

def iterate(initial, limit, stepsNo):

   cast = roundIfNeeded(type(initial))

   for i in range(stepsNo):
      x = initial + 1. * i / (stepsNo - 1) * (limit - initial) if stepsNo > 1 else initial
      yield cast(x)     

class Iterable(object):

   def __init__(self, label, setter, startValue, stopValue, interationsNo):
      self.keys = list(iterate(startValue, stopValue, interationsNo))
      self.name = label
      self.setter = setter
      self.stepsNo = interationsNo

_timeout = 30

def wrap_exc(F, q, *args):
   try:
      F(q, *args)
   except Exception, exc:
      exc_type, exc_value, exc_traceback = sys.exc_info()
      res = [("Exception", [str(exc)]+traceback.format_tb(exc_traceback))]
      q.put(res)
      return 

def _compute_impl(q, (mod, opt, method_obj)):
   begin = time()
   res = method_obj(opt, mod)
   end = time()
   res.append(("Time", end - begin))
   q.put(res)

def _fetch_scalar(queue):
   return queue.get(timeout=_timeout)

def _fetch_iterations(iterations):
   def inner(queue):
      data = []
      begin = time()
      for i in range(iterations):
         q = queue.get(timeout=begin + _timeout - time())
         if q[0][0] == "Exception":
            raise Exception, q[0][1][0]
         if len(data) == 0:
            for k,v in q:
               data.append((k, [v]))
         else:
            idx = 0
            for k,v in q:
                assert data[idx][0] == k
                data[idx][1].append(v)
                idx += 1
      return data
   return inner

def _spawn_process(G, keys, indata, fetcher):
   try: 
      queue = Queue()
      process = Process(target = wrap_exc, args = (G, queue, indata))
      try:            
         process.start()
         data = fetcher(queue)

      except Empty, exc:
         process.terminate()
         raise Exception("Method has worked more than " + str(_timeout) + "s. Please try another parameter combination")
                  
      process.join(timeout=_timeout)
      return keys+data
   except Exception, exc:
      return [("Exception", [str(exc)]+traceback.format_tb(sys.exc_info()[2]))]

def _compute_scalar(indata):
   return _spawn_process(_compute_impl, [], indata, _fetch_scalar)

def _compute_iteration_1d(indata, iteration):

   def G(q, indata):
      for x in iteration.keys:
         iteration.setter(x)
         _compute_impl(q, indata)
  
   keys = [
      (iteration.name, iteration.keys)
   ]

   return _spawn_process(G, keys, indata, _fetch_iterations(iteration.stepsNo))

def _compute_iteration_2d(indata, iteration_1, iteration_2):

   def G(q, indata):
      for x_1 in iteration_1.keys:
         for x_2 in iteration_2.keys:
            iteration_2.setter(x_2)
            iteration_1.setter(x_1)
            _compute_impl(q, indata)

   keys = [
      (iteration_1.name, iteration_1.keys),
      (iteration_2.name, iteration_2.keys)
   ]
   
   return _spawn_process(G, keys, indata, _fetch_iterations(iteration_1.stepsNo * iteration_2.stepsNo))

def compute(*args, **kwargs):
   [asset, model, model_params], [family, option, option_params], [method, method_params] = _parse(kwargs)

   iterables = []

   def append_iterables(label, setter, startValue, stopValue, interationsNo):
      iterables.append(Iterable(label, setter, startValue, stopValue, interationsNo))

   indata = (_lookupModel(model).create(model_params, append_iterables),
             _lookupOption(family, option).create(option_params, append_iterables),
             _lookupMethod(model, family, method).create(method_params, append_iterables))
   
   result = _compute_scalar(indata) if iterables == [] else \
            _compute_iteration_1d(indata, iterables[0]) if len(iterables) == 1 else \
            _compute_iteration_2d(indata, iterables[0], iterables[1]) if len(iterables) == 2 else \
            []

   _return(result)

def adjust(*args, **kwargs):

   [[[asset, model, model_params], [family, option, option_params], [method, method_params]], [path, value]] = _parse(kwargs)

   def nope(*args): pass

   indata = {'model' : _lookupModel(model).create(model_params, nope), 
             'option' : _lookupOption(family, option).create(option_params, nope), 
             'method' : _lookupMethod(model, family, method).create(method_params, nope) }

   obj = indata[path[0]]
   setattr(obj, path[1], value)
   
   _return([indata["model"].meta(), indata["option"].meta(), indata["method"].meta()])




