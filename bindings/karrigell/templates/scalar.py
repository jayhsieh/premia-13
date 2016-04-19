from HTMLTags import *
from kspremia.field_base import *

class Scalar (FieldBase):
  
   def __init__(self, propertyName,
                     friendlyName,
                     fullName,
                     constraint,
                     onChange,
                     isIterable,
                     converter):
      super(Scalar,self).__init__(propertyName, friendlyName, fullName)
      self.constraint = constraint
      self.onChange = onChange
      self.iterable = isIterable
      self.fromString = converter

   def load(self, v):
      v.get(self.fullName, 
         lambda x: setattr(v.entity, self.propertyName, self.fromString(x)))
          
   def render(self, v):
      mc = INPUT(name=self.fullName,onchange=self.onChange,value=getattr(v.entity, self.propertyName))  
      return [v.rowEx(self.friendlyName, mc, self.constraint)] if self.fromString != str else []

   def renderHistoryEx(self, v):
      return [v.rowEx(self.friendlyName, getattr(v.entity, self.propertyName))] if self.fromString != str else []

   def getIterables(self, v):
      ctx = v.ctx

      if self.iterable:
         ctx._iterables.append(self.friendlyName)
         ctx._iterables_corr.append(self.fullName)
         ctx._iterables_getter.append(lambda: getattr(v.entity, self.propertyName))
         ctx._iterables_setter.append(lambda x: setattr(v.entity, self.propertyName, x))

def f(): pass