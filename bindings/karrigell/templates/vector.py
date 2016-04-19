from HTMLTags import *
from kspremia.field_base import *

def brackets(s, i):
  return s + '[' + str(i) + ']'

class Vector(FieldBase):

  def __init__(self,
               propertyName,
               friendlyName,
               fullName):
    super(Vector,self).__init__(propertyName,friendlyName,fullName)

  def render(self, v):

    pmem = getattr(v.entity, self.propertyName)

    def mc(idx):
      return INPUT(name=brackets(self.fullName, idx),value=pmem[idx])

    return [v.spannedRowsEx(self.friendlyName, map(mc, range(len(pmem))), '&#8477;')]

  def renderHistoryEx(self, v):

    pmem = getattr(v.entity, self.propertyName)
    return [v.spannedRowsEx(self.friendlyName, pmem)]

  def load(self, v):

      pmem = getattr(v.entity, self.propertyName)
      for i in range(len(pmem)):
        def setItem(x): pmem[i] = float(x)
        v.get(brackets(self.fullName, i), setItem)

  def getIterables(self, v):

    ctx = v.ctx

    pmem = getattr(v.entity, self.propertyName)

    for i in range(len(pmem)): 
      ctx._iterables.append(brackets(self.friendlyName, i))
      ctx._iterables_corr.append(self.fullName + str(i))

    ctx._iterables_getter.extend(map(lambda x: (lambda: x), pmem))
    ctx._iterables_setter.extend(map(lambda i: (lambda z: pmem.__setitem__(i, z)), range(len(pmem))))

def f(): pass