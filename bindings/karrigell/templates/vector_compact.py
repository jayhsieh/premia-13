from HTMLTags import *
from kspremia.field_base import *

class VectorCompact(FieldBase):

  def __init__(self,
               propertyName,
               friendlyName,
               fullName):

    super(VectorCompact, self).__init__(propertyName, friendlyName, fullName)
    self.varnameType = self.fullName + '__t'
    self.varnameConst = self.fullName + '__c'
    self.propnameMode = '__mode_'+self.propertyName

  def varnameIdx(self, i):
      return brackets(self.fullName, i)

  def load(self, v):

       pmem = getattr(v.entity, self.propertyName)

       isConstant = v.REQUEST[self.varnameType] == '0' \
                      if self.varnameType in v.REQUEST \
                      else all(map(lambda x: x == pmem[0], pmem))

       def fillPmem(val):
          for i in range(len(pmem)):
            pmem[i] = float(val)

      
       setattr(v.entity, self.propnameMode, isConstant)

       if not v.get(self.varnameConst, fillPmem):
          for i in range(len(pmem)):
            def setitem(x): pmem[i] = float(x)
            v.get(self.varnameIdx(i), setitem)

  def render(self, v):

    pmem = getattr(v.entity, self.propertyName)
    isConstant = getattr(v.entity, self.propnameMode)

    L = SELECT(name = self.varnameType, onchange='submit();').from_list(['Constant','Array'])

    if isConstant:
      L.select(value=0)
      mc = INPUT(name=self.varnameConst,value=pmem[0])
      return [v.spannedRowsEx(self.friendlyName, [L,mc],'&#8477;')]
    else:
      def mc(i): 
          return INPUT(name=self.varnameIdx(i),value=pmem[i])
      L.select(value=1)
      return [v.spannedRowsEx(self.friendlyName, [L]+map(mc, range(len(pmem))), '&#8477;')] 

  def renderHistoryEx(self, v):

    pmem = getattr(v.entity, self.propertyName)
    isConstant = getattr(v.entity, self.propnameMode)
     
    if isConstant:
      return [v.rowEx(self.friendlyName, pmem[0])]
    else:
      return [v.spannedRowsEx(self.friendlyName, pmem)]

  def getIterables(self, v):

     ctx = v.ctx

     pmem = getattr(v.entity, self.propertyName)

     isConstant = getattr(v.entity, self.propnameMode)
        
     if isConstant:
        ctx._iterables.append(self.friendlyName)
        ctx._iterables_corr.append(self.varnameConst)
        ctx._iterables_getter.append(lambda: pmem[0])
        def F(x):
           for i in range(len(pmem)): pmem[i] = x          
        ctx._iterables_setter.append(F)
     else:
        for i in range(len(pmem)): 
           ctx._iterables.append(self.varnameIdx(i))
           ctx._iterables_corr.append(self.fullName + str(i))
        ctx._iterables_getter.extend(map(lambda x: (lambda: x), pmem))
        ctx._iterables_setter.extend(map(lambda i: (lambda z: pmem.__setitem__(i, z)), range(len(pmem))))

def f(): pass
