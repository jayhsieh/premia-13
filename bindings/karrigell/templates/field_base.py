def brackets(s, i):
  return s + '[' + str(i) + ']'

class FieldBase(object):

  def __init__(self,
               propertyName,
               friendlyName,
               fullName):

    self.propertyName = propertyName
    self.friendlyName = friendlyName
    self.fullName = fullName

  def process(self, v): return v.process(self)
