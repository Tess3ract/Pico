class ECPoint:

  def __init__(self,x,y,i):
      self.x=x
      self.y=y
      self.i=i
  

  def __str__(self):
    return "x:{}  y:{}  identity:{}".format(self.x,self.y,self.i)
