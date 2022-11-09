class EllipticCurve:
# E: y^2 = x^3 + Ax + B mod p
  def __init__(self, A, B, p):
    self.A = A
    self.B = B
    self.p = p

  def __str__(self):
    return "E: y^2 = x^3 + {}x + {} mod {}".format(self.A,self.B,self.p)

  