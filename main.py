# main.py -- put your code here!

#import machine
#import utime


##########################################################################
class ECPoint:

  def __init__(self,x,y,i):
    if(i):
      self.x=None
      self.y=None
      self.i=True
    else:
        self.x=x
        self.y=y
        self.i=False

  def __str__(self):
    return "x:{}  y:{}  identity:{}".format(self.x,self.y,self.i)
##########################################################################


##########################################################################
class EllipticCurve:
# E: y^2 = x^3 + Ax + B mod p
    def __init__(self, A, B, p):
        self.A = A
        self.B = B
        self.p = p

    def __str__(self):
        return "E: y^2 = x^3 + {}x + {} mod {}".format(self.A,self.B,self.p)

    def addPoints(self,p1,p2):
        prime=self.p
        if(p1.i==True):
            return ECPoint(p2.x,p2.y,False)
        elif (p2.i == True):
            return ECPoint(p1.x,p1.y,False)
        elif(p1.x != p2.x):
            m = (((p2.y-p1.y)%prime) * modInv(p2.x-p1.x,prime)%prime) %prime
            x3 = (pow(m,2,prime) -p1.x-p2.x) %prime
            y3 = (m*(p1.x-x3)-p1.y) %prime
            return ECPoint (x3,y3,False)
        elif(p1.y != p2.y):
            return ECPoint (None,None,True)
        elif(p1.y!=0):
            m = (((3*pow(p1.x,2,prime) +self.A)%prime) * modInv(2*p1.y,prime)) %prime
            x3 = (pow(m,2,prime) -2*p1.x) %prime
            y3 = (m*(p1.x-x3) -p1.y) %prime
        else:
            return ECPoint (None,None,True)
        
    
    def onCurve (self, point):
        if point.i==True: return True
        return ( pow(point.y,2,self.p) == (pow(point.x,3,self.p) + (self.A*point.x)%self.p + self.B)%self.p )
        
##########################################################################
#extended gcd
def extendedGcd(a, n):
    if a == 0:
        return (n, 0, 1)
    else:
        gcd, x, y = extendedGcdIterative(n % a, a)
        return (gcd, y - (n // a) * x, x)

def extendedGcdIterative(a,n):
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        (q, a), n = divmod(n, a), a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return (n, x0, y0)




def modInv(a,n):
    while(a<0): a+=n
    gcd, x, y = extendedGcd(a,n)
    while(x<0):
        x=x+n
    return x
#################################################################

#led = machine.Pin(25, machine.Pin.OUT)
#while True:
    #led.value(1)
    #utime.sleep(1)
    #led.value(0)
    #utime.sleep(1)
el=EllipticCurve(5,10,17)
p1=ECPoint(15,3,False)
p2=ECPoint(12,9,False)
print(el.addPoints(p1,p2))
print(el.onCurve(ECPoint(333,222,True)))
print(el.onCurve(p2))