# main.py -- put your code here!

#import machine
#import utime
import time


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
        if(p1.i==True and p2.i==True):
            return ECPoint(None,None,True)
        elif(p1.i==True):
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
            return ECPoint (x3,y3,False)
        else:
            return ECPoint (None,None,True)
        
    
    def onCurve (self, point):
        if point.i==True: return True
        return ( pow(point.y,2,self.p) == (pow(point.x,3,self.p) + (self.A*point.x)%self.p + self.B)%self.p )
    
    def mult(self,k,point):
        if(point.i==True):
            return ECPoint(None,None,True)
        result = ECPoint(point.x,point.y,False)
        for i in range(2,k+1):
            result = self.addPoints(result,point)
        return result

    def multEfficient(self,k,point):
        #to binary
        binaryK = []
        while(k>0):
            binaryK.insert(0,k&1)
            k = k >> 1
        
        #square abd multiply
        result = ECPoint(None,None,True)
        for bit in binaryK:
            result = self.addPoints(result,result)
            if(bit&1):
                result = self.addPoints(result,point)
        return result
            

        
##########################################################################
#erweiterter euklidischer Algorithmus
def EEA(a, n):
    if a == 0:
        return (n, 0, 1)
    else:
        ggt, x, y = EEAIterativ(n % a, a)
        return (ggt, y - (n // a) * x, x)

def EEAIterativ(a,n):
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        (q, a), n = divmod(n, a), a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return (n, x0, y0)




def modInv(a,n):
    while(a<0): a+=n
    gcd, x, y = EEA(a,n)
    while(x<0):
        x=x+n
    return x
#################################################################
def performanceTest(k):
    print("Starting test...")
    el = EllipticCurve(15792089237316195423570985008687907853269984665640564039457583998564650230,115792089237316195423570985008687907853269984665640564039457583998564650230197,115792089237316195423570985008687907853269984665640564039457584007913129640233)
    p = ECPoint(104643561111080986962651636923349329984053149284385423934626192020695774567758,112131754190385317238677253571605670478896761682179953703038465984768056609039,False)
    start_time = time.time()
    pmult=(el.mult(k,p))
    delta_time = time.time() - start_time
    print("{} additions took {} seconds".format(k,delta_time))
    print("The average addition took {} ms".format(delta_time/k*1000))
    print(pmult.x == 47038396355578291256360099332357863034825867246118772816397986217504313036723 and pmult.y == 24508769116294836097001848252570143147127863609696015815968961242500353703447)
    print("Test ended.")



#################################################################

def standardTest():
    el = EllipticCurve(2,10,449)
    points1 = [
        ECPoint(48,306,False),
        ECPoint(50,89,False),
        ECPoint(50,360,False),
        ECPoint(51,204,False),
        ECPoint(51,245,False),
        ECPoint(58,99,False),
        ECPoint(58,350,False),
        ECPoint(60,126,False),
        ECPoint(60,323,False),
        ECPoint(61,154,False),
        ECPoint(61,295,False),
        ECPoint(63,156,False),
        ECPoint(63,293,False),
        ECPoint(65,110,False),
        ECPoint(None,None,True)
    ]
    points2 = [
        ECPoint(17,359,False),
        ECPoint(18,85,False),
        ECPoint(18,364,False),
        ECPoint(20,42,False),
        ECPoint(20,407,False),
        ECPoint(23,10,False),
        ECPoint(23,439,False),
        ECPoint(26,24,False),
        ECPoint(26,425,False),
        ECPoint(27,201,False),
        ECPoint(27,248,False),
        ECPoint(31,45,False),
        ECPoint(31,404,False),
        ECPoint(33,175,False),
        ECPoint(33,274,False),
        ECPoint(37,0,False),
        ECPoint(None,None,True)
    ]

    for p1 in points1:
        for p2 in points2:
            print(p1," + ",p2, " = ", el.addPoints(p1,p2))

#################################################################

def squareAndMultiplyPerformanceTest(k):
    el = EllipticCurve(15792089237316195423570985008687907853269984665640564039457583998564650230,226240555579959135501798302826772606856863713236549343666662115354037011922517,231584178474632390847141970017375815706539969331281128078915168015826259280027)
    p = ECPoint(84643561111080986662651636923349329984053149284385423934626192020695774567758,16942054190385317238677253571605670478896761682179953703038465984768056609039,False)
    start_time = time.time_ns()
    print(el.multEfficient(k,p))
    delta_time_efficient = time.time_ns() - start_time
    print("Efficient Multiplication of {}*point took {} milli seconds".format(k,delta_time_efficient/1000000))

    start_time = time.time()
    print(el.mult(k,p))
    delta_time_slow = time.time_ns() - start_time
    print("Inefficient Multiplication of {}*point took {} milli seconds".format(k,delta_time_slow/1000000))
    print("That is a speedup of {}".format(delta_time_slow/delta_time_efficient))
    


#led = machine.Pin(25, machine.Pin.OUT)
#while True:
    #led.value(1)
    #utime.sleep(1)
    #led.value(0)
    #utime.sleep(1)
#performanceTest(1000)
squareAndMultiplyPerformanceTest(10000)