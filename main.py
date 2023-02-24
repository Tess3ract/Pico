# main.py -- put your code here!

import time
import random
import machine


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

    def __str__(self) -> str:
        return "[x:{}  y:{}  identity:{}]".format(self.x,self.y,self.i)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, ECPoint):
            if(self.i==True and other.i==True):
                return True
            else:
                return ((self.i == other.i) and (self.x == other.x) and (self.y == other.y))
        return False
##########################################################################




##########################################################################
class EllipticCurve:
# E: y^2 = x^3 + Ax + B mod p
    def __init__(self, A:int, B:int, p:int):
        self.A = A
        self.B = B
        self.p = p

    def __str__(self) -> str:
        return "[E: y^2 = x^3 + {}x + {} mod {}]".format(self.A,self.B,self.p)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, EllipticCurve):
            return ((self.A == other.A) and (self.B == other.B) and (self.p == other.p))
        return False

    def addPointsProjective(self, p1:ECPoint, p2:ECPoint) -> ECPoint:
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

    def doubleProjective(self, p:ECPoint) -> ECPoint:
        return self.addPointsProjective(p,p)    
    
    def onCurve (self, point:ECPoint) -> bool:
        if point.i==True: return True
        return ( pow(point.y,2,self.p) == (pow(point.x,3,self.p) + (self.A*point.x)%self.p + self.B)%self.p )
    
    def mult(self, k:int, point:ECPoint) -> ECPoint:
        if(point.i==True):
            return ECPoint(None,None,True)
        result = ECPoint(point.x,point.y,False)
        for i in range(2,k+1):
            result = self.addPointsProjective(result,point)
        return result

    def multEfficient(self, k:int, point:ECPoint) -> ECPoint:
        #to binary
        binaryK = []
        while(k>0):
            binaryK.insert(0,k&1)
            k = k >> 1
        
        #square and multiply
        result = ECPoint(None,None,True)
        for bit in binaryK:
            result = self.doubleProjective(result)
            if(bit&1):
                result = self.addPointsProjective(result,point)
        return result
            

        
##########################################################################
class Curve25519 (EllipticCurve):
    # E: B*y^2 = x^3 + A x^2 + x mod p
    def __init__(self):
        EllipticCurve.__init__(self,486662,1, pow(2,255) - 19)
    
    def onCurve (self, point:ECPoint) -> bool:
        if point.i==True: return True
        return (self.B*pow(point.y,2,self.p) % self.p)  == ((pow(point.x,3,self.p) + self.A* pow(point.x,2,self.p) + point.x)%self.p )

    
    def addPointsAffine(self,  p1:ECPoint, p2:ECPoint) -> ECPoint:
        if(p1.i==True and p2.i == True):
            return ECPoint(None,None,True)
        elif (p1.i == True):
            return ECPoint(p2.x,p2.y,False)
        elif (p2.i == True):
            return ECPoint(p1.x,p1.y,False)
        elif(p1.x == p2.x):
            if (p1.y == p2.y):   #p1 = p2
                return self.doubleAffine(p1)  
            else:  # p1 = p2^-1
                return ECPoint(None,None,True)


        inverseDeltaX = modInv(p2.x -p1.x, self.p)
        x3 = (pow(p2.y-p1.y,2,self.p) * pow(inverseDeltaX,2,self.p) - self.A - p1.x -p2.x) %self.p
        y3 = ((2*p1.x+p2.x+self.A)*(p2.y-p1.y)* inverseDeltaX - pow(p2.y-p1.y,3,self.p)* pow(inverseDeltaX,3,self.p) -p1.y) %self.p
 

        return ECPoint(x3,y3,False)

    
    def addPointsProjective(self, p1:ECPoint, p2:ECPoint) -> ECPoint:
        #in case some point is point at infinity
        if(p1.i==True and p2.i == True):
            return ECPoint(None,None,True)
        elif (p1.i == True):
            return ECPoint(p2.x,p2.y,False)
        elif (p2.i == True):
            return ECPoint(p1.x,p1.y,False)

        #affine to projective
        X1 = p1.x
        Y1 = p1.y
        Z = 1
        X2 = p2.x
        Y2 = p2.y

        #calculation
        U1 = Y2*Z % self.p
        U2 = Y1*Z % self.p
        V1 = X2*Z % self.p
        V2 = X1*Z % self.p
        if(V1==V2):
            if(U1 != U2):
                return ECPoint(None,None,True)
            else:
                return self.doubleProjective(p1)
        U = U1 - U2 % self.p
        V = V1 - V2 % self.p
        W = Z*Z % self.p
        A = (pow(U,2,self.p)*W - pow(V,3,self.p) - 2* pow(V,2,self.p) * V2) % self.p
        X3 = V*A % self.p
        Y3 = (U*(pow(V,2,self.p)*V -A) - pow(V,3,self.p) * U2) % self.p
        Z3 = (pow(V,3,self.p)*W) % self.p

        #projective to affine
        Z_inverse = modInv(Z3,self.p)
        x_new = (X3 * Z_inverse - self.A) % self.p
        y_new_squared = (pow(x_new,3,self.p) + self.A* pow(x_new,2,self.p) + x_new)%self.p
        y_new = modular_sqrt(y_new_squared,self.p) %self.p

        #correct y_new? maybe we got -(y_new) here

        #gradient of line (correct m)
        m = (p2.y - p1.y) * modInv(p2.x - p1.x,self.p)
        #gradient with new point
        m_new = (p2.y - y_new) * modInv(p2.x - x_new,self.p)
        #if m != m_new we need to recalculate y_new
        if(m != m_new):
            y_new = self.p - y_new

        return ECPoint(x_new,y_new,False)


        
    def doubleAffine(self, point:ECPoint) -> ECPoint:
        if(point.i == True):
            return ECPoint(None,None,True)
        
        c = ((3 * pow(point.x,2,self.p) + 2 * self.A * point.x +1) * modInv(2 * self.B * point.y, self.p)) %self.p
        

        x_new = (self.B * pow(c,2,self.p) -self.A - 2 * point.x) %self.p

        y_new = ((2 * point.x + point.x + self.A)*c - pow(c,3,self.p) - point.y) %self.p

        

        return ECPoint(x_new,y_new,False)

    def doubleProjective(self, point:ECPoint) -> ECPoint:
        if(point.i == True):
            return ECPoint(None,None,True)

        X1 = point.x
        Z = 1
        XX1 = pow(X1,2,self.p)
        X3 = pow(XX1-1,2,self.p)
        Z3 = 4*X1*(XX1 + self.A*X1 +1) % self.p

        Z_inverse = modInv(Z3, self.p)
        x_new = X3 * Z_inverse % self.p

        #y_new = ((pow(point.x,2,self.p)+1)*(2*point.x+2*self.A) - 2*self.A) * modInv(2*point.y,self.p) % self.p
        y_new_squared = (pow(x_new,3,self.p) + self.A* pow(x_new,2,self.p) + x_new)%self.p
        y_new = modular_sqrt(y_new_squared,self.p) %self.p

        #correct y_new? maybe we got -(y_new) here

        #gradient of line (correct m)
        #m = (p2.y - p1.y) * modInv(p2.x - p1.x,self.p)
        #gradient with new point
        # = (p2.y - y_new) * modInv(p2.x - x_new,self.p)
        #if m != m_new we need to recalculate y_new
        #if(m != m_new):
            #y_new = self.p - y_new


        return ECPoint(x_new,y_new,False)

    def mult(self, k:int, point:ECPoint) -> ECPoint:
        if(point.i==True):
            return ECPoint(None,None,True)
        result = ECPoint(point.x,point.y,False)
        for i in range(2,k+1):
            print("p*", i-1, " : ", result)
            result = self.addPointsAffine(result,point)
        return result

    def multEfficient(self, k:int, point:ECPoint) -> ECPoint:
        #to binary
        binaryK = []
        while(k>0):
            binaryK.insert(0,k&1)
            k = k >> 1
        
        #square and multiply
        result = ECPoint(None,None,True)
        for bit in binaryK:
            result = self.doubleAffine(result)
            #result = self.double(result)
            if(bit&1):
                result = self.addPointsAffine(result,point)
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

def modular_sqrt(a, p):
    """ Find a quadratic residue (mod p) of 'a'. p
        must be an odd prime.

        Solve the congruence of the form:
            x^2 = a (mod p)
        And returns x. Note that p - x is also a root.

        0 is returned is no square root exists for
        these a and p.

        The Tonelli-Shanks algorithm is used (except
        for some simple cases in which the solution
        is known from an identity). This algorithm
        runs in polynomial time (unless the
        generalized Riemann hypothesis is false).
    """
    # Simple cases
    #
    if legendre_symbol(a, p) != 1:
        return 0
    elif a == 0:
        return 0
    elif p == 2:
        return 0
    elif p % 4 == 3:
        var2 = (p+1) >> 2
        return pow(a, var2, p)

    # Partition p-1 to s * 2^e for an odd s (i.e.
    # reduce all the powers of 2 from p-1)
    #
    s = p - 1
    e = 0
    while s % 2 == 0:
        s = s >> 1
        e += 1

    # Find some 'n' with a legendre symbol n|p = -1.
    # Shouldn't take long.
    #
    n = 2
    while legendre_symbol(n, p) != -1:
        n += 1

    # Here be dragons!
    # Read the paper "Square roots from 1; 24, 51,
    # 10 to Dan Shanks" by Ezra Brown for more
    # information
    #

    # x is a guess of the square root that gets better
    # with each iteration.
    # b is the "fudge factor" - by how much we're off
    # with the guess. The invariant x^2 = ab (mod p)
    # is maintained throughout the loop.
    # g is used for successive powers of n to update
    # both a and b
    # r is the exponent - decreases with each update
    #
    var1 = (s+1) >> 1
    x = pow(a, var1, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e

    while True:
        t = b
        m = 0
        for m in range(r):
            if t == 1:
                break
            t = pow(t, 2, p)

        if m == 0:
            return x

        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m


def legendre_symbol(a, p):
    """ Compute the Legendre symbol a|p using
        Euler's criterion. p is a prime, a is
        relatively prime to p (if p divides
        a, then a|p = 0)

        Returns 1 if a has a square root modulo
        p, -1 otherwise.
    """
    p_half = (p-1) >> 1

    ls = pow(a, p_half, p)
    return -1 if ls == p - 1 else ls

#################################################################

#Elliptic Curve curve, generator g
class DiffieHellman:
    def __init__(self, curve:EllipticCurve, g:ECPoint):
        self.curve = curve
        self.g = g
    
    def __str__(self) -> str:
        return "Diffie-Hellman Protocol with prime {} and generator {}".format(self.curve.p,self.g)
    
    def generatePrivateKey(self) -> int:
        return random.randint(2
        ,self.curve.p-1) #this is not secure. It's only pseudo-random generator. But it's not really relevant either

    #private key a
    def generatePublicKey(self, a:int) -> ECPoint:
        return self.curve.multEfficient(a,self.g)
    
    def generateSharedKey(self, a:int, gb:ECPoint) -> ECPoint:
        return self.curve.multEfficient(a,gb)
    
    def demonstration(self, a:int, b:int):
        #A sends g^a to B
        ga = self.generatePublicKey(a)
        print("A sends {} to B.".format(ga))
        #B sends g^b to A
        gb = self.generatePublicKey(b)
        print("B sends {} to A.".format(gb))
        #A calculates shared key
        ka = self.generateSharedKey(a,gb)
        print("A calculates shared key: ",ka)
        #B calculates shared key
        kb = self.generateSharedKey(b,ga)
        print("B calculates shared key: ",kb)
        led = machine.Pin(25, machine.Pin.OUT)
        if(ka == kb):
            print("Key exchange was successful!")
            for i in range(8):
                led.value(1)
                time.sleep(1)
                led.value(0)
                time.sleep(1)
        else:
            print("Error! \t Key exchange was not successful!")
            for i in range(40):
                led.value(1)
                time.sleep(0.2)
                led.value(0)
                time.sleep(0.2)





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
            print(p1," + ",p2, " = ", el.addPointsProjective(p1,p2))

#################################################################

def squareAndMultiplyPerformanceTest(k):
    el = EllipticCurve(15792089237316195423570985008687907853269984665640564039457583998564650230,226240555579959135501798302826772606856863713236549343666662115354037011922517,231584178474632390847141970017375815706539969331281128078915168015826259280027)
    p = ECPoint(84643561111080986662651636923349329984053149284385423934626192020695774567758,16942054190385317238677253571605670478896761682179953703038465984768056609039,False)
    start_time = time.ticks_us()
    print(el.multEfficient(k,p))
    delta_time_efficient = time.ticks_us() - start_time
    print("Efficient Multiplication of {}*point took {} micro seconds".format(k,delta_time_efficient))

    start_time = time.ticks_us()
    print(el.mult(k,p))
    delta_time_slow = time.ticks_us() - start_time
    print("Inefficient Multiplication of {}*point took {} micro seconds".format(k,delta_time_slow))
    print("That is a speedup of {}".format(delta_time_slow/delta_time_efficient))
    
###############################################################
def DiffieHellmanTest():
    el = EllipticCurve(15792089237316195423570985008687907853269984665640564039457583998564650230,115792089237316195423570985008687907853269984665640564039457583998564650230197,115792089237316195423570985008687907853269984665640564039457584007913129640233)
    g = ECPoint(104643561111080986962651636923349329984053149284385423934626192020695774567758,112131754190385317238677253571605670478896761682179953703038465984768056609039,False)
    dh =DiffieHellman(el,g)
    dh.demonstration(57896044618658097711785492504343953926512770110598059797506569779613311366406, 47606440156792287940606943253909558533971493099538253817755912803560908337955)

def DiffieHellmanTestCurve25519():
    
    g = ECPoint(9, 14781619447589544791020593568409986887264606134616475288964881837755586237401,False)
    curve = Curve25519()
    dh = DiffieHellman(curve,g)
    dh.demonstration(57896044618658097711785492504343953926512770110598059797506569779613311366406,47606440156792287940606943253909558533971493099538253817755912803560908337955)

###############################################################
#a^(-1) equals a^(p-2)
def inverseMultRatioTest(number:int, epoch:int):
    led = machine.Pin(25, machine.Pin.OUT)
    p = 115792089237316195423570985008687907853269984665640564039457584007913129640233
    time_mult=0
    time_inverse=0
    for i in range(0,epoch):
        start_time = time.ticks_us()
        pow(number,2,p)
        delta_time = time.ticks_us() - start_time
        time_mult += delta_time
        led.value(1)

        start_time = time.ticks_us()
        modInv(number,p)
        delta_time = time.ticks_us() - start_time
        time_inverse += delta_time
        led.value(0)

    print("Average inverse in micro seconds: ", time_inverse/epoch)
    print("Average mult in micro seconds: ", time_mult/epoch)
    print("Ratio is ", time_inverse/time_mult)
###############################################################

def Curve25519Test():
    curve = Curve25519()
    g = ECPoint(9, 14781619447589544791020593568409986887264606134616475288964881837755586237401,False)
    p = ECPoint(9, 43114425171068552920764898935933967039370386198203806730763910166200978582548, False)
    q = ECPoint(14847277145635483483963372537557091634710985132825781088887140890597596352251,48981431527428949880507557032295310859754924433568441600873610210018059225738,False)
    print("On Curve? -> True, True, False")
    print(curve.onCurve(g))
    print(curve.onCurve(q))
    print(curve.onCurve(ECPoint(1,4434243,False)))
    pq = curve.addPointsProjective(p,q)
    p2 = curve.doubleProjective(g)
    print("Addition ------------------------------------")
    print("P2 -------------------------------------")
    print(p2)
    print(curve.doubleAffine(g))
    print("P3 -------------------------------------")
    print(curve.addPointsAffine(p2,g))
    print("P4 -------------------------------------")
    print(curve.doubleProjective(p2))
    print(curve.doubleAffine(p2))
    print("P8 -------------------------------------")
    print(curve.doubleProjective(curve.doubleProjective(p2)))
    print(curve.doubleAffine(curve.doubleAffine(p2)))
    print("P + Q -------------------------------------")
    print(pq)
    print("P + Q affine -----------------------------------------" )
    print(curve.addPointsAffine(p,q))
    #print("P4 -------------------------------------")
    #print(curve.double(p2))
    #p3 = curve.addPoints(p1,p2)
    #print(curve.addPoints(p3,p1))
    #print(curve.addPoints(p2,p2))
    print("100 * g: ----------------------------------------------------------")
    print(curve.mult(1000,g))


def AffineVsProjectiveTest():
    curve = Curve25519()
    a_start = ECPoint(9, 14781619447589544791020593568409986887264606134616475288964881837755586237401,False)
    b_start = ECPoint(9, 43114425171068552920764898935933967039370386198203806730763910166200978582548, False)
    c_start = ECPoint(14847277145635483483963372537557091634710985132825781088887140890597596352251,48981431527428949880507557032295310859754924433568441600873610210018059225738,False)
    d_start = ECPoint(12697861248284385512127539163427099897745340918349830473877503196793995869202,39113539887452079713994524130201898724087778094240617142109147539155741236674,False)

    #affine
    start_time = time.ticks_us()
    a = a_start
    b = b_start
    c = c_start
    d = d_start
    for x in range (0,50):
        a = curve.addPointsAffine(a,b)
        b = curve.addPointsAffine(d,c)
        c = curve.addPointsAffine(a,d)
        d = curve.addPointsAffine(b,c)
    delta_time_affine_addition = time.ticks_us() - start_time
    print("affine addition: ", delta_time_affine_addition / 200, " micro seconds")

    start_time = time.ticks_us()
    for x in range (0,200):
        a = curve.doubleAffine(a)
    delta_time_affine_double = time.ticks_us() - start_time
    print("affine doubling: ", delta_time_affine_double / 200, " micro seconds")

    #projective
    start_time = time.ticks_us()
    a = a_start
    b = b_start
    c = c_start
    d = d_start
    for x in range (0,50):
        a = curve.addPointsProjective(a,b)
        b = curve.addPointsProjective(d,c)
        c = curve.addPointsProjective(a,d)
        d = curve.addPointsProjective(b,c)
    delta_time_projective_addition = time.ticks_us() - start_time

    print("projective addition: ", delta_time_projective_addition / 200, " micro seconds")

    start_time = time.ticks_us()
    for x in range (0,200):
        a = curve.doubleProjective(a)
    delta_time_projective_double = time.ticks_us() - start_time
    print("projective doubling: ", delta_time_projective_double / 200 , " micro seconds")

    print("projective to affine ratio:")
    print("addition: ", delta_time_projective_addition / delta_time_affine_addition)
    print("doubling ", delta_time_projective_double / delta_time_affine_double)



###############################################################################################
def WeierstrassVsCurve25519Test():
    #initialize Weierstraß Curve
    weierstrass = EllipticCurve(15792089237316195423570985008687907853269984665640564039457583998564650230,115792089237316195423570985008687907853269984665640564039457583998564650230197,115792089237316195423570985008687907853269984665640564039457584007913129640233)
    g_w = ECPoint(104643561111080986962651636923349329984053149284385423934626192020695774567758,112131754190385317238677253571605670478896761682179953703038465984768056609039,False)
    #initialize Curve25519
    curve = Curve25519()
    g_c = ECPoint(14847277145635483483963372537557091634710985132825781088887140890597596352251,48981431527428949880507557032295310859754924433568441600873610210018059225738,False)

    start_time = time.ticks_us()
    weierstrass.multEfficient(24623137305878816398398911160594559162372272573335115162610566893532021787955,g_w)
    time_w = time.ticks_us() - start_time

    start_time = time.ticks_us()
    curve.multEfficient(24623137305878816398398911160594559162372272573335115162610566893532021787955,g_c)
    time_c = time.ticks_us() - start_time

    print("Multiplication on Weierstrass curve took {} micro seconds".format(time_w))
    print("Multiplication on Curve25519 took {} micro seconds".format(time_c))
    print("Curve25519 is {} times faster than Weierstrass Curve".format(time_w/time_c))





###############################################################

#performanceTest(1000)
#squareAndMultiplyPerformanceTest(100000)
#DiffieHellmanTestCurve25519()
#DiffieHellmanTest()
#inverseMultRatioTest(1046435611110809869626516369233493299840531492843854239346266420534567758,1000)
#Curve25519Test()
#AffineVsProjective()
WeierstrassVsCurve25519Test()