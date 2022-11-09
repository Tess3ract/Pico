# Python program for the extended Euclidean algorithm

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
    gcd, x, y = extendedGcd(a,n)
    while(x<0):
        x=x+n
    return x
