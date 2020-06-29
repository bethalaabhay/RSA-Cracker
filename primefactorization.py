from random import randint
from math import gcd]

def primes(limit):
    b = [True] * (limit + 1)
    primenumberlist= []
    for p in range(2, limit  + 1):
        if b[p]:
            primenumberlist.append(p)
            for i in range(p, limit + 1, p):
                b[i] = False
    #print(primenumberlist)            
    return primenumberlist

# Returns modulo inverse, unused helper and gcd
def modular_inv(a, b):
    if b == 0:
        return 1, 0, a
    q, r = divmod(a, b) #returns a/b and a%b
    x, y, g = modular_inv(b, r)
    return y, x - q * y, g

# Addition in Elliptic curve modulo m space
def elliptic_add(p, q, a, b, n):
    # If one point is infinity, return other one
    if p[2] == 0: return q
    if q[2] == 0: return p
    if p[0] == q[0]:
        if (p[1] + q[1]) % n== 0:
            return 0, 1, 0  # Infinity
        num = (3 * p[0] * p[0] + a) % n
        denom = (2 * p[1]) % n
    else:
        num = (q[1] - p[1]) % n
        denom = (q[0] - p[0]) % n
    inv, _, g = modular_inv(denom, n)
   
   
   # inv=modInverse(denom,n)
    #g=GCD(denom,n)
    # Unable to find inverse, arithmetic breaks
    if g > 1:
        return 0, 0, denom  # Failure
    
    lam=num*inv
    x3 = (lam*lam - p[0] - q[0]) % n
    y3 = (lam*(p[0]-x3)-p[1])%n
    
    return x3,y3, 1

# Multiplication (repeated addition and doubling)
def elliptic_mul(p, q, a, b, n):
    r = (0, 1, 0)  # Infinity
    while p > 0:
        # q is failure, return it
        if q[2] > 1:
            return q
        if p % 2 == 1:
            r = elliptic_add(q, r, a, b, n)
        p = p // 2
        q= elliptic_add(q, q, a, b, n)
    #print(r)    
    return r

# Limit specifies the amount of work permitted
def lenstra(n, limit):
    g = n
    while g == n:
        # Randomized x and y
        q = randint(0, n - 1), randint(0, n - 1), 1
        # Randomized curve coefficient a, computed b
        a = randint(0, n - 1)
        b = (q[1] * q[1] - q[0] * q[0] * q[0] - a * q[0]) % n
        g = gcd(4 * a * a * a + 27 * b * b, n)  # singularity check
    # If we got lucky, return lucky factor
    if g > 1:
        return g
    # increase k step by step until lcm(1, ..., limit)
    primenumberlist=primes(limit)
    for p in primenumberlist:
        pp = p
        while pp < limit:
            q = elliptic_mul(p, q, a, b, n)
            # Elliptic arithmetic breaks
            #print(q[2])
            if q[2] > 1:
                return gcd(q[2], n)
            pp = p * pp
    return False


def main():
    n=int(input("Enter the Number to Factorize: "))
    limit=int(input("Enter the Search Limit:  "))


    res=lenstra(n,limit)

    if(res==False):
        print("The Number Does not have a prime factor in the given B-Smooth Limit")

    else:
        print(f"The Prime Factors are {res} and {n//res}  ")
if __name__ == '__main__':
    main()