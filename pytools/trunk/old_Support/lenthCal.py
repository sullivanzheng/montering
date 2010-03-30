from  scipy.integrate import quad
from math import sqrt
def cal_1(a,n,l):
    print "Converted a=%f, n=%d, l=%f"%(a,n,l)
    return sqrt(n*l**2*((1+a)/(1-a)-2*a*(1-a**n)/(n*(1-a)**2)))
def cal_2(g,n,l):
    pi=3.1415926536
    a1=quad(lambda x: sin(x)*cos(x)*exp(-g*x*x),0,pi)
    a2=quad(lambda x: sin(x)*exp(-g*x*x),0,pi)
    a=a1/a2
    return cal_1(a,n,l)
def cal_3(m,n,l):
    return cal_1((m-1.0)/(m+1.0),n,l)
    
print "=======Worm-like chain calculator========"
print "Choose one of the parameter sets:"
print "1. <Cos theta>, n and l"
print "2. g, n and l"
print "3. m(segments per Kuhn length),n(number of small segments, and l"
choice=int(input("Your choice: "))
if choice==1:
    cos_mean=float(raw_input("<COS theta>:"))
    n=int(raw_input("Number of small segments-n:"))
    l=float(raw_input("length of each segments-l:"))
    print "The length sqr<r**2>(a,n,l)= %f"%cal_1(a,n,l)
elif choice==2:
    g=float(raw_input("regidity:"))
    n=int(raw_input("Number of small segments-n:"))
    l=float(raw_input("length of each segments-l:"))
    print "the length sqr<r**2>(g,n,l)= %f"%cal_2(g,n,l)
elif choice==3:
    m=int(raw_input("Number of small segments/ Kuhn length:"))
    n=int(raw_input("Number of small segments-n:"))
    l=float(raw_input("length of each segments-l:"))
    print "The length sqr<r**2>(m,n,l)= %f"%cal_3(m,n,l)
else:
    print "error!"
