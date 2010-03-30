from scipy.integrate import quad
from scipy.optimize import fsolve
from math import *
pi=3.14159
def numerator(g):
    result=quad(lambda theta: cos(theta)*sin(theta)*exp(-g*theta*theta),0,pi)
    #print "NUMRT:IN:g=%12.6f OUT=%12.6f ERR=%12.6f"%(g,result[0],result[1])
    return result[0]

def denominator(g):
    result= quad(lambda theta: sin(theta)*exp(-g*theta*theta),0,pi)
    #print "DENMT:IN:g=%12.6f OUT=%12.6f ERR=%12.6f"%(g,result[0],result[1])
    return result[0]

for k in [28.57142857,]:
#k=5.0#282.686 #Kuhn length in the unit of segments.
    m=(k-1.0)/(k+1.0)
    print "%3f\t%12.9f"%(k,fsolve(lambda g:numerator(g)/denominator(g)-m,4))
