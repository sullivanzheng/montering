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

for k in [282./21,282./42,282./84]:#[282.686/21,]:#[28.57142857,]:
#k=5.0#282.686 #Kuhn length in the unit of segments.
    m=(k-1.0)/(k+1.0)
    g=fsolve(lambda g:numerator(g)/denominator(g)-m,4);
    print "k/2=%3f\tg=%12.9f\tfabs(g*2/(k/2)-1)=%12.9f"%(k/2.0,g,fabs(g*2/(k/2)-1))


