from math import sqrt
from random import gauss, sample
def stat(s):
    tot=sum(s)
    sqtot=sum([i*i for i in s])
    avg=tot/len(s)
    std=sqrt((sqtot-tot*tot/len(s))/(len(s)-1))
    return avg,std

l=[str(gauss(0,100))+'\n' for i in xrange(1000000)]
file('randomdata.txt','w').writelines(l)

