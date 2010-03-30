from math import sqrt
from random import gauss,sample
debug=0
firstcut=00
segnum=3
#f=file('data.txt','r')
#ls=f.readlines()
#f.close()
#l=[float(i) for i in ls][firstcut:]
l=[gauss(0,1) for i in xrange(10000)]
print "Initial number of subsegments: ", len(l)+firstcut
print "Cut preequilibrium subsegments: 1~%d "% firstcut
print "Subsegments remains: %d~%d"%(firstcut+1,len(l)+firstcut)
for segnum in [2,3,5,7,10,15,20,30,50,70,100,150,200]:
    batch=int(len(l)/segnum)
    print "Subsegments involved: %d~%d, divided into %d segs (%d/seg)"% \
          (firstcut+1,int(len(l)/batch)*batch+firstcut,int(len(l)/batch),batch)
    #MAIN
    s=[]
    for j in range(100):
        s.append(sum(sample(l,batch))/batch)
    if debug:
        for i in s:
            print i
    tot=sum(s)
    sqtot=sum([i*i for i in s])
    avg=tot/len(s)
    std=sqrt((sqtot-tot*tot/len(s))/(len(s)-1)/len(s))
    print "%.6f+/-%.6f"%( avg,std)
