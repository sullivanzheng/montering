from math import sqrt
debug=0
firstcut=5000
segnum=3
#f=file('subseg_data.txt','r')
f=file('subseg_data.txt','r')
ls=f.readlines()
f.close()
l=[float(i) for i in ls][firstcut:]
print "Initial number of subsegments: ", len(l)+firstcut
print "Cut preequilibrium subsegments: 1~%d "% firstcut
print "Subsegments remains: %d~%d"%(firstcut+1,len(l)+firstcut)
for segnum in [5,7,10,15,20,30,50,70,100,150,200,500,1000]:
    batch=int(len(l)/segnum)
    print "Subsegments involved: %d~%d, divided into %d segs (%d/seg)"% \
          (firstcut+1,int(len(l)/batch)*batch+firstcut,int(len(l)/batch),batch)
    #MAIN
    s=[]
    for j in range(len(l)/batch):
        if debug:
            print "Summing from %d to %d"%(j*batch,(j+1)*batch-1)
        s.append(sum(l[j*batch:(j+1)*batch])/batch)
    if debug:
        for i in s:
            print i
    tot=sum(s)
    sqtot=sum([i*i for i in s])
    avg=tot/len(s)
    std=sqrt((sqtot-tot*tot/len(s))/((len(s)-1)*len(s)))
    print "%.6f+/-%.6f"%( avg,std)
