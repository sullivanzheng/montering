from math import sqrt

def show_distribution(l):
    d={}.fromkeys(l,())
    s=d.keys()
    for i in s:
        d[i]=l.count(i)
        print i,d[i]
    return d
        
def process(s):
    s=s.replace('=',' ')
    t=s.split()
    step,kmin,kmax,cenmark,e=[None]*5
    for j in xrange(len(t)):
        if t[j].find("Enow")!=-1:
            e=float(t[j+1])
        if t[j].find("Knot")!=-1:
            kmin=int(t[j+3])
            kmax=int(t[j+4])
        if t[j].find("Centralization")!=-1:
            cenmark="Y"
        if t[j].find("Angle_Adapter::Finished")!=-1:
            step=int(t[j+2])
    return [step,kmin,kmax,e,cenmark]

def processfile(infilename):
    f=file(infilename,'r')
    l=f.readlines()
    f.close()
    s=""
    r=[]
    for i in l:
        if i.find("Angle_")!=-1:
            r.append(process(s))
            s=i
        else:
            s=s+i
            
    energy_correct=1
    accu_E=0.0
    i=0
    if energy_correct==1:
        while i<len(r):
            if r[i][-2]!=None:
                r[i][-2]-=accu_E
            if r[i][-1]=="Y":
                print accu_E
                accu_E+=r[i][-2]
            i+=1
            
    r=[(i[0],i[2]-i[1],i[1],i[2],i[3],i[4]) for i in r if i[1]!=None]
    clear_r=[]
    i=0
    equ=100
    cen_cnt=0
    centers=0
    while i<len(r):
        if cen_cnt>0:
            cen_cnt-=1
        else:
            clear_r.append(r[i])
        if r[i][-1]=='Y':
            cen_cnt=equ
            centers+=1
        i+=1

    print "$ Number of raw records",len(r)
    print "$ Number of effective records",len(clear_r)
    print "$ Centralization time:",centers
#    show_distribution([i[0] for i in r])
    return clear_r,r

def subseg(l):
    firstcut=5000
    segnum=3
    l=l[firstcut:]
    print "Initial number of subsegments: ", len(l)+firstcut
    print "Cut preequilibrium subsegments: 1~%d "% firstcut
    print "Subsegments remains: %d~%d"%(firstcut+1,len(l)+firstcut)
    for segnum in [7,10,15,20,30,70,150]:
        batch=int(len(l)/segnum)
        print "Subsegments involved: %d~%d, divided into %d segs (%d/seg)"% \
              (firstcut+1,int(len(l)/batch)*batch+firstcut,int(len(l)/batch),batch)
        #MAIN
        s=[]
        for j in range(len(l)/batch):
            s.append(sum(l[j*batch:(j+1)*batch])/batch)
        tot=sum(s)
        sqtot=sum([i*i for i in s])
        avg=tot/len(s)
        std=sqrt((sqtot-tot*tot/len(s))/((len(s)-1)*len(s)))
        print "%.6f \t %.6f"%(avg,std)
    return None

def writecsv(outfile,l):
    f=file(outfile,'w')
    for i in l:
        s=[str(j) for j in i]
        s=','.join(s)
        print >>f,s
    f.close()
    return None

cnt=1
for i in ['./f35/f35.log']:#['./f35/f35.log','./f70/f70.log','f130/f130.log','f200/f200.log']:
    print "--------------------------"
    print "$ Processing file:",i
    clean,raw=processfile(i)
    subseg([float(p[1]) for p in clean])
    writecsv(str(cnt)+'_raw.csv',[[j[0],j[-2]] for j in raw])
    writecsv(str(cnt)+'_clean.csv',[[j[0],j[-2]] for j in clean])
    cnt+=1

