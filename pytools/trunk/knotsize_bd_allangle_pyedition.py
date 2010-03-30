from math import sin,cos,pi
import cProfile,os
def readJmol(name,delta=0):
    '''readJmol(name,delta=0), read Jmol xyz file with filename name
and when delta !=0, returns a list of vectors of each bond, otherwise, by
default, returns coordinates.'''
    f=file(name,'r')
    l=f.readlines()
    j=int(l[0].split()[0])
    #print "File: '%s' has %d vertices"%(name,j-1)
    l=[i.split() for i in l[1:j]]
    x=[float(i[2]) for i in l]
    y=[float(i[3]) for i in l]
    z=[float(i[4]) for i in l]
    #print "x,y,z loaded, length:",len(x),len(y),len(z)
    if delta!=0:
        x=[x[i+1]-x[i] for i in range(len(x)-1)]
        y=[y[i+1]-y[i] for i in range(len(y)-1)]
        z=[z[i+1]-z[i] for i in range(len(z)-1)]
    f.close()
    return x,y,z

def intsecxdx(x1,y1,dx1,dy1,
              x2,y2,dx2,dy2,l=1.0):
    '''Judge if (x1,y1)->dx1,dy2 and (x2,y2)->dx2,dy2 intersects.
if yes, return 0, if no, return 1.'''
    if abs(x1-x2)>2*l:return 1
    if abs(y1-y2)>2*l:return 1
    d=dx2*dy1-dx1*dy2
    if abs(d)<1E-9:return 1
    rx=(-dx1*(x2*dy2-y2*dx2)+dx2*(x1*dy1-y1*dx1))/d
    if dx1>0:
        head,tail=x1,x1+dx1
    else:
        head,tail=x1+dx1,x1
    if rx<head or rx>tail: return 1
    if dx2>0:
        head,tail=x2,x2+dx2
    else:
        head,tail=x2+dx2,x2
    if rx<head or rx>tail: return 1
    return 0

def intsecxx(x11,y11,x12,y12,
             x21,y21,x22,y22,l=1.0):
    '''Judge if (x11,y11)->(x12,y12) and (x21,y21)-(x22,y22) intersects.
if yes, return 0, if no, return 1.'''
    return intsecxdx(x11,y11,(x12-x11),(y12-y11),
                     x21,y21,(x22-x21),(y22-y21))

def knotlist(x,y):
    '''Return a list of tuples that indicate the sequence number of
intersected bonds. Sequence number started from 1, which is compatiable
with Jmol and xyz and fortran standard.'''
    if len(x)<3: return []
    secls=[]
    for i in range(len(x)-3):
        for j in range(i+2,len(x)-1):
            if intsecxx(x[i],y[i],x[i+1],y[i+1],
                        x[j],y[j],x[j+1],y[j+1])==0:
                secls.append((i+1,j+1))
                secls.sort()
    return secls

def knotrange(x,y):
    maxa=len(x)-1 #Max index of vectex
    head,tail=-1,-1
    for i in xrange(maxa-2):
        for j in xrange(i+2,maxa):
            if intsecxx(x[i],y[i],x[i+1],y[i+1],
                        x[j],y[j],x[j+1],y[j+1])==0:
                head=i
                break
        if head!=-1: break
        
    for i in xrange(maxa-1,1,-1):
        for j in xrange(i-2,-1,-1):
            if intsecxx(x[i],y[i],x[i+1],y[i+1],
                        x[j],y[j],x[j+1],y[j+1])==0:
                tail=i
                break
        if tail!=-1:break
        
    if (head,tail)==(-1,-1): return 0,0
    return head+1,tail+1
    
def Jmolfile_list(dir_name='.'):
    import os
    out=[]
    l=os.listdir(dir_name)
    for i in l:
        try:
            temp=int(i)
            out.append(str(i))
        except ValueError:
            continue
    return out

def main(sampleFreq):
    ls=[i for i in os.listdir('.') if 'dif_' in i]
    print "File listing complete, %d file(s) found."%len(ls)
    fo=file('knotsize_list.csv','w')
    proj_per=8 #how many perspectives will be examined
    proj_vec=[(cos(theta*pi/proj_per), sin(theta*pi/proj_per)) \
          for theta in range(proj_per)]
    count=0
    for f in ls[::sampleFreq][:]:
        print '>',f
        count+=1
        x,y,z=[i[100:300] for i in readJmol(f)]
        yzlist=zip(y,z)
        knot_range_list=[knotrange(x, [proj1*yi+proj2*zi \
                              for yi,zi in yzlist]) \
                        for proj1,proj2 in proj_vec]
        knot_size_list=[abs(last-first) for last,first in knot_range_list]
        print >>fo,','.join((f[f.find('_')+1:],str(sum(knot_size_list)/(len(knot_size_list))),
                             str(max(knot_size_list)),str(min(knot_size_list))))
        if count%10==0:
            print "%d"%count
        if count%100==0:
            print
    fo.close()

if __name__=='__main__':
    main(20)
#    cProfile.run('main()')



