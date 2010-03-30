from convertFor2Jmol import write_Jmol_file
from numpy import matrix,concatenate
from scipy.linalg import det,norm
from math import sin,cos,pi,acos


def readJmol(name,delta=0):
    '''readJmol(name,delta=0), read Jmol xyz file with filename name
and when delta !=0, returns a list of vectors of each bond, otherwise, by
default, returns coordinates.'''
    f=file(name,'r')
    l=f.readlines()
    j=int(l[0].split()[0])
    print "File: '%s' has %d vectors"%(name,j-1)
    l=[i.split() for i in l[1:j+1]]
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

def x2dx(x):
    return [x[i+1]-x[i] for i in range(len(x)-1)]

def dx2x(dx,base=0):
    """Input: dx (N-vector)
    Return: x (N+1 vector) with base of (by default) 0"""
    x=[base]
    for j in dx:
        x.append(x[-1]+j)
    return x

def dx2xN(dx):
    return dx2x(dx,0)[1:]

def write_Jmol_file_dx(filename,dx,dy,dz):
    write_Jmol_file(filename,dx2x(dx),dx2x(dy),dx2x(dz))

def makeCircle(n):
    dx=[cos(2*pi/n*i) for i in range(n)]
    dy=[sin(2*pi/n*i) for i in range(n)]
    dz=[0]*n
    return dx,dy,dz

def makeCircle_File(filename,n):
    dx,dy,dz=makeCircle(n)
    write_Jmol_file_dx(filename,dx,dy,dz)

def rotdx(dx,dy,dz,rotMat):
    return rotMat*matrix([dx,dy,dz])

def Xpt(v1,v2):
    """cross product"""
    return(v1[1]*v2[2]-v1[2]*v2[1],
           -v1[0]*v2[2]+v1[2]*v2[0],
           v1[0]*v2[1]-v1[1]*v2[0])

def betav1v2(mv1,mv2):
    return acos(mv1*mv2.T/norm(mv1)/norm(mv2))

def MatRotVec1toVec2(v1,v2):
    """return the rotation matrix that rotates from
direction 1 (vector 1, not necessarily normalized) to
direction 2."""
    n=mXpt(v1)*v2.T
    #Euler's Angle/Axis rotation form
    b=betav1v2(v1,v2)
    #Rodrigues' rotation formula
    #crossproduct matrix
    return MatRotAngleAxis(b,n)

def mXpt(n):
    if type(n)!=type(matrix([1,2])):
        raise Exception("Type Error: Not a matrix")
    if len(n.tolist())>1:
        n=n.T
    if len(n.tolist())>1:
        raise Exception("Type Error: Not a vector matrix")
      
    return matrix([[      0, -n[0,2],  n[0,1] ],
              [ n[0,2],       0,      -n[0,0] ],
              [-n[0,1],  n[0,0],           0] ])

def MatRotAngleAxis(b,n):
    n=n/norm(n)
    n=mXpt(n)
    #Eye
    I=matrix("""[1. 0 0;
                 0  1 0;
                 0  0 1]""")
    return I+sin(b)*n+(1-cos(b))*n*n

def rotDxVec1toVec2(dx,dy,dz,v1,v2):
    DX=rotdx(dx,dy,dz,MatRotVec1toVec2(v1,v2))
    return [DX[0].tolist()[0],
            DX[1].tolist()[0],
            DX[2].tolist()[0]]

def integrateToCircular(dx2,dy2,dz2,
                        dx1,dy1,dz1,h,m):
    """insert dX2 to dX1 and circle up.
m2,n2 is the starting and ending segment of dX2
m is the point you want to insert to dX1
h is the hinge on dX1.
dX1 circle will be cut open at the starting point of m,
segments dX[h:m] will be rotated to creat an opening for dX2[m2:n2]
Therefore, 0<=h<<m<=dX[-1].index()"""
    if not( h>=0 and m>h and len(dx1)>=m ):
        raise StandardError("parameter error! at integrateToCircular")
    #inserted vector
    dXmain=matrix([dx1,dy1,dz1])
    print "dXmain\n",dXmain
    dXins=matrix([dx2,dy2,dz2])
    print "dXins\n",dXins
    vmain=matrix([dx2x(dx1[h:m])[-1],
            dx2x(dy1[h:m])[-1],
            dx2x(dz1[h:m])[-1]])
    print "vmain",vmain
    
    vins=matrix([dx2x(dx2)[-1],
          dx2x(dy2)[-1],
          dx2x(dz2)[-1]])
    print "vins",vins
    #insertion vector
    l=norm(matrix(vmain))
    s=norm(matrix(vins)) #length of insertion vector
    print "l,s-->",l,s

    if l<s/2:
        raise StandardError("parameter error! The circle has too small segment selected"+ \
       "Which is in sufficent to accomodate the insertion chain [at integrateToCircular]")
    
    beta=acos((2*l*l-s*s)/(2*l*l))
    print "beta",beta/pi*180
    n=mXpt(vmain)*vins.T
    n=n.T
    print "axis rot:",n
    rm=MatRotAngleAxis(beta,n)
    print "rotation matrix of mainChain\n",rm
    print "rotation sample:[1,0,0]\n",rm*matrix([1,0,0]).T
    vmain2=rm*vmain.T
    vmain2=vmain2.T
    print "vmain2",vmain2
    print "before rotation dXmain[:,h:m]\n",dXmain[:,h:m]
    dXmain[:,h:m]=rm*dXmain[:,h:m]
    print "rotated dXmain[:,h:m]\n",dXmain[:,h:m]
    dXins=MatRotVec1toVec2(vins,vmain-vmain2)*dXins
    print "dXins rotated\n",dXins
    dXfinal=concatenate((concatenate((dXmain[:,:m],dXins),1),
                        dXmain[:,m:]),1)
    return [dXfinal[0].tolist()[0],
            dXfinal[1].tolist()[0],
            dXfinal[2].tolist()[0]]

def main_case1():
    #case 1. 2* 4seg/res site, 192 seg circle =200*10.5bp/seg=2.1kbase
    dx1,dy1,dz1=readJmol("main.txt",1)
    dx2,dy2,dz2=readJmol("res.txt",1)
    dx1,dy1,dz1=integrateToCircular(dx2,dy2,dz2,
             dx1,dy1,dz1,1,103)
    dx1,dy1,dz1=integrateToCircular(dx2,dy2,dz2,
             dx1,dy1,dz1,1,7)
    write_Jmol_file_dx("main.txt",dx1,dy1,dz1)

makeCircle_File("main.txt",192)
main()
