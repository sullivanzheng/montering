def readVMCfile(filename):
#Read all coordinates in the file and return the X containing all coordinates X
#expressed by a 3-tuple of lists of x coord., y coor. and z coor..
    lines=file(filename).readlines()[:-1]
    dx=[float(j) for l in lines for j in l.split()]
    if len(dx)%3 !=0:
        print "file error, the x y z do not have the same number of coordinates."
        return []
    dX=[dx[:len(dx)/3],dx[len(dx)/3:2*len(dx)/3],dx[2*len(dx)/3:]]
    return dX
    
def write_Jmol_file(filename,x,y,z):
    num_of_atoms=len(x)
    fp=file(filename,'w')
    print >>fp,"    %d   In the main body"%num_of_atoms
    num=range(1,num_of_atoms+1)
    atom=["c"]*3+["C"]*(num_of_atoms-3)
    typeatom=["1"]*num_of_atoms
    pos1=["2"]+[str(i) for i in range(1,num_of_atoms)]
    pos2=[" "]+[str(i) for i in range(3,num_of_atoms+1)]+[" "]
    for i in range(num_of_atoms):
        #print len(num),len(atom),len(x),len(y),len(z)
        print >>fp,"%6d%4s%12.6f%12.6f%12.6f%6s%6s%6s"%\
          (num[i],atom[i],x[i],y[i],z[i],typeatom[i],pos1[i],pos2[i])
    fp.close()

def integrate(l):
    """Numerical integration of a list.
x[i]=sum(l[:i+1] i.e l(0~i)."""
    x=[l[0]]
    for j in l[1:]:
        x.append(x[-1]+j)
    return x

def convertf2j(filename):
    dX=readVMCfile(filename)
    X=[[0]+integrate(i) for i in dX]
    write_Jmol_file(filename+'.txt',X[0],X[1],X[2])
                   
if __name__=="__main__":
    print "input filename"
    filename=raw_input()
    #filename='init_51.data'
    dX=readVMCfile(filename)
    X=[[0]+integrate(i) for i in dX]
    write_Jmol_file(filename+'.txt',X[0],X[1],X[2])
 
