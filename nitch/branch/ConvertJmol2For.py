import os

def writeVMCfile_from_dxdydz(filename,dx,dy,dz):
    f=file(filename,'w')
    j=1
    for i in dx:
        f.write("%12.6f"%i)
        if j!=5 and dx.index(i)!=len(dx)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1

    if j!=1:f.write("\n")

    j=1
    for i in dy:
        f.write("%12.6f"%i)
        if j!=5 and dy.index(i)!=len(dy)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1
    if j!=1:f.write("\n")

    j=1
    for i in dz:
        f.write("%12.6f"%i)
        if j!=5 and dz.index(i)!=len(dz)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1
    if j!=1:f.write("\n")
    f.write("1111111")

    f.close()

def writeVMCfile_from_xyz(filename,x,y,z):
    dx=[x[i+1]-x[i] for i in range(len(x)-1)]
    dy=[y[i+1]-y[i] for i in range(len(y)-1)]
    dz=[z[i+1]-z[i] for i in range(len(z)-1)]    
    f=file(filename,'w')
    j=1
    for i in dx:
        f.write("%12.6f"%i)
        if j!=5 and dx.index(i)!=len(dx)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1

    if j!=1:f.write("\n")

    j=1
    for i in dy:
        f.write("%12.6f"%i)
        if j!=5 and dy.index(i)!=len(dy)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1
    if j!=1:f.write("\n")

    j=1
    for i in dz:
        f.write("%12.6f"%i)
        if j!=5 and dz.index(i)!=len(dz)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1
    if j!=1:f.write("\n")
    f.write("1111111")

    f.close()
    
def convertJ2F(name):
    f=file(name,'r')
    l=f.readlines()
    j=int(l[0].split()[0])
    print j
    l=[i.split() for i in l[1:j+1]]
    x=[float(i[2]) for i in l]
    dx=[x[i+1]-x[i] for i in range(len(x)-1)]
    y=[float(i[3]) for i in l]
    dy=[y[i+1]-y[i] for i in range(len(y)-1)]
    z=[float(i[4]) for i in l]
    dz=[z[i+1]-z[i] for i in range(len(z)-1)]
    f.close()
    writeVMCfile_from_dxdydz(name+'.vmc',dx,dy,dz)
    print "DONE"
    
if __name__=="__main__":
    name=''
    while name.strip()=='' or not os.path.isfile(name):
        print "Name the input file"
        name=raw_input()
    convertJ2F(name)

