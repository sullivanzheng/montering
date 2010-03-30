import os

def writeVMCfile(filename,x,y,z):
    f=file(filename,'w')
    j=1
    for i in x:
        f.write("%12.6f"%i)
        if j!=5 and x.index(i)!=len(x)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1

    if j!=1:f.write("\n")

    j=1
    for i in y:
        f.write("%12.6f"%i)
        if j!=5 and y.index(i)!=len(y)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1
    if j!=1:f.write("\n")

    j=1
    for i in z:
        f.write("%12.6f"%i)
        if j!=5 and z.index(i)!=len(z)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1
    if j!=1:f.write("\n")
    f.write("1111111")

    f.close()

def convert(name,outname=''):
    f=file(name,'r')
    l=f.readlines()
    j=int(l[0].split()[0])
    print j
    l=[i.split() for i in l[1:j+1]]
    x=[float(i[2]) for i in l]
    x=[x[i+1]-x[i] for i in range(len(x)-1)]
    y=[float(i[3]) for i in l]
    y=[y[i+1]-y[i] for i in range(len(y)-1)]
    z=[float(i[4]) for i in l]
    z=[z[i+1]-z[i] for i in range(len(z)-1)]
    f.close()
    if outname=='':
        outname==name+'.vmc'
    writeVMCfile(outname,x,y,z)
    print "DONE"
    
if __name__=="__main__":
    name=''
    while name.strip()=='' or not os.path.isfile(name):
        print "Name the input file"
        name=raw_input()
    convert(name)
        

