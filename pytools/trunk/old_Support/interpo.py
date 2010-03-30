
def conv(filename,outname,num):
    print filename,outname,num
    ln=file(filename).readlines()
    loosed=[i.split() for i in file(filename).readlines()[0:50]]
    dx=[float(i) for j in loosed for i in j ]
    loosed=[i.split() for i in file(filename).readlines()[50:100]]
    dy=[float(i) for j in loosed for i in j ]
    loosed=[i.split() for i in file(filename).readlines()[100:150]]
    dz=[float(i) for j in loosed for i in j ]
    x=[[i]*num for i in dx]
    x=[i for j in x for i in j]
    y=[[i]*num for i in dy]
    y=[i for j in y for i in j]
    z=[[i]*num for i in dz]
    z=[i for j in z for i in j]
    f=file(outname,'w')
    j=1
    for i in x:
        f.write("%12.6f"%i)
        if j!=5 and x.index(i)!=len(x)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1

    if j!=1:
        f.write("\n")

    j=1
    for i in y:
        f.write("%12.6f"%i)
        if j!=5 and y.index(i)!=len(y)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1
    if j!=1:
        f.write("\n")

    j=1
    for i in z:
        f.write("%12.6f"%i)
        if j!=5 and z.index(i)!=len(z)-1:
            f.write(" ")
            j+=1
        else:
            f.write("\n")
            j=1
    if j!=1:
        f.write("\n")
    f.write("1111111")
    f.close()




filename=raw_input("input filename: ")
outname=raw_input("output filename: ")
num=int(raw_input("fold to expand: "))
conv(filename,outname,num)
