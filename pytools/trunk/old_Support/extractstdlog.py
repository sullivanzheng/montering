f=file("stdout.log",'r')
l=f.readlines()
f.close()
# step,Gs2,Ls2, Gsst, Lsst,Gx,Lx,Gxst,Lxst
s=""

def process(s):
    t=s.split()
    step,Gs2,Ls2, Gsst, Lsst,Gx,Lx,Gxst,Lxst,e=[None]*10
    for j in range(len(t)):
        if t[j].find("Angle_Adapter::Finished")!=-1:
            step=int(t[j+2])
        if t[j].find("G<LEN^2>")!=-1:
            Gs2=float(t[j+1])
            Gsst=float(t[j+3])
            Ls2=float(t[j+5])
            Lsst=float(t[j+7])
        if t[j].find("G<x>")!=-1:
            Gx=float(t[j+1])
            Gxst=float(t[j+3])
            Lx=float(t[j+5])
            Lxst=float(t[j+7])
        if t[j].find("Enow")!=-1:
            e=float(t[j+1])
    return (step,Gs2,Gsst,Ls2,Lsst,Gx,Gxst,Lx,Lxst,e)

r=[]
for i in l:
    if i.find("Angle_")!=-1:
        r.append(process(s))
        s=i
    else:
        s=s+i
f=file("out.txt",'w')
r=[str(i)[1:-1]+'\n' for i in r]
f.writelines(r)
f.close()

