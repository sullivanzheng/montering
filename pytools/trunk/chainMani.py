from ConvertFor2Jmol import readVMCfile
from ConvertJmol2For import writeVMCfile_from_dxdydz
from math import *

if __name__=="__main__":
    print "input filename"
    filename=raw_input()
    dX=readVMCfile(filename)
    print "[File \"%s\" read. Number of segments: %d]"%(filename,len(dX[0]))
    
    print """========MAIN MENU=========
Select function:
1.Elongate the chain with straight segments on X
2.Cut and paste
3.Cut chain
4.Ligate 2 chains
5.Enlarge by interpolation
6.Normalize the segments
7.Shrink by delete every other segment
?
"""
    sel=int(raw_input())

    #Elongation
    if sel==1:
        n=int(raw_input("Number of segments to elongate:"))
        seglen=float(raw_input("segment length:"))
        elong=[[seglen]*n,[0.0]*n,[0.0]*n]
        dX=[dX[i]+elong[i] for i in range(3)]
        print "Number of segments after elongation:%d"%len(dX[0])
    #Cut and paste
    elif sel==2:
        print """Select:
1.Cut from left and paste to the right
2.Cut from right and paste to the left
?"""
        direction=int(raw_input())
        print "Number of segments to cut:"
        n=int(raw_input())
        if direction==1:
            dX=[j[n:]+j[:n] for j in dX]
        elif direction==2:
            n=len(dX[0])-n
            dX=[j[n:]+j[:n] for j in dX]
        else:
            print "invalid direction"
    #Cut the chain
    elif sel==3:
        m=int(raw_input("Cut at [1]left or [2]right:"))
        n=int(raw_input("Number of segments to cut:"))
        if m==1:
            print "Cut on the left side"
            dX=[i[n:] for i in dX]
        elif m==2:
            print "Cut on the right side"
            dX=[i[:-n] for i in dX]
        print "Number of segments after cut:%d"%len(dX[0])
    elif sel==4:
        file2=raw_input("Select the chain file you want to paste to current chain:")
        dX2=readVMCfile(file2)
        print "[File \"%s\" read. Number of segments: %d]"%(file2,len(dX2[0]))
        m=int(raw_input("load from ?th segment (1st segment is 0):"))
        n=int(raw_input("load to ?th(excluded) segment:"))
        side=int(raw_input("ligate to the left[1] or right[2]:"))
        if side==1:
            dX=[dX2[j][m:n]+dX[j] for j in range(3)]
        elif side==2:
            dX=[dX[j]+dX2[j][m:n] for j in range(3)]
    elif sel==5:
        zoom=int(raw_input("Enter the points to interpolate per segment:"))
        dX_temp=[]
        for j in dX:
            temp=[]
            for k in j:
                temp+=[k]*zoom
            dX_temp.append(temp)
        dX=dX_temp
    elif sel==6:
        segl=float(raw_input("Enter the segment length you need to normalized with:"))
        for i in range(len(dX[0])):
            segl_old=sqrt(dX[0][i]**2+dX[1][i]**2+dX[2][i]**2)
            dX[0][i],dX[1][i],dX[2][i]= [dX[0][i]/segl_old*segl,
                                         dX[1][i]/segl_old*segl,
                                         dX[2][i]/segl_old*segl]
    elif sel==7:
        x=[dX[0][i] for i in range(len(dX[0])) if i%2==0]
        y=[dX[1][i] for i in range(len(dX[1])) if i%2==0]
        z=[dX[2][i] for i in range(len(dX[2])) if i%2==0]
        dX=[x,y,z]
    writeVMCfile_from_dxdydz(filename,dX[0],dX[1],dX[2])
    
    
