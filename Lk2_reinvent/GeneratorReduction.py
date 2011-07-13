import os
def generatorplot(g,mfile,plotname):
    g=g[:-1]
    fp=file("xyz.txt")#file("rotated.txt")
    L1=111
    L2=111
    X=[ [float(i) for i in l.split()]+[1] for l in fp.readlines() ]
    fp.close()

    g=[(int(i.split()[4]),float(i.split()[5]),i.split()[0][2:]) for i in g.split('\n')]
    X_=X #keep an original copy of X for annotation.

    LL1=L1
    LL2=L2
    ra=0.95
    for i in range(len(g)-1,-1,-1):
        print g[i]
        if g[i][0]>L1-1:
            x0=X[g[i][0]  ]
            x1=X[g[i][0]+1]
            xx1=[x0[j]+ (x1[j]-x0[j])*   g[i][1] *ra for j in range(3)]
            xx2=[x1[j]- (x1[j]-x0[j])*(1-g[i][1])*ra for j in range(3)]
            print xx1
            print xx2
            X=X[:g[i][0]+1]+[xx1+[0],xx2+[2]]+X[g[i][0]+1:]
            LL2+=2
        else:
            x0=X[g[i][0]  ]
            x1=X[g[i][0]+1]
            xx1=[x0[j]+ (x1[j]-x0[j])*   g[i][1] *ra for j in range(3)]
            xx2=[x1[j]- (x1[j]-x0[j])*(1-g[i][1])*ra for j in range(3)]
            print xx1
            print xx2
            X=X[:g[i][0]+1]+[xx1+[0],xx2+[2]]+X[g[i][0]+1:]
            LL1+=2
            
    fp=file(mfile,'a')
    print >>fp,"close;figure;set(gcf,'Position',[0,50,1000,1000]);set(gca,'Position',[0.05,0.05,0.9,0.9]);"
    for i in range(len(X)-1):
        if i==LL1:
            continue
        if X[i][3]>0:
            if i>LL1:
                if X[i+1][3]==0:
                    print >>fp,"line([%f,%f],[%f,%f],'Color','%s','LineWidth',5)"%(X[i][0],X[i+1][0],X[i][1],X[i+1][1],'g')
                elif X[i][3]==2:
                    print >>fp,"line([%f,%f],[%f,%f],'Color','%s','LineWidth',5)"%(X[i][0],X[i+1][0],X[i][1],X[i+1][1],'g')
                else:
                    print >>fp,"line([%f,%f],[%f,%f],'Color','%s','LineWidth',2)"%(X[i][0],X[i+1][0],X[i][1],X[i+1][1],'r')
            else:
                if X[i+1][3]==0:
                    print >>fp,"line([%f,%f],[%f,%f],'Color','%s','LineWidth',5)"%(X[i][0],X[i+1][0],X[i][1],X[i+1][1],'m')
                elif X[i][3]==2:
                    print >>fp,"line([%f,%f],[%f,%f],'Color','%s','LineWidth',5)"%(X[i][0],X[i+1][0],X[i][1],X[i+1][1],'m')
                else:
                    print >>fp,"line([%f,%f],[%f,%f],'Color','%s','LineWidth',2)"%(X[i][0],X[i+1][0],X[i][1],X[i+1][1],'b')
    for gi in g:
           print >>fp,"text(%f,%f,'%s');"%(X_[gi[0]][0],X_[gi[0]][1],gi[2])
    print >>fp,"saveas(gcf,'%s','bmp');"%plotname
                
    fp.close()

if __name__=='__main__':
    mfile='generatorplot.m'
    if os.path.exists(mfile):
        os.remove(mfile)
        
    log='generators.txt'
    fp=file(log)
    recording=0
    s=''
    t=0

    for i in fp.xreadlines():
        if i.find('###')==0:
            recording=0
            t+=1
            plotfname='plot_in_directory_synapsis%02d_.bmp'%t
            generatorplot(s,mfile,plotfname)
            s=''
        print recording and "*" or " ",i[:-1]
        if recording==1:
            s+=i
        if i.find('---Gen')==0:
            recording=1
    
    
