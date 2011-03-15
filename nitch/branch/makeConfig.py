from vec3d import makeCircle_File  #import makeCircleFile
from ConvertJmol2For import convertJ2F
import os,shutil
import math
from math import exp,log,sin,cos,tan

def prepare_a_JF_simualtion(ds,ss,SPECIAL_ANGLE,BREAKANGLE=16,seeding=1,shift=4):
    if ds % 2 != 0:
        print "ds should be an even number!"
        return
    dsg=72.0
    ssg=1.007078
    hinge=1.9862
    l=[hinge]+[dsg]*(ds-1)+[hinge]+[ssg]*(ss-1)
    l=l[shift:]+l[:shift]
    fp=file('_rigid.cfg','w')
    print >>fp, "specify_rigidity"
    for i in range(len(l)):
        print >>fp, ">",i,l[i]
    print >>fp,"end_specify_rigidity"
    fp.close()

    fp=file('_config.cfg')
    lines=fp.readlines()
    for i in range(len(lines)):
        if lines[i].find('seeding')==0:
            lines[i]='seeding %d\n'% seeding
            break
    for i in range(len(lines)):
        if lines[i].find('SPECIAL_ANGLE')==0:
            lines[i]='SPECIAL_ANGLE %d\n'% ( SPECIAL_ANGLE - shift )
            break
    for i in range(len(lines)):
        if lines[i].find('totsegnum')==0:
            lines[i]='totsegnum %d\n' % ( ds + ss )
            break
    for i in range(len(lines)):
        if lines[i].find('BREAKANGLE')==0:
            lines[i]='BREAKANGLE %f  #%5.3f/180*PI \n' % ( float(BREAKANGLE)/180*3.14159,BREAKANGLE )
            break
    fp.close()

    fp=file('_config.cfg','w')
    fp.writelines(lines)
    fp.close()

    makeCircle_File('ligate', ds + ss )
    convertJ2F('ligate')
    print "Made a circle of %d segments"%(ds+ss)

    print "Config files and coordinate files successfully composed."

def do_a_JF_simulation(fname):
    os.mkdir(fname)
    shutil.copy('_config.cfg',os.path.join(fname,'_config.cfg'))
    shutil.copy('_rigid.cfg',os.path.join(fname,'_rigid.cfg'))
    shutil.copy(r'D:\_rch\GradProg\nitch\release\nitch.exe', os.path.join(fname,'nitch.exe'))
    shutil.copy('ligate.vmc',os.path.join(fname,'ligate.vmc'))
    print 'start "%s" /LOW /affinity 1 /D"%s" %s'%(fname,os.path.join('.',fname),'nitch.exe')
    os.system('start "%s" /LOW /affinity 1 /D"%s" %s'%(fname,os.path.join('.',fname),'nitch.exe'))
    
def get_JF(fname):
    f=file(os.path.join(fname,'ligate_log.txt'))
    s=f.readlines()[-1]
    st=len('>>>>>>>>>>>>>>> J factor:')
    f.close()
    try:
        s=float(s[st:])
    except ValueError:
        s=1
    return s

#------------------------------------
#ds=18, ss varies
#hinge 17xds hinge ss-1*ss
#0 1-8 9-17
#special angle is 9 (0+18/2).

def Shimada_Yamakawa(bp):
    pi=3.14159
    Na=6.022e23
    b=100e-9
    Lbp=0.34e-9
    L=Lbp*bp
    jf=4*pi**3 * b**3/(Na*L**6) * math.exp(-pi**2*b/L + 0.514*L/b)
    return jf

def getCoeff(dis,ang):
    PI=3.14159
    return 1/6.022141 * (3.0/ (4.0 * PI * (pow(0.34*dis,3))) ) *  2.0/(1-cos(ang)) * 1e4

def test():
    for i in range(10):
        prepare_a_JF_simualtion(100,0,200,360,12305981+i)
        do_a_JF_simulation(r'.\test%d'%i)
        
def testreport():
    for i in range(10):
        print log(get_JF(r'.\test%d'%i),10)

def theory():
    print "="*50
    print "Shimada Yamakawa equation gives JF"
    print " Jf=",Shimada_Yamakawa(100)," Log10(JF)=",log(Shimada_Yamakawa(100),10)
    print "Inside my program, P is linked with JF with coeffcient."
    print " Coef=%8.5e"%getCoeff(0.1,5.0/180.0*3.14159)
    print "Therefore the P should be (backed by S-Y)"
    print " P=",Shimada_Yamakawa(100)/getCoeff(0.1,5.0/180.0*3.14159),


def run():
    for breakangle in [6,10,14,18,22,26]:
        for i in range(12,34,3):
            for j in range(3):
                prepare_a_JF_simualtion(18,i,9,breakangle,j+10000)
                do_a_JF_simulation(r'.\Bkangle%2d_DS%2dSS%2d_run%d'%(breakangle,18,i,j))
        
def report(output):    
    r=[]
    fp=file(output,'w')
    print >>fp,"a={};"
    print >>fp,"ba=[];"
    print >>fp,"""experi=[12	9.8000 
        15	9.6000 
        18	9.3000 
        21	9.1000 
        24	8.6500 
        26	7.9000 
        27	7.8500 
        30	6.9000 
        33	5.9500];"""
    ii=0
    breakanglelist=[6,10,14,18,22,26]
    for breakangle in breakanglelist:
        ii=ii+1
        print >>fp,"ba(%d)=%f;"%(ii,breakangle)
        print >>fp,"a{%d}=["%ii
        for i in range(12,34,3):
            print >>fp,i,
            for j in range(3):
                print >>fp,-math.log(get_JF(r'.\Bkangle%2d_DS%2dSS%2d_run%d'%(breakangle,18,i,j))),
            print >>fp,' '
        print >>fp,"];"
    print >>fp, """figure;hold on;
color={'b','g','k','r','m','y'};
for i=1:6
    plot(a{i}(:,1),log(mean(exp(a{i}(:,2:end)),2)),[color{i} '+-'],'LineWidth',2)
end
plot(experi(:,1),experi(:,2),'k','LineWidth',3);"""
    print >>fp, "legend('%s','experi');"%"','".join([str(i) for i in breakanglelist])
    

#run()
report(r'D:\_rch\GradProg\nitch\__project\20110309\ana.m')
