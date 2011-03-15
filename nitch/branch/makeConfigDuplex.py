from vec3d import makeCircle_File  #import makeCircleFile
from ConvertJmol2For import convertJ2F
import os,shutil
import math

def prepare_a_JF_simualtion(ds,ss,SPECIAL_ANGLE,seeding=1,shift=4):
    if ds % 2 != 0:
        print "ds should be an even number!"
        return
    dsg=72.0
    ssg=1.007078
    hinge=1.9862
    l=[hinge]+[dsg]*(ds-1)+[hinge]+[ssg]*(ss-1)
    l=l*2
    l=l[shift:]+l[:shift]
    print l
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
            lines[i]='totsegnum %d\n' % (( ds + ss )*2)
            break
    fp.close()

    fp=file('_config.cfg','w')
    fp.writelines(lines)
    fp.close()

    makeCircle_File('ligate', (ds + ss)*2 )
    convertJ2F('ligate')
    print "Made a circle of %d segments"% (ds+ss) * 2

    print "Config files and coordinate files successfully composed."

def do_a_JF_simulation(fname):
    os.mkdir(fname)
    shutil.copy('_config.cfg',os.path.join(fname,'_config.cfg'))
    shutil.copy('_rigid.cfg',os.path.join(fname,'_rigid.cfg'))
    shutil.copy(r'.\release\nitch.exe', os.path.join(fname,'nitch.exe'))
    shutil.copy('ligate.vmc',os.path.join(fname,'ligate.vmc'))
    print 'start "%s" /D"%s" /LOW %s'%(fname,os.path.join('.',fname),'nitch.exe')
    os.system('start "%s" /D"%s" /LOW %s'%(fname,os.path.join('.',fname),'nitch.exe'))
    
def get_JF(fname):
    f=file(os.path.join(fname,'ligate_log.txt'))
    s=f.readlines()[-1]
    st=s.find('Ligation P:')+len('Ligation P:')
    se=s.find('===',st)
    return float(s[st:se])
    f.close()
#------------------------------------
#ds=18, ss varies
#hinge 17xds hinge ss-1*ss
#0 1-8 9-17
#special angle is 9 (0+18/2).

def run():
    for i in range(12,34,3):
        for j in range(1):
            prepare_a_JF_simualtion(18,i,9,j)
            do_a_JF_simulation(r'.\__project\20110303\Dup_DS%2dSS%2drun%d'%(18,i,j))
        
def report():    
    r=[]
    for i in range(12,34,3):
        print i,
        for j in range(1):
            print math.log(get_JF(r'.\__project\20110303\Dup_DS%2dSS%2drun%d'%(18,i,j))),
        print 
#run()
report()
