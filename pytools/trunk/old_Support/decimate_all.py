from os import *

def decimate_dir(dirname):
    l=[int(i) for i in listdir(dirname) if i.isdigit()]
    print 'Total number of snapshot files;%d'%len(l)
    l.sort()
    del l[9::10]
    print 'Deleting %d files'%len(l)
    ##l=range(10000,600000001,10000)
    ##del l[9::10]
    ##print l[:20]

    for i in l:
        try:
            remove(path.join(dirname,str(i)))
        except WindowsError:
            print '.'
            pass

if __name__=='__main__':
    r,d,f=walk('.').next()
    for j in d:
        print 'deleting in directory:%s'%j
        decimate_dir(j)

    
