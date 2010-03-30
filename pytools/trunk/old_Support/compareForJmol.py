ini=file('ini.txt').readlines()[1:193]
ini=[float(i.split()[2]) for i in ini]
loosed=[i.split() for i in file('loosed.vmc').readlines()[0:39]]
loosed=[float(i) for j in loosed for i in j ]

print len(ini)
print len(loosed)
for i in range(191):
    if (abs(ini[i+1]-ini[i]-loosed[i])>1e-6):
        print 'ini-> %f, %f<-loosed  (%d)'% \
              (ini[i+1]-ini[i],loosed[i],i+1)
    

            
