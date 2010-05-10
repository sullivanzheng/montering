import re

f=file('main_log.txt').read()
b=re.findall(r'Branch=[0-9]+',f)
b=[int(re.search(r'[0-9]',i).group()) for i in b]
hist=[0]*8
for i in b:
    hist[i-1]+=1
fo=file('branchHistory.txt','w')
print 'Total Samples = %d'%len(b)
for i,j in zip(range(1,9),hist):
    print 'Branch[%d]=%5.3f'%(i,float(j)/sum(hist))
for i in b:
    print >>fo,i
fo.close()
