from os import *
l=[i for i in listdir('.') if ('.cpp' in i ) or ('.h' in i)]
print l
f=file('allcode.txt','w')
u=[];
for i in l:
    print >>f ,'======='+i
    f.write(file(i).read())
    print >>f,'\n------------\n\n\n\n'

f.close()
    
