from knotsize import *
x,y,z=readJmol('deformed.txt')
f.close()
print knotlist(x,y)

