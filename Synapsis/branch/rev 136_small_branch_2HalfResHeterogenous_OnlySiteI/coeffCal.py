from math import *
PI=3.14159
DELTA_TW_K=0.0
def G_t(Lk, L_bp,kink_num):
    C = 2.4e-28 / (1.3806503e-23 * 300.) / 3.4e-10
    E = + 2 * PI * PI * C / L_bp * \
        +(Lk - L_bp / 10.5 + DELTA_TW_K * kink_num)**2
    return E;

def P_t_over_2t(L_bp,kinkNum):
    Lk=L_bp/10.5
    P=0.
    C = 2.4e-28 / (1.3806503e-23 * 300.) / 3.4e-10
    denominator=sqrt(L_bp/(2*PI*C))
    for i in range(int(Lk)-5,int(Lk)+5):
        P+=exp(-G_t(i,L_bp,kinkNum))
    return P/denominator

def noTorCoef(r,theta):
    return 10/6.022*1.5/(PI*(0.34*r)**3.0)/(1-cos(theta/180.0*PI))

def ECondP(r,theta,totR,segnum):
    return (r/float(totR))**(3.0/segnum)*(theta/180.)**(2.0/segnum)

def listCoef(r,theta,totR,segnum,DELTA_TW_K=0):
    print 'NoTorCoef=',noTorCoef(r,theta)
    print 'Density=',P_t_over_2t(totR,DELTA_TW_K)
    print 'Expected CondP=',ECondP(r,theta,totR,segnum)

listCoef(4.2,5.7,84,15)
