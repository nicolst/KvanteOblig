# -*- coding: utf-8 -*-
#Python-program som finner energi-egenverdier og tilhørende bølgefunksjoner
#for en endimensjonal harmonisk oscillator. Potensialet er kvadratisk over
#en bredde på ca 200 Å, og har verdien 1 eV i posisjon +- 100 Å.
#På hver side av dette er det et område med konstant
#potensial og bredde 400 Å. Utenfor dette igjen
#er potensialet uendelig. Partikkelen har masse lik elektronmassen.
#Problemet løses ved å diskretisere den 2. deriverte, som
#inngår i operatoren for kinetisk energi. Den resulterende Hamiltonmatrisen H
#diagonaliseres med bruk av numpy-funksjonen np.linalg.eigh(H), som regner ut
#og returnerer samtlige egenverdier og tilhørende egenvektorer.
#Programmet plotter sannsynlighetstettheten til de 10 laveste egentilstandene.
import matplotlib.pyplot as plt
import numpy as np

#m=elektronmassen, V0 = 1 eV
hbar=1.05E-34
m=9.11E-31
V0=1.6E-19
#N = 100 = antall posisjonsverdier i halve det harmoniske området
N=100
#dz = 1 Å = intervallbredde
dz=1E-10
#V = liste med potensialverdier
V = [V0]*4*N + [V0*((n-N)/(N*1.0))**2 for n in range(2*N+1)] + [V0]*4*N
#z = liste med posisjonsverdier (i nm)
#Med len(V) angis antall elementer i vektoren V
z = [dz*1E9*n for n in range(len(V))]
#d = liste med diagonalelementer i Hamiltonmatrisen H
d = [v + hbar**2/(m*dz**2) for v in V]
#e = verdi til ikke-diagonale elementer i H, dvs H(i,i+-1) = -e
e = - hbar**2/(2*m*dz**2)
#Initialisering av matrisen H: Legger inn verdi 0 i samtlige elementer
H = [[0]*(len(V)) for n in range(len(V))]
#Dobbel for-løkke som lager den tridiagonale matrisen H
for i in range(len(V)):
    for j in range(len(V)):
        if i==j:
            H[i][j]=d[i]
        if abs(i-j)==1:
            H[i][j]=e
#Finner w = egenverdiene og v = egenvektorene til matrisen H
w,v = np.linalg.eigh(H)
#evalues = liste med egenverdier i enheten eV, med nullpunkt for
#energien i bunnen av det harmoniske potensialet:
evalues = w/1.6E-19
#Hvis ønskelig skriver neste linje ut de 4 laveste energiegenverdiene
#print(evalues[0],evalues[1],evalues[2],evalues[3])
#psi0 = bølgefunksjonen til grunntilstanden, psi1 = 1. eksiterte tilstand osv
psi0 = v[:,0]
psi1 = v[:,1]
psi2 = v[:,2]
psi3 = v[:,3]
psi4 = v[:,4]
psi5 = v[:,5]
psi6 = v[:,6]
psi7 = v[:,7]
psi8 = v[:,8]
psi9 = v[:,9]
#Lager tabell VineV med potensialet i enheten eV
VineV = [x/1.6E-19 for x in V]
#Lager sannsynlighetstettheten til de 10 laveste tilstandene, med
#null-nivå forskjøvet med energiverdien, for å gi en bra figur
rho0 = np.abs(psi0)**2+evalues[0]
rho1 = np.abs(psi1)**2+evalues[1]
rho2 = np.abs(psi2)**2+evalues[2]
rho3 = np.abs(psi3)**2+evalues[3]
rho4 = np.abs(psi4)**2+evalues[4]
rho5 = np.abs(psi5)**2+evalues[5]
rho6 = np.abs(psi6)**2+evalues[6]
rho7 = np.abs(psi7)**2+evalues[7]
rho8 = np.abs(psi8)**2+evalues[8]
rho9 = np.abs(psi9)**2+evalues[9]
#Plotter sannsynlighetstettheten som tilsvarer de 10 laveste egenverdiene
#sammen med potensialet:
plt.figure('Bølgefunksjoner')
plt.plot(z,psi0,z,psi1,z,psi2,z,psi3,z,psi4)
plt.show()
plt.xlim(40,60)
plt.figure('Sannsynlighetstettheter')
plt.plot(z,rho0,z,rho1,z,rho2,z,rho3,z,rho4,z,rho5,z,rho6,z,rho7,z,rho8,z,rho9,z,VineV)
plt.show()
plt.title('Navn: Jon Andreas Støvneng',fontsize=20)
plt.xlabel('$z$ (nm)',fontsize=20)
plt.ylabel('$|\psi(z)|^2$',fontsize=20)
plt.xlim(40,60)
plt.ylim(0.0,0.4)