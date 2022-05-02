# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:15:43 2022

@author: ETUDIANT
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

r=0.45
K=10
a=0.5

Qn=0.5
En=1
Qp=0.5

e=0.4
m=0.4
E=1

def isoP_N():
        
    NL=np.linspace(0,200,1000)
    PPlot=(r/a)*(1 - (NL/K)) - (Qn*En*E)/a

    A= (m+Qp*(1-En)*E)/(e*a)
    plt.xlim([0,K])
    plt.ylim([0,3])
    plt.plot(NL,PPlot,'r',label='Isocline de P')
    plt.legend()
    plt.plot([A,A],[0,A],'g',label='Isocline de N')
    plt.legend()
    print(A)
    
    
def sys_schaefer(t,Y):
# calcule les membres de droite des deux équations différentielles
    N,P = Y
    dN= N * (r*(1-(N/K))-a*P - Qn*En*E)
    dP= P * (e*a*N - m -Qp*(1-En)*E)
    return [dN,dP]
    
def plot_traj(Yini):
# calcule et affiche la trajectoire obtenue en partant de la condition initiale Yini
    Nini,Pini = Yini
    tmax=500
    tL=np.linspace(0,tmax,100000)
    sol= solve_ivp(sys_schaefer,[0,tmax],[Nini,Pini],t_eval=tL)
    plt.plot(sol.y[0],sol.y[1],'b',label="Portrait de phase")
    plt.legend()
    plt.xlabel('N')
    plt.ylabel('P')
    plt.title('Trajectoire du système')

def plot_dynamique(Yini):
    Nini,Pini=Yini
    tmax=500
    tL=np.linspace(0,tmax,100000)
    sol= solve_ivp(sys_schaefer,[0,tmax],[Nini,Pini],t_eval=tL)
    plt.xlim([0,100])
    plt.plot(tL,sol.y[1],'r',label="Prédateurs")
    plt.plot(tL,sol.y[0],'b',label="Proies")
    plt.title('Evolution de la population avec le temps')
    plt.xlabel('temps')
    plt.ylabel('Population')
    plt.legend()
    
        
def equil(Yini):
    # calcule l'équilibre vers lequel on converge en partant de la condition initial Yini
    Nini,Pini=Yini
    tmax=500
    sol=solve_ivp(sys_schaefer,[0,tmax],[Nini,Pini])
    return[sol.y[0][-1] , sol.y[1][-1]]

def jacobienne(Y):
    # calcule la jacobienne du système d'équations différentielles au point Y (=[N,P])
    epsilon=1e-5
    N,P=Y
    FN= lambda N,P : N * (r*(1-(N/K))-a*P - Qn*En*E)
    FP= lambda N,P : P * (e*a*N - m -Qp*(1-En)*E)
    j11= (FN(N+epsilon,P) - FN(N,P)) / epsilon
    j12= (FN(N,P+epsilon) - FN(N,P)) / epsilon
    j21= (FP(N+epsilon,P) - FN(N,P)) / epsilon
    j22= (FP(N,P+epsilon) - FN(N,P)) / epsilon   
    Det= (j11*j22)-(j21*j12)
    print(Det)
    return np.array([[j11, j12],
                     [j21, j22]])

def stab(Y):
    # calcule les valeurs propres de la jacobienne au point Y (=[H,P])
    # permet de déteminer la stabilité des équilibre
    J= jacobienne(Y)
    return np.linalg.eigvals(J)

#évolution du max des valeurs propres en fonction de l'effort de pêche
# =============================================================================
# def DF(Y,EN):
#     epsilon=1e-5
#     N,P=Y
#     ENL=np.arange(0,EN,100)
#     FN= lambda N,P : N * (r*(1-(N/K))-a*P - Qn*ENL*E)
#     FP= lambda N,P : P * (e*a*N - m -Qp*(1-ENL)*E)
#         
#     j11= (FN(N+epsilon,P) - FN(N,P)) / epsilon
#     j12= (FN(N,P+epsilon) - FN(N,P)) / epsilon
#     j21= (FP(N+epsilon,P) - FN(N,P)) / epsilon
#     j22= (FP(N,P+epsilon) - FN(N,P)) / epsilon
#         
#     JL=DF(Y)
#     MVP=max(np.linalg.eigvals(DF))
#     plt.plot(ENL,MVP)
#     
#     
# =============================================================================
    
 

    


       


    
        
        
        
        
        
        
        