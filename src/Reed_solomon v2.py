# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 21:38:51 2017

@author: Aurelien
"""

import numpy as np
import random as rdm
import time as tm

#un polynome irreductible primitif de degre 8 ie irred et tq X est gen de F2[X]/P8(X)
#ce sera le polynome générateur de nos 'octets'
P8 = np.array([1,0,0,0,1,1,1,0,1])

def mod2(A):
    #applique "mod 2" à tous les coeff de A
    for i in range(len(A)):
        A[i] = A[i] % 2
    return A

def deg(A):
    #renvoie le degré de A
    d = len(A)-1
    while d >= 0 and A[len(A)-1-d] == 0:
        d -= 1
    return d

def somme_pol(A,B):
    l1 = len(A)
    l2 = len(B)
    if l1 == l2:
        return mod2(A+B)
    if l1 < l2:
        h = l2 - l1
        for k in range(h):
            A = np.insert(A,0,0)
        return mod2(A+B)
    if l1 > l2:
        return somme_pol(B,A)

def mult_pol(A,B):
    #renvoie la multiplication des polynomes A et B
    l1 = len(A)
    l2 = len(B)
    produit = np.zeros(l1+l2-1, int)
    for i in range(l1+l2-1):
        #on utilise la def du produit de polynome ie le produit de Cauchy
        coeff = 0
        for j in range(i+1):
            if j < l1 and i-j < l2:
                coeff += A[j]*B[i-j]
        produit[i] = coeff
    return mod2(produit)

def reste_div_eucl(A,B):
    #renvoie le reste de la division euclidienne de A par B
    D = np.zeros(0, int)
    R = A
    while deg(R) >= deg(B):
        h = deg(R) - deg(B)
        #on créé le polynome X^h
        P = np.zeros(h+1, int)
        P[0] = 1
        #que l'on ajoute à D
        D = somme_pol(P,D)
        #et on met R à jour
        R = somme_pol(R,mult_pol(B,P))
    return R

def quot_div_eucl(A,B):
    #renvoie le reste de la division euclidienne de A par B
    D = np.zeros(0, int)
    R = A
    while deg(R) >= deg(B):
        h = deg(R) - deg(B)
        #on créé le polynome X^h
        P = np.zeros(h+1, int)
        P[0] = 1
        #que l'on ajoute à D
        D = somme_pol(P,D)
        #et on met R à jour
        R = somme_pol(R,mult_pol(B,P))
    return D

def mult_modP8(A,B):
    M = reste_div_eucl(mult_pol(A,B), P8)
    return M[-8:]

def inv_pol(A):
    #calcul de l'inverse de A grace à l'algorithme d'euclide etendu
    r1,r2=P8,A
    u1,u2 = np.array([0]), np.array([1])
    while deg(r2) >= 0:
        q = quot_div_eucl(r1,r2)
        r1,r2,u1,u2 = r2, somme_pol(r1, mult_pol(q,r2)), u2, somme_pol(u1, mult_pol(q,u2))
    return reste_div_eucl(u1,P8)

def convertit(n):
    a = []
    for k in range(8):
        a.append(n%2)
        n //= 2
    return np.flipud(np.array(a))
    
def convert(a):
    i = 0
    n = len(a)
    for k in range(n):
        i += a[k]*2**(n-k-1)
    return str(i)

def charge_inv():
    f = open("C:/Users/Aurelien/Desktop/Programme TIPE/table_inv.txt",'w')
    for k in range(255):
        f.write(str(convert(inv_pol(convertit(k+1))))+'\n')
    f.close()

#crée la liste des 255 elts non nuls de F2[X]/P8(X)
def calcul_L():
    L = [np.array([0,0,0,0,0,0,1,0])]
    for k in range(254):
        L.append(mult_modP8(L[k],np.array([1,0])))
    return L
L = calcul_L()


        
Lp = []

def charge_L():
    for k in range(255):
        Lp.append(convert(L[k]))

def charge_log():
    f = open("C:/Users/Aurelien/Desktop/Programme TIPE/table_log.txt",'w')
    for k in range(255):
        f.write(str(Lp.index(str(k+1))+1)+'\n')
    f.close()


def mod2_oct(M):
    for i in range(len(M)):
        M[i] = mod2(M[i])
    return M

def emonde(M):
    d = deg_oct(M)
    return M[-d-1:]

def deg_oct(A):
    #renvoie le degré de A
    d = len(A)-1
    while d >= 0 and (A[len(A)-1-d] == 0).all():
        d -= 1
    return d

def somme_pol_oct(A,B):
    l1 = len(A)
    l2 = len(B)
    if l1 == l2:
        return mod2_oct(A+B)
    if l1 < l2:
        h = l2 - l1
        for k in range(h):
            #insere un coefficient 0 en tete du ploynome
            A = np.insert(A,0,0,axis=0)
        return emonde(mod2_oct(A+B))
    if l1 > l2:
        return somme_pol_oct(B,A)

def mult_pol_oct(A,B):
    #renvoie la multiplication des polynomes A et B
    l1 = len(A)
    l2 = len(B)
    produit = np.zeros((l1+l2-1,8), int)
    for i in range(l1+l2-1):
        #on utilise la def du produit de polynome ie le produit de Cauchy
        coeff = np.zeros(8,int)
        for j in range(i+1):
            if j < l1 and i-j < l2:
                coeff = somme_pol(coeff, mult_modP8(A[j],B[i-j]))
        produit[i] = coeff
    return mod2_oct(produit)

def dom(A):
    return A[len(A)-deg_oct(A)-1]

def X_puiss(n):
    P = np.zeros((n+1,8),int)
    P[0] = np.array([0,0,0,0,0,0,0,1])
    return P

def reste_DE_oct(A,B):
    #renvoie le reste de la division euclidienne de A par B
#    D = np.zeros((0,8), int)
    R = A
    while deg_oct(R) >= deg_oct(B):
        h = deg_oct(R) - deg_oct(B)
        #on créé le polynome en X^h
        P = np.zeros((h+1,8), int)
        P[0] = mult_modP8(dom(R),inv_pol(dom(B)))
        #que l'on ajoute à D
#        D = somme_pol_oct(P,D)
        #et on met R à jour
        R = somme_pol_oct(R,mult_pol_oct(B,P))
    return R

def eval_pol_oct(M,x):
    y = np.array([0])
    for i in range(len(M)):
        y = somme_pol(M[i], mult_modP8(x,y))
    return y

def mult_oct_par_pol(P,M):
    N = np.zeros_like(M)
    for k in range(len(M)):
        N[k] = mult_modP8(P,M[k])
    return N

def renverse(M):
    M = emonde(M)
    return np.flipud(M)

#calcul g, le polynome generateur de notre code de Reed-Solomon
def calcul_g():
    g = np.array([[0,0,0,0,0,0,0,1]])
    for k in range(16):
        g = mult_pol_oct(g, np.array([[0,0,0,0,0,0,0,1],L[k]]))
    return g

g = calcul_g()






def code_RS(M):
    t1 = tm.clock()
    #P est le pol X^(n-k)
    P = X_puiss(16)
    #Q est M * X^(n-k)
    Q=mult_pol_oct(M,P)
    #CK est le controle de parite
    CK=reste_DE_oct(Q,g)
    print("codage : ",tm.clock()-t1)
    return somme_pol_oct(Q,CK)
    




def syndrome(R):
    S = [[0,0,0,0,0,0,0,0]]
    for i in range(16):
        S.append(eval_pol_oct(R,L[i]))
    return np.flipud(np.array(S))

def berlekamp_massey(R,S):
    sigma = X_puiss(0)
    omega = X_puiss(0)
    tau = X_puiss(0)
    gamma = np.zeros((0,8))
    D = 0
    B = 0
    for k in range(16):
        print(reste_DE_oct(mult_pol_oct(sigma, somme_pol_oct(S,X_puiss(0))),X_puiss(k+1))==omega)
        C = somme_pol_oct(reste_DE_oct(mult_pol_oct(sigma,somme_pol_oct(X_puiss(0),S)),X_puiss(k+2)),omega)
        print(C)
        if deg_oct(C) < 0:
            delta = np.array([0,0,0,0,0,0,0,0])
        else:
            delta = dom(C)
        print(delta)
        if (delta == 0).all() or D > (k+1)/2:
            tau = mult_pol_oct(tau,X_puiss(1))
            gamma = mult_pol_oct(gamma,X_puiss(1))
            print('A')
        elif D < (k+1)/2:
            d = inv_pol(delta)
            tau = mult_oct_par_pol(d,sigma)
            gamma = mult_oct_par_pol(d,omega)
            D = k+1-D
            B = 1-B
            print('B')
        elif B == 0:
            tau = mult_pol_oct(tau,X_puiss(1))
            gamma = mult_pol_oct(gamma,X_puiss(1))
            print('C')
        else:
            d = inv_pol(delta)
            tau = mult_oct_par_pol(d,sigma)
            gamma = mult_oct_par_pol(d,omega)
            B = 1-B
            print('D')
        sigma = somme_pol_oct(sigma,mult_pol_oct(np.array([delta,[0,0,0,0,0,0,0,0]]),tau))
        omega = somme_pol_oct(omega,mult_pol_oct(np.array([delta,[0,0,0,0,0,0,0,0]]),gamma))
    return sigma,omega

def chien_search(sigma):
    r = []
    for k in L:
        if deg(eval_pol_oct(sigma,k)) < 0:
            r.append(inv_pol(k))
    return r

def forney(r,omega):
    w = []
    omegat = renverse(omega)
    for a in r:
        P = a
        for k in r:
            if (k != a).any():
                P = mult_modP8(P,somme_pol(a,k))
        P = inv_pol(P)
        w.append(mult_modP8(P,eval_pol_oct(omegat,a)))
    return w

def idx(a,L):
    i = 0
    while (a != L[i]).all():
        i += 1
    return i
        
def decode_RS(R):
    t1 = tm.clock()
    S = syndrome(R)
    t2 = tm.clock()
    sigma,omega = berlekamp_massey(R,S)
    t3 = tm.clock()
    r = chien_search(sigma)
    t4 = tm.clock()
    w = forney(r, omega)
    t5 = tm.clock()
    l = []
    n = len(r)
    for k in range(n):
        while len(r[k]) != 8:
            r[k] = np.insert(r[k],0,0)
        l.append(idx(r[k],L))
    print(r,l,w)
    E = np.zeros_like(R)
    for k in range(n):
        E = somme_pol_oct(E,mult_oct_par_pol(w[k],X_puiss(l[k])))
    print(E)
    M = somme_pol_oct(R,E)
    t6 = tm.clock()
    print("syndromes:", t2-t1)
    print("berlekamp-massey:", t3-t2)
    print("chien search:", t4-t3)
    print("forney:", t5-t4)
    print("calcul de l'erreur:", t6 - t5)
    print("total decodage:", t6-t1)
    return M[:-16]








def canal_bin_sym(M):
    N = M.copy()
    for k in range(len(M)):
        for l in range(8):
            r = rdm.random()
            if r < 0.00245:
                N[k,l] = (N[k,l] + 1) % 2
    return N

def gen_message(n):
    M=np.zeros((n,8), int)
    for i in range(n):
        for j in range(8):
            M[i,j]=rdm.randint(0,1)
    return M

def nb_differences(M,N):
    n = 0
    for k in range(len(M)):
        if not (M[k] == N[k]).all():
            n += 1
    return n

def test():
    M = gen_message(239)
    C = code_RS(M)
    R = canal_bin_sym(C)
    E = R - C
    n1 = nb_differences(R,C)
    Mp= decode_RS(R)
    n2 =nb_differences(M,Mp)
    print((Mp==M).all(),n1,n2)
    return E