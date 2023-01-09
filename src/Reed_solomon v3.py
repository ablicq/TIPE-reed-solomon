# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 19:30:48 2017

@author: Aurelien
"""

import numpy as np
import random as rdm



#15 -> 76 : fonctions pour manipuler les octets

table_exp = []
table_log = []
table_inv = []

def gen_tables():
    #charge les fichiers contenant les tables log exp et inv dans des listes
    f = open("table_log.txt",'r')
    for ligne in f:
        ligne = ligne[:-1]
        table_log.append(np.uint8(ligne))
    f.close()
    f = open("table_exp.txt",'r')
    for ligne in f:
        ligne = ligne[:-1]
        table_exp.append(np.uint8(ligne))
    f.close()
    f = open("table_inv.txt",'r')
    for ligne in f:
        ligne = ligne[:-1]
        table_inv.append(np.uint8(ligne))
    f.close()

gen_tables()

#définitions de exp, log et inv à partir des tables
def exp(x):
    return table_exp[x-1]

def log(x):
    return table_log[x-1]

def inv(x):
    return table_inv[x-1]


#somme et produit dans GF(256)
def somme_oct(x1,x2):
    return np.bitwise_xor(x1,x2)

def pdt_oct(x1,x2):
    if x1 == 0 or x2 == 0:
        return np.uint8(0)
    else:
        return exp((int(log(x1))+int(log(x2))) % 255)


def puiss(a,b):
    #exponentiation rapide dans GF(256)
    if b == 0:
        return 1
    if b == 1:
        return a
    x = puiss(a,b//2)
    if b%2 == 0:
        return pdt_oct(x,x)
    return pdt_oct(a,pdt_oct(x,x))

def pdt_itere(i,a):
    #a+...+a , i fois dans GF(256), un corps de caractéristique 2 (ie a+a = 0)
    if i%2 == 0:
        return 0
    else:
        return a

#81 -> 186 fonctions de manipulation des messages (polynomes sur GF(256))

def deg(P):
    #calcul le degre du polynome P
    d = len(P)-1
    while d >= 0 and P[len(P)-1-d] == 0:
        d -= 1
    return d

def emonde(P):
    #retire tous les zeros superflus en tête des polynomes
    d = deg(P)
    if d >= 0:
        return P[-deg(P)-1:]
    else:
        return np.array([0])

def dom(P):
    #coefficient dominant de P
    return P[-deg(P)-1]

def X_puiss(n):
    #renvoie le monome X^n
    P = np.zeros(n+1,np.uint8)
    P[0] = 1
    return P

def eval_pol(P,x):
    #calcule P(x) par la methode de Horner
    P=emonde(P)
    y = np.array(0,np.uint8)
    for i in range(len(P)):
        y = somme_oct(P[i], pdt_oct(x,y))
    return y

def derive(P):
    #dérivée formelle du polynome P
    P = emonde(P)
    Pprime = np.zeros_like(P)
    for i in range(1,len(P)):
        Pprime[-i] = pdt_itere(i,P[-i-1])
    return emonde(Pprime)

def somme_pol(P,Q):
    #somme terme à terme les polynomes P et Q en faisant attention à la taille 
    #des polynomes
    P = emonde(P)
    Q = emonde(Q)
    l1 = len(P)
    l2 = len(Q)
    if l1 <= l2:
        R = np.zeros(l2, dtype=np.uint8)
        for k in range(l1):
            R[l2-k-1]=somme_oct(P[l1-k-1],Q[l2-k-1])
        for k in range(l1,l2):
            R[l2-k-1]=Q[l2-k-1]
        return emonde(R)
    else:
        return somme_pol(Q,P)

def pdt_par_scal(alpha,P):
    #renvoie alpha.P
    P = emonde(P)
    Q = np.zeros_like(P)
    for k in range(len(P)):
        Q[k] = pdt_oct(P[k],alpha)
    return emonde(Q)

def pdt_pol(P,Q):
    #calcule de P*Q par la definition du produit de polynomes
    P = emonde(P)
    Q = emonde(Q)
    l1 = len(P)
    l2 = len(Q)
    produit = np.zeros(l1+l2-1, np.uint8)
    for i in range(l1+l2-1):
        coeff = np.uint8(0)
        for j in range(i+1):
            if j < l1 and i-j < l2:
                coeff = somme_oct(coeff,pdt_oct(P[j],Q[i-j]))
        produit[i] = coeff
    return emonde(produit)

def reste_DE(P,Q):
    #renvoie le reste de la division euclidienne de P par Q (algorithme usuel)
    P = emonde(P)
    Q = emonde(Q)
    R = P
    while deg(R) >= deg(Q):
        h = deg(R) - deg(Q)
        A = np.zeros(h+1, np.uint8)
        A[0] = pdt_oct(dom(R), inv(dom(Q)))
        R = somme_pol(R,pdt_pol(Q,A))
    return emonde(R)

def quotient_DE(P,Q):
    #renvoie le quotient de la division euclidienne de P par Q
    P = emonde(P)
    Q = emonde(Q)
    R = P
    quotient = np.array([0])
    while deg(R) >= deg(Q):
        h = deg(R) - deg(Q)
        A = np.zeros(h+1, np.uint8)
        A[0] = pdt_oct(dom(R), inv(dom(Q)))
        quotient = somme_pol(quotient,A)
        R = somme_pol(R,pdt_pol(Q,A))
    return emonde(quotient)

#192 -> 206 : codage d'un mot

def calcul_g():
    #renvoie le polynome generateur du code de Reed-Solomon
    g = np.array([1],np.uint8)
    for k in range(1,17):
        g = pdt_pol(g, np.array([1,exp(k)]))
    return g

g = calcul_g()

def code_RS(M):
    #renvoie le mot de code associe à M
    P = X_puiss(16)
    Q = pdt_pol(M,P)
    CK = reste_DE(Q,g)
    return somme_pol(Q,CK)


#211 -> 255 : decodage d'un mot

def syndromes(R):
    #calcul des syndromes
    S = []
    for k in range(1,17):
        S.append(eval_pol(R,exp(k)))
    return S

def algorithme_euclidien(S):
    #determination de sigma et omega par l'algorithme d'euclide etendu
    s = X_puiss(16)
    t = np.flipud(np.array(S))
    A = np.array([[np.array([1]),np.array([0])],[np.array([0]),np.array([1])]])
    while deg(t)>= 8:
        Q = quotient_DE(s,t)
        s,t = t,somme_pol(s,pdt_pol(Q,t))
        A1 = A[1,0]
        A2 = A[1,1]
        A3 = somme_pol(A[0,0],pdt_pol(Q,A[1,0]))
        A4 = somme_pol(A[0,1],pdt_pol(Q,A[1,1]))
        A = np.array([[A1,A2],[A3,A4]])
    delta = A[1,1][-1]
    sigma = pdt_par_scal(inv(delta),A[1,1])
    omega = pdt_par_scal(inv(delta),t)
    return sigma,omega

def chien_search(sigma):
    #recherhce des racines de sigma dans GF(256) et donc des position des erreurs
    r = []
    for k in range(1,256):
        if eval_pol(sigma,k) == 0:
            r.append(inv(k))
    return r

def forney(sigma,omega,r):
    #utilisation de l'algorithme de forney pour déterminer la valeur des erreurs
    y = []
    sigmap = derive(sigma)
    for i in range(len(r)):
        a = eval_pol(omega,inv(r[i]))
        b = eval_pol(sigmap,inv(r[i]))
        y.append(pdt_oct(a,inv(b)))
    return y

def decode_RS(R):
    #correction du mot R a partir des position et des valeurs des erreus trouvees
    S = syndromes(R)
    sigma,omega = algorithme_euclidien(S)
    r = chien_search(sigma)
    y = forney(sigma,omega,r)
    E = np.array([0])
    for i in range(len(r)):
        coeff = pdt_par_scal(y[i],X_puiss(log(r[i])))
        E = somme_pol(E, coeff)
    C = somme_pol(R,E)
    return C[:-16]



def canal_bin_sym(M,p):
    #a une probabilite p de modifier chaque bit de M
    N = M.copy()
    for k in range(len(M)):
        for l in range(8):
            r = rdm.random()
            if r < p:
                N[k] = somme_oct(N[k],exp(l))
    return N

def gen_message(n):
    #genere un message aleatoire
    M = np.zeros(n, np.uint8)
    for k in range(n):
        M[k]=rdm.randint(0,255)
    return M


#289 -> 336 : procedures de test

def nb_differences(M,N):
    #compte le nombre de composantes en lesquels les messages M et N different
    S = somme_pol(M,N)
    n=0
    for k in range(len(S)):
        if S[k] != 0:
            n+=1
    return n

def test():
    #procedure de test de la correction d'erreur en dessous de 8 composantes
    M = gen_message(239)
    C = code_RS(M)
    R = canal_bin_sym(C, 0.00306)
    n1 = nb_differences(C,R)
    Mp = decode_RS(R)
    n2 = nb_differences(M,Mp)
    print("\nerreurs introduites:",n1)
    print("erreurs apres decodage:",n2)

def test_erreurs(M,n):
    #procedure de test de la corrction d'erreur au dela de 9 composantes
    C = code_RS(M)
    corrige = 0
    partiel = 0
    ajout = 0
    i = 0
    for k in range(n):
        R = canal_bin_sym(C,0.00613)
        E = somme_pol(R,C)
        Cp = decode_RS(R)
        Ep = somme_pol(C,Cp)
        j = nb_differences(E,np.array([0]))
        if j > 8:
            i+=1
            if list(Ep)==[0]:
                corrige += 1
            part = 0
            for k in range(len(E)):
                if E[k] != 0 and k < len(Ep) and Ep[k] == 0:
                    part = 1
            partiel += part
            aj = 0
            for k in range(len(Ep)):
                if Ep[k] != 0 and k < len(E) and E[k] == 0:
                    aj = 1
            ajout += aj
    return corrige/i, partiel/i, ajout/i

if __name__ == "__main__":
    test()