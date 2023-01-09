# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 15:45:37 2016

@author: Aurelien
"""

import numpy as np
import random as rdm



def code743(M):
    G = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]])
    A = G*M
    for x in A:
        x[0] %= 2
    return A

def perturbe(A):
    for x in A:
        r = rdm.random()
        if r < .05:
            x[0] = (x[0] + 1) % 2
    return A

def decode743(R):
    H = np.matrix([[1,0,1,1,1,0,0],[1,1,0,1,0,1,0],[1,1,1,0,0,0,1]])
    M = H * R
    for x in M:
        x[0] %= 2
    for i in range(7):
        if list(M) == list(H[:,i]):
            R[i,0]= (R[i,0] + 1) % 2
    return R[0:4,0]

print(decode743(perturbe(code743(np.array([[1],[0],[1],[1]])))))
