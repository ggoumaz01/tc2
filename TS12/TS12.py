#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 10:59:14 2023

@author: gonzalo
"""

# ---- Modulos y funciones necesarias ---- #

import sympy as sp
from pytc2.general import print_latex, a_equal_b_latex_s
import pytc2.cuadripolos as tc2


# ---------------------------------------- #

s = sp.symbols('s', complex=True) # Defino la variable compleja S

# ---------------------------------------- #

# ------------- Ejercicio 1 -------------- #

# # Defino los símbolos necesarios

R1 = sp.nsimplify(5)
R2 = sp.nsimplify(2)
R3 = sp.nsimplify(10)
RL = sp.nsimplify(1)
C2 = sp.nsimplify(1/8)
C3 = sp.nsimplify(1/10)

# # Armo la MAI

#             ---C2---     ---C3---             #
#             -      -     -      -             #
# 0------------      -------      ---------  2  #
#        -    -      -     -      -     -       #
#        -    ---R2---     ---R3---     -       #
#        R1                             RL      #
#        -                              -       #
#        -                              -       #
# 1----------------------------------------  1  #

Ya = s*C2 + 1/R2 
Yb = s*C3 + 1/R3
Za = 1/Ya
Zb = 1/Yb

Ymai = sp.Matrix([
        [1/R1 + 1/(Za+Zb), -1/R1, -1/(Za+Zb)],
        [-1/R1, 1/R1 + 1/RL, -1/RL],
        [-1/(Za+Zb), -1/RL, 1/(Za+Zb) + 1/RL]
        
    ])

Tz = tc2.calc_MAI_ztransf_ij_mn(Ymai, ii=2, jj=1, mm=0, nn=1)

# Luego I1 = V2/RL -> Como RL es 1, I1 será de la misma forma que V2

Ti = Tz/RL 
print("La transferencia de corrientes de la red está dada por:")
print_latex(a_equal_b_latex_s('T(S)', sp.simplify(sp.expand(Ti))))


# ---------------------------------------- #

# ------------- Ejercicio 2 -------------- #

RL = sp.nsimplify(1)
C1 = sp.nsimplify(27/17)
C2 = sp.nsimplify(27/289)
C3 = sp.nsimplify(7/17)
L2 = sp.nsimplify(289/243)

# # Armo la MAI

#             ---L2---                     #
#             -      -                     #
# 0------------      -------------------2  #
#        -    -      -       -       -     #
#        -    ---C2---       -       -     #
#        C1                  C3      RL    #
#        -                   -       -     #
#        -                   -       -     #
# 1-------------------------------------1  #

Ymai = sp.Matrix([
        [s*C1 + s*C2 + 1/(s*L2), -s*C1, -(s*C2 + 1/(s*L2))],
        [-s*C1, s*C1 + s*C3 + 1/RL, -(s*C3 + 1/RL)],
        [-(s*C2 + 1/(s*L2)), -(s*C3 + 1/RL), s*C2 + 1/(s*L2) + s*C3 + 1/RL]   
    ])

Tz = tc2.calc_MAI_ztransf_ij_mn(Ymai, ii=2, jj=1, mm=0, nn=1)

print("La transimpedancia del circuito está dada por:")
print_latex(a_equal_b_latex_s('T(S)', sp.simplify(sp.expand(Tz))))

# ---------------------------------------- #