#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 20:14:15 2023

@author: gonzalo
"""

# ---- Modulos y funciones necesarias ---- #

import sympy as sp
from pytc2.general import print_latex, a_equal_b_latex_s
from pytc2.remociones import modsq2mod_s
import pytc2.cuadripolos as tc2


# ---------------------------------------- #

s = sp.symbols('s', complex=True) # Defino la variable compleja S

# ---------------------------------------- #

# ------------- Calculos utilizados para desarrollo -------------- #

s21 = sp.Rational('15')/(s**3+ sp.Rational('6')*s**2 + sp.Rational('15')*s + sp.Rational('15'))
module_s21_square = s21 * s21.subs(s, -s) 

module_s11_square = sp.factor(1 - module_s21_square)
s11 = modsq2mod_s(module_s11_square)

print_latex(a_equal_b_latex_s('s11(s)', s11))

# ---------Impedancia de entrada---------- #

# Armo la MAI
#                     2 
# 0----Rg----L1---------------L2-----------  3  #
#                     -                 -       #
#                     -                 -       #
#                     C                 RL      #
#                     -                 -       #
#                     -                 -       #
# 1----------------------------------------  1  #

Rg = 1
RL = 1
C  = 0.552
L1 = 1.2547
L2 = 0.1925



Ymai = sp.Matrix([
        [1/(Rg + s*L1), 0, -1/(Rg+ s*L1), 0],
        [0, s*C + 1/RL, -s*C, -1/RL],
        [-1/(Rg + s*L1), -s*C, s*C + 1/(Rg + s*L1) + 1/(s*L2), -1/(s*L2)],
        [0, -1/RL, -1/(s*L2), 1/RL + 1/(s*L2)]
        
    ])
 
Zin = tc2.calc_MAI_impedance_ij(Ymai, ii=0, jj=1)

print("La impedancia de entrada de la red está dada por:")
print_latex(a_equal_b_latex_s("Z_in(s)", sp.factor(Zin)));

# ---------------------------------------- #