# ---- Modulos y funciones necesarias ---- #

import sympy as sp
from pytc2.dibujar import (display, dibujar_puerto_entrada, dibujar_puerto_salida, dibujar_funcion_exc_abajo, 
                           dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_tanque_RC_derivacion,
                           dibujar_espacio_derivacion, Capacitor, Resistor, ResistorIEC) 
from pytc2.general import print_latex, a_equal_b_latex_s
import pytc2.cuadripolos as tc2


# ---------------------------------------- #

s = sp.symbols('s', complex=True) # Defino la variable compleja S

# ---------------------------------------- #

# ------------- Ejercicio 1 -------------- #

# # Defino los símbolos necesarios

# C1 = sp.nsimplify(1)
# C2 = sp.nsimplify(2)
# C3 = sp.nsimplify(1/3)
# L2 = sp.nsimplify(1/2)
# L3 = sp.nsimplify(1)

# # Armo la MAI

# #         2                 #
# # 0---C1-----L3--C3-------3 #
# #         -                 #
# #         L2                #
# #         -                 #
# #         C2                #
# #         -                 #
# # 1-----------------------1 #

# Ymai = sp.Matrix([
#         [s*C1 ,0, -s*C1, 0],
#         [0, 1/(s*L2 + 1/(s*C2)), -1/(s*L2 + 1/(s*C2)), 0],
#         [-s*C1,  -1/(s*L2 + 1/(s*C2)) ,s*C1 + 1/(s*L2 + 1/(s*C2)) + 1/(s*L3 + 1/(s*C3)), -1/(s*L3 + 1/(s*C3))],
#         [0, 0, -1/(s*L3 + 1/(s*C3)), 1/(s*L3 + 1/(s*C3))]
        
#     ])



# ---------------------------------------- #

# ------------- Ejercicio 2 -------------- #

R1 = sp.nsimplify(1)
R2 = sp.nsimplify(1/4)
C1 = sp.nsimplify(2/5)
C2 = sp.nsimplify(4)
C3 = sp.nsimplify(2)

# Armo la MAI

#                              #
# 0---R1---C1---------------2  #
#               -      -       #
#               C2     -       #
#               -      C3      #
#               R2     -       #
#               -      -       #
# 1-------------------------1  #

Ymai = sp.Matrix([
        [1/(R1 + 1/(s*C1)), 0, -1/(R1 + 1/(s*C1))],
        [0, s*C3 + 1/(R2 + 1/(s*C2)), -(s*C3 + 1/(R2 + 1/(s*C2)))],
        [-1/(R1 + 1/(s*C1)), -(s*C3 + 1/(R2 + 1/(s*C2))), 1/(R1 + 1/(s*C1)) + s*C3 + 1/(R2 + 1/(s*C2)) ]   
    ])

Tf = tc2.calc_MAI_vtransf_ij_mn(Ymai, ii=2, jj=1, mm=0, nn=1, verbose=False)

print("La transferencia del circuito es")
print_latex(a_equal_b_latex_s('T(S)', sp.factor(Tf)))
# ---------------------------------------- #

# ------------- Ejercicio 2 - Por admitancias -------------- #

R1 = sp.nsimplify(2/35)
R2 = sp.nsimplify(4/5)
R3 = sp.nsimplify(2/3)
C1 = sp.nsimplify(5)
C2 = sp.nsimplify(5/4)

# Armo la MAI

#               ---C2---           #
#          2    -      -           #
# 0---R1---------      ---------3  #
#          -    -      -   -       #
#          -    ---R2---   -       #
#          C1              R3      #
#          -               -       #
#          -               -       #
# 1-----------------------------1  #

Ymai = sp.Matrix([
        [1/R1, 0, -1/R1, 0],
        [0, s*C1 + 1/R3, -s*C1, -1/R3],
        [-1/R1, -s*C1, 1/R1 + s*C1 + s*C2 + 1/R2, -(s*C2 + 1/R2)],
        [0, -1/R3, -(s*C2 + 1/R2), s*C2 + 1/R2 + 1/R3]   
    ])

Tf = tc2.calc_MAI_vtransf_ij_mn(Ymai, ii=3, jj=1, mm=0, nn=1, verbose=False)

print("La transferencia del circuito es")
print_latex(a_equal_b_latex_s('T(S)', sp.factor(Tf)))

# --------------------------------------------------------- #

