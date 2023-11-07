# ---- Modulos y funciones necesarias ---- #

import sympy as sp
from pytc2.dibujar import (display, dibujar_puerto_entrada, dibujar_puerto_salida, dibujar_funcion_exc_abajo, 
                           dibujar_elemento_serie, dibujar_elemento_derivacion, dibujar_tanque_RC_derivacion,
                           dibujar_espacio_derivacion, Capacitor, Resistor, ResistorIEC) 
import pytc2.cuadripolos as tc2


# ---------------------------------------- #

s = sp.symbols('s', complex=True) # Defino la variable compleja S

# ---------------------------------------- #

# ------------- Ejercicio 1 -------------- #

