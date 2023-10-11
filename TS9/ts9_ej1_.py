# ---- M贸dulos y funciones necesarias ---- #

import sympy as sp
from pytc2.sintesis_dipolo import foster, cauer_LC, remover_polo_dc
from pytc2.dibujar import dibujar_foster_serie, dibujar_foster_derivacion, dibujar_cauer_LC
from pytc2.general import print_latex, print_subtitle, a_equal_b_latex_s



# ---------------------------------------- #

s = sp.symbols('s', complex=True) # Defino la variable compleja S

# ------------- Ejercicio 1 -------------- #

print('Se tiene la siguiente funci贸n de excitaci贸n, correspondiente a una Impedancia')

Z = ((s**2 + 3)*(s**2 + 1))/(s*(s**2 + 2))
Y = 1/Z

## --- Foster Serie --- ##

str_aux = a_equal_b_latex_s('Z(S)', Z)

Z = (s**4 + 4*s**2 + 3)/(s**3 + 2*s)
print_latex(a_equal_b_latex_s(str_aux[1:-1], Z))

print('Si se expande la funci贸n seg煤n Foster, Z(S) puede expresarse como')

k0, koo, ki_wi , _, Z_foster = foster(Z);

print_latex(a_equal_b_latex_s('Z(S)', Z_foster))

print("Con lo cual el circuito Foster serie resulta")
dibujar_foster_serie(k0, koo, ki_wi, z_exc=Z)

## --- Foster Paralelo --- ##

print("Ahora se trabaja con la admitancia")

str_aux = a_equal_b_latex_s('Y(S)', Y)
Y = 1/Z
print_latex(a_equal_b_latex_s(str_aux[1:-1], Y))

print("La cual expandida por Foster se corresponde con")

k0, koo, ki_wi , _, Y_foster = foster(Y);

print_latex(a_equal_b_latex_s('Y(S)', Y_foster))

print("Con lo cual el circuito Foster paralelo resulta")
dibujar_foster_derivacion(k0, koo, ki_wi, y_exc=Y)

## --- Cauer I --- ##

print("Se realizan remociones de los residuos en infinito")

koo, Z_cauer_oo, rem = cauer_LC(Z, remover_en_inf=True)
print_latex(a_equal_b_latex_s(a_equal_b_latex_s('Z(S)', Z)[1:-1], Z_cauer_oo))
dibujar_cauer_LC(koo, z_exc=Z_cauer_oo)

## --- Cauer II --- ##

print("Se realizan remociones de los residuos en cero")

k0, Z_cauer_oo, rem = cauer_LC(Z, remover_en_inf=False)
print_latex(a_equal_b_latex_s(a_equal_b_latex_s('Z(S)', Z)[1:-1], Z_cauer_oo))
dibujar_cauer_LC(k0, z_exc=Z_cauer_oo)

# ---------------------------------------- #

# ------ Funciones para dibujar Ej2 ------ #

from pytc2.dibujar import (dibujar_puerto_entrada, dibujar_funcion_exc_abajo, dibujar_espacio_derivacion,
                           dibujar_tanque_derivacion, dibujar_elemento_serie )

from schemdraw import Drawing
from schemdraw.elements import  Capacitor

from IPython.display import display


def dibujar_ej2(k0_w0=None, ki=None, y_exc=None):
    
## Fuente: Libreria pytc2
    
    d = Drawing(unit=4)
    
    bComponenteDibujado = False

    
    d = dibujar_puerto_entrada(d,
                                   voltage_lbl = ('+', '$V$', '-'), 
                                   current_lbl = '$I$')
    
    d, _ = dibujar_funcion_exc_abajo(d, 
                                              'Y',  
                                              y_exc, 
                                              hacia_salida = True,
                                              k_gap_width = 0.5)

    if not(k0_w0 is None):
    
        d = dibujar_elemento_serie(d, Capacitor, 1/k0_w0)


    if not(ki is None):

        for un_tanque in ki:

            if bComponenteDibujado:
                
                dibujar_espacio_derivacion(d)
             
            d = dibujar_tanque_derivacion(d, inductor_lbl = un_tanque[1], capacitor_lbl = 1/un_tanque[0])
            bComponenteDibujado = True

    display(d)


# ------------- Ejercicio 2 -------------- #

omega_resonancia = 1

Y = sp.nsimplify((3*s*(s**2 + 7/3))/((s**2 + 2)*(s**2 + 5)))
Z = 1/Y

print("Se tiene una funci髇 de admitancia dada por")
print_latex(a_equal_b_latex_s("Y(S)", Y))

print(f"Para poder fijar un cero en {omega_resonancia}rad/s " +
      "es necesario realizar una remoci髇 parcial del polo en cero de la impendacia Z = 1/Y")

Z_w0, k0_w0 = remover_polo_dc(Z, omega_zero=omega_resonancia)

print("Resulta entonces que")
print_latex(a_equal_b_latex_s(a_equal_b_latex_s("Z(S)",Z)[1:-1], Z_w0+k0_w0))
k0_w0 *= s ## Por como lo devuelve la funcion remover_polo_dc (solo me interesa el coeficiente)

print("Con la impedancia que resulta de realizar la remoci髇, se realiza un foster derivaci髇 para llegar al circuito dado en la consigna")

k0, koo, ki_wi, _, Y_foster = foster(1/Z_w0)
print_latex(a_equal_b_latex_s('Y_2(S)',Y_foster ))

dibujar_ej2(k0_w0, ki_wi, Y)

# ---------------------------------------- #