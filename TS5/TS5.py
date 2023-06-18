#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 21:24:32 2023

@author: gonzalo
"""

import scipy.signal as sig
import matplotlib.pyplot as plt
import pytc2.sistemas_lineales as tc2

#Requisitos de plantilla

fc = 300       # [Hz]
fz = 100
n  = 3         # [Hz]

# Normalización de frecuencia
wc = 1
wz = fz/fc

# Plantilla del pasabajos prototipo
Omega_c = 1
Omega_z = 1/wz

# La frecuencia de corte se define como la frecuencia en la cual la transferencia cae 3dB respecto de la banda de paso
# Entonces diseño como Butter para esa frecuencia.

z_lp, p_lp, k_lp = sig.buttap(n)

num_lp, den_lp = sig.zpk2tf(z_lp, p_lp, k_lp)

# Agrego el cero de transmision

num_zero = [1/(Omega_z**2), 0, 1]

# Armo la transferencia del filtro prototipo
tf_proto = sig.TransferFunction(num_zero, den_lp)
print("Transferencia del filtro prototipo")
tc2.pretty_print_lti(num_zero, den_lp)

# factorizo la expresion
sos = tc2.tf2sos_analog(num_zero, den_lp)
print("\nTransferencia del filtro prototipo en su forma factorizada:\n")
tc2.pretty_print_SOS(sos)


# Aplico el núcleo de transformación PasaBajo-PasaAlto
num_hp, den_hp = sig.lp2hp(num_zero, den_lp)
tf_hp = sig.TransferFunction(num_hp, den_hp)

print("Transferencia del filtro pasa altos normalizado")
tc2.pretty_print_lti(num_hp, den_hp)

# factorizo la expresion
sos = tc2.tf2sos_analog(num_hp, den_hp)
print("\nTransferencia del filtro pasa altos normalizado en su forma factorizada:\n")
tc2.pretty_print_SOS(sos)

# Analizo las transferencias obtenidas
tc2.analyze_sys(tf_hp, sys_name='Pasa Altos')
tc2.analyze_sys(tf_proto, sys_name='Prototipo')


