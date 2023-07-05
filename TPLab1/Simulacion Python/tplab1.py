#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 21:24:32 2023

@author: gonzalo
"""

import scipy.signal as sig
import matplotlib.pyplot as plt
import pytc2.sistemas_lineales as tc2
import math as m

# Requisitos de plantilla
f0 = 50         # [Hz]
BW = 10         # [Hz]              
n  = 2          # [Hz]

w0 = 2*m.pi*f0
Q = f0/BW

# Armo la funcion transferencia normalizada
den = [1, 1/Q, 1]
num = [1, 0, 1]

tf_notch = sig.TransferFunction(num, den)

print("Transferencia del filtro notch normalizado")
tc2.pretty_print_lti(num, den)

# Analizo la transferencia
tc2.analyze_sys(tf_notch, sys_name='Filtro Notch')


