import sympy as sp
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import pytc2.sistemas_lineales as tc2
from pytc2.general import print_latex

##  Resolución simbólica  ##

S, Z = sp.symbols("S Z", complex=True)
k = sp.symbols("k", real=True, positive=True)

Hs = 1/(S**2+S*sp.sqrt(2)+1)

f_bi = k * (Z-1)/(Z+1)

Hz = sp.collect(sp.simplify(sp.expand(Hs.subs(S, f_bi))),Z)

print("Función de transferencia del filtro analógico butterworth de segundo orden:")
print_latex('H(S)='+ sp.latex(Hs))
print("Función de transferencia del filtro digital asociado:")
print_latex('H(Z)='+ sp.latex(Hz))

## ---------------------- ##

## Defino función que resuelve el filtro para las frecuencias indicadas

def ejercicio_2(fc, fs):
    
    ## Normalizo sobre la frecuencia de corte
    norma = fc
    fc = fc/norma
    fs = fs/norma

    ## Genero el filtro analógico

    z,p,k = sig.buttap(2)
    z,p,k = sig.lp2lp_zpk(z, p, k, wo=2*np.pi*fc)
    num_s, den_s = sig.zpk2tf(z, p, k)
    tf_s = sig.TransferFunction(num_s, den_s)

    ## Genero la transferencia del filtro digital

    num_z, den_z = sig.bilinear(num_s, den_s ,fs=fs)
    tf_z = sig.TransferFunction(num_z, den_z, dt=1/fs)

    ## Analizo las funciones

    ## Para poder plotear ##
    z,p,k = sig.lp2lp_zpk(z, p, k, wo=fc/fs)
    num_s, den_s = sig.zpk2tf(z, p, k)
    tf_s = sig.TransferFunction(num_s, den_s)
    ## ------------------ ##

    plt.close('all')

    wrad_z, hh_z = sig.freqz(num_z, den_z)
    ww_z = wrad_z / np.pi
    
    _, hh_s = sig.freqresp(tf_s, w=ww_z)
    
    plt.figure(1)
    
    ## Plot de respuesta de módulo
    plt.subplot(2, 1, 1)
    plt.grid(visible=True)
    plt.title("Module response")
    plt.ylabel("Module [dB]")
    plt.plot(ww_z, 20* np.log10(abs(hh_z)), color ='r', label='Second Order Digital Butterworth')

    axes_hdl = plt.gca()
    axes_hdl.legend()

    ## Plot de respuesta de fase
    plt.subplot(2, 1, 2)
    plt.grid(visible=True)
    plt.title("Phase response")
    plt.ylabel("Phase [deg]")
    plt.xlabel("Normalized Frequency [fs/2]")
    plt.plot(ww_z, np.angle(hh_z, deg=True) , color ='r', label='Second Order Digital Butterworth')

    axes_hdl = plt.gca()
    axes_hdl.legend()
    
    ## Mapa de polos y ceros
    tc2.pzmap(tf_z, annotations=True ,filter_description="Second Order Digital Butterworth",fig_id=2 , digital=True)

    ## --------------------- ##


## -Resolución numérica- ##

# Defino datos de consigna

fc = 1e3
fs = 100e3

ejercicio_2(fc, fs)

