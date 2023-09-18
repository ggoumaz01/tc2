import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import pytc2.sistemas_lineales as tc2

## --- Parámetros del filtro ---- ##

fs = 10e3
fs /= fs/2 #normalizada
N = [2, 4]
alpha = 0.8
i = 0
## ------- Implementación ------- ##

plt.close('all')

while i<2:
    
    bn = np.zeros_like(np.arange(N[i] + 1) , dtype=float)
    an = np.zeros_like(np.arange(N[i] + 1) , dtype=float)
    
    bn[0]  = alpha
    bn[-1] = 1
    
    an[0]  = 1
    an[-1] = alpha
    
    tf = sig.TransferFunction(bn, an)
    WW_Z, HH_Z = sig.freqz(bn, an, fs=fs)
    
    # Grafico de respuesta de módulo
    
    plt.figure(i)
    
    plt.subplot(2, 1, 1)
    plt.grid(visible=True)
    plt.title("Module response")
    plt.ylabel("Module [dB]")
    plt.plot(WW_Z, 20* np.log10(abs(HH_Z)), label=f'order {N[i]} digital all-pass filter')
    plt.ylim(-1, 1)
    axes_hdl = plt.gca()
    axes_hdl.legend()
    
    # Grafico de respuesta de fase
    
    plt.subplot(2, 1, 2)
    plt.grid(visible=True)
    plt.title("Phase response")
    plt.ylabel("Phase [deg]")
    plt.xlabel("Normalized frequency [fs/2]")
    plt.plot(WW_Z, np.angle(HH_Z, deg=True), label=f'order {N[i]} digital all-pass filter')
    axes_hdl = plt.gca()
    axes_hdl.legend()
    
    # Mapa de polos y ceros
    tc2.pzmap(tf, annotations=True ,filter_description=f'order {N[i]} digital all-pass filter',fig_id = i+2 , digital=True)

    i+=1
## ------------------------------ ## 

