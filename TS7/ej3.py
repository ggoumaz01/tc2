import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import pytc2.sistemas_lineales as tc2
    
## Frecuencia de sampling y normalización

fs = 1e3
norma = fs/2
fs = fs/norma

## -- Filtros de media móvil -- ##

# Coeficientes para h1

bn1 = np.array([1, 1])
an1 = np.array([1, 0])

# Coeficientes para h2

bn2 = np.array([1, 1, 1])
an2 = np.array([1, 0, 0])

# Armo la transferencia

H1 = sig.TransferFunction(bn1, an1, dt=1/fs)
H2 = sig.TransferFunction(bn2, an2, dt=1/fs)

WW_Z, HH_1 = sig.freqz(bn1, an1, fs=fs)
_, HH_2  = sig.freqz(bn2, an2, fs=fs)

# Grafico respuesta de módulo
plt.close('all')
plt.figure(1)

plt.subplot(2, 1, 1)
plt.grid(visible=True)
plt.title("Module response")
plt.ylabel("Module [dB]")
plt.plot(WW_Z, 20* np.log10(abs(HH_1)), color ='r', label='H1 filter')
plt.plot(WW_Z, 20* np.log10(abs(HH_2)), color ='b', label='H2 filter')

axes_hdl = plt.gca()
axes_hdl.legend()

# Grafico respuesta de fase

plt.subplot(2, 1, 2)
plt.grid(visible=True)
plt.title("Phase response")
plt.ylabel("Phase [deg]")
plt.xlabel("Normalized frequency [fs/2]")
plt.plot(WW_Z, np.angle(HH_1, deg=True) , color ='r', label='H1 filter')
plt.plot(WW_Z, np.angle(HH_2, deg=True) , color ='b', label='H2 filter')

axes_hdl = plt.gca()
axes_hdl.legend()

# Mapa de polos y ceros
tc2.pzmap(H1, annotations=True ,filter_description="H1 poles and zeros",fig_id=2 , digital=True)
tc2.pzmap(H2, annotations=True ,filter_description="H2 poles and zeros",fig_id=2, digital=True)

## ---------------------------- ##

## - Filtros diferenciadores -  ##

# Coeficientes para h1

bn1 = np.array([1, -1])
an1 = np.array([1, 0])

# Coeficientes para h2

bn2 = np.array([1, 0 ,-1])
an2 = np.array([1, 0, 0])

# Armo la tarnsferencia

H1 = sig.TransferFunction(bn1, an1, dt=1/fs)
H2 = sig.TransferFunction(bn2, an2, dt=1/fs)

WW_Z, HH_1 = sig.freqz(bn1, an1, fs=fs)
_, HH_2  = sig.freqz(bn2, an2, fs=fs)

# Grafico respuesta de módulo
plt.figure(3)

plt.subplot(2, 1, 1)
plt.grid(visible=True)
plt.title("Module response")
plt.ylabel("Module [dB]")
plt.plot(WW_Z, 20* np.log10(abs(HH_1)), color ='r', label='H1 filter')
plt.plot(WW_Z, 20* np.log10(abs(HH_2)), color ='b', label='H2 filter')

axes_hdl = plt.gca()
axes_hdl.legend()

# Grafico respuesta de fase

plt.subplot(2, 1, 2)
plt.grid(visible=True)
plt.title("Phase response")
plt.ylabel("Phase [deg]")
plt.xlabel("Normalized frequency [fs/2]")
plt.plot(WW_Z, np.angle(HH_1, deg=True) , color ='r', label='H1 filter')
plt.plot(WW_Z, np.angle(HH_2, deg=True) , color ='b', label='H2 filter')

axes_hdl = plt.gca()
axes_hdl.legend()

# Mapa de polos y ceros
tc2.pzmap(H1, annotations=True ,filter_description="H1 poles and zeros",fig_id=4 , digital=True)
tc2.pzmap(H2, annotations=True ,filter_description="H2 poles and zeros",fig_id=4, digital=True)
## ---------------------------- ##
