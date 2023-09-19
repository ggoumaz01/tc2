# -------- MÃ³dulos a utilizar -------- # 

import scipy.signal as sig
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.io as sio
from pytc2.sistemas_lineales import plot_plantilla

import warnings
warnings.filterwarnings('ignore')

# ------------------------------------ #

# Para fijar el estilo de gráficos en el notebook

fig_sz_x = 13
fig_sz_y = 10
fig_dpi = 80
fig_font_size = 12

plt.rcParams['figure.figsize'] = (fig_sz_x, fig_sz_y)
plt.rcParams['figure.dpi'] = (fig_dpi)
plt.rcParams.update({'font.size':fig_font_size})

# ----- Requisitos de plantilla ------ #     

fs = 1e3 #Hz
f_nyq = fs/2

fs1 = 1.0 #Hz
fp1 = 3.0 #Hz
fp2 = 25.0 #Hz
fs2 = 35.0 #Hz
ripple = 0.5
att = 40

# ------------------------------------ #

# ---- DiseÃ±o de filtro tipo IIR ----- #
worN=4096

f_pass = [fp1, fp2]
f_stop = [fs1, fs2]

bp_IIR_sos = sig.iirdesign(f_pass, f_stop, ripple, att, analog=False ,ftype='butter',output='sos', fs=fs)

num_IIR, den_IIR = sig.sos2tf(bp_IIR_sos)

ww_z, hh_z = sig.sosfreqz(bp_IIR_sos, worN=worN ,fs=fs)
#ww_z = ww_z * fs/(2*np.pi)

plt.close('all')
plt.figure(1)

## Plot de respuesta de mÃ³dulo
plt.subplot(2, 1, 1)
plt.grid(visible=True)
plt.title("IIR Module response")
plt.ylabel("[dB]")
plt.plot(ww_z, 20* np.log10(abs(hh_z)), color ='purple')
plt.axis([0, fs/20, -50, 10])

## Plot de plantilla
plot_plantilla(filter_type = 'bandpass', fpass = [fp1, fp2] , ripple = ripple , fstop = [fs1, fs2], attenuation = att, fs = fs)

## Plot de respuesta de fase
plt.subplot(2, 1, 2)
plt.grid(visible=True)
plt.title("IIR Phase response")
plt.ylabel("[rad]")
plt.xlabel('[Hz]')
plt.plot(ww_z, np.angle(hh_z) , color ='purple')
plt.axis([0, fs/20, -2*np.pi, 2*np.pi])

## Retardo de grupo
ww_gd, gd_IIR = sig.group_delay((num_IIR, den_IIR), w=worN ,fs=fs)

plt.figure(2)
plt.grid(visible=True)
plt.title("IIR Group Delay")
plt.ylabel("[samples]")
plt.xlabel('[Hz]')
plt.plot(ww_gd, gd_IIR, color ='purple')
plt.axis([0, fs/20, -2, 6])

## Respuesta al impulso

samples = 2000

impulse = sig.unit_impulse(samples)
IIR_h = sig.sosfilt(bp_IIR_sos, impulse)

plt.figure(3)
plt.grid(visible=True)
plt.title("IIR Impulse response x10")
plt.ylabel("[Amplitude]")
plt.xlabel('[samples]')
plt.plot(np.arange(0, samples), 10*IIR_h , color ='purple', label='IIR impulse response')
plt.plot(np.arange(0, samples), impulse , color ='red', label='impulse')

axes_hdl = plt.gca()
axes_hdl.legend()

# ------------------------------------ #

# ---- Diseño de filtro tipo FIR ----- #

n_taps = 5001
freq = np.array([0 , fs1, fp1, fp2, fs2 ,f_nyq])/f_nyq
gains = 10**(np.array([-att, -att, -ripple, -ripple, -att, -att])/20)

num_win = sig.firwin2(n_taps, freq, gains, window='hamming')
den = np.zeros(n_taps)
den[0] = 1.0

w, hh_win = sig.freqz(num_win, den, worN=n_taps)

# renormalizo el eje de frecuencia
w = w / np.pi * f_nyq

plt.figure(4)
## Plot de respuesta de modulo 
plt.subplot(2, 1, 1)
plt.grid(visible=True)
plt.title("FIR Module response")
plt.ylabel("[dB]")
plt.plot(w, 20 * np.log10(abs(hh_win)),color='purple' ,label='FIR-Win {:d}'.format(num_win.shape[0]))
plt.axis([0, f_nyq/10, -60, 5 ]);

axes_hdl = plt.gca()
axes_hdl.legend()

## Plot de plantilla
plot_plantilla(filter_type = 'bandpass', fpass = [fp1, fp2] , ripple = ripple , fstop = [fs1, fs2], attenuation = att, fs = fs)

## Plot de respuesta de fase
phase = np.angle(hh_win)

plt.subplot(2, 1, 2)
plt.grid(visible=True)
plt.title("FIR Phase response")
plt.ylabel("[rad]")
plt.xlabel('[Hz]')
plt.plot(w, phase , color ='purple', label='FIR-Win {:d}'.format(num_win.shape[0]))
plt.axis([0, f_nyq/10, -2*np.pi, 2*np.pi])

## Retardo de grupo
ww_gd, gd_FIR = sig.group_delay((num_win, den), w=n_taps ,fs=fs)

plt.figure(5)
plt.grid(visible=True)
plt.title("FIR Group Delay")
plt.ylabel("[samples]")
plt.xlabel('[Hz]')
plt.plot(ww_gd, gd_FIR, color ='purple')
plt.axis([0, f_nyq/10, 0, n_taps-1])

## Respuesta al impulso

samples = 3000

impulse = sig.unit_impulse(samples)
FIR_h = sig.lfilter(num_win, den, impulse)

plt.figure(6)
plt.grid(visible=True)
plt.title("FIR Impulse response (X10)")
plt.ylabel("[Amplitude]")
plt.xlabel('[samples]')
plt.plot(np.arange(0, samples), 10*FIR_h , color ='purple', label='FIR impulse response')
plt.plot(np.arange(0, samples), impulse , color ='red', label='impulse')

axes_hdl = plt.gca()
axes_hdl.legend()

# ----------------------------- #

# --------- Parte 2 ----------- #

## - Señal de ECG registrada a 1 kHz, con contaminación de diversos orígenes - ##

# para listar las variables que hay en el archivo
#io.whosmat('ecg.mat')
mat_struct = sio.loadmat('ecg.mat')

ecg_one_lead = mat_struct['ecg_lead']
ecg_one_lead = ecg_one_lead.flatten()
cant_muestras = len(ecg_one_lead)

# ---- Uso del filtro FIR ----  #

ECG_f_win = sig.lfilter(num_win, den, ecg_one_lead)

## Calculo la demora introducida

demora = int(np.round(gd_FIR[500])) # FIR tiene gd cte


# Segmentos de interés con ALTA contaminación

regs_interes = ( 
        np.array([5, 5.2]) *60*fs, # minutos a muestras
        np.array([12, 12.4]) *60*fs, # minutos a muestras
        np.array([15, 15.2]) *60*fs, # minutos a muestras
        )

for ii in regs_interes:
    
    # intervalo limitado de 0 a cant_muestras
    zoom_region = np.arange(np.max([0, ii[0]]), np.min([cant_muestras, ii[1]]), dtype='uint')
    
    plt.figure(figsize=(fig_sz_x, fig_sz_y), dpi= fig_dpi, facecolor='w', edgecolor='k')
    plt.plot(zoom_region, ecg_one_lead[zoom_region], label='ECG', linewidth=2)
    #plt.plot(zoom_region, ECG_f_butt[zoom_region], label='Butter')
    plt.plot(zoom_region, ECG_f_win[zoom_region + demora], label='Win')
    
    plt.title('ECG filtering example from ' + str(ii[0]) + ' to ' + str(ii[1]) )
    plt.ylabel('Adimensional')
    plt.xlabel('Muestras (#)')
    
    axes_hdl = plt.gca()
    axes_hdl.legend()
    axes_hdl.set_yticks(())
            
    plt.show()
    
# Segmentos de interés con BAJA contaminación
    
regs_interes = ( 
        [4000, 5500], # muestras
        [10e3, 11e3], # muestras
        )

for ii in regs_interes:
    
    # intervalo limitado de 0 a cant_muestras
    zoom_region = np.arange(np.max([0, ii[0]]), np.min([cant_muestras, ii[1]]), dtype='uint')
    
    plt.figure(figsize=(fig_sz_x, fig_sz_y), dpi= fig_dpi, facecolor='w', edgecolor='k')
    plt.plot(zoom_region, ecg_one_lead[zoom_region], label='ECG', linewidth=2)
    #plt.plot(zoom_region, ECG_f_butt[zoom_region], label='Butter')
    plt.plot(zoom_region, ECG_f_win[zoom_region + demora], label='Win')
    
    plt.title('ECG filtering example from ' + str(ii[0]) + ' to ' + str(ii[1]) )
    plt.ylabel('Adimensional')
    plt.xlabel('Muestras (#)')
    
    axes_hdl = plt.gca()
    axes_hdl.legend()
    axes_hdl.set_yticks(())
            
    plt.show()
# ----------------------------- #

# ---- Uso del filtro IIR ----  #

ECG_f_butt = sig.sosfilt(bp_IIR_sos, ecg_one_lead)

# Segmentos de interés con ALTA contaminación

regs_interes = ( 
        np.array([5, 5.2]) *60*fs, # minutos a muestras
        np.array([12, 12.4]) *60*fs, # minutos a muestras
        np.array([15, 15.2]) *60*fs, # minutos a muestras
        )

for ii in regs_interes:
    
    # intervalo limitado de 0 a cant_muestras
    zoom_region = np.arange(np.max([0, ii[0]]), np.min([cant_muestras, ii[1]]), dtype='uint')
    
    plt.figure(figsize=(fig_sz_x, fig_sz_y), dpi= fig_dpi, facecolor='w', edgecolor='k')
    plt.plot(zoom_region, ecg_one_lead[zoom_region], label='ECG', linewidth=2)
    plt.plot(zoom_region, ECG_f_butt[zoom_region], label='Butter')
    #plt.plot(zoom_region, ECG_f_win[zoom_region + demora], label='Win')
    
    plt.title('ECG filtering example from ' + str(ii[0]) + ' to ' + str(ii[1]) )
    plt.ylabel('Adimensional')
    plt.xlabel('Muestras (#)')
    
    axes_hdl = plt.gca()
    axes_hdl.legend()
    axes_hdl.set_yticks(())
            
    plt.show()
    
# Segmentos de interés con BAJA contaminación
    
regs_interes = ( 
        [4000, 5500], # muestras
        [10e3, 11e3], # muestras
        )

for ii in regs_interes:
    
    # intervalo limitado de 0 a cant_muestras
    zoom_region = np.arange(np.max([0, ii[0]]), np.min([cant_muestras, ii[1]]), dtype='uint')
    
    plt.figure(figsize=(fig_sz_x, fig_sz_y), dpi= fig_dpi, facecolor='w', edgecolor='k')
    plt.plot(zoom_region, ecg_one_lead[zoom_region], label='ECG', linewidth=2)
    plt.plot(zoom_region, ECG_f_butt[zoom_region], label='Butter')
    #plt.plot(zoom_region, ECG_f_win[zoom_region + demora], label='Win')
    
    plt.title('ECG filtering example from ' + str(ii[0]) + ' to ' + str(ii[1]) )
    plt.ylabel('Adimensional')
    plt.xlabel('Muestras (#)')
    
    axes_hdl = plt.gca()
    axes_hdl.legend()
    axes_hdl.set_yticks(())
            
    plt.show()


# ----------------------------- #

# ---- Análisis mediante filtrado bidireccional ----  #


## --- Segmentos de interés con alta contaminación --- ##

# Procedemos al filtrado
ECG_f_butt = sig.sosfiltfilt(bp_IIR_sos, ecg_one_lead)
ECG_f_win = sig.filtfilt(num_win, den, ecg_one_lead)

# Segmentos de interés
regs_interes = ( 
        np.array([5, 5.2]) *60*fs, # minutos a muestras
        np.array([12, 12.4]) *60*fs, # minutos a muestras
        np.array([15, 15.2]) *60*fs, # minutos a muestras
        )

for ii in regs_interes:
    
    # intervalo limitado de 0 a cant_muestras
    zoom_region = np.arange(np.max([0, ii[0]]), np.min([cant_muestras, ii[1]]), dtype='uint')
    
    plt.figure(figsize=(fig_sz_x, fig_sz_y), dpi= fig_dpi, facecolor='w', edgecolor='k')
    plt.plot(zoom_region, ecg_one_lead[zoom_region], label='ECG', lw=2)
    plt.plot(zoom_region, ECG_f_butt[zoom_region], label='Butter')
    plt.plot(zoom_region, ECG_f_win[zoom_region], label='Win')
    
    plt.title('ECG filtering example from ' + str(ii[0]) + ' to ' + str(ii[1]) )
    plt.ylabel('Adimensional')
    plt.xlabel('Muestras (#)')
    
    axes_hdl = plt.gca()
    axes_hdl.legend()
    axes_hdl.set_yticks(())
            
    plt.show()

## --- Segmentos de interés con baja contaminación --- #

regs_interes = ( 
        [4000, 5500], # muestras
        [10e3, 11e3], # muestras
        )

for ii in regs_interes:
    
    # intervalo limitado de 0 a cant_muestras
    zoom_region = np.arange(np.max([0, ii[0]]), np.min([cant_muestras, ii[1]]), dtype='uint')
    
    plt.figure(figsize=(fig_sz_x, fig_sz_y), dpi= fig_dpi, facecolor='w', edgecolor='k')
    plt.plot(zoom_region, ecg_one_lead[zoom_region], label='ECG', lw=2)
    plt.plot(zoom_region, ECG_f_butt[zoom_region], label='Butter')
    plt.plot(zoom_region, ECG_f_win[zoom_region], label='Win')
    
    plt.title('ECG filtering example from ' + str(ii[0]) + ' to ' + str(ii[1]) )
    plt.ylabel('Adimensional')
    plt.xlabel('Muestras (#)')
    
    axes_hdl = plt.gca()
    axes_hdl.legend()
    axes_hdl.set_yticks(())
            
    plt.show()

# ----------------------------- #