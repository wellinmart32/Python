'
Curso: Procesamiento Digital de Senales

Este es un conjunto de funciones utiles para ejecutar codigo en sus
tareas y proyectos de curso

Como usar?
Para activar las funciones en sus notebook jupyter deben incluir la 
siguiente sentencia antes de ingresar su codigo.

import dsp as dsp

Si quieren usar una funcion como por ejemplo la funcion 'step':

dsp.step(0, 1, 10)

'

import numpy as np
import scipy.signal as signal
from scipy.io import wavfile
import matplotlib.pyplot as plt
from matplotlib import mlab 


def unit_sample(n0,n1,n2):
    'Generate an unit sample sequence of length n2-n1+1, with a shift at n0'
    n = np.arange(n1,n2,1) # secuencia "tiempo" time-index
    x = np.zeros(n2-n1) # secuencia de ceros
    x[n0] = 1 # impulso
    return n,x


def step(n0,n1,n2):
    n = np.arange(n1, n2+1, 1) # time sequence
    x = np.zeros(n2+1-n1)
    x[np.where(n-n0 >= 0)] = 1
    return n, x



def exp_cpl(n, a):
    'Entrada es un arreglo o sequencia 1D llamada n, y la base exponencial es el numero a'
    return a**n




def sigshift(x,m,k):
    ' Implementa el desplazamiento y(n) = x(n-k)'
    n = m + k
    return n,x

def senal_sum(x1,n1,x2,n2):
    'Sumar dos secuencias cuyo resultado incluye los tiempos definidos por n1 y n2'
    n_min = np.min(np.concatenate((n1,n2))) # encontrar el valor minimo del vector concatenado n1,n2
    n_max = np.max(np.concatenate((n1,n2))) # igual pero el maximo
    n = np.arange(n_min, n_max+1,1) # duracion total de la secuencia suma y(n)  
    y1 = np.zeros(np.size(n)) # inicializacion de senales
    y2 = np.zeros(np.size(n))
    # Asignacion
    y1[np.nonzero(np.logical_and((n>=np.min(n1)) ,(n<=np.max(n1))))] = x1
    y2[np.nonzero(np.logical_and((n>=np.min(n2)) ,(n<=np.max(n2))))] = x2
    # Sumar las dos secuencias
    return n, y1 + y2


def upsample():
    ' Lab 3 tarea a resolver'
    return

def downsample():
    ' Lab 3 tarea a resolver'
    return

# Practica 3

def prin_alias(f_in,fs):
    """
    The code is from module ssd.py which was developed to accompany the book
    Signals and Systems for Dummies published by Wiley Publishing.

    Copyright 2012, 2013 Mark Wickert, mwickert@uccs.edu
    
    Calculate the principle alias frequencies.
    
    Given an array of input frequencies the function returns an
    array of principle alias frequencies.
    
    Parameters
    ----------
    f_in : ndarray of input frequencies
    fs : sampling frequency
    
    Returns
    -------
    f_out : ndarray of principle alias frequencies
    
    Examples
    --------
    >>> # Linear frequency sweep from 0 to 50 Hz
    >>> f_in = arange(0,50,0.1)
    >>> # Calculate principle alias with fs = 10 Hz
    >>> f_out = prin_alias(f_in,10)
    """
    return abs(np.rint(f_in/fs)*fs - f_in)

    """
    Principle alias via recursion
    f_out = np.copy(f_in)
    for k in range(len(f_out)):
        while f_out[k] > fs/2.:
            f_out[k] = abs(f_out[k] - fs)
    return f_out
    """

# Practica 4    

def from_wav(filename):
    """
    Read a wave file.

    A wrapper function for scipy.io.wavfile.read
    that also includes int16 to float [-1,1] scaling.

    Parameters
    ----------
    filename : file name string

    Returns
    -------
    fs : sampling frequency in Hz
    x : ndarray of normalized to 1 signal samples

    Examples
    --------
    >>> fs,x = from_wav('test_file.wav')
    """
    fs, x = wavfile.read(filename)
    return fs, x/32767.


def to_wav(filename,rate,x):
    """
    Write a wave file.

    A wrapper function for scipy.io.wavfile.write
    that also includes int16 scaling and conversion.
    Assume input x is [-1,1] values.

    Parameters
    ----------
    filename : file name string
    rate : sampling frequency in Hz

    Returns
    -------
    Nothing : writes only the *.wav file

    Examples
    --------
    >>> to_wav('test_file.wav', 8000, x)

    """
    x16 = np.int16(x*32767)
    wavfile.write(filename, rate, x16)

# Practica 5

def dft(xn,N):
    'Calcula la DFT de una secuencia'
    n = np.arange(0,N,1) # indice en el 'tiempo'
    m = np.arange(0,N,1) # indice en 'frecuencia'
    WN = np.exp(-1j*2*np.pi/N)
    # Conversion a matriz
    xn = np.mat(xn)
    n = np.mat(n)
    m = np.mat(m)
    # Creacion de la matriz WN
    nm = n.T * m
    WNnm = np.power(WN,nm)
    # Ejecutar la operacion DFT en forma matricial
    Xm = xn * WNnm
    Xm = np.squeeze(np.asarray(Xm)) # convertir a numpy array
    return Xm

def idft(Xm,N):
    'Calcula la IDFT de una secuencia'
    m = np.arange(0,N,1)
    n = np.arange(0,N,1)
    WN = np.exp(-1j*2*np.pi/N)
    # Conversion a matriz
    Xm = np.mat(Xm)
    n = np.mat(n)
    m = np.mat(m)
    # Creacion de la matriz WN
    nm = n.T * m
    WNnm = np.power(WN,-nm)
    # Ejecutar la operacion IDFT en forma matricial
    xn =  (Xm * WNnm)/N
    xn = np.squeeze(np.asarray(xn)) # convertir a numpy array
    return xn

# Practica 6

def ideal_lpf(wc, M):
    'wc: es la frecuencia de corte \     M: es la longtud del filtro'
    alpha = (M-1)/2
    n = np.arange(0,M)
    m = n - alpha
    fc = wc / np.pi
    h = fc * np.sinc(fc*m) # respuesta al impulso que es la transformada inversa de fourier de un filtro ideal
    return h


def freqz_mod(b,a):
    [w, H] = signal.freqz(b,a,1000,'True') # mil frecuencias espaciadas equidistantemente
    H = np.transpose(H[0:501])
    w = np.transpose(w[0:501])
    magH = np.abs(H)
    dB = 20 * np.log10(magH/magH.max())  # magnitud normalizada a escala logaritmica
    phase = np.angle(H)                  # fase de la respuesta de frecuencia
    grpdelay = signal.group_delay((b,a)) # retardo de grupo usa funcion de libreria scipy.signal
    return dB, phase, grpdelay, w

# Practica 7

def zplane(b,a,auto_scale=True,size=2,detect_mult=True,tol=0.001):
    """
    The code module ssd.py was developed to accompany the book
    Signals and Systems for Dummies published by Wiley Publishing.

    Copyright 2012, 2013 Mark Wickert, mwickert@uccs.edu

    Create an z-plane pole-zero plot.

    Create an z-plane pole-zero plot using the numerator
    and denominator z-domain system function coefficient
    ndarrays b and a respectively. Assume descending powers of z.

    Parameters
    ----------
    b : ndarray of the numerator coefficients
    a : ndarray of the denominator coefficients
    auto_scale : bool (default True)
    size : plot radius maximum when scale = False

    Returns
    -------
    (M,N) : tuple of zero and pole counts + plot window
    
    Notes
    -----
    This function tries to identify repeated poles and zeros and will 
    place the multiplicity number above and to the right of the pole or zero.
    The difficulty is setting the tolerance for this detection. Currently it
    is set at 1e-3 via the function signal.unique_roots.

    Examples
    --------
    >>> # Here the plot is generated using auto_scale
    >>> zplane(b,a)
    >>> # Here the plot is generated using manual scaling
    >>> zplane(b,a,False,1.5)
    """
    M = len(b) - 1
    N = len(a) - 1
    # Plot labels if multiplicity greater than 1
    x_scale = 1.5*size
    y_scale = 1.5*size   
    x_off = 0.02
    y_off = 0.01
    #N_roots = np.array([1.0])
    if M > 0:
        N_roots = np.roots(b)
    #D_roots = np.array([1.0])
    if N > 0:
        D_roots = np.roots(a)
    if auto_scale:
        if M > 0 and N > 0:
            size = max(np.max(np.abs(N_roots)),np.max(np.abs(D_roots)))+.1
        elif M > 0:
            size = max(np.max(np.abs(N_roots)),1.0)+.1
        elif N > 0:
            size = max(1.0,np.max(np.abs(D_roots)))+.1
        else:
            size = 1.1
    plt.figure(figsize=(5,5))
    plt.axis('equal')
    r = np.linspace(0,2*np.pi,200)
    plt.plot(np.cos(r),np.sin(r),'r--')
    plt.plot([-size,size],[0,0],'k-.')
    plt.plot([0,0],[-size,size],'k-.')
    if M > 0:
        if detect_mult == True:
            N_uniq, N_mult = unique_cpx_roots(N_roots,tol=tol)
            plt.plot(np.real(N_uniq),np.imag(N_uniq),'ko',mfc='None',ms=8,mew=2)
            #idx_N_mult = mlab.find(N_mult>1)
            idx_N_mult = find_indices(N_mult, lambda x: x>1) 
            for k in range(len(idx_N_mult)):
                x_loc = np.real(N_uniq[idx_N_mult[k]]) + x_off*x_scale
                y_loc =np.imag(N_uniq[idx_N_mult[k]]) + y_off*y_scale
                plt.text(x_loc,y_loc,str(N_mult[idx_N_mult[k]]),ha='center',va='bottom',fontsize=10)
        else:
            plt.plot(np.real(N_roots),np.imag(N_roots),'ko',mfc='None',ms=8,mew=2)                
    if N > 0:
        if detect_mult == True:
            D_uniq, D_mult=unique_cpx_roots(D_roots,tol=tol)
            plt.plot(np.real(D_uniq),np.imag(D_uniq),'kx',ms=8,mew=2)
            #idx_D_mult = mlab.find(D_mult>1)
            idx_D_mult = find_indices(D_mult, lambda x: x>1) 
            for k in range(len(idx_D_mult)):
                x_loc = np.real(D_uniq[idx_D_mult[k]]) + x_off*x_scale
                y_loc =np.imag(D_uniq[idx_D_mult[k]]) + y_off*y_scale
                plt.text(x_loc,y_loc,str(D_mult[idx_D_mult[k]]),ha='center',va='bottom',fontsize=10)            
        else:
            plt.plot(np.real(D_roots),np.imag(D_roots),'kx',ms=8,mew=2)                
    if M - N < 0:
        plt.plot(0.0,0.0,'bo',mfc='None',ms=8)
    elif M - N > 0:
        plt.plot(0.0,0.0,'kx',ms=8)
    if abs(M - N) > 1:
        plt.text(x_off*x_scale,y_off*y_scale,str(abs(M-N)),ha='center',va='bottom',fontsize=10)        
    plt.xlabel('Real Part')
    plt.ylabel('Imaginary Part')
    plt.title('Pole-Zero Plot')
    #plt.grid()
    plt.axis([-size,size,-size,size])
    return M,N

def unique_cpx_roots(rlist,tol = 0.001):
    """
    
    The average of the root values is used when multiplicity 
    is greater than one.

    Mark Wickert October 2016
    """
    uniq = [rlist[0]]
    mult = [1]
    for k in range(1,len(rlist)):
        N_uniq = len(uniq)
        for m in range(N_uniq):
            if abs(rlist[k]-uniq[m]) <= tol:
                mult[m] += 1
                uniq[m] = (uniq[m]*(mult[m]-1) + rlist[k])/float(mult[m])
                break
        uniq = np.hstack((uniq,rlist[k]))
        mult = np.hstack((mult,[1]))
    return np.array(uniq), np.array(mult)

def find_indices(a, func):
    """
    Replace mlab.find which was deprecated since matplotlib 2.2
    """
    return [i for (i, val) in enumerate(a) if func(val)]



def plot_freq_domain(w1, dB1, phase1, fs):
    fig,ax = subplots(2, sharex=True)
    freq = fs*w1/(2*np.pi) # convierte frecuencias digitales a analogicas
    ax[0].plot(freq,dB1)
    ax[0].set_ylabel(r"$ |H(f)| [dB]$",fontsize=20)
    ymin,ymax = -50 ,5+max(dB1)
    ax[0].axis(ymin = ymin,ymax=ymax)
    ax[0].grid()
    ax[1].plot(freq,phase1*180/np.pi)
    ax[1].set_xlabel("frequency (Hz)",fontsize=16)
    ax[1].set_ylabel(r"$ Phase [deg]$",fontsize=20)
    ax[1].grid()
    

# Practica 8

def u_buttap(N,omega_c):
    """
    Prototipo de filtro analogico Butterworth No-normalizado
    N:      orden del filtro
    omegac: frecuencia de corte en rad/seg
    """

    [z,p,k] = signal.buttap(N)
    p = p*omega_c
    k = k*omega_c**N
    B = np.real(np.poly(z))
    b0 = k
    b = k*B
    a = np.real(np.poly(p))
    return b,a



def afd_butt(wp,ws,rp,As):
    """
    Diseno de Filtros Butterworth bajo especificaciones generales
    SALIDA:
    b: coeficientes del numerador de H(s)
    a: coeficientes del denominador de H(s)
    ENTRADAS:
    wp: frecuencia limite en pasobanda (rad/seg)
    ws: frecuencia limite en banda de rechazo (stopband) (rad/seg)
    rp: sobreoscilacion (ripple) en pasobanda (dB)
    As: Atenuacion en banda de rechazo (dB)
    """
        
    if wp <= 0:
        print ("error la frecuencia pasobanda debe ser mayor que cero")
        return
    if ws <= wp:
        print ("error banda de rechazo debe ser mayor que la banda de paso")
        return
    if (rp <=0) and (As<0):
        print ("Ripple de banda de paso y/o la atenuacion debe ser mayor a cero")
        return
    
    N = np.ceil((np.log10((10**(rp/10.)-1)/(10**(As/10.)-1)))/(2*np.log10(wp/ws)))
    print (N)
    print ("Filtro Butterworth de orden N=", N)
    omega_c = wp/((10**(rp/10.)-1)**(1/(2*N)))
    b,a = u_buttap(N,omega_c)
    if b.size == 1:
        b = b.reshape(1)
    if a.size == 1:
        a = a.reshape(1)
    return b,a



def freqs_m(b,a,wmax):
    w = np.arange(0,500)*wmax/50.
    w,H = signal.freqs(b,a,worN=w)
    mag = np.abs(H)
    eps = np.finfo(float).eps
    db = 20*np.log10((mag+eps)/np.max(mag))
    phase = np.angle(H)
    return db, mag, phase,w


def plot_freq_domain(w1, dB1, phase1, fs):
    fig,ax = subplots(2, sharex=True)
    freq = fs*w1/(2*np.pi) # convierte frecuencias digitales a analogicas
    ax[0].plot(freq,dB1)
    ax[0].set_ylabel(r"$ |H(j\Omega)| [dB]$",fontsize=20)
    ymin,ymax = -80,5
    ax[0].axis(ymin = ymin,ymax=ymax)
    ax[0].axis(xmin = w1[0],xmax=1.0)
    ax[0].grid()
    ax[1].plot(freq,phase1*180/np.pi)
    ax[1].set_xlabel("frequency $\Omega / \pi$",fontsize=16)
    ax[1].set_ylabel(r"$ Phase [deg]$",fontsize=20)
    ax[1].grid()


def u_chb1ap(N,Rp,omegac):
    """
    Prototipo de filtro analogico Chebyshev I No-normalizado
    
    b:      coeficientes del polinomio numerador
    a:      coeficientes del polinomio denominador
    N:      orden del filtro
    Rp:     Ripple de banda de paso (dB)
    omegac: frecuencia de corte en rad/seg
    """
    [z,p,k] = signal.cheb1ap(N,Rp)
    a = np.real(np.poly(p))
    aNn = a[N]
    p = p*omegac
    a = np.real(np.poly(p))
    aNu = a[N]
    k = k*aNu/aNn
    b0 = k
    B = np.real(np.poly(z))
    b = k*B
    return b,a


def afd_chb1(wp,ws,rp,As):
    """
    Diseno de Filtros Chebyshev tipo I bajo especificaciones generales
    SALIDA:
    b:    coeficientes del numerador de H(s)
    a:    coeficientes del denominador de H(s)
    ENTRADAS:
    wp:   frecuencia limite en pasobanda (rad/seg)
    ws:   frecuencia limite en banda de rechazo (stopband) (rad/seg)
    rp:   sobreoscilacion (ripple) en pasobanda (dB)
    As:   Atenuacion en banda de rechazo (dB)
    """
    
    if wp <= 0:
        print ("error la frecuencia pasobanda debe ser mayor que cero")
        return
    if ws <= wp:
        print ("error banda de rechazo debe ser mayor que la banda de paso")
        return
    if (rp <=0) and (As<0):
        print ("Ripple de banda de paso y/o la atenuacion debe ser mayor a cero")
        return
    
    ep = np.sqrt(10**(rp/10.)-1.)
    A = 10**(As/20.)
    omega_c = wp
    omega_r = ws/wp
    g = np.sqrt(A*A-1.)/ep
    
    N = int(np.ceil(np.log10(g+np.sqrt(g*g-1.0))/np.log10(omega_r+np.sqrt(omega_r*omega_r-1.0))))
    
    print ("Filtro Chebyshev tipo I de orden N=", N)

    b,a = u_chb1ap(N,rp,omega_c)
    if b.size == 1:
        b = b.reshape(1)
    if a.size == 1:
        a = a.reshape(1)
    return b,a



def imp_invar(c,d,T):
    """
    Implementa el metodo 2 de Invarianza al impulso
    a: coeficientes del denominador de funcion transferencia dominio-z
    b: coeficientes del numerador de funcion transferencia dominio-z
    c: coeficientes del numerador de funcion H(s)
    d: coeficientes del denominador de funcion H(s)
    T: parametro de tiempo de muestreo
    """
    R,p,k = signal.residue(c,d) # Calcula la expansion de fracciones parciales de H(s)
    p = np.exp(p*T) # conversion de polos en plano-s al plano-z
    b,a = signal.invresz(R,p,k) # Calcula los coeficientes H(z) a partir de fracciones parciales
    b = np.real(b) # garantiza que los valores de coeficientes sean reales
    a = np.real(a)
    return b,a

