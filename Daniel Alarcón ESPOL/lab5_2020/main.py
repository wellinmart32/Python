# Serie de Fourier, con n coeficientes
# Ref.Chapra 5ed Ejemplo 19.2 p548/pdf572
import numpy as np
# import sympy as sym
# import matplotlib.pyplot as plt
#
# # INGRESO
T = 2*np.pi
#
# t = sym.Symbol('t')
# ft = sym.Piecewise((-1,t <-T/2),
#                    (-1,t <-T/4),
#                    ( 1,t < T/4),
#                    (-1,t < T/2),
#                    (-1,True),)
# # intervalo
# a = -T/2
# b = T/2
#
# # número de coeficientes
# n = 4
#
# # PROCEDIMIENTO
# k = sym.Symbol('k')
# w0 = 2*np.pi/T
#
# # Términos ak para coseno
# enintegral = ft*sym.cos(k*w0*t)
# yaintegrado = sym.integrate(enintegral,(t,a,b))
# ak = (2/T)*yaintegrado
# ak = sym.simplify(ak)
#
# # Términos bk para seno
# enintegral = ft*sym.sin(k*w0*t)
# yaintegrado = sym.integrate(enintegral,(t,a,b))
# bk = (2/T)*yaintegrado
# bk = sym.simplify(bk)
#
# print(' expresión ak:')
# sym.pprint(ak)
# print('\n ak formato latex')
# print(sym.latex(ak))
#
# print('\n expresión bk:')
# sym.pprint(bk)
# print('\n bk formato latex')
# print(sym.latex(bk))