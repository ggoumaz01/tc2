import sympy as sp
import pytc2.cuadripolos as tc2
from pytc2.general import print_latex


'''           1
   0-----SL1-------SL3---------2
              -            -
            1/SC2          R
              -            -
   3---------------------------
'''

# Defino los simbolos necesarios 

S = sp.symbols('S', complex=True)
Z1, Z3 = sp.symbols('L1 L2', complex=True)
Y1 = 1/(S*Z1)
Y2 = S*sp.symbols('C2', complex=True)
Y3 = 1/(S*Z3)
R = sp.symbols('R', real=True, positive = True)
G = 1/R

# Armo la MAI

Ymai = sp.Matrix([
        [Y1 , -Y1     ,  0  , 0   ],
        [-Y1, Y1+Y2+Y3, -Y3 , -Y2 ],
        [ 0 ,   -Y3   , Y3+G, -G  ],
        [ 0 ,   -Y2   ,  -G , Y2+G]
    ])

#con_detalles = True
con_detalles = False

# Calculo la funcion de transferencia a partir de la MAI

print("Funcion Transferencia de tension:")
tf = tc2.calc_MAI_vtransf_ij_mn(Ymai, 2, 3, 0 ,3, verbose=con_detalles)
print_latex('T(S)' + '=' + sp.latex(tf))

print("Si reemplazo por los valores de los componentes, obtengo la Transferencia de tensiones normalizada")
tf_n = sp.simplify(tf.subs(Z1, (3/2)))
tf_n = sp.simplify(tf_n.subs(Z3, (1/2)))
tf_n = sp.simplify(tf_n.subs(Y2, S*(4/3)))
tf_n = sp.simplify(tf_n.subs(R, 1))

print_latex('T(S)' + '=' + sp.latex(tf_n))
# Calculo la impedancia de entrada

print("Impedancia de entrada:")
Zin = tc2.calc_MAI_impedance_ij(Ymai, 0, 3, verbose=con_detalles)
print_latex('Z_{in}(S)' +'='+ sp.latex(Zin))
