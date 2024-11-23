#-----------------------------------------------------------------------------#
#=============================================================================#
#----------------\

# PME-3480 - Motores de Combustão Interna
# 1D Otto cycle simulator - 2024
# Implementation 2 - Group 4

#-----------------------------------------------------------------------------#
#=============================================================================#
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Importando bibliotecas
import math
import matplotlib.pyplot as plt
import numpy as np
import OttoCycle as oc 
from scipy.optimize import curve_fit
import os
os.system('cls')

#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
#=============================================================================#
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Parâmetros do motor
n = 3500/60  # Velocidade rotacional [rev/s]
x = 2        # Número de tempos por ciclo
Zc = 1       # Número de cilindros
PCI = 50000  # Poder calorífico inferior do combustível [kJ/kg]
#-----------------------------------------------------------------------------#
#=============================================================================#
#-----------------------------------------------------------------------------#

# Parâmetros do grupo:
p = [[9.00, 528.0, 1.431, 2.148],
     [10.00, 518.0, 1.436, 2.227],
     [11.00, 510.0, 1.442, 2.300]]  
# rv, Texh [°C], mpF [kg/h], Mt [kgfm]

# Dados da função de Wiebe
W = np.array([[-15.0, 0.0, 0.0, 0.0],
              [-14.0, 0.0, 0.0, 0.0],
              [-13.0, 0.0, 0.0, 0.0],
              [-12.0, 0.0, 0.0, 0.0],
              [-11.0, 0.0, 0.0, 0.0],
              [-10.0, 0.0, 0.0, 0.0],
              [-9.0, 1.249219075419194525e-03, 0.0, 0.0],
              [-8.0, 9.950166250832004344e-03, 0.0, 0.0],
              [-7.0, 3.318682227845726196e-02, 0.0, 0.0],
              [-6.0, 7.688365361336424453e-02, 0.0, 0.0],
              [-5.0, 1.446546726925775905e-01, 0.0, 0.0],
              [-4.0, 2.366205056631469628e-01, 2.958577720375998865e-03, 0.0],
              [-3.0, 3.486772604862288238e-01, 2.342497754313854763e-02, 0.0],
              [-2.0, 4.727075759569516755e-01, 7.688365361336435555e-02, 0.0],
              [-1.0, 5.979786169053452616e-01, 1.727345281701816448e-01, 0.0],
              [0.0, 7.134952031398098526e-01, 3.095214495228908458e-01, 0.0],
              [1.0, 8.105727058287343079e-01, 4.727075759569517865e-01, 9.950166250832004344e-03],
              [2.0, 8.846748789619376385e-01, 6.380670467389586431e-01, 7.688365361336457759e-02],
              [3.0, 9.358319592586000768e-01, 7.806391168261050950e-01, 2.366205056631470738e-01],
              [4.0, 9.676130592270930642e-01, 8.846748789619377495e-01, 4.727075759569520086e-01],
              [5.0, 9.852829707013648353e-01, 9.483343931211155597e-01, 7.134952031398100747e-01],
              [6.0, 9.940239771049941275e-01, 9.806236824544479758e-01, 8.846748789619377495e-01],
              [7.0, 9.978477683533923948e-01, 9.940239771049941275e-01, 9.676130592270929531e-01],
              [8.0, 9.993176719472436353e-01, 9.985109688860511756e-01, 9.940239771049941275e-01],
              [9.0, 9.998109974403737166e-01, 9.997055408783267483e-01, 9.993176719472436353e-01],
              [10.0, 1.0, 1.0, 1.0],
              [11.0, 1.0, 1.0, 1.0],
              [12.0, 1.0, 1.0, 1.0],
              [13.0, 1.0, 1.0, 1.0],
              [14.0, 1.0, 1.0, 1.0]])
# primeira coluna - CAD (em graus)
# segunda coluna - fração de massa queimada xb1
# terceira coluna - fração de massa queimada xb2
# quarta coluna - fração de massa queimada xb3

# partindo dos dados experimentais, podemos observar:
theta0 = [-9, -4, 1]
dtheta = [12, 9, 6]
#Valores de theta0 e dtheta para os casos 1, 2 e 3


#-----------------------------------------------------------------------------#
# Loop sobre os casos
#-----------------------------------------------------------------------------#
for z in [1, 2, 3]:  # Diferentes valores de xb

    # Definição dos parâmetros através dos dados experimentais
    theta = W[:, 0]  # CAD em graus
    xb = W[:, z]

    # Selecionando os valores de dtheta e theta0
    theta0caso = theta0[z -1]
    dthetacaso = dtheta[z -1]

    # Definição da função de Wiebe
    def Wiebe(theta, a, m):
        theta_adj = np.maximum((theta - theta0caso) / dthetacaso, 0)
        return 1 - np.exp(-a * theta_adj**m)

    initial = [5, 4]  # Chute inicial para os parâmetros
    params, covariance = curve_fit(Wiebe, theta, xb, p0=initial, maxfev=100000)
    a, m = params

    theta_fit = np.linspace(min(theta), max(theta), 500)
    xb_fit = Wiebe(theta_fit, a, m)

    # Plotando a função de Wiebe com os parametros calculados
    texto_parametros = f"Parâmetros:\na = {a:.2f}\nm = {m:.2f}"
    plt.annotate(texto_parametros, xy=(0.05, 0.85), xycoords='axes fraction', fontsize=10, 
             bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))
    plt.plot(theta_fit, xb_fit)
    plt.xlabel("Ângulo do Virabrequim (Radianos)")
    plt.ylabel("Fração Mássica Queimada (xb)")
    plt.title(f"Função de Wiebe para xb{z}")
    plt.show()

    # Plotando o gráfico de xb em funcção do CAD experimental, para comparação
    plt.plot(theta, xb)
    plt.xlabel("Ângulo do Virabrequim (Radianos)")
    plt.ylabel("Fração Mássica Queimada (xb)")
    plt.title(f"Curva da fração mássica experimental do caso {z}")
    plt.show()

    # Calculando a derivada da função de Wiebe
    dx = np.gradient(theta_fit)
    dy = np.gradient(xb_fit)
    derivative = dy / dx

    # Plotando o gráfico da derivada da função de Wiebe
    plt.figure()
    plt.plot(theta_fit, derivative)
    plt.xlabel("Ângulo de Manivela (Radianos)")
    plt.ylabel("Derivada de xb em relação a θ")
    plt.title(f"Derivada da Função de Wiebe para xb{z}")
    plt.legend()
    plt.show()

    for i in [1]:  # Para S2, apenas o caso com rv=10 será simulado
        #-----------------------------------------------------------------------------#
        # Parâmetros do motor
        #-----------------------------------------------------------------------------#
        case1 = 'fired'

        #-----------------------------------------------------------------------------#
        # Cilindro
        #-----------------------------------------------------------------------------#
        b = 90/1000  # Diâmetro do cilindro [m]
        s = 90/1000  # Curso [m]
        l = 120/1000  # Comprimento da biela [m]
        rv = p[i][0]  # Taxa de compressão

        #-----------------------------------------------------------------------------#
        # Cronograma das válvulas
        #-----------------------------------------------------------------------------#
        ThIVO = np.radians(+360.0)  # Abertura da válvula de admissão [radianos]
        ThIVC = np.radians(-150.0)  # Fechamento da válvula de admissão [radianos]
        ThEVO = np.radians(+150.0)  # Abertura da válvula de escape [radianos]
        ThEVC = np.radians(-360.0)  # Fechamento da válvula de escape [radianos]

        #-----------------------------------------------------------------------------#
        # Parâmetros da função de Wiebe
        #-----------------------------------------------------------------------------#
        ThSOC = -10.*(np.pi/180)  # Start of combustion [radians]
        ThEOC = +10.*(np.pi/180)  # End of combustion [radians]
        aWF = a     # Fator de eficiência de Wiebe
        mWF = m     # Fator de forma de Wiebe

        #-----------------------------------------------------------------------------#
        # Parâmetros de contorno
        #-----------------------------------------------------------------------------#
        pint = 100e3     # Pressão de admissão [Pa]
        Tint = 273.15 + 25  # Temperatura de admissão [K]
        pexh = 100e3     # Pressão de escape [Pa]
        Texh1 = p[i][1] + 273.15  # Temperatura de escape [K]
        phi = 1.0        # Razão de equivalência
        fuel = 'CH4'     # Tipo de combustível

        #-----------------------------------------------------------------------------#
        # Parâmetros do simulador
        #-----------------------------------------------------------------------------#
        pars1 = (case1, b, s, l, rv, n, ThIVO, ThIVC, ThEVO, ThEVC,
                 ThSOC, ThEOC, aWF, mWF, pint, Tint, pexh, Texh1, phi, fuel)

        #-----------------------------------------------------------------------------#
        # Configuração do ângulo de manivela
        #-----------------------------------------------------------------------------#
        Th0 = np.radians(-360.0)  # Ângulo inicial [radianos]
        Th1 = np.radians(+360.0)  # Ângulo final [radianos]
        Ths = np.radians(1.0)     # Passo do ângulo [radianos]
        Thn = int(((Th1 - Th0)/Ths) + 1)  # Número de pontos
        Th = np.linspace(start=Th0, stop=Th1, num=Thn, endpoint=True)  # Ângulos de manivela [radianos]

        #-----------------------------------------------------------------------------#
        # Simulação do Ciclo Otto
        #-----------------------------------------------------------------------------#
        v1, m1, t1, p1 = oc.ottoCycle(Th, pars1)  # Executa a simulação do ciclo Otto

        plt.plot(Th, p1)
        plt.title(f"Pressão em função do ângulo do virabrequim no caso xb{z}")
        plt.xlabel('Ângulo do Virabrequim (Radianos)')
        plt.ylabel('Pressão [Pa]')
        plt.show()

        plt.plot(Th, t1)
        plt.title(f"Temperatura em função do ângulo do virabrequim no caso xb{z}")
        plt.xlabel('Ângulo do Virabrequim (Radianos)')
        plt.ylabel('Temperatura [K]')
        plt.show()