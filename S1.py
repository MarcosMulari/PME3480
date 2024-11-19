#-----------------------------------------------------------------------------#
#=============================================================================#
#----------------\

# PME-3480 - Motores de Combustão Interna
# 1D Otto cycle simulator - 2024
# Implementation 1 - Group 4

#-----------------------------------------------------------------------------#
#=============================================================================#
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Importing libraries
import math
import matplotlib.pyplot as plt
import OttoCycle as oc 
import os
os.system('cls')

#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
#=============================================================================#
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# engine parameters
n = 3500/60  # Rotational speed [rev/s]
x = 2  # Number of strokes per cycle
Zc = 1  # Number of cylinders
PCI = 50000  # Fuel Lower Heating Value [kJ/kg]
#-----------------------------------------------------------------------------#
#=============================================================================#
#-----------------------------------------------------------------------------#

# PME-3480 - Motores de Combustão Interna
# 1D Otto cycle simulator - 2024
# Example file

#-----------------------------------------------------------------------------#
#=============================================================================#
#-----------------------------------------------------------------------------#

#-----------------------------------------------------------------------------#
# Importing libraries
#-----------------------------------------------------------------------------#
import numpy as np # NumPy is a library for numerical operations with multi-dimensional arrays and matrices, and provides mathematical functions.
import OttoCycle as oc  # Imports the OttoCycle module and assigns it the alias 'oc' for easier reference within the code.

#-----------------------------------------------------------------------------#
# Data storage
Ne_l = []  # Effective power [W]
Mt_l = []  # Torque [kgfm]
Nt_l = []  # Thermal power [W]
mpF_l = []  # Fuel mass flow rate [kg/h]
Ce_l = []  # Specific fuel consumption [kg/W·s]
nt_l = []  # Thermal efficiency
j = []  # Iteration index
Ni_l = []  # Indicated power [W]
nm_l = []  # Mechanical efficiency
Na_l = []  # Friction power [W]
ng_l = []  # Global efficiency

#Group parameters:
p = [[9.00, 528.0, 1.431, 2.148],
     [10.00, 518.0, 1.436, 2.227],
     [11.00, 510.0, 1.442, 2.300]]  # rv, Texh [°C], mpF [kg/h], Mt [kgfm]

#-----------------------------------------------------------------------------#
# Loop over cases
#-----------------------------------------------------------------------------#
for i in [0, 1, 2]:
    #-----------------------------------------------------------------------------#
    # engine parameters
    #-----------------------------------------------------------------------------#
    case1 = 'fired'
    
    #-----------------------------------------------------------------------------#
    # Cylinder
    #-----------------------------------------------------------------------------#
    b = 90/1000  # Bore [m]
    s = 90/1000  # Stroke [m]
    l = 120/1000  # Connecting rod length [m]
    rv = p[i][0]  # Compression ratio

    #-----------------------------------------------------------------------------#
    # Valve timing
    #-----------------------------------------------------------------------------#
    ThIVO = +360.*(np.pi/180.)  # Intake valve opening [radians]
    ThIVC = -150.*(np.pi/180.)  # Intake valve closing [radians]
    ThEVO = +150.*(np.pi/180)  # Exhaust valve opening [radians]
    ThEVC = -360.*(np.pi/180)  # Exhaust valve closing [radians]

    #-----------------------------------------------------------------------------#
    # Wiebe function parameters
    #-----------------------------------------------------------------------------#
    ThSOC = -10.*(np.pi/180)  # Start of combustion [radians]
    ThEOC = +10.*(np.pi/180)  # End of combustion [radians]
    aWF = 5  # Wiebe efficiency factor
    mWF = 2  # Wiebe form factor

    #-----------------------------------------------------------------------------#
    # Boundary parameters
    #-----------------------------------------------------------------------------#
    pint = 100e3  # Intake pressure [Pa]
    Tint = 273.15 + 25  # Intake temperature [K]
    pexh = 100e3  # Exhaust pressure [Pa]
    Texh1 = p[i][1] + 273.15  # Exhaust temperature [K]
    phi = 1.0  # Equivalence ratio
    fuel = 'CH4'  # Fuel type

    #-----------------------------------------------------------------------------#
    # Simulator parameters
    #-----------------------------------------------------------------------------#
    pars1 = (case1, b, s, l, rv, n, ThIVO, ThIVC, ThEVO, ThEVC,
             ThSOC, ThEOC, aWF, mWF, pint, Tint, pexh, Texh1, phi, fuel)

    #-----------------------------------------------------------------------------#
    # Crank angle setup
    #-----------------------------------------------------------------------------#
    Th0 = -360.*(np.pi/180)  # Start angle [radians]
    Th1 = +360.*(np.pi/180)  # End angle [radians]
    Ths = 1*(np.pi/180)  # Step angle [radians]
    Thn = int(((Th1 - Th0)/Ths) + 1)  # Number of points
    Th = np.linspace(start=Th0, stop=Th1, num=Thn, endpoint=True)  # Crank angles [radians]
    
    cad = Th*(180/np.pi)  # Crank Angle Degrees [°]

    #-----------------------------------------------------------------------------#
    # Otto Cycle simulation
    #-----------------------------------------------------------------------------#
    v1, m1, t1, p1 = oc.ottoCycle(Th, pars1)  # Run Otto cycle simulation
    
    # Work calculation [J]
    Wliq = np.trapz(p1, v1)  

    #-----------------------------------------------------------------------------#
    # Performance calculations
    #-----------------------------------------------------------------------------#

    # Effective power [W]
    Ne = p[i][3] * 2 * math.pi * n * 9.81
    Ne_l.append(Ne)
    Mt_l.append(p[i][3])

    # Indicated power [W]
    Ni = Zc * (n / x) * Wliq 
    Ni_l.append(Ni)

    # Mechanical efficiency
    nm = Ne / Ni
    nm_l.append(nm)

    # Thermal power [W]
    Nt = (p[i][2] / 3600) * PCI * 1000  # Fuel mass flow rate [kg/s] and PCI [J/kg]
    Nt_l.append(Ne)
    mpF_l.append(p[i][2])

    # Thermal efficiency
    nt = Ni / Nt
    nt_l.append(nt)

    # Specific fuel consumption [kg/W·s]
    Ce = (p[i][2] / 3600) / Ne
    Ce_l.append(Ce)

    # Global efficiency
    ng = Ne/Nt
    ng_l.append(ng)

    # Rendimento volumétrico (fluxo mássico de ar seria necessário para cálculo)

    j.append(i + 1)

#-----------------------------------------------------------------------------#

# Converte Ce_l de kg/kW para g/kWh

Ce_l_g_kWh = Ce_l * 1000 * 3600  # Multiplicando por 1000 e por 3600
Ce_l_g_kWh = []
for i in Ce_l:
    Ce_l_g_kWh.append(i*1000*3600)

output_dir = './Files/S1'

# Gráfico 1: Potência efetiva em função do torque
plt.plot(Mt_l, Ne_l)
plt.title("Potência efetiva em função do torque")
plt.xlabel('Torque (kgfm)')
plt.ylabel('Potência efetiva (kW)')
plt.savefig(os.path.join(output_dir, 'potencia_efetiva.png'))
plt.close()

# Gráfico 2: Potência térmica em função da vazão mássica de combustível
plt.plot(mpF_l, Nt_l)
plt.title("Potência térmica em função da vazão mássica de combustível")
plt.xlabel('Vazão (kg/s)')
plt.ylabel('Potência térmica (kW)')
plt.savefig(os.path.join(output_dir, 'potencia_termica.png'))
plt.close()

# Gráfico 3: Consumo específico em função da vazão mássica de combustível
plt.plot(mpF_l, Ce_l_g_kWh) 
plt.title("Consumo específico em função da vazão mássica de combustível")
plt.xlabel('Vazão (kg/s)')
plt.ylabel('Consumo específico (g/kWh)')  # Atualizando a unidade
plt.savefig(os.path.join(output_dir, 'consumo_especifico.png'))
plt.close()

# Gráfico 4: Rendimento térmico para os três estados
plt.plot(j, [nt * 100 for nt in nt_l])
plt.title("Rendimento térmico para os três estados")
plt.xlabel('Estado')  
plt.ylabel('Rendimento térmico [%]')
plt.savefig(os.path.join(output_dir, 'rendimento_termico.png'))
plt.close()

# Gráfico 5: Potência Indicada para os Três Estados
plt.plot(j, Ni_l)
plt.title("Potência Indicada para os Três Estados")
plt.xlabel('Estado')  
plt.ylabel('Potência Indicada')
plt.savefig(os.path.join(output_dir, 'potencia_indicada.png'))
plt.close()

# Gráfico 6: Rendimento Mecânico
plt.plot(j, [nm * 100 for nm in nm_l])
plt.title("Rendimento Mecânico")
plt.xlabel('Estado')  # Alterado para "Estado"
plt.ylabel('Rendimento Mecânico [%]')
plt.savefig(os.path.join(output_dir, 'rendimento_mecanico.png'))
plt.close()

# Gráfico 7: Rendimento global
plt.plot(j, [ng * 100 for ng in ng_l])
plt.title("Rendimento global")
plt.xlabel('Estado')  # Alterado para "Estado"
plt.ylabel('Rendimento global [%]')
plt.savefig(os.path.join(output_dir, 'rendimento_global.png'))
plt.close()