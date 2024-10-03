import numpy as np
import cmath
import random

#Valores das distâncias das Linhas de transmissão (km)
dist_LT1 = 80
dist_LT2 = 80
dist_LT3 = 80
dist_LT4 = 120
dist_LT5 = 120
dist_LT6 = 100

# Impedancia serie Thévenin
Rf = 4
Xf = 0.38
Zf = Rf + 1j * Xf

w = 2 * np.pi * 60

# Cargas
Z1 = 8000 + 1j * w * 42
Z2 = 1350.55 + 1j * w * 7.83
Z3 = 649 + 1j * w * 3.2


def Linha_Transmissao(Distancia):
    w = 2 * np.pi * 60  # rad/s

    # R(ohm), L(H), Cap(F)
    R = 0.172 * Distancia
    L = 2.18E-3 * Distancia
    Cap = 0.0136E-6 * Distancia
    # Parametros da Matriz de Transmissão da Linha
    A = 1 + (1j * w * Cap * R - (w ** 2) * L * Cap) / 2
    B = R + 1j * w * L
    C = 1j * w * Cap + ((R + 1j * w * L) * ((w * Cap) ** 2)) / 4
    D = 1 + (1j * w * Cap * R - (w ** 2) * L * Cap) / 2

    matriz = np.array(
        [
            [A, B],
            [C, D]
        ]
    )
    return matriz


def Transformador(Trafo):
    R1 = 7.6E-3
    X1 = 3.8E-3
    R2 = 33.9E-3
    X2 = 0.85E-3
    v1 = Trafo[0]
    v2 = Trafo[1]
    Rm = Trafo[2]
    Xm = Trafo[3]

    Z1 = R1 + 1j * X1
    Z2 = R2 + 1j * X2
    Y = (Rm + 1j * Xm) / (1j * Rm * Xm)
    N = v1 / v2
    # Parametros da matriz de trasmissão do Traformador
    A = N * (1 + Y * Z1)
    B = (1 / N) * (Z1 + Z2 + Y * Z1 * Z2)
    C = N * Y
    D = (1 / N) * (1 + Y * Z2)

    matriz = np.array(
        [
            [A, B],
            [C, D]
        ]
    )
    return matriz


def Carga(Z):
    # Parametros da matriz de trasmissão da Carga
    A = 1
    B = Z
    C = 0
    D = 1

    matriz = np.array(
        [
            [A, B],
            [C, D]
        ]
    )

    return matriz


def Carga_Derivada(Z):
    # Parametros da matriz de trasmissão da Carga Derivada
    A = 1
    B = 0
    C = 1 / Z
    D = 1

    matriz = np.array(
        [
            [A, B],
            [C, D]
        ]
    )

    return matriz


def Matrizes_Cascata(m1, m2):
    # Parametros para matriz de transmissão da associção em cascata
    matriz = np.dot(m1, m2)
    return matriz


def Matrizes_Paralelo(m1, m2):
    A1 = m1[0][0]
    B1 = m1[0][1]
    C1 = m1[1][0]
    D1 = m1[1][1]
    A2 = m2[0][0]
    B2 = m2[0][1]
    C2 = m2[1][0]
    D2 = m2[1][1]
    # Parametros para matriz de transmissão da associção em paralelo
    A = (A1 * B2 + A2 * B1) / (B1 + B2)
    B = (B1 * B2) / (B1 + B2)
    C = C1 + C2 + ((A1 - A2) * (D2 - D1)) / (B1 + B2)
    D = (B2 * D1 + B1 * D2) / (B1 + B2)

    matriz = np.array(
        [
            [A, B],
            [C, D]
        ]
    )

    return matriz


resposta = input("Quer calcular com ajuste de TAP ou Compensão ( adicionar Reatores) para tensões em Z1, Z2 e Z3 seja igual, respectivamente, a 500kV, 230kV e 69kV?"
                 "\nAjuste de TAP - digite 1 \nAjuste com Reatores - digite 2 \nNenhum - digite qualquer outro número\n ")

# Converte a resposta para inteiro
try:
    resposta = int(resposta)
except ValueError:
    resposta = None  # Caso não digite um número válido

# Exibe a resposta capturada
if resposta == 1:
    print("Vamos calcular o ajuste de TAP.")
    vaVz1 = 0
    vaVz2 = 0
    vaVz3 = 0
    # O programa entra em loop até que seja obtido os valores de tensões nos intervalos desejados
    while not (500000 < vaVz1 < 500300 and 230000 < vaVz2 < 230300 and vaVz3 == 69000):
        # Trasformadores = [v1,v2,Rm,Xm], com v1 e v2 em kV, com Rm e Xm em Ohm
        # COM AJUSTE DE TAP
        Trafo1 = [69000, random.randint(490000, 510000), 4320, 5050]
        Trafo2 = [500000, random.randint(220000, 240000), 432000, 505000]
        Trafo3 = [230000, random.randint(67000, 70000), 402000, 607000]

        # Matrizes de Transmissão para cada componente do sistema
        zth = Carga(Zf)

        T1 = Transformador(Trafo1)
        LT1 = Linha_Transmissao(dist_LT1)
        LT2 = Linha_Transmissao(dist_LT2)

        z1 = Carga_Derivada(Z1)
        LT3 = Linha_Transmissao(dist_LT3)

        T2 = Transformador(Trafo2)
        z2 = Carga_Derivada(Z2)
        LT4 = Linha_Transmissao(dist_LT4)

        T3 = Transformador(Trafo3)
        z3 = Carga_Derivada(Z3)

        # Modelo do Sistema
        matriz1 = Matrizes_Cascata(zth, T1)
        matriz2 = Matrizes_Cascata(matriz1, Matrizes_Paralelo(LT1, LT2))
        matriz3 = Matrizes_Cascata(matriz2, z1)
        matriz4 = Matrizes_Cascata(matriz3, LT3)
        matriz5 = Matrizes_Cascata(matriz4, T2)
        matriz6 = Matrizes_Cascata(matriz5, z2)
        matriz7 = Matrizes_Cascata(matriz6, LT4)
        matriz8 = Matrizes_Cascata(matriz7, T3)

        # Parte 1: Encontrar V,I,P e Q na Entrada do circuito
        Vz3 = 69E3
        Iz3 = Vz3 / Z3
        Sz3 = Vz3 * (Iz3.conjugate())
        Pz3 = Sz3.real
        Qz3 = Sz3.imag
        absIz3, angIz3 = cmath.polar(Iz3)
        angIz3 = np.rad2deg(angIz3)
        matriz9 = [Vz3, Iz3]
        matrizentrada = Matrizes_Cascata(matriz8, matriz9)
        Vac = matrizentrada[0]
        Iac = matrizentrada[1]
        Sac = Vac * (Iac.conjugate())
        Pac = Sac.real
        Qac = Sac.imag
        absVac, angVac = cmath.polar(Vac)
        absIac, angIac = cmath.polar(Iac)
        angVac = np.rad2deg(angVac)
        angIac = np.rad2deg(angIac)

        # Parte 2: Encontrar V,I,P e Q em Z1
        Vz1 = (Vac * Z1) / (matriz2[0][0] * Z1 + matriz2[0][1])
        Iz1 = Vz1 / Z1
        Sz1 = Vz1 * (Iz1.conjugate())
        Pz1 = Sz1.real
        Qz1 = Sz1.imag
        absVz1, angVz1 = cmath.polar(Vz1)
        absIz1, angIz1 = cmath.polar(Iz1)
        angVz1 = np.rad2deg(angVz1)
        angIz1 = np.rad2deg(angIz1)

        # Parte 3: Encontrar V,I,P e Q em Z2
        Vz2 = (Vac * Z2) / (matriz5[0][0] * Z2 + matriz5[0][1])
        Iz2 = Vz2 / Z2
        Sz2 = Vz2 * (Iz2.conjugate())
        Pz2 = Sz2.real
        Qz2 = Sz2.imag
        absVz2, angVz2 = cmath.polar(Vz2)
        absIz2, angIz2 = cmath.polar(Iz2)
        angVz2 = np.rad2deg(angVz2)
        angIz2 = np.rad2deg(angIz2)

        # Calculo de Perdas de Potência Ativa e Reativa
        S_totalcargas = Sz2 + Sz1 + Sz3
        S_perdas = Sac - S_totalcargas

        # Variaveis para verificar no While
        vaVz1 = absVz1
        vaVz2 = absVz2
        vaVz3 = Vz3



    print(f"T = {matriz8}")
    print()
    print(f"Vz3 = {Vz3} V")
    print(f"Iz3 = {Iz3} A")
    print(f"Iz3 = {absIz3:.2f} ∠ {angIz3:.2f}° V")
    print(f"Potência Ativa em Z3 = {Pz3:.2f}W\nPotência Reativa em Z3 = {Qz3:.2f}VAR")
    print("Tensão e corrente na entrada:")
    print(f"Vac = {Vac} V")
    print(f"Iac = {Iac} A")
    print(f"Vac = {absVac:.2f} ∠ {angVac:.2f}° V")
    print(f"Iac = {absIac:.2f} ∠ {angIac:.2f}° A")
    print(f"Potência Ativa em Vac = {Pac:.2f}W\nPotência Reativa em Vac = {Qac:.2f}VAR")
    print("Tensão e corrente na carga Z1:")
    print(f"Vz1 = {Vz1} V")
    print(f"Iz1 = {Iz1} A")
    print(f"Vz1 = {absVz1:.2f}∠{angVz1:.2f}° V")
    print(f"Iz1 = {absIz1:.2f}∠{angIz1:.2f}° A")
    print(f"Potência Ativa em Z1 = {Pz1:.2f}W\nPotência Reativa em Z1 = {Qz1:.2f}VAR")
    print("Tensão e corrente na carga Z2:")
    print(f"Vz2 = {Vz2} V")
    print(f"Iz2 = {Iz2} A")
    print(f"Vz2 = {absVz2:.2f}∠{angVz2:.2f}° V")
    print(f"Iz2 = {absIz2:.2f}∠{angIz2:.2f}° A")
    print(f"Potência Ativa em Z2 = {Pz2:.2f}W\nPotência Reativa em Z2 = {Qz2:.2f}VAR")
    print(f"Potência Ativa de Perdas = {S_perdas.real:.2f}W\nPotência Reativa de Perdas = {S_perdas.imag:.2f}VAR")

    print()
    print(f"v2 do ajuste de TAP do primeiro Transformador {Trafo1[1]}V")
    print(f"v2 do ajuste de TAP do primeiro Transformador {Trafo2[1]}V")
    print(f"v2 do ajuste de TAP do primeiro Transformador {Trafo3[1]}V")

elif resposta == 2:
    print("Vamos calcular o ajuste adicionando os Reatores.")
    vaVz1 = 0
    vaVz2 = 0
    vaVz3 = 0
    # O programa entra em loop até que seja obtido os valores de tensões nos intervalos desejados
    while not (475000 < vaVz1 < 525000 and 218500 < vaVz2 < 241500 and vaVz3 == 69000):
        # Trasformadores = [v1,v2,Rm,Xm], com v1 e v2 em kV, com Rm e Xm em Ohm
        # COM AJUSTE ADICIONANDO OS REATORES
        Trafo1 = [69000, 500000, 4320, 5050]
        Trafo2 = [500000, 230000, 432000, 505000]
        Trafo3 = [230000, 69000, 402000, 607000]

        L1 = random.randint(1, 50)
        L2 = random.randint(1, 50)
        L3 = random.randint(1, 50)

        reator1 = 1j * w * (L1)
        reator2 = 1j * w * (L2)
        reator3 = 1j * w * (L3)

        # Matrizes de Transmissão para cada componente do sistema
        zth = Carga(Zf)

        zreator1 = Carga_Derivada(reator1)
        zreator2 = Carga_Derivada(reator2)
        zreator3 = Carga_Derivada(reator3)

        T1 = Transformador(Trafo1)
        LT1 = Linha_Transmissao(dist_LT1)
        LT2 = Linha_Transmissao(dist_LT2)

        z1 = Carga_Derivada(Z1)
        LT3 = Linha_Transmissao(dist_LT3)

        T2 = Transformador(Trafo2)
        z2 = Carga_Derivada(Z2)
        LT4 = Linha_Transmissao(dist_LT4)

        T3 = Transformador(Trafo3)
        z3 = Carga_Derivada(Z3)

        # Modelo do Sistema
        matriz1 = Matrizes_Cascata(zth, T1)
        matriz2 = Matrizes_Cascata(matriz1, Matrizes_Paralelo(LT1, LT2))
        matrizreator1 = Matrizes_Cascata(matriz2, zreator1)
        matriz3 = Matrizes_Cascata(matrizreator1, z1)
        matriz4 = Matrizes_Cascata(matriz3, LT3)
        matriz5 = Matrizes_Cascata(matriz4, T2)
        matrizreator2 = Matrizes_Cascata(matriz5, zreator2)
        matriz6 = Matrizes_Cascata(matrizreator2, z2)
        matriz7 = Matrizes_Cascata(matriz6, LT4)
        matriz8 = Matrizes_Cascata(matriz7, T3)
        matrizreator3 = Matrizes_Cascata(matriz8, zreator3)

        # Parte 1: Encontrar V,I,P e Q na Entrada do circuito
        Vz3 = 69E3
        Iz3 = Vz3 / Z3
        Sz3 = Vz3 * (Iz3.conjugate())
        Pz3 = Sz3.real
        Qz3 = Sz3.imag
        absIz3, angIz3 = cmath.polar(Iz3)
        angIz3 = np.rad2deg(angIz3)
        matriz9 = [Vz3, Iz3]
        matrizentrada = Matrizes_Cascata(matrizreator3, matriz9)
        Vac = matrizentrada[0]
        Iac = matrizentrada[1]
        Sac = Vac * (Iac.conjugate())
        Pac = Sac.real
        Qac = Sac.imag
        absVac, angVac = cmath.polar(Vac)
        absIac, angIac = cmath.polar(Iac)
        angVac = np.rad2deg(angVac)
        angIac = np.rad2deg(angIac)

        # Parte 2: Encontrar V,I,P e Q em Z1
        Vz1 = (Vac * Z1) / (matrizreator1[0][0] * Z1 + matrizreator1[0][1])
        Iz1 = Vz1 / Z1
        Sz1 = Vz1 * (Iz1.conjugate())
        Pz1 = Sz1.real
        Qz1 = Sz1.imag
        absVz1, angVz1 = cmath.polar(Vz1)
        absIz1, angIz1 = cmath.polar(Iz1)
        angVz1 = np.rad2deg(angVz1)
        angIz1 = np.rad2deg(angIz1)

        # Parte 3: Encontrar V,I,P e Q em Z2
        Vz2 = (Vac * Z2) / (matrizreator2[0][0] * Z2 + matrizreator2[0][1])
        Iz2 = Vz2 / Z2
        Sz2 = Vz2 * (Iz2.conjugate())
        Pz2 = Sz2.real
        Qz2 = Sz2.imag
        absVz2, angVz2 = cmath.polar(Vz2)
        absIz2, angIz2 = cmath.polar(Iz2)
        angVz2 = np.rad2deg(angVz2)
        angIz2 = np.rad2deg(angIz2)
        # Calculo de Perdas de Potência Ativa e Reativa
        S_totalcargas = Sz2 + Sz1 + Sz3
        S_perdas = Sac - S_totalcargas

        # Variaveis para verificar no While
        vaVz1 = absVz1
        vaVz2 = absVz2
        vaVz3 = Vz3

    print(f"T = {matriz8}")
    print()
    print(f"Vz3 = {Vz3} V")
    print(f"Iz3 = {Iz3} A")
    print(f"Iz3 = {absIz3:.2f} ∠ {angIz3:.2f}° V")
    print(f"Potência Ativa em Z3 = {Pz3:.2f}W\nPotência Reativa em Z3 = {Qz3:.2f}VAR")
    print("Tensão e corrente na entrada:")
    print(f"Vac = {Vac} V")
    print(f"Iac = {Iac} A")
    print(f"Vac = {absVac:.2f} ∠ {angVac:.2f}° V")
    print(f"Iac = {absIac:.2f} ∠ {angIac:.2f}° A")
    print(f"Potência Ativa em Vac = {Pac:.2f}W\nPotência Reativa em Vac = {Qac:.2f}VAR")
    print("Tensão e corrente na carga Z1:")
    print(f"Vz1 = {Vz1} V")
    print(f"Iz1 = {Iz1} A")
    print(f"Vz1 = {absVz1:.2f}∠{angVz1:.2f}° V")
    print(f"Iz1 = {absIz1:.2f}∠{angIz1:.2f}° A")
    print(f"Potência Ativa em Z1 = {Pz1:.2f}W\nPotência Reativa em Z1 = {Qz1:.2f}VAR")
    print("Tensão e corrente na carga Z2:")
    print(f"Vz2 = {Vz2} V")
    print(f"Iz2 = {Iz2} A")
    print(f"Vz2 = {absVz2:.2f}∠{angVz2:.2f}° V")
    print(f"Iz2 = {absIz2:.2f}∠{angIz2:.2f}° A")
    print(f"Potência Ativa em Z2 = {Pz2:.2f}W\nPotência Reativa em Z2 = {Qz2:.2f}VAR")
    print(f"Potência Ativa de Perdas = {S_perdas.real:.2f}W\nPotência Reativa de Perdas = {S_perdas.imag:.2f}VAR")
    print()
    print(f"Valor para o Reator acoplado em Z1: Zreator1 = {reator1}Ohm e o valor da sua indutâcia L = {L1}H")
    print(f"Valor para o Reator acoplado em Z2: Zreator2 = {reator2}Ohm e o valor da sua indutâcia L = {L2}H")
    print(f"Valor para o Reator acoplado em Z3: Zreator3 = {reator3}Ohm e o valor da sua indutâcia L = {L3}H")


else:
    print("Você escolheu NÃO, cálculo de TAP não será realizado e nem o cálculo com Reator.")
    # Trasformadores = [v1,v2,Rm,Xm], com v1 e v2 em kV, com Rm e Xm em Ohm
    # SEM AJUSTE DE TAP
    Trafo1 = [69000, 500000, 4320, 5050]
    Trafo2 = [500000, 230000, 432000, 505000]
    Trafo3 = [230000, 69000, 402000, 607000]
    # Matrizes de Transmissão para cada componente do sistema
    zth = Carga(Zf)

    T1 = Transformador(Trafo1)
    LT1 = Linha_Transmissao(dist_LT1)
    LT2 = Linha_Transmissao(dist_LT2)

    z1 = Carga_Derivada(Z1)
    LT3 = Linha_Transmissao(dist_LT3)

    T2 = Transformador(Trafo2)
    z2 = Carga_Derivada(Z2)
    LT4 = Linha_Transmissao(dist_LT4)

    T3 = Transformador(Trafo3)
    z3 = Carga_Derivada(Z3)

    # Modelo do Sistema
    matriz1 = Matrizes_Cascata(zth, T1)
    matriz2 = Matrizes_Cascata(matriz1, Matrizes_Paralelo(LT1, LT2))
    matriz3 = Matrizes_Cascata(matriz2, z1)
    matriz4 = Matrizes_Cascata(matriz3, LT3)
    matriz5 = Matrizes_Cascata(matriz4, T2)
    matriz6 = Matrizes_Cascata(matriz5, z2)
    matriz7 = Matrizes_Cascata(matriz6, LT4)
    matriz8 = Matrizes_Cascata(matriz7, T3)
    print(f"T = {matriz8}")
    print()

    # Parte 1: Encontrar V,I,P e Q na Entrada do circuito
    print("Tensão e corrente na carga Z3:")
    Vz3 = 69E3
    print(f"Vz3 = {Vz3} V")
    Iz3 = Vz3 / Z3
    Sz3 = Vz3 * (Iz3.conjugate())
    Pz3 = Sz3.real
    Qz3 = Sz3.imag
    absIz3, angIz3 = cmath.polar(Iz3)
    angIz3 = np.rad2deg(angIz3)
    print(f"Iz3 = {absIz3:.2f} ∠ {angIz3:.2f}° V")
    matriz9 = [Vz3, Iz3]
    print(f"Iz3 = {Iz3} A")
    print(f"Potência Ativa em Z3 = {Pz3:.2f}W\nPotência Reativa em Z3 = {Qz3:.2f}VAR\n")
    print("Tensão e corrente na entrada:")
    matrizentrada = Matrizes_Cascata(matriz8, matriz9)
    Vac = matrizentrada[0]
    Iac = matrizentrada[1]
    Sac = Vac * (Iac.conjugate())
    Pac = Sac.real
    Qac = Sac.imag
    print(f"Vac = {Vac} V")
    print(f"Iac = {Iac} A")
    absVac, angVac = cmath.polar(Vac)
    absIac, angIac = cmath.polar(Iac)
    angVac = np.rad2deg(angVac)
    angIac = np.rad2deg(angIac)
    print(f"Vac = {absVac:.2f} ∠ {angVac:.2f}° V")
    print(f"Iac = {absIac:.2f} ∠ {angIac:.2f}° A")
    print(f"Potência Ativa em Vac = {Pac:.2f}W\nPotência Reativa em Vac = {Qac:.2f}VAR\n")

    # Parte 2: Encontrar V,I,P e Q em Z1
    print("Tensão e corrente na carga Z1:")
    Vz1 = (Vac * Z1) / (matriz2[0][0] * Z1 + matriz2[0][1])
    print(f"Vz1 = {Vz1} V")
    Iz1 = Vz1 / Z1
    Sz1 = Vz1 * (Iz1.conjugate())
    Pz1 = Sz1.real
    Qz1 = Sz1.imag
    print(f"Iz1 = {Iz1} A")
    absVz1, angVz1 = cmath.polar(Vz1)
    absIz1, angIz1 = cmath.polar(Iz1)
    angVz1 = np.rad2deg(angVz1)
    angIz1 = np.rad2deg(angIz1)
    print(f"Vz1 = {absVz1:.2f}∠{angVz1:.2f}° V")
    print(f"Iz1 = {absIz1:.2f}∠{angIz1:.2f}° A")
    print(f"Potência Ativa em Z1 = {Pz1:.2f}W\nPotência Reativa em Z1 = {Qz1:.2f}VAR\n")

    # Parte 3: Encontrar V,I,P e Q em Z2
    print("Tensão e corrente na carga Z2:")
    Vz2 = (Vac * Z2) / (matriz5[0][0] * Z2 + matriz5[0][1])
    print(f"Vz2 = {Vz2} V")
    Iz2 = Vz2 / Z2
    print(f"Iz2 = {Iz2} A")
    Sz2 = Vz2*(Iz2.conjugate())
    Pz2 = Sz2.real
    Qz2 = Sz2.imag
    absVz2, angVz2 = cmath.polar(Vz2)
    absIz2, angIz2 = cmath.polar(Iz2)
    angVz2 = np.rad2deg(angVz2)
    angIz2 = np.rad2deg(angIz2)
    print(f"Vz2 = {absVz2:.2f}∠{angVz2:.2f}° V")
    print(f"Iz2 = {absIz2:.2f}∠{angIz2:.2f}° A")
    print(f"Potência Ativa em Z2 = {Pz2:.2f}W\nPotência Reativa em Z2 = {Qz2:.2f}VAR\n")
    # Calculo de Perdas de Potência Ativa e Reativa
    S_totalcargas = Sz2 + Sz1 + Sz3
    S_perdas = Sac - S_totalcargas
    print(f"Potência Ativa de Perdas = {S_perdas.real:.2f}W\nPotência Reativa de Perdas = {S_perdas.imag:.2f}VAR")

