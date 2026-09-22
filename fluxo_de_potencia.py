# ============================================================
# FLUXO DE POTÊNCIA - VARREDURA BACKWARD / FORWARD
# (método da soma de correntes, Shirmohammadi et al.)
# ============================================================
#
# Convenção adotada:
#   - S_i  = potência da CARGA da barra i (positiva = consumo)
#   - I_i  = (S_i / V_i)*          corrente nodal (Eq. 2 do artigo)
#   - J_l  = I_j + soma(J das linhas que saem de j)   (backward)
#   - V_j  = V_i - Z_ij * J_ij                        (forward)
#
# Como as cargas são consumo, o sinal negativo do artigo
# (Eq. 3) não é necessário aqui: as correntes de linha já
# saem positivas no sentido SE -> cargas.
#
# Sistema em pu (por unidade), fase equivalente monofásica:
#   S_base (kVA) = parâmetro (padrão 10 MVA)
#   V_base (kV)  = tensão do arquivo (linha-linha)
#   Z_base (ohm) = V_base² / S_base
# ============================================================

import math


def fluxo_de_potencia(
    dados,
    Ni_ordenado,
    Nj_ordenado,
    id_ordenado,
    S_base_kVA=10000.0,
    V_SE_pu=1.0,
    tol=1e-5,
    max_iter=100,
):
    """
    Executa o fluxo de potência radial pelo método backward/forward.

    Parâmetros
    ----------
    dados        : dicionário retornado por ler_arquivo()
    Ni_ordenado, Nj_ordenado, id_ordenado : saída de ordenar_linhas()
    S_base_kVA   : potência base em kVA
    V_SE_pu      : tensão da subestação em pu (módulo, ângulo 0)
    tol          : tolerância do erro de potência (pu)
    max_iter     : número máximo de iterações

    Retorna
    -------
    Dicionário com tensões, correntes, perdas e informações de convergência.
    """

    # --------------------------------------------------------
    # 1. Bases do sistema
    # --------------------------------------------------------

    kV_base = dados["voltage_kv"]

    # kV² / MVA = ohm
    Z_base = (kV_base ** 2) / (S_base_kVA / 1000.0)

    # kVA / kV = A   (corrente base do sistema monofásico)
    I_base_A = S_base_kVA / ( kV_base)


    # --------------------------------------------------------
    # 2. Encontrar a subestação e montar S (pu) de cada barra
    # --------------------------------------------------------

    SE = None
    S = {}

    for barra in dados["barras"]:

        if barra["type"] == "SE" and SE is None:
            SE = barra["id"]

        # Potência da carga em pu (kW + j kvar) / S_base
        S[barra["id"]] = complex(barra["SaReal"], barra["SaImag"]) / S_base_kVA

    if SE is None:
        raise ValueError("Nenhuma barra do tipo SE foi encontrada.")

    # A SE é a barra de referência (slack): não tem carga no cálculo
    S[SE] = 0j


    # --------------------------------------------------------
    # 3. Impedâncias das linhas em pu (indexadas pelo id)
    # --------------------------------------------------------

    Z = {}

    for linha in dados["linhas"]:
        Z[linha["id"]] = complex(linha["ZaaReal"], linha["ZaaImag"]) / Z_base


    # --------------------------------------------------------
    # 4. Inicialização
    # --------------------------------------------------------

    n_linhas = len(Ni_ordenado)

    # Tensão inicial de todas as barras = tensão da SE (flat start)
    V = {barra["id"]: complex(V_SE_pu, 0.0) for barra in dados["barras"]}

    # Correntes nas linhas (mesma ordem das listas ordenadas)
    J = [0j] * n_linhas

    convergiu = False
    erro_max = float("inf")
    iteracao = 0


    # --------------------------------------------------------
    # 5. Processo iterativo
    # --------------------------------------------------------

    for iteracao in range(1, max_iter + 1):

        # ----------------------------------------------------
        # PASSO 1: Corrente nodal  I_i = (S_i / V_i)*
        # ----------------------------------------------------

        I_no = {}

        for b in V:
            if b == SE:
                I_no[b] = 0j
            else:
                I_no[b] = (S[b] / V[b]).conjugate()


        # ----------------------------------------------------
        # PASSO 2: Backward - soma das correntes (folhas -> SE)
        # Percorre a ordem de trás para frente: os filhos de um
        # nó sempre aparecem depois dele na ordem BFS, então
        # quando chegamos na linha k, tudo que está a jusante
        # de Nj já foi somado em soma_jusante[Nj].
        # ----------------------------------------------------

        soma_jusante = {b: 0j for b in V}

        for k in reversed(range(n_linhas)):

            i = Ni_ordenado[k]
            j = Nj_ordenado[k]

            # Corrente na linha = corrente da carga em j + tudo a jusante
            J[k] = I_no[j] + soma_jusante[j]

            # Repassa essa corrente para o nó de montante
            soma_jusante[i] += J[k]


        # ----------------------------------------------------
        # PASSO 3: Forward - atualização das tensões (SE -> folhas)
        #   V_j = V_i - Z_ij * J_ij
        # ----------------------------------------------------

        V[SE] = complex(V_SE_pu, 0.0)

        for k in range(n_linhas):

            i = Ni_ordenado[k]
            j = Nj_ordenado[k]

            V[j] = V[i] - Z[id_ordenado[k]] * J[k]


        # ----------------------------------------------------
        # PASSO 4: Erro de potência (Eq. 5 do artigo)
        #   dS_i = V_i(k) * (I_i(k))* - S_i
        # Critério: parte real e imaginária <= tol (pu)
        # ----------------------------------------------------

        erro_max = 0.0

        for b in V:

            if b == SE:
                continue

            dS = V[b] * I_no[b].conjugate() - S[b]

            erro_max = max(erro_max, abs(dS.real), abs(dS.imag))

        if erro_max <= tol:
            convergiu = True
            break


    # --------------------------------------------------------
    # 6. Pós-processamento: perdas
    # --------------------------------------------------------

    perdas_pu = 0j

    for k in range(n_linhas):
        perdas_pu += Z[id_ordenado[k]] * abs(J[k]) ** 2

    # Perdas em kW + j kvar
    perdas_kVA = perdas_pu * S_base_kVA

    # Potência total fornecida pela SE: soma das potências que saem dela
    S_SE = 0j
    for k in range(n_linhas):
        if Ni_ordenado[k] == SE:
            S_SE += V[SE] * J[k].conjugate()
    S_SE_kVA = S_SE * S_base_kVA


    # --------------------------------------------------------
    # 7. Resultado
    # --------------------------------------------------------

    return {
        "convergiu": convergiu,
        "iteracoes": iteracao,
        "erro_max": erro_max,
        "SE": SE,
        "V": V,                           # tensões complexas em pu
        "J": J,                           # correntes de linha em pu (ordem ordenada)
        "Ni": Ni_ordenado,
        "Nj": Nj_ordenado,
        "id": id_ordenado,
        "perdas_kW": perdas_kVA.real,
        "perdas_kvar": perdas_kVA.imag,
        "S_SE_kW": S_SE_kVA.real,
        "S_SE_kvar": S_SE_kVA.imag,
        "kV_base": kV_base,
        "I_base_A": I_base_A,
        "S_base_kVA": S_base_kVA,
    }


# ============================================================
# IMPRESSÃO DOS RESULTADOS
# ============================================================

def imprimir_resultados(res):

    print("\n" + "=" * 60)
    print("RESULTADO DO FLUXO DE POTÊNCIA")
    print("=" * 60)

    if res["convergiu"]:
        print(f"Convergiu em {res['iteracoes']} iterações")
    else:
        print(f"NÃO convergiu após {res['iteracoes']} iterações")

    print(f"Erro máximo de potência : {res['erro_max']:.3e} pu")
    print(f"Potência base           : {res['S_base_kVA']:.0f} kVA")
    print(f"Tensão base             : {res['kV_base']} kV")


    # ---------------- Tensões ----------------

    print("\n" + "=" * 60)
    print("TENSÕES NAS BARRAS")
    print("=" * 60)
    print(f"{'barra':>6} {'|V| (pu)':>10} {'ang (graus)':>12} {'|V| (kV)':>10}")

    for b in sorted(res["V"]):
        v = res["V"][b]
        print(
            f"{b:>6} "
            f"{abs(v):>10.5f} "
            f"{math.degrees(math.atan2(v.imag, v.real)):>12.4f} "
            f"{abs(v) * res['kV_base']:>10.4f}"
        )

    # Menor tensão do sistema
    barra_min = min(res["V"], key=lambda b: abs(res["V"][b]))
    print(
        f"\nMenor tensão: barra {barra_min} "
        f"= {abs(res['V'][barra_min]):.5f} pu"
    )


    # ---------------- Correntes ----------------

    print("\n" + "=" * 60)
    print("CORRENTES NAS LINHAS")
    print("=" * 60)
    print(f"{'id':>4} {'Ni':>4} {'Nj':>4} {'|J| (pu)':>12} {'|J| (A)':>10}")

    for k in range(len(res["J"])):
        j_pu = abs(res["J"][k])
        print(
            f"{res['id'][k]:>4} "
            f"{res['Ni'][k]:>4} "
            f"{res['Nj'][k]:>4} "
            f"{j_pu:>12.6f} "
            f"{j_pu * res['I_base_A']:>10.3f}"
        )


    # ---------------- Perdas ----------------

    print("\n" + "=" * 60)
    print("BALANÇO DE POTÊNCIA")
    print("=" * 60)
    print(f"Potência da SE : {res['S_SE_kW']:.3f} kW + j{res['S_SE_kvar']:.3f} kvar")
    print(f"Perdas totais  : {res['perdas_kW']:.3f} kW + j{res['perdas_kvar']:.3f} kvar")
    