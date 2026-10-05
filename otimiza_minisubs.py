# ============================================================
# ALOCAÇÃO DE MINISUBESTAÇÕES - SIMULATED ANNEALING
# ============================================================
#
# Problema MESTRE  : escolher em quais barras instalar as minisubestações
# Problema ESCRAVO : fluxo de potência backward/forward (calcula as perdas)
#
# Função objetivo:
#   custo = N * preco_instalacao + preco_kW * perdas_kW  (+ penalidade de tensão)
#
# Solução = conjunto de N barras distintas (a SE não é candidata).
# Vizinhança = trocar UMA barra escolhida por outra livre.
# ============================================================

import math
import random
from fluxo_de_potencia import fluxo_de_potencia


# ------------------------------------------------------------
# Entrada de dados do usuário (chamada pelo main.py)
# ------------------------------------------------------------

def _ler(msg, tipo=float, minimo=None):
    while True:
        try:
            v = tipo(input(msg).replace(",", "."))
            if minimo is not None and v < minimo:
                print(f"  Valor deve ser >= {minimo}")
                continue
            return v
        except ValueError:
            print("  Entrada inválida, tente novamente.")


def pedir_parametros():
    print("=" * 60)
    print("PARÂMETROS DAS MINISUBESTAÇÕES")
    print("=" * 60)
    n = _ler("Quantas minisubestações deseja inserir? ", int, 1)
    preco_inst = _ler("Preço de instalação de CADA uma (R$): ", float, 0)
    p_kw = _ler("Potência ativa de cada uma (kW): ", float, 0)
    q_kvar = _ler("Potência reativa de cada uma (kvar): ", float)
    preco_kw = _ler("Preço por kW perdido (R$/kW): ", float, 0)
    return n, preco_inst, p_kw, q_kvar, preco_kw


# ------------------------------------------------------------
# Função objetivo
# ------------------------------------------------------------

def avaliar(solucao, ctx, cache):
    """Retorna (custo_total, resultado_fluxo). Usa cache p/ não recalcular."""
    chave = frozenset(solucao)
    if chave in cache:
        return cache[chave]

    injecoes = {b: ctx["S_mini"] for b in solucao}
    res = fluxo_de_potencia(
        ctx["dados"], ctx["Ni"], ctx["Nj"], ctx["id"],
        S_base_kVA=ctx["S_base"], injecoes=injecoes,
    )

    custo = len(solucao) * ctx["preco_inst"] + ctx["preco_kw"] * res["perdas_kW"]

    # Penalidade opcional por violar a faixa de tensão (0,95 - 1,05 pu)
    if ctx["v_min"] is not None:
        for b, v in res["V"].items():
            viol = max(ctx["v_min"] - abs(v), abs(v) - ctx["v_max"], 0.0)
            custo += ctx["peso_tensao"] * viol

    if not res["convergiu"]:
        custo += 1e12  # solução inviável

    cache[chave] = (custo, res)
    return cache[chave]


# ------------------------------------------------------------
# Simulated Annealing
# ------------------------------------------------------------

def simulated_annealing(ctx, n, candidatas, T0=None, alfa=0.95,
                        iter_por_temp=20, T_final=1e-3, semente=None):
    rng = random.Random(semente)
    cache = {}

    atual = rng.sample(candidatas, n)
    custo_atual, _ = avaliar(atual, ctx, cache)
    melhor, custo_melhor = list(atual), custo_atual

    # Temperatura inicial proporcional ao custo (aceita ~piora de 10%)
    T = T0 if T0 else 0.1 * custo_atual
    T_final = T_final * custo_atual / 1000.0 if T_final < 1 else T_final

    historico = [custo_melhor]

    while T > T_final:
        for _ in range(iter_por_temp):
            vizinho = list(atual)
            livres = [b for b in candidatas if b not in atual]
            if not livres:
                break
            vizinho[rng.randrange(n)] = rng.choice(livres)

            custo_viz, _ = avaliar(vizinho, ctx, cache)
            delta = custo_viz - custo_atual

            if delta < 0 or rng.random() < math.exp(-delta / T):
                atual, custo_atual = vizinho, custo_viz
                if custo_atual < custo_melhor:
                    melhor, custo_melhor = list(atual), custo_atual

        historico.append(custo_melhor)
        T *= alfa

    return sorted(melhor), custo_melhor, historico, cache


# ------------------------------------------------------------
# Execução completa (reinícios) + relatório de ONDE instalar
# ------------------------------------------------------------

def otimizar_minisubestacoes(dados, Ni, Nj, ids, n, preco_inst, p_kw,
                             q_kvar, preco_kw, S_base=10000.0,
                             v_min=0.95, v_max=1.05, reinicios=5):
    """Retorna dicionário com as barras escolhidas e comparação antes/depois."""

    ctx = {
        "dados": dados, "Ni": Ni, "Nj": Nj, "id": ids, "S_base": S_base,
        "S_mini": complex(p_kw, q_kvar),
        "preco_inst": preco_inst, "preco_kw": preco_kw,
        "v_min": v_min, "v_max": v_max,       # v_min=None desativa a restrição
        "peso_tensao": 1e6 * max(preco_kw, 1.0),
    }

    SE = next(b["id"] for b in dados["barras"] if b["type"] == "SE")
    candidatas = [b["id"] for b in dados["barras"] if b["id"] != SE]
    if n > len(candidatas):
        raise ValueError("Mais minisubestações do que barras disponíveis.")

    base = fluxo_de_potencia(dados, Ni, Nj, ids, S_base_kVA=S_base)

    melhor, custo_melhor, cache_melhor = None, float("inf"), None
    for r in range(reinicios):
        sol, c, _, cache = simulated_annealing(ctx, n, candidatas, semente=r)
        if c < custo_melhor:
            melhor, custo_melhor, cache_melhor = sol, c, cache

    _, res = avaliar(melhor, ctx, cache_melhor)

    return {
        "barras": melhor, "base": base, "res": res, "ctx": ctx,
        "custo_total": custo_melhor,
        "custo_instalacao": n * preco_inst,
        "custo_perdas_antes": preco_kw * base["perdas_kW"],
        "custo_perdas_depois": preco_kw * res["perdas_kW"],
    }


def imprimir_otimizacao(o):
    base, res, ctx = o["base"], o["res"], o["ctx"]
    S = ctx["S_mini"]

    print("\n" + "=" * 60)
    print("ONDE INSTALAR AS MINISUBESTAÇÕES")
    print("=" * 60)
    print(f"{'barra':>6} {'P (kW)':>10} {'Q (kvar)':>10} "
          f"{'|V| antes':>10} {'|V| depois':>11}")
    for b in o["barras"]:
        print(f"{b:>6} {S.real:>10.2f} {S.imag:>10.2f} "
              f"{abs(base['V'][b]):>10.5f} {abs(res['V'][b]):>11.5f}")
    print(f"\n>>> Instalar nas barras: {o['barras']}")

    print("\n" + "=" * 60)
    print("COMPARAÇÃO ANTES x DEPOIS")
    print("=" * 60)
    red = 100 * (1 - res["perdas_kW"] / base["perdas_kW"])
    print(f"Perdas ativas          : {base['perdas_kW']:.3f} kW -> "
          f"{res['perdas_kW']:.3f} kW ({red:.1f}% de redução)")
    print(f"Perdas reativas        : {base['perdas_kvar']:.3f} kvar -> "
          f"{res['perdas_kvar']:.3f} kvar")
    print(f"Menor tensão           : "
          f"{min(abs(v) for v in base['V'].values()):.5f} pu -> "
          f"{min(abs(v) for v in res['V'].values()):.5f} pu")
    print(f"Custo das perdas antes : R$ {o['custo_perdas_antes']:,.2f}")
    print(f"Custo das perdas depois: R$ {o['custo_perdas_depois']:,.2f}")
    print(f"Custo de instalação    : R$ {o['custo_instalacao']:,.2f}")
    print(f"CUSTO TOTAL            : R$ {o['custo_total']:,.2f}")
