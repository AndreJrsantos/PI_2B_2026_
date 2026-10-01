# ============================================================
# ALOCAÇÃO DE SUBESTAÇÕES - SISTEMA IMUNOLÓGICO ARTIFICIAL (CLONALG)
# ============================================================
#
# Problema MESTRE : escolher em quais barras instalar as N subestações.
# Problema ESCRAVO: fluxo de potência backward/forward (fluxo_de_potencia.py)
#
# Função objetivo (minimizar):
#
#   custo_total = SOMA(preço das subestações)
#               + custo_por_kW * perdas_kW        (potência dissipada)
#               + penalidades (não convergência / tensão fora da faixa)
#
# Cada subestação é modelada como uma INJEÇÃO de potência (P + jQ)
# na barra escolhida (carga negativa). A SE original continua sendo
# a barra de referência.
#
# Estratégia pedida:
#   - a 1ª subestação é testada PONTO A PONTO (todas as barras);
#   - para cada ponto fixo da 1ª, as demais começam ALEATÓRIAS e são
#     refinadas pelo ciclo imunológico (clonagem + hipermutação + seleção);
#   - no final, um .txt informa onde instalar cada subestação.
# ============================================================

import math
import random
import sys
import time

from ler_dados import ler_arquivo
from ordena import ordenar_linhas
from fluxo_de_potencia import fluxo_de_potencia


# ------------------------------------------------------------
# Avaliação (afinidade) de uma configuração
# ------------------------------------------------------------

class Avaliador:
    """Roda o fluxo de potência e calcula o custo total (com cache)."""

    def __init__(self, dados, subestacoes, custo_perda,
                 v_min=0.95, v_max=1.05, penal_tensao=1e7, S_base_kVA=10000.0):
        self.dados = dados
        self.subs = subestacoes            # lista de dicts: preco, P, Q
        self.custo_perda = custo_perda     # R$ por kW dissipado
        self.v_min, self.v_max = v_min, v_max
        self.penal_tensao = penal_tensao
        self.S_base = S_base_kVA
        self.Ni, self.Nj, self.ids = ordenar_linhas(dados)
        self.preco_total = sum(s["preco"] for s in subestacoes)
        self.cache = {}
        self.n_fluxos = 0

    def avaliar(self, posicoes):
        """posicoes[k] = barra onde está a subestação k. Retorna dict de resultados."""
        chave = tuple(posicoes)
        if chave in self.cache:
            return self.cache[chave]

        inj = {}
        for s, barra in zip(self.subs, posicoes):
            inj[barra] = inj.get(barra, 0j) + complex(s["P"], s["Q"])

        res = fluxo_de_potencia(self.dados, self.Ni, self.Nj, self.ids,
                                S_base_kVA=self.S_base, injecoes=inj)
        self.n_fluxos += 1

        if not res["convergiu"]:
            out = {"custo": 1e15, "perdas_kW": float("nan"), "violacao": float("nan"),
                   "convergiu": False, "res": res}
        else:
            mod = [abs(v) for v in res["V"].values()]
            viol = sum(max(0.0, self.v_min - m) + max(0.0, m - self.v_max) for m in mod)
            custo = (self.preco_total
                     + self.custo_perda * res["perdas_kW"]
                     + self.penal_tensao * viol)
            out = {"custo": custo, "perdas_kW": res["perdas_kW"], "violacao": viol,
                   "convergiu": True, "res": res}

        self.cache[chave] = out
        return out


# ------------------------------------------------------------
# Operadores imunológicos
# ------------------------------------------------------------

def anticorpo_aleatorio(fixa, candidatas, n_sub, rng):
    """1ª subestação = 'fixa'; as demais sorteadas, sem repetir barra."""
    outras = [b for b in candidatas if b != fixa]
    return [fixa] + rng.sample(outras, n_sub - 1)


def hipermutar(ab, fixa, candidatas, taxa, rng):
    """
    Hipermutação: cada gene (menos o 1º, que é fixo) troca de barra com
    probabilidade 'taxa' (maior para anticorpos de baixa afinidade).
    Garante pelo menos 1 mutação.
    """
    novo = list(ab)
    livres = list(range(1, len(novo)))
    mutar = [k for k in livres if rng.random() < taxa] or [rng.choice(livres)]
    for k in mutar:
        usadas = set(novo)
        opcoes = [b for b in candidatas if b not in usadas]
        if opcoes:
            novo[k] = rng.choice(opcoes)
    return novo


def ciclo_imunologico(av, fixa, candidatas, n_sub, rng,
                      pop=10, n_sel=5, beta=2.0, rho=3.0,
                      n_novos=2, geracoes=15, paciencia=6):
    """CLONALG para as subestações 2..N, com a 1ª fixada em 'fixa'."""

    # N == 1: só existe um ponto a avaliar
    if n_sub == 1:
        r = av.avaliar([fixa])
        return [fixa], r

    # População inicial ALEATÓRIA
    populacao = [anticorpo_aleatorio(fixa, candidatas, n_sub, rng) for _ in range(pop)]

    melhor_ab, melhor_r, sem_melhora = None, None, 0

    for _ in range(geracoes):

        # 1. Afinidade = 1 / custo  (menor custo -> maior afinidade)
        aval = [(ab, av.avaliar(ab)) for ab in populacao]
        aval.sort(key=lambda t: t[1]["custo"])

        if melhor_r is None or aval[0][1]["custo"] < melhor_r["custo"] - 1e-9:
            melhor_ab, melhor_r, sem_melhora = aval[0][0], aval[0][1], 0
        else:
            sem_melhora += 1
            if sem_melhora >= paciencia:
                break

        # afinidade normalizada em [0,1] para controlar a mutação
        aff = [1.0 / t[1]["custo"] for t in aval]
        a_max, a_min = max(aff), min(aff)
        faixa = (a_max - a_min) or 1.0

        # 2. Seleção dos n_sel melhores, 3. clonagem e 4. hipermutação
        clones = []
        for rank, (ab, r) in enumerate(aval[:n_sel], start=1):
            a_norm = (aff[rank - 1] - a_min) / faixa
            taxa = math.exp(-rho * a_norm)              # melhor afinidade -> menos mutação
            qtd = max(1, round(beta * pop / rank))      # melhor afinidade -> mais clones
            for _ in range(qtd):
                clones.append(hipermutar(ab, fixa, candidatas, taxa, rng))

        # 5. Re-seleção: melhores entre população atual + clones
        todos = {tuple(ab): ab for ab in populacao + clones}
        ordenados = sorted(todos.values(), key=lambda ab: av.avaliar(ab)["custo"])
        populacao = ordenados[:pop]

        # 6. Edição de receptores: troca os piores por anticorpos novos (aleatórios)
        for k in range(1, n_novos + 1):
            populacao[-k] = anticorpo_aleatorio(fixa, candidatas, n_sub, rng)

    # avaliação final da população remanescente
    for ab in populacao:
        r = av.avaliar(ab)
        if melhor_r is None or r["custo"] < melhor_r["custo"]:
            melhor_ab, melhor_r = ab, r
    return melhor_ab, melhor_r


# ------------------------------------------------------------
# Algoritmo completo: varredura ponto a ponto da 1ª subestação
# ------------------------------------------------------------

def otimizar(dados, subestacoes, custo_perda, seed=None, verbose=True, **params):
    rng = random.Random(seed)
    av = Avaliador(dados, subestacoes, custo_perda)
    n_sub = len(subestacoes)

    SE = next(b["id"] for b in dados["barras"] if b["type"] == "SE")
    candidatas = sorted(b["id"] for b in dados["barras"] if b["id"] != SE)

    if n_sub > len(candidatas):
        raise ValueError("Há mais subestações do que barras disponíveis.")

    # caso base (sem subestações) para comparação
    base = fluxo_de_potencia(dados, av.Ni, av.Nj, av.ids)
    mod_base = [abs(v) for v in base["V"].values()]
    custo_base = custo_perda * base["perdas_kW"]

    varredura = []      # (barra da 1ª, melhor custo, melhor configuração)
    melhor_global = None

    t0 = time.time()
    for idx, barra1 in enumerate(candidatas, start=1):
        ab, r = ciclo_imunologico(av, barra1, candidatas, n_sub, rng, **params)
        varredura.append((barra1, r["custo"], ab))
        if melhor_global is None or r["custo"] < melhor_global[1]["custo"]:
            melhor_global = (ab, r)
        if verbose:
            print(f"[{idx:>3}/{len(candidatas)}] 1ª subestação na barra {barra1:>3} "
                  f"-> melhor custo {r['custo']:>14,.2f}  "
                  f"(global: {melhor_global[1]['custo']:,.2f})")

    return {
        "melhor_ab": melhor_global[0], "melhor_r": melhor_global[1],
        "varredura": varredura, "base": base, "custo_base": custo_base,
        "vmin_base": min(mod_base), "preco_total": av.preco_total,
        "n_fluxos": av.n_fluxos, "tempo_s": time.time() - t0,
        "custo_perda": custo_perda, "subs": subestacoes, "SE": SE,
    }


# ------------------------------------------------------------
# Relatório .txt
# ------------------------------------------------------------

def gerar_txt(out, caminho="resultado_subestacoes.txt"):
    r = out["melhor_r"]
    res = r["res"]
    pos = out["melhor_ab"]
    cp = out["custo_perda"]
    V = res["V"]
    barra_min = min(V, key=lambda b: abs(V[b]))

    L = []
    L.append("=" * 64)
    L.append("LOCALIZACAO OTIMA DAS SUBESTACOES - ALGORITMO IMUNOLOGICO")
    L.append("=" * 64)
    L.append("")
    L.append("ONDE INSTALAR:")
    L.append(f"{'Subest.':>8} {'Barra':>6} {'P (kW)':>12} {'Q (kvar)':>12} {'Preco':>16}")
    for k, (s, b) in enumerate(zip(out["subs"], pos), start=1):
        L.append(f"{k:>8} {b:>6} {s['P']:>12.2f} {s['Q']:>12.2f} {s['preco']:>16,.2f}")
    L.append("")
    L.append("CUSTOS:")
    L.append(f"  Custo das subestacoes        : {out['preco_total']:>16,.2f}")
    L.append(f"  Perdas (kW)                  : {r['perdas_kW']:>16.3f}")
    L.append(f"  Custo por kW perdido         : {cp:>16,.2f}")
    L.append(f"  Custo das perdas             : {cp * r['perdas_kW']:>16,.2f}")
    L.append(f"  CUSTO TOTAL                  : {r['custo']:>16,.2f}")
    L.append("")
    L.append("COMPARACAO COM O SISTEMA SEM NOVAS SUBESTACOES:")
    L.append(f"  Perdas (kW)                  : {out['base']['perdas_kW']:>16.3f}")
    L.append(f"  Custo das perdas             : {out['custo_base']:>16,.2f}")
    L.append(f"  Menor tensao (pu)            : {out['vmin_base']:>16.5f}")
    L.append("")
    L.append("QUALIDADE DA SOLUCAO:")
    L.append(f"  Menor tensao (pu)            : {abs(V[barra_min]):>16.5f}  (barra {barra_min})")
    L.append(f"  Violacao de tensao (soma pu) : {r['violacao']:>16.5f}")
    L.append(f"  Fluxo convergiu em           : {res['iteracoes']} iteracoes")
    L.append("")
    # ---------------- Fluxo de potência da solução ----------------
    base = out["base"]
    L.append("=" * 64)
    L.append("FLUXO DE POTENCIA DA SOLUCAO (com as subestacoes instaladas)")
    L.append("=" * 64)
    L.append("")
    L.append("TENSOES NAS BARRAS (S = barra com subestacao nova):")
    L.append(f"{'Barra':>6} {'|V| (pu)':>10} {'ang (graus)':>12} {'|V| (kV)':>10} {'|V| base (pu)':>14}")
    for b in sorted(V):
        v = V[b]
        marca = " S" if b in pos else ""
        L.append(f"{b:>6} {abs(v):>10.5f} {math.degrees(math.atan2(v.imag, v.real)):>12.4f} "
                 f"{abs(v) * res['kV_base']:>10.4f} {abs(base['V'][b]):>14.5f}{marca}")
    L.append("")
    L.append("CORRENTES NAS LINHAS (modulo):")
    L.append(f"{'id':>4} {'Ni':>4} {'Nj':>4} {'|J| (pu)':>12} {'|J| (A)':>10} {'|J| base (A)':>13}")
    base_J = {base["id"][k]: abs(base["J"][k]) * base["I_base_A"] for k in range(len(base["J"]))}
    for k in range(len(res["J"])):
        j_pu = abs(res["J"][k])
        L.append(f"{res['id'][k]:>4} {res['Ni'][k]:>4} {res['Nj'][k]:>4} {j_pu:>12.6f} "
                 f"{j_pu * res['I_base_A']:>10.3f} {base_J[res['id'][k]]:>13.3f}")
    L.append("")
    P_inj = sum(s["P"] for s in out["subs"])
    Q_inj = sum(s["Q"] for s in out["subs"])
    L.append("BALANCO DE POTENCIA:")
    L.append(f"  {'':28} {'COM subestacoes':>20} {'SEM subestacoes':>20}")
    L.append(f"  {'Perdas totais':28} "
             f"{res['perdas_kW']:>10.3f} + j{res['perdas_kvar']:<8.3f} "
             f"{base['perdas_kW']:>10.3f} + j{base['perdas_kvar']:<8.3f}")
    L.append("")
    L.append("EXECUCAO:")
    L.append(f"  Fluxos de potencia calculados: {out['n_fluxos']}")
    L.append(f"  Tempo (s)                    : {out['tempo_s']:.1f}")
    L.append("")
    L.append("VARREDURA PONTO A PONTO DA 1a SUBESTACAO (melhor custo por barra):")
    L.append(f"{'Barra 1a':>9} {'Custo total':>16}   Configuracao completa")
    for barra1, custo, ab in sorted(out["varredura"], key=lambda t: t[1]):
        L.append(f"{barra1:>9} {custo:>16,.2f}   {ab}")

    with open(caminho, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    return caminho


# ------------------------------------------------------------
# Entrada do usuário
# ------------------------------------------------------------

def _num(msg, tipo=float, minimo=None):
    while True:
        try:
            v = tipo(input(msg).replace(",", "."))
            if minimo is not None and v < minimo:
                print(f"  Valor deve ser >= {minimo}.")
                continue
            return v
        except ValueError:
            print("  Valor inválido, tente de novo.")


def pedir_dados_usuario():
    n = _num("Quantas subestações podem ser instaladas? ", int, 1)
    subs = []
    for k in range(1, n + 1):
        print(f"\n--- Subestação {k} ---")
        preco = _num("  Preço (R$): ", float, 0)
        P = _num("  Potência real P (kW): ", float)
        Q = _num("  Potência imaginária Q (kvar): ", float)
        subs.append({"preco": preco, "P": P, "Q": Q})
    print()
    custo_perda = _num("Custo por kW de potência desperdiçada/dissipada (R$/kW): ", float, 0)
    return subs, custo_perda


if __name__ == "__main__":
    caminho = sys.argv[1] if len(sys.argv) > 1 else r"C:\Users\55149\Downloads\069BusSystem.txt"
    dados = ler_arquivo(caminho)

    subs, custo_perda = pedir_dados_usuario()

    saida = otimizar(dados, subs, custo_perda, seed=None)
    arq = gerar_txt(saida)

    print("\nMelhor solução:", saida["melhor_ab"], f"| custo total = {saida['melhor_r']['custo']:,.2f}")
    print(f"Arquivo gerado: {arq}")