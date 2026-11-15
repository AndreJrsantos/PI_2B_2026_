# ============================================================
# ALOCAÇÃO DE SUBESTAÇÕES E/OU CAPACITORES
# ALGORITMO IMUNOLÓGICO ARTIFICIAL (CLONALG)
# ============================================================
#
# Problema MESTRE : escolher em quais barras instalar os equipamentos.
# Problema ESCRAVO: fluxo de potência backward/forward (fluxo_de_potencia.py)
#
# EQUIPAMENTOS
#   - Subestação (SE) : injeta potência ativa e reativa (P + jQ) na barra.
#   - Capacitor  (CAP): injeta apenas potência reativa (jQ, em kvar).
#                       É o mesmo que uma carga reativa negativa.
#
# MODOS DE EXECUÇÃO
#   1 = só subestações
#   2 = só capacitores
#   3 = subestações e capacitores juntos
#
# FUNÇÃO OBJETIVO (minimizar):
#
#   custo_total = soma dos preços dos equipamentos instalados
#               + custo_por_kW * perdas_kW
#               + penalidades (fluxo não convergiu / tensão fora da faixa)
#
# COMO O ALGORITMO FUNCIONA (resumo)
#   Um "anticorpo" é uma lista de barras. A posição k da lista diz
#   em qual barra está o equipamento k.
#       Exemplo: [61, 62, 22]  ->  equip. 1 na barra 61, equip. 2 na 62...
#
#   1) O 1º equipamento é testado PONTO A PONTO (todas as barras).
#   2) Para cada barra do 1º equipamento, os outros começam ALEATÓRIOS
#      e são melhorados pelo ciclo imunológico:
#        avaliar -> selecionar melhores -> clonar -> hipermutar -> re-selecionar
#   3) No final, um .txt informa onde instalar cada equipamento.
# ============================================================

import math
import os
import random
import sys
import time
from dataclasses import dataclass

from ler_dados import ler_arquivo, resolver_caminho
from ordena import ordenar_linhas
from fluxo_de_potencia import fluxo_de_potencia




# ============================================================
# 1. ESTRUTURAS DE DADOS
# ============================================================

MODO_SO_SE = 1
MODO_SO_CAP = 2
MODO_AMBOS = 3


@dataclass
class Equipamento:
    """Um equipamento a ser instalado (subestação ou capacitor)."""
    tipo: str          # "SE" ou "CAP"
    preco: float       # preço de instalação (R$)
    P: float = 0.0     # potência ativa injetada (kW)   -> capacitor = 0
    Q: float = 0.0     # potência reativa injetada (kvar)

    @property
    def nome(self):
        return "Subestacao" if self.tipo == "SE" else "Capacitor"


@dataclass
class ParametrosImunologico:
    """Parâmetros do algoritmo (mexa aqui para ajustar a busca)."""
    pop: int = 10            # tamanho da população de anticorpos
    n_sel: int = 5           # quantos melhores são clonados por geração
    beta: float = 2.0        # fator de clonagem (mais alto = mais clones)
    rho: float = 3.0         # decaimento da mutação (mais alto = muta menos os bons)
    n_novos: int = 2         # anticorpos aleatórios inseridos por geração
    geracoes: int = 15       # máximo de gerações por ponto da varredura
    paciencia: int = 6       # para se ficar N gerações sem melhorar
    melhoria_min: float = 0.02   # só adiciona mais um equipamento se o custo total
                                 # cair pelo menos 2% (0.02) em relação ao anterior


# ============================================================
# 2. AVALIAÇÃO (AFINIDADE) DE UMA CONFIGURAÇÃO
# ============================================================

class Avaliador:
    """Roda o fluxo de potência e calcula o custo total (com cache)."""

    def __init__(self, dados, equipamentos, custo_perda,
                 v_min=0.95, v_max=1.05, penal_tensao=1e7, S_base_kVA=10000.0,
                 despachar=True, penal_refluxo=1e4):
        self.dados = dados
        self.equipamentos = equipamentos
        self.custo_perda = custo_perda          # R$ por kW dissipado
        self.v_min, self.v_max = v_min, v_max
        self.penal_tensao = penal_tensao        # R$ por pu de violação de tensão
        self.S_base = S_base_kVA
        self.Ni, self.Nj, self.ids = ordenar_linhas(dados)
        self.preco_total = sum(e.preco for e in equipamentos)
        # Despacho das SEs: P e Q informados são a CAPACIDADE da SE. Ela só
        # injeta o que a demanda do sistema comporta (sem refluxo para a SE
        # principal). Com despachar=False, injeta sempre P e Q cheios.
        self.despachar = despachar
        self.penal_refluxo = penal_refluxo      # R$ por kW devolvido à SE principal
        self.carga_kW = sum(b["SaReal"] for b in dados["barras"])
        self.cache = {}                         # evita recalcular a mesma configuração
        self.n_fluxos = 0

    def fluxo_sem_equipamentos(self):
        """Caso base: rede original, nada instalado."""
        return fluxo_de_potencia(self.dados, self.Ni, self.Nj, self.ids,
                                 S_base_kVA=self.S_base)

    def _violacao(self, res):
        return sum(max(0.0, self.v_min - abs(v)) + max(0.0, abs(v) - self.v_max)
                   for v in res["V"].values())

    def custo_sem_equipamentos(self):
        """Custo total do sistema original (perdas + penalidade de tensão)."""
        res = self.fluxo_sem_equipamentos()
        viol = self._violacao(res)
        return {"custo": self.custo_perda * res["perdas_kW"] + self.penal_tensao * viol,
                "violacao": viol, "res": res}

    def avaliar(self, posicoes):
        """posicoes[k] = barra do equipamento k. Retorna dict com custo e detalhes."""
        chave = tuple(posicoes)
        if chave in self.cache:
            return self.cache[chave]

        # Soma as injeções de cada barra (pode haver SE e CAP na mesma barra)
        cap_se = sum(e.P for e in self.equipamentos if e.tipo == "SE")
        alfa = min(1.0, self.carga_kW / cap_se) if (self.despachar and cap_se > 0) else 1.0

        injecoes = {}
        for equip, barra in zip(self.equipamentos, posicoes):
            f = alfa if equip.tipo == "SE" else 1.0     # capacitor injeta sempre o nominal
            injecoes[barra] = injecoes.get(barra, 0j) + complex(equip.P, equip.Q) * f

        res = fluxo_de_potencia(self.dados, self.Ni, self.Nj, self.ids,
                                S_base_kVA=self.S_base, injecoes=injecoes)
        self.n_fluxos += 1

        if not res["convergiu"]:
            resultado = {"custo": 1e15, "perdas_kW": float("nan"), "alfa": alfa,
                         "violacao": float("nan"), "convergiu": False, "res": res}
        else:
            violacao = self._violacao(res)
            refluxo = max(0.0, -res["S_SE_kW"])          # SE principal recebendo potência
            custo = (self.preco_total
                     + self.custo_perda * res["perdas_kW"]
                     + self.penal_tensao * violacao
                     + self.penal_refluxo * refluxo)
            resultado = {"custo": custo, "perdas_kW": res["perdas_kW"], "alfa": alfa,
                         "violacao": violacao, "convergiu": True, "res": res}

        self.cache[chave] = resultado
        return resultado


# ============================================================
# 3. OPERADORES IMUNOLÓGICOS
# ============================================================
# Regra de posicionamento: dois equipamentos do MESMO tipo não
# ficam na mesma barra. Uma SE e um capacitor podem dividir a barra.

def barras_proibidas(equipamentos, anticorpo, k):
    """Barras já ocupadas por outros equipamentos do mesmo tipo que o k."""
    return {anticorpo[j] for j in range(len(anticorpo))
            if j != k and anticorpo[j] is not None
            and equipamentos[j].tipo == equipamentos[k].tipo}


def anticorpo_aleatorio(equipamentos, fixa, candidatas, rng):
    """Cria um anticorpo: o 1º equipamento fica em 'fixa', os outros são sorteados."""
    anticorpo = [None] * len(equipamentos)
    anticorpo[0] = fixa
    for k in range(1, len(equipamentos)):
        opcoes = [b for b in candidatas
                  if b not in barras_proibidas(equipamentos, anticorpo, k)]
        anticorpo[k] = rng.choice(opcoes)
    return anticorpo


def hipermutar(equipamentos, anticorpo, candidatas, taxa, rng):
    """
    Hipermutação: cada equipamento (menos o 1º, que é fixo) troca de barra
    com probabilidade 'taxa'. Anticorpos ruins mutam mais. Sempre muta pelo menos 1.
    """
    novo = list(anticorpo)
    moveis = list(range(1, len(novo)))
    a_mutar = [k for k in moveis if rng.random() < taxa] or [rng.choice(moveis)]

    for k in a_mutar:
        proibidas = barras_proibidas(equipamentos, novo, k) | {novo[k]}
        opcoes = [b for b in candidatas if b not in proibidas]
        if opcoes:
            novo[k] = rng.choice(opcoes)
    return novo


def ciclo_imunologico(av, fixa, candidatas, rng, p):
    """CLONALG: otimiza os equipamentos 2..N com o 1º fixado na barra 'fixa'."""
    equipamentos = av.equipamentos

    # Só um equipamento: não há o que otimizar, só avaliar
    if len(equipamentos) == 1:
        return [fixa], av.avaliar([fixa])

    # População inicial ALEATÓRIA
    populacao = [anticorpo_aleatorio(equipamentos, fixa, candidatas, rng)
                 for _ in range(p.pop)]

    melhor_ab, melhor_r, geracoes_sem_melhora = None, None, 0

    for _ in range(p.geracoes):

        # --- 1. Avaliar e ordenar (menor custo = melhor afinidade) ---
        avaliados = sorted(((ab, av.avaliar(ab)) for ab in populacao),
                           key=lambda t: t[1]["custo"])

        # --- 2. Guardar o melhor / critério de parada ---
        if melhor_r is None or avaliados[0][1]["custo"] < melhor_r["custo"] - 1e-9:
            melhor_ab, melhor_r = avaliados[0]
            geracoes_sem_melhora = 0
        else:
            geracoes_sem_melhora += 1
            if geracoes_sem_melhora >= p.paciencia:
                break

        # Afinidade normalizada em [0,1] (usada para definir a taxa de mutação)
        afin = [1.0 / r["custo"] for _, r in avaliados]
        a_max, a_min = max(afin), min(afin)
        faixa = (a_max - a_min) or 1.0

        # --- 3. Selecionar os melhores, clonar e hipermutar ---
        clones = []
        for posto, (ab, _) in enumerate(avaliados[:p.n_sel], start=1):
            a_norm = (afin[posto - 1] - a_min) / faixa
            taxa = math.exp(-p.rho * a_norm)               # bom -> muta pouco
            n_clones = max(1, round(p.beta * p.pop / posto))  # bom -> mais clones
            for _ in range(n_clones):
                clones.append(hipermutar(equipamentos, ab, candidatas, taxa, rng))

        # --- 4. Re-seleção: os melhores entre população atual + clones ---
        unicos = {tuple(ab): ab for ab in populacao + clones}
        populacao = sorted(unicos.values(),
                           key=lambda ab: av.avaliar(ab)["custo"])[:p.pop]

        # --- 5. Edição de receptores: troca os piores por anticorpos novos ---
        for k in range(1, min(p.n_novos, len(populacao) - 1) + 1):
            populacao[-k] = anticorpo_aleatorio(equipamentos, fixa, candidatas, rng)

    # Avaliação final da população remanescente
    for ab in populacao:
        r = av.avaliar(ab)
        if melhor_r is None or r["custo"] < melhor_r["custo"]:
            melhor_ab, melhor_r = ab, r
    return melhor_ab, melhor_r


# ============================================================
# 4. ALGORITMO COMPLETO (varredura ponto a ponto do 1º equipamento)
# ============================================================

def ordem_de_adicao(equipamentos):
    """
    Define a ordem em que os equipamentos são "adicionados" ao testar 1, 2, 3...
    Com SE e CAP juntos, intercala (SE1, CAP1, SE2, CAP2...) para que os dois
    tipos sejam considerados antes de se pensar em um segundo do mesmo tipo.
    """
    ses = [e for e in equipamentos if e.tipo == "SE"]
    caps = [e for e in equipamentos if e.tipo == "CAP"]
    ordem = []
    for i in range(max(len(ses), len(caps))):
        if i < len(ses):
            ordem.append(ses[i])
        if i < len(caps):
            ordem.append(caps[i])
    return ordem


def _varrer(av, candidatas, rng, params, n_execucoes, verbose, rotulo):
    """Varredura ponto a ponto do 1º equipamento, repetida n_execucoes vezes."""
    melhor_por_barra = {}
    melhor_global = None
    custos_execucoes = []

    for execucao in range(1, n_execucoes + 1):
        melhor_exec = None

        for idx, barra1 in enumerate(candidatas, start=1):
            ab, r = ciclo_imunologico(av, barra1, candidatas, rng, params)

            if barra1 not in melhor_por_barra or r["custo"] < melhor_por_barra[barra1][0]:
                melhor_por_barra[barra1] = (r["custo"], ab)
            if melhor_exec is None or r["custo"] < melhor_exec[1]["custo"]:
                melhor_exec = (ab, r)
            if melhor_global is None or r["custo"] < melhor_global[1]["custo"]:
                melhor_global = (ab, r)

            if verbose >= 2:
                print(f"[{rotulo}] [exec {execucao}/{n_execucoes}] [{idx:>3}/{len(candidatas)}] "
                      f"1º equipamento na barra {barra1:>3} "
                      f"-> melhor custo {r['custo']:>14,.2f}  "
                      f"(global: {melhor_global[1]['custo']:,.2f})")

        custos_execucoes.append(melhor_exec[1]["custo"])
        if verbose:
            print(f">>> [{rotulo}] Execução {execucao} terminou: melhor custo = "
                  f"{melhor_exec[1]['custo']:,.2f}  configuração {melhor_exec[0]}")

    varredura = [(b, c, ab) for b, (c, ab) in melhor_por_barra.items()]
    return melhor_global, varredura, custos_execucoes


def otimizar(dados, equipamentos, custo_perda, seed=None, verbose=True,
             params=None, n_execucoes=1):
    """
    Otimiza a localização dos equipamentos pedidos, mas só usa a quantidade
    que realmente vale a pena:

      - testa com 1 equipamento, depois 2, depois 3... (até o máximo pedido);
      - adicionar o equipamento k só é aceito se o CUSTO TOTAL cair pelo menos
        'params.melhoria_min' (padrão 2%) em relação à melhor solução com k-1;
      - na primeira vez que a melhoria fica abaixo disso, o algoritmo para e
        mantém a solução com k-1 equipamentos.

    Ex.: pediu 4, mas de 3 para 4 o custo caiu só 0,8% -> instala apenas 3.
    """
    rng = random.Random(seed)
    params = params or ParametrosImunologico()

    ordem = ordem_de_adicao(equipamentos)

    SE = next(b["id"] for b in dados["barras"] if b["type"] == "SE")
    candidatas = sorted(b["id"] for b in dados["barras"] if b["id"] != SE)

    for tipo in ("SE", "CAP"):
        qtd = sum(1 for e in ordem if e.tipo == tipo)
        if qtd > len(candidatas):
            raise ValueError("Há mais equipamentos do mesmo tipo do que barras disponíveis.")

    # Caso base (nada instalado), para comparação
    av0 = Avaliador(dados, ordem, custo_perda)
    base = av0.fluxo_sem_equipamentos()
    modulos_base = [abs(v) for v in base["V"].values()]
    cb = av0.custo_sem_equipamentos()

    historico = []        # um item por quantidade testada
    escolhido = None      # solução da quantidade aceita
    n_fluxos = 0

    t0 = time.time()
    for k in range(1, len(ordem) + 1):
        av = Avaliador(dados, ordem[:k], custo_perda)
        # Com 1 equipamento a varredura é determinística: 1 execução basta.
        n_exec_k = 1 if k == 1 else n_execucoes
        melhor_global, varredura, custos_exec = _varrer(
            av, candidatas, rng, params, n_exec_k, verbose, f"{k} equip.")
        n_fluxos += av.n_fluxos

        custo_k = melhor_global[1]["custo"]
        sol = {"k": k, "av": av, "melhor_ab": melhor_global[0],
               "melhor_r": melhor_global[1], "varredura": varredura,
               "custos_execucoes": custos_exec}

        if escolhido is None:
            ganho, aceito = None, True
        else:
            custo_ant = escolhido["melhor_r"]["custo"]
            ganho = (custo_ant - custo_k) / custo_ant
            aceito = ganho >= params.melhoria_min

        historico.append({"k": k, "custo": custo_k, "ganho": ganho, "aceito": aceito,
                          "config": melhor_global[0]})

        if verbose:
            if ganho is None:
                print(f"\n### {k} equipamento(s): custo {custo_k:,.2f}\n")
            else:
                print(f"\n### {k} equipamento(s): custo {custo_k:,.2f} | "
                      f"melhoria vs {k-1}: {ganho*100:.2f}% "
                      f"({'ACEITO' if aceito else 'REJEITADO'}; mínimo {params.melhoria_min*100:.1f}%)\n")

        if not aceito:
            break
        escolhido = sol

    av = escolhido["av"]
    return {
        "melhor_ab": escolhido["melhor_ab"], "melhor_r": escolhido["melhor_r"],
        "varredura": escolhido["varredura"], "base": base,
        "custo_base": custo_perda * base["perdas_kW"],
        "custo_base_total": cb["custo"], "violacao_base": cb["violacao"],
        "vmin_base": min(modulos_base), "preco_total": av.preco_total,
        "n_fluxos": n_fluxos, "tempo_s": time.time() - t0,
        "custo_perda": custo_perda, "equipamentos": av.equipamentos, "SE": SE,
        "params": params, "n_execucoes": n_execucoes,
        "custos_execucoes": escolhido["custos_execucoes"],
        "n_pedidos": len(ordem), "n_usados": escolhido["k"], "historico": historico,
    }


# ============================================================
# 5. RELATÓRIO .txt
# ============================================================

def gerar_txt(out, caminho="resultado_subestacoes.txt"):
    r = out["melhor_r"]
    res = r["res"]
    base = out["base"]
    pos = out["melhor_ab"]
    cp = out["custo_perda"]
    V = res["V"]
    barra_min = min(V, key=lambda b: abs(V[b]))

    L = []
    L.append("=" * 72)
    L.append("LOCALIZACAO OTIMA DOS EQUIPAMENTOS - ALGORITMO IMUNOLOGICO")
    L.append("=" * 72)
    L.append("")
    L.append(f"QUANTIDADE: pedidos {out['n_pedidos']}, instalados {out['n_usados']} "
             f"(melhoria minima para adicionar mais um: {out['params'].melhoria_min*100:.1f}%)")
    L.append(f"{'Qtd':>4} {'Custo total':>16} {'Melhoria':>10}  Decisao")
    for h in out["historico"]:
        mel = "-" if h["ganho"] is None else f"{h['ganho']*100:.2f}%"
        L.append(f"{h['k']:>4} {h['custo']:>16,.2f} {mel:>10}  {'aceito' if h['aceito'] else 'rejeitado'}")
    L.append("")
    L.append("ONDE INSTALAR:")
    L.append(f"{'N':>3} {'Tipo':>11} {'Barra':>6} {'P (kW)':>12} {'Q (kvar)':>12} {'Preco':>16}")
    for k, (e, b) in enumerate(zip(out["equipamentos"], pos), start=1):
        L.append(f"{k:>3} {e.nome:>11} {b:>6} {e.P:>12.2f} {e.Q:>12.2f} {e.preco:>16,.2f}")
    L.append("")
    L.append("CUSTOS:")
    L.append(f"  Custo dos equipamentos       : {out['preco_total']:>16,.2f}")
    L.append(f"  Perdas (kW)                  : {r['perdas_kW']:>16.3f}")
    L.append(f"  Custo por kW perdido         : {cp:>16,.2f}")
    L.append(f"  Custo das perdas             : {cp * r['perdas_kW']:>16,.2f}")
    L.append(f"  Penalidades (tensao/refluxo) : {r['custo'] - out['preco_total'] - cp * r['perdas_kW']:>16,.2f}")
    L.append(f"  CUSTO TOTAL                  : {r['custo']:>16,.2f}")
    L.append("")
    L.append("DESPACHO DAS SUBESTACOES NOVAS (P e Q informados = capacidade):")
    L.append(f"  Fator de despacho            : {r['alfa']:>16.4f}")
    L.append(f"  Potencia da SE principal (kW): {res['S_SE_kW']:>16.3f}")
    L.append("")
    L.append("COMPARACAO COM O SISTEMA SEM NOVOS EQUIPAMENTOS:")
    L.append(f"  Perdas (kW)                  : {base['perdas_kW']:>16.3f}")
    L.append(f"  Custo das perdas             : {out['custo_base']:>16,.2f}")
    L.append(f"  Violacao de tensao (soma pu) : {out['violacao_base']:>16.5f}")
    L.append(f"  CUSTO TOTAL sem equipamentos : {out['custo_base_total']:>16,.2f}")
    L.append(f"  Menor tensao (pu)            : {out['vmin_base']:>16.5f}")
    L.append("")
    L.append("QUALIDADE DA SOLUCAO:")
    L.append(f"  Menor tensao (pu)            : {abs(V[barra_min]):>16.5f}  (barra {barra_min})")
    L.append(f"  Violacao de tensao (soma pu) : {r['violacao']:>16.5f}")
    L.append(f"  Fluxo convergiu em           : {res['iteracoes']} iteracoes")
    L.append("")

    # ---------------- Fluxo de potência da solução ----------------
    L.append("=" * 72)
    L.append("FLUXO DE POTENCIA DA SOLUCAO (com os equipamentos instalados)")
    L.append("=" * 72)
    L.append("")
    L.append("TENSOES NAS BARRAS (S = subestacao nova, C = capacitor):")
    L.append(f"{'Barra':>6} {'|V| (pu)':>10} {'ang (graus)':>12} {'|V| (kV)':>10} {'|V| base (pu)':>14}")
    for b in sorted(V):
        v = V[b]
        marca = ""
        for e, barra in zip(out["equipamentos"], pos):
            if barra == b:
                marca += " S" if e.tipo == "SE" else " C"
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
    L.append("BALANCO DE POTENCIA:")
    L.append(f"  {'':28} {'COM equipamentos':>20} {'SEM equipamentos':>20}")
    L.append(f"  {'Perdas totais':28} "
             f"{res['perdas_kW']:>10.3f} + j{res['perdas_kvar']:<8.3f} "
             f"{base['perdas_kW']:>10.3f} + j{base['perdas_kvar']:<8.3f}")
    L.append("")
    L.append("EXECUCAO:")
    L.append(f"  Tamanho da populacao         : {out['params'].pop}")
    L.append(f"  Quantidade de execucoes      : {out['n_execucoes']}")
    L.append(f"  Melhor custo de cada execucao: "
             + ", ".join(f"{c:,.2f}" for c in out["custos_execucoes"]))
    L.append(f"  Fluxos de potencia calculados: {out['n_fluxos']}")
    L.append(f"  Tempo (s)                    : {out['tempo_s']:.1f}")
    L.append("")
    L.append("VARREDURA PONTO A PONTO DO 1o EQUIPAMENTO (melhor custo por barra):")
    L.append(f"{'Barra 1o':>9} {'Custo total':>16}   Configuracao completa")
    for barra1, custo, ab in sorted(out["varredura"], key=lambda t: t[1]):
        L.append(f"{barra1:>9} {custo:>16,.2f}   {ab}")

    with open(caminho, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    return caminho


# ============================================================
# 6. ENTRADA DO USUÁRIO
# ============================================================

def _num(msg, tipo=float, minimo=None):
    """Lê um número do teclado, repetindo a pergunta até ser válido."""
    while True:
        try:
            v = tipo(input(msg).replace(",", "."))
            if minimo is not None and v < minimo:
                print(f"  Valor deve ser >= {minimo}.")
                continue
            return v
        except ValueError:
            print("  Valor inválido, tente de novo.")


def _pedir_subestacoes():
    n = _num("Quantas subestações podem ser instaladas? ", int, 1)
    lista = []
    for k in range(1, n + 1):
        print(f"\n--- Subestação {k} ---")
        preco = _num("  Preço de instalação (R$): ", float, 0)
        P = _num("  Potência real P (kW): ", float)
        Q = _num("  Potência imaginária Q (kvar): ", float)
        lista.append(Equipamento("SE", preco, P, Q))
    return lista


def _pedir_capacitores():
    n = _num("Quantos capacitores podem ser instalados? ", int, 1)
    lista = []
    for k in range(1, n + 1):
        print(f"\n--- Capacitor {k} ---")
        preco = _num("  Preço de instalação (R$): ", float, 0)
        Q = _num("  Capacitância / potência reativa do banco (kvar): ", float, 0)
        lista.append(Equipamento("CAP", preco, 0.0, Q))
    return lista


def pedir_dados_usuario():
    """Pergunta o modo de execução e os dados. Retorna (equipamentos, custo_perda)."""
    print("\nO que deseja otimizar?")
    print("  1 - Somente subestações")
    print("  2 - Somente capacitores")
    print("  3 - Subestações e capacitores juntos")
    modo = _num("Opção: ", int, 1)
    while modo not in (MODO_SO_SE, MODO_SO_CAP, MODO_AMBOS):
        print("  Escolha 1, 2 ou 3.")
        modo = _num("Opção: ", int, 1)

    equipamentos = []
    if modo in (MODO_SO_SE, MODO_AMBOS):
        print("\n=== SUBESTAÇÕES ===")
        equipamentos += _pedir_subestacoes()
    if modo in (MODO_SO_CAP, MODO_AMBOS):
        print("\n=== CAPACITORES ===")
        equipamentos += _pedir_capacitores()

    print()
    custo_perda = _num("Custo por kW de potência desperdiçada/dissipada (R$/kW): ", float, 0)
    return equipamentos, custo_perda


def pedir_parametros():
    """Pergunta o tamanho da população e a quantidade de execuções."""
    print("\n=== PARÂMETROS DO ALGORITMO ===")
    pop = _num("Tamanho da população de anticorpos (ex.: 10): ", int, 2)
    n_execucoes = _num("Quantidade de execuções do algoritmo (ex.: 3): ", int, 1)
    params = ParametrosImunologico(pop=pop, n_sel=max(1, min(5, pop // 2)))
    return params, n_execucoes


if __name__ == "__main__":
    nome = sys.argv[1] if len(sys.argv) > 1 else "069BusSystem.txt"
    dados = ler_arquivo(resolver_caminho(nome))

    equipamentos, custo_perda = pedir_dados_usuario()
    params, n_execucoes = pedir_parametros()

    saida = otimizar(dados, equipamentos, custo_perda,
                     params=params, n_execucoes=n_execucoes)
    arq = gerar_txt(saida)

    print("\nMelhor solução:", saida["melhor_ab"],
          f"| custo total = {saida['melhor_r']['custo']:,.2f}")
    print(f"Arquivo gerado: {arq}")