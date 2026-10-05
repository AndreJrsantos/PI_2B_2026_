import os
import sys

from ler_dados import ler_arquivo, imprimir_dados
from ordena import ordenar_linhas, imprimir_linhas_ordenadas
from fluxo_de_potencia import fluxo_de_potencia, imprimir_resultados
from otimiza_minisubs import (
    pedir_parametros,
    otimizar_minisubestacoes,
    imprimir_otimizacao,
)


# ============================================================
# PROGRAMA PRINCIPAL
# ============================================================

if __name__ == "__main__":

    S_BASE_KVA = 10000.0

    # --------------------------------------------------------
    # 1. Caminho do arquivo (argumento > pasta do script > Downloads)
    # --------------------------------------------------------

    if len(sys.argv) > 1:
        caminho_arquivo = sys.argv[1]
    else:
        local = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "069BusSystem.txt")
        caminho_arquivo = local if os.path.exists(local) else \
            r"C:\Users\55149\Downloads\069BusSystem.txt"

    # --------------------------------------------------------
    # 2. Leitura e ordenação da rede
    # --------------------------------------------------------

    dados = ler_arquivo(caminho_arquivo)
    imprimir_dados(dados)

    Ni_ordenado, Nj_ordenado, id_ordenado = ordenar_linhas(dados)
    imprimir_linhas_ordenadas(Ni_ordenado, Nj_ordenado, id_ordenado)

    # --------------------------------------------------------
    # 3. Fluxo de potência do caso base (sem minisubestações)
    # --------------------------------------------------------

    resultado = fluxo_de_potencia(
        dados, Ni_ordenado, Nj_ordenado, id_ordenado,
        S_base_kVA=S_BASE_KVA, tol=1e-5,
    )
    print("\n##### CASO BASE (SEM MINISUBESTAÇÕES) #####")
    imprimir_resultados(resultado)

    # --------------------------------------------------------
    # 4. Dados digitados pelo usuário
    # --------------------------------------------------------

    n, preco_inst, p_kw, q_kvar, preco_kw = pedir_parametros()

    # --------------------------------------------------------
    # 5. Otimização (Simulated Annealing)
    # --------------------------------------------------------

    print("\nOtimizando, aguarde...")

    otim = otimizar_minisubestacoes(
        dados, Ni_ordenado, Nj_ordenado, id_ordenado,
        n, preco_inst, p_kw, q_kvar, preco_kw,
        S_base=S_BASE_KVA,
    )

    # --------------------------------------------------------
    # 6. Resultado: ONDE instalar + fluxo de potência final
    # --------------------------------------------------------

    imprimir_otimizacao(otim)

    print("\n##### FLUXO DE POTÊNCIA COM AS MINISUBESTAÇÕES #####")
    imprimir_resultados(otim["res"])
