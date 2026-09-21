import sys

from ler_dados import ler_arquivo, imprimir_dados
from ordena import ordenar_linhas, imprimir_linhas_ordenadas
from fluxo_de_potencia import fluxo_de_potencia, imprimir_resultados


# ============================================================
# PROGRAMA PRINCIPAL
# ============================================================

if __name__ == "__main__":

    # --------------------------------------------------------
    # 1. Definir o caminho do arquivo
    # --------------------------------------------------------
 
    if len(sys.argv) > 1:

        caminho_arquivo = sys.argv[1]

    else:

        caminho_arquivo = r"C:\Users\55149\Downloads\069BusSystem.txt"

    # --------------------------------------------------------
    # 2. Ler o arquivo
    # --------------------------------------------------------

    dados = ler_arquivo(caminho_arquivo)

    # --------------------------------------------------------
    # 3. Imprimir os dados originais
    # --------------------------------------------------------

    imprimir_dados(dados)

    # --------------------------------------------------------
    # 4. Ordenar as linhas
    # --------------------------------------------------------

    Ni_ordenado, Nj_ordenado, id_ordenado = ordenar_linhas(dados)

    # --------------------------------------------------------
    # 5. Imprimir as linhas ordenadas
    # --------------------------------------------------------

    imprimir_linhas_ordenadas(
        Ni_ordenado,
        Nj_ordenado,
        id_ordenado
    )

    # --------------------------------------------------------
    # 6. Fluxo de potência (backward / forward)
    # --------------------------------------------------------

    resultado = fluxo_de_potencia(
        dados,
        Ni_ordenado,
        Nj_ordenado,
        id_ordenado,
        S_base_kVA=10000.0,   # base de potência (kVA)
        tol=1e-5              # critério de convergência (pu)
    )

    # --------------------------------------------------------
    # 7. Imprimir os resultados do fluxo de potência
    # --------------------------------------------------------

    imprimir_resultados(resultado)
    