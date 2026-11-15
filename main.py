import os
import sys

from ler_dados import ler_arquivo, imprimir_dados, resolver_caminho
from ordena import ordenar_linhas, imprimir_linhas_ordenadas
from fluxo_de_potencia import fluxo_de_potencia, imprimir_resultados
from imunologico import pedir_dados_usuario, pedir_parametros, otimizar, gerar_txt


# Arquivo usado quando nenhum nome é informado
ARQUIVO_PADRAO = "069BusSystem.txt"

# ============================================================
# PROGRAMA PRINCIPAL
# ============================================================

if __name__ == "__main__":

    # --------------------------------------------------------
    # 1. Definir o caminho do arquivo
    # --------------------------------------------------------
 
    # Basta o NOME do arquivo (na mesma pasta do programa). Também aceita
    # caminho completo. Sem argumento, usa ARQUIVO_PADRAO.
    nome_arquivo = sys.argv[1] if len(sys.argv) > 1 else ARQUIVO_PADRAO
    caminho_arquivo = resolver_caminho(nome_arquivo)
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

    # --------------------------------------------------------
    # 8. Alocação de subestações (algoritmo imunológico)
    # --------------------------------------------------------

    resp = input("\nDeseja otimizar a localização de subestações e/ou capacitores? (s/n): ").strip().lower()

    if resp == "s":

        # usuário escolhe o modo (só SE / só capacitor / ambos) e digita os dados
        equipamentos, custo_perda = pedir_dados_usuario()

        # usuário digita o tamanho da população e a quantidade de execuções
        params, n_execucoes = pedir_parametros()

        saida = otimizar(dados, equipamentos, custo_perda,
                         params=params, n_execucoes=n_execucoes)

        arquivo = gerar_txt(saida, "resultado_subestacoes.txt")

        # fluxo de potência da solução encontrada (com as subestações instaladas)
        print("\n>>> FLUXO DE POTÊNCIA COM OS EQUIPAMENTOS INSTALADOS EM:", saida["melhor_ab"])
        imprimir_resultados(saida["melhor_r"]["res"])

        print(f"\nEquipamentos pedidos: {saida['n_pedidos']} | realmente instalados: {saida['n_usados']}")
        print("Melhor solução (barra de cada equipamento instalado):", saida["melhor_ab"])
        print(f"Custo total: {saida['melhor_r']['custo']:,.2f}")
        print(f"Custo total SEM equipamentos: {saida['custo_base_total']:,.2f}")
        if saida["melhor_r"]["custo"] >= saida["custo_base_total"]:
            print("ATENÇÃO: instalar os equipamentos NÃO compensa frente ao sistema original.")
        print(f"Arquivo gerado: {arquivo}")