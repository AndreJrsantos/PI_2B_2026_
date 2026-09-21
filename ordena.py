# ============================================================
# ORDENAÇÃO DAS LINHAS DA REDE
# ============================================================

def ordenar_linhas(dados):


    # --------------------------------------------------------
    # 1. Pegar as linhas do arquivo
    # --------------------------------------------------------

    linhas = dados["linhas"]


    # --------------------------------------------------------
    # 2. Encontrar a subestação
    # --------------------------------------------------------

    SE = None

    # Percorre todas as barras
    for barra in dados["barras"]:

        # Verifica se a barra é do tipo SE
        if barra["type"] == "SE":

            # Guarda o ID da barra da subestação
            SE = barra["id"]

            # Como encontramos a SE, podemos parar
            break


    # Se nenhuma SE foi encontrada, gera um erro
    if SE is None:
        raise ValueError("Nenhuma barra do tipo SE foi encontrada.")


    # --------------------------------------------------------
    # 3. Criar uma estrutura de conexões
    # --------------------------------------------------------

    conexoes = {}


    # Percorre todas as linhas originais
    for linha in linhas:

        Ni = linha["Ni"]

        if Ni not in conexoes:
            conexoes[Ni] = []

        conexoes[Ni].append(linha)

    # --------------------------------------------------------
    # 4. Criar as listas de resultado
    # --------------------------------------------------------

    # Lista dos Ni na nova ordem
    Ni_ordenado = []

    # Lista dos Nj na nova ordem
    Nj_ordenado = []

    # Lista dos IDs ORIGINAIS na nova ordem
    id_ordenado = []


    # --------------------------------------------------------
    # 5. Criar a fila dos nós que ainda serão analisados
    # --------------------------------------------------------

    fila = [SE]

    # --------------------------------------------------------
    # 6. Controle dos nós que já foram analisados
    # --------------------------------------------------------

    nos_visitados = set()

    # --------------------------------------------------------
    # 7. Percorrer a rede
    # --------------------------------------------------------

    while fila:

        no_atual = fila.pop(0)


        if no_atual in nos_visitados:
            continue


        nos_visitados.add(no_atual)


        # ----------------------------------------------------
        # 8. Procurar TODAS as linhas que saem do nó atual
        # ----------------------------------------------------

        if no_atual not in conexoes:

            continue

        for linha in conexoes[no_atual]:

            # Pega o Ni da linha
            Ni = linha["Ni"]

            # Pega o Nj da linha
            Nj = linha["Nj"]

            # Pega o ID ORIGINAL da linha
            id_original = linha["id"]


            # ------------------------------------------------
            # 9. Guardar a conexão na nova ordem
            # ------------------------------------------------

            # Guarda o Ni
            Ni_ordenado.append(Ni)

            # Guarda o Nj
            Nj_ordenado.append(Nj)

            id_ordenado.append(id_original)


            # ------------------------------------------------
            # 10. Colocar o Nj na fila
            # ------------------------------------------------

            fila.append(Nj)


    # --------------------------------------------------------
    # 11. Verificar se todas as linhas foram encontradas
    # --------------------------------------------------------

    if len(Ni_ordenado) != len(linhas):

        print("\nATENÇÃO:")

        print(
            f"{len(linhas) - len(Ni_ordenado)} "
            "linha(s) não foram alcançadas a partir da SE."
        )


    # --------------------------------------------------------
    # 12. Retornar os resultados
    # --------------------------------------------------------

    return Ni_ordenado, Nj_ordenado, id_ordenado


# ============================================================
# IMPRESSÃO DAS LINHAS ORDENADAS
# ============================================================

def imprimir_linhas_ordenadas(
    Ni_ordenado,
    Nj_ordenado,
    id_ordenado
):

    # Título
    print("\n" + "=" * 60)
    print("LINHAS ORDENADAS A PARTIR DA SUBESTAÇÃO")
    print("=" * 60)


    # Cabeçalho
    print(
        f"{'Ni_ord':>10} "
        f"{'Nj_ord':>10} "
        f"{'id_ord':>10}"
    )


    # Percorre simultaneamente as três listas
    for Ni, Nj, id_linha in zip(
        Ni_ordenado,
        Nj_ordenado,
        id_ordenado
    ):

        # Imprime uma linha da tabela
        print(
            f"{Ni:>10} "
            f"{Nj:>10} "
            f"{id_linha:>10}"
        )
