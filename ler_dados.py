import sys


def ler_arquivo(caminho):
    """Lê o arquivo e retorna um dicionário com todas as informações extraídas."""

    with open(caminho, "r", encoding="utf-8") as f:
        linhas = [linha.strip() for linha in f.readlines()]

    dados = {
        "number_bars": None,
        "number_lines": None,
        "voltage_kv": None,
        "barras": [],
        "linhas": [],
    }

    i = 0
    total_linhas = len(linhas)

    while i < total_linhas:
        linha = linhas[i]

        if linha == "NumberBars":
            i += 1
            dados["number_bars"] = int(linhas[i])

        elif linha == "NumberLines":
            i += 1
            dados["number_lines"] = int(linhas[i])

        elif linha.startswith("Voltage"):
            i += 1
            dados["voltage_kv"] = float(linhas[i])

        elif linha.startswith("Bars_data"):
            # avança até encontrar a linha de cabeçalho com "id"
            i += 1
            while i < total_linhas and "id" not in linhas[i]:
                i += 1
            i += 1  # pula o cabeçalho "id type SaReal SaImag Ca"

            # lê as linhas de dados das barras
            while i < total_linhas and linhas[i] != "" and not linhas[i].startswith("Lines_Data"):
                partes = linhas[i].split()
                if len(partes) >= 5:
                    barra = {
                        "id": int(partes[0]),
                        "type": partes[1],
                        "SaReal": float(partes[2]),
                        "SaImag": float(partes[3]),
                        "Ca": float(partes[4]),
                    }
                    dados["barras"].append(barra)
                i += 1
            continue

        elif linha.startswith("Lines_Data"):
            i += 1
            # avança até encontrar o cabeçalho com "id"
            while i < total_linhas and "Ni" not in linhas[i]:
                i += 1
            i += 1  # pula o cabeçalho "id Ni Nj ZaaReal ZaaImag"

            while i < total_linhas and linhas[i] != "":
                partes = linhas[i].split()
                if len(partes) >= 5:
                    linha_dados = {
                        "id": int(partes[0]),
                        "Ni": int(partes[1]),
                        "Nj": int(partes[2]),
                        "ZaaReal": float(partes[3]),
                        "ZaaImag": float(partes[4]),
                    }
                    dados["linhas"].append(linha_dados)
                i += 1
            continue

        i += 1

    return dados


def imprimir_dados(dados):
    """Imprime os dados lidos de forma organizada no terminal."""

    print("=" * 60)
    print("DADOS GERAIS DO SISTEMA")
    print("=" * 60)
    print(f"Número de barras : {dados['number_bars']}")
    print(f"Número de linhas : {dados['number_lines']}")
    print(f"Tensão (kV)      : {dados['voltage_kv']}")

    print("\n" + "=" * 60)
    print(f"DADOS DAS BARRAS ({len(dados['barras'])} registros)")
    print("=" * 60)
    print(f"{'id':>4} {'tipo':>6} {'SaReal':>10} {'SaImag':>10} {'Ca':>6}")
    for b in dados["barras"]:
        print(f"{b['id']:>4} {b['type']:>6} {b['SaReal']:>10.2f} {b['SaImag']:>10.2f} {b['Ca']:>6.0f}")

    print("\n" + "=" * 60)
    print(f"DADOS DAS LINHAS ({len(dados['linhas'])} registros)")
    print("=" * 60)
    print(f"{'id':>4} {'Ni':>4} {'Nj':>4} {'ZaaReal':>10} {'ZaaImag':>10}")
    for l in dados["linhas"]:
        print(f"{l['id']:>4} {l['Ni']:>4} {l['Nj']:>4} {l['ZaaReal']:>10.4f} {l['ZaaImag']:>10.4f}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        caminho_arquivo = sys.argv[1]
    else:
      caminho_arquivo = r"C:\Users\55149\Downloads\069BusSystem.txt"

    dados = ler_arquivo(caminho_arquivo)
    imprimir_dados(dados)

