"""
Simulador probabilístico de risco de inundação - Rio Rolante (RS)
Mini TCC - IFRS

Método: simulação de Monte Carlo.
A cada rodada, as variáveis incertas são sorteadas dentro de uma faixa de
incerteza e o modelo calcula o nível máximo do rio. Depois de milhares de
rodadas, a probabilidade de transbordar é a fração das rodadas que passaram
da cota.

Cadeia do modelo (hora a hora):
  chuva prevista --> hietograma (chuva ao longo do tempo)
  umidade do solo + chuva antecedente + radiação --> umidade efetiva
  umidade efetiva --> Curve Number (CN)           [método SCS-CN]
  chuva acumulada + CN --> escoamento superficial [método SCS-CN]
  escoamento --> vazão do rio (reservatório linear, gera a defasagem)
  vazão --> nível do rio (curva-chave)

ATENÇÃO: todo valor marcado com [CALIBRAR] é um PALPITE inicial e precisa
ser trocado por dado real (ANA/Hidroweb, INMET, CEMADEN, SGB) ou ajustado
comparando o modelo com eventos históricos (função validar). Modelo
educativo: NÃO substitui os alertas da Defesa Civil.

Usa SOMENTE módulos que já vêm com o Python (math, random, sys).
Nada precisa ser instalado ou baixado, e tudo roda direto no terminal.

Uso:
  python simulador_enchente_rolante.py             # simula um cenário
  python simulador_enchente_rolante.py --validar   # compara com eventos reais
"""

import math
import random
import sys

# =====================================================================
# 1. PARÂMETROS (tudo que precisa de fonte ou calibração fica aqui)
# =====================================================================
PARAMETROS = {
    # --- Bacia ---
    "area_bacia_km2": 300.0,         # [CALIBRAR] área da bacia até o ponto de medição
    "armazenamento_solo_mm": 100.0,  # [CALIBRAR] quanto o solo "guarda" de água

    # --- Escoamento superficial (SCS Curve Number) ---
    "cn_solo_seco": 55,              # [CALIBRAR] CN com solo seco
    "cn_solo_umido": 90,             # [CALIBRAR] CN com solo saturado
    "incerteza_cn": 4.0,             # desvio padrão do CN sorteado

    # --- Resposta da bacia (defasagem) ---
    "k_reservatorio_h": (3.0, 8.0),  # [CALIBRAR] faixa do tempo de resposta (horas)
    "horas_pos_chuva": 36,           # horas simuladas depois que a chuva acaba

    # --- Curva-chave: nível = h0 + a * vazão^b ---
    "h0": 0.5,                       # [CALIBRAR] com dados da estação fluviométrica
    "a": 0.12,                       # [CALIBRAR]
    "b": 0.55,                       # [CALIBRAR]
    "incerteza_a": 0.10,             # incerteza relativa (10%) da curva-chave

    # --- Umidade do solo ---
    "peso_sensor_umidade": 0.6,      # [CALIBRAR] peso do sensor x chuva antecedente
    "chuva_5d_saturacao_mm": 80.0,   # [CALIBRAR] chuva de 5 dias que satura o solo
    "incerteza_umidade": 0.05,       # erro da medida de umidade (fração, 0 a 1)

    # --- Chuva e radiação ---
    "incerteza_previsao": 0.25,      # erro relativo da previsão de chuva (25%)
    "incerteza_radiacao": 0.15,      # erro relativo da radiação (15%)
}

# Cotas de referência do rio, em metros.  [CALIBRAR] buscar fonte oficial
COTAS = {
    "atencao": 3.5,
    "transbordamento": 4.5,
    "severa": 5.5,
}

# Eventos históricos para validação. Preencha com dados reais, ex.:
# {"nome": "Evento X", "chuva_mm": 150, "duracao_h": 24, "umidade_solo": 0.8,
#  "chuva_5d_mm": 90, "radiacao": 8, "nivel_inicial_m": 1.5,
#  "nivel_observado_m": 5.2},
EVENTOS_HISTORICOS = []


# =====================================================================
# 2. FUNÇÕES AUXILIARES
# =====================================================================
def limitar(x, minimo, maximo):
    return max(minimo, min(maximo, x))


def evapotranspiracao_mm_dia(radiacao_mj_m2_dia):
    """ET de referência por radiação (Makkink simplificado), em mm/dia.
    Mais radiação solar = mais água evaporada do solo."""
    return max(0.0, 0.149 * radiacao_mj_m2_dia - 0.12)


def umidade_efetiva(umidade_sensor, chuva_5d_mm, radiacao, duracao_h, ruido):
    """Combina sensor de umidade, chuva dos últimos 5 dias e radiação.
    Retorna a umidade do solo no começo da chuva (0 = seco, 1 = saturado)."""
    p = PARAMETROS
    umid_chuva = limitar(chuva_5d_mm / p["chuva_5d_saturacao_mm"], 0.0, 1.0)
    umid = (p["peso_sensor_umidade"] * umidade_sensor
            + (1 - p["peso_sensor_umidade"]) * umid_chuva)
    # a radiação seca o solo durante o evento
    perda_mm = evapotranspiracao_mm_dia(radiacao) / 24.0 * duracao_h
    umid -= perda_mm / p["armazenamento_solo_mm"]
    return limitar(umid + ruido, 0.0, 1.0)


def hietograma(chuva_total_mm, duracao_h):
    """Distribui a chuva total ao longo das horas com um pico em posição
    aleatória (chuva real não cai uniforme)."""
    pico = random.uniform(0, duracao_h)
    sigma = max(0.5, duracao_h * random.uniform(0.15, 0.40))
    pesos = [math.exp(-0.5 * (((h + 0.5) - pico) / sigma) ** 2)
             for h in range(duracao_h)]
    soma = sum(pesos)
    return [chuva_total_mm * w / soma for w in pesos]


def escoamento_scs(chuva_acumulada_mm, cn):
    """Escoamento superficial acumulado (mm) pelo método SCS-CN."""
    s = 25400.0 / cn - 254.0   # retenção potencial
    ia = 0.2 * s               # perdas iniciais (infiltração, interceptação)
    if chuva_acumulada_mm <= ia:
        return 0.0
    return (chuva_acumulada_mm - ia) ** 2 / (chuva_acumulada_mm - ia + s)


def vazao_de_nivel(nivel_m, a):
    p = PARAMETROS
    if nivel_m <= p["h0"]:
        return 0.0
    return ((nivel_m - p["h0"]) / a) ** (1.0 / p["b"])


def nivel_de_vazao(vazao_m3s, a):
    p = PARAMETROS
    return p["h0"] + a * vazao_m3s ** p["b"]


# =====================================================================
# 3. UMA RODADA DO MODELO
# =====================================================================
def uma_rodada(cenario):
    """Roda o modelo uma vez com valores sorteados.
    Retorna (nível máximo em metros, hora em que o pico ocorreu)."""
    p = PARAMETROS
    duracao = cenario["duracao_h"]

    # --- sorteio das variáveis incertas ---
    s = p["incerteza_previsao"]
    fator_previsao = random.lognormvariate(-0.5 * s * s, s)  # média 1
    chuva_total = cenario["chuva_mm"] * fator_previsao

    radiacao = max(0.0, cenario["radiacao"]
                   * random.gauss(1.0, p["incerteza_radiacao"]))

    umid = umidade_efetiva(
        cenario["umidade_solo"], cenario["chuva_5d_mm"], radiacao, duracao,
        random.gauss(0.0, p["incerteza_umidade"]))

    cn_base = p["cn_solo_seco"] + (p["cn_solo_umido"] - p["cn_solo_seco"]) * umid
    cn = limitar(random.gauss(cn_base, p["incerteza_cn"]), 30, 98)

    k = random.uniform(*p["k_reservatorio_h"])
    a = p["a"] * random.gauss(1.0, p["incerteza_a"])
    a = max(a, 0.01)

    # --- estado inicial do rio ---
    vazao_base = vazao_de_nivel(cenario["nivel_inicial_m"], a)
    vazao = vazao_base
    alfa = 1.0 - math.exp(-1.0 / k)

    chuva_h = hietograma(chuva_total, duracao) + [0.0] * p["horas_pos_chuva"]

    # --- simulação hora a hora ---
    acumulada = 0.0
    escoamento_anterior = 0.0
    nivel_max = cenario["nivel_inicial_m"]
    hora_pico = 0

    for hora, chuva in enumerate(chuva_h, start=1):
        acumulada += chuva
        escoamento_total = escoamento_scs(acumulada, cn)
        escoamento_hora = escoamento_total - escoamento_anterior  # mm nesta hora
        escoamento_anterior = escoamento_total

        # mm/h sobre a bacia (km²) -> m³/s
        entrada_m3s = escoamento_hora * p["area_bacia_km2"] / 3.6

        # reservatório linear: a vazão "persegue" a entrada com atraso k
        vazao += alfa * ((vazao_base + entrada_m3s) - vazao)

        nivel = nivel_de_vazao(vazao, a)
        if nivel > nivel_max:
            nivel_max = nivel
            hora_pico = hora

    return nivel_max, hora_pico


# =====================================================================
# 4. MONTE CARLO
# =====================================================================
def percentil(ordenada, fracao):
    return ordenada[int(fracao * (len(ordenada) - 1))]


def monte_carlo(cenario, n_rodadas=10000, semente=None):
    if semente is not None:
        random.seed(semente)
    resultados = [uma_rodada(cenario) for _ in range(n_rodadas)]
    picos = sorted(r[0] for r in resultados)
    horas = [r[1] for r in resultados]

    return {
        "n": n_rodadas,
        "picos": picos,
        "prob_atencao": sum(x >= COTAS["atencao"] for x in picos) / n_rodadas,
        "prob_transbordo": sum(x >= COTAS["transbordamento"] for x in picos) / n_rodadas,
        "prob_severa": sum(x >= COTAS["severa"] for x in picos) / n_rodadas,
        "p5": percentil(picos, 0.05),
        "p50": percentil(picos, 0.50),
        "p95": percentil(picos, 0.95),
        "media": sum(picos) / n_rodadas,
        "hora_pico_media": sum(horas) / n_rodadas,
    }


# =====================================================================
# 5. SAÍDAS (texto e histograma no terminal)
# =====================================================================
def histograma_ascii(picos, n_faixas=12, largura=40):
    minimo, maximo = picos[0], picos[-1]
    if maximo - minimo < 1e-9:
        print(f"  Todos os resultados ficaram em {minimo:.2f} m")
        return
    passo = (maximo - minimo) / n_faixas
    contagem = [0] * n_faixas
    for x in picos:
        i = min(int((x - minimo) / passo), n_faixas - 1)
        contagem[i] += 1
    maior = max(contagem)
    for i, c in enumerate(contagem):
        ini, fim = minimo + i * passo, minimo + (i + 1) * passo
        barra = "█" * int(c / maior * largura)
        marca = " <- transbordo" if ini <= COTAS["transbordamento"] < fim else ""
        print(f"  {ini:5.2f}-{fim:5.2f} m | {barra} {c / len(picos):.1%}{marca}")


def mostrar_resultado(res):
    print("\n=== RESULTADO (Monte Carlo, "
          f"{res['n']} rodadas) ===")
    print(f"Nível máximo mais provável (mediana): {res['p50']:.2f} m")
    print(f"Faixa de 90% dos casos (P5 a P95):    "
          f"{res['p5']:.2f} m a {res['p95']:.2f} m")
    print(f"Hora média do pico após o início:     "
          f"{res['hora_pico_media']:.1f} h")
    print("\nProbabilidades:")
    print(f"  Atenção (>= {COTAS['atencao']} m):         {res['prob_atencao']:.1%}")
    print(f"  Transbordamento (>= {COTAS['transbordamento']} m): "
          f"{res['prob_transbordo']:.1%}")
    print(f"  Inundação severa (>= {COTAS['severa']} m):  {res['prob_severa']:.1%}")
    print("\nDistribuição do nível máximo:")
    histograma_ascii(res["picos"])
    print("\nAviso: modelo educativo. Siga a Defesa Civil e o CEMADEN.")


# =====================================================================
# 6. VALIDAÇÃO COM EVENTOS REAIS
# =====================================================================
def validar(n_rodadas=5000):
    if not EVENTOS_HISTORICOS:
        print("Preencha a lista EVENTOS_HISTORICOS no começo do arquivo "
              "com dados reais para validar o modelo.")
        return
    dentro = 0
    print("\n=== VALIDAÇÃO ===")
    for ev in EVENTOS_HISTORICOS:
        res = monte_carlo(ev, n_rodadas, semente=42)
        obs = ev["nivel_observado_m"]
        ok = res["p5"] <= obs <= res["p95"]
        dentro += ok
        print(f"{ev['nome']}: observado {obs:.2f} m | "
              f"modelo P5-P95 {res['p5']:.2f}-{res['p95']:.2f} m "
              f"(mediana {res['p50']:.2f}) -> {'DENTRO' if ok else 'FORA'}")
    print(f"\n{dentro} de {len(EVENTOS_HISTORICOS)} eventos dentro da faixa P5-P95.")
    print("Se muitos ficarem FORA, ajuste os parâmetros [CALIBRAR] e repita.")


# =====================================================================
# 7. ENTRADA DE DADOS
# =====================================================================
def ler_float(msg, minimo=0.0, maximo=None):
    while True:
        try:
            valor = float(input(msg).replace(",", "."))
        except ValueError:
            print("  Valor inválido. Digite um número.")
            continue
        if valor < minimo or (maximo is not None and valor > maximo):
            faixa = f"entre {minimo} e {maximo}" if maximo is not None \
                else f"maior ou igual a {minimo}"
            print(f"  Digite um valor {faixa}.")
            continue
        return valor


def ler_int(msg, minimo=1):
    while True:
        try:
            valor = int(input(msg))
        except ValueError:
            print("  Valor inválido. Digite um número inteiro.")
            continue
        if valor < minimo:
            print(f"  Digite um valor maior ou igual a {minimo}.")
            continue
        return valor


def coletar_cenario():
    print("=== SIMULADOR PROBABILÍSTICO DE INUNDAÇÃO - ROLANTE ===\n")
    return {
        "chuva_mm": ler_float("Chuva total prevista (mm): "),
        "duracao_h": ler_int("Duração estimada da chuva (horas): "),
        "chuva_5d_mm": ler_float("Chuva acumulada nos últimos 5 dias (mm): "),
        "umidade_solo": ler_float(
            "Umidade do solo medida na estação (0 a 100%): ", 0, 100) / 100.0,
        "radiacao": ler_float(
            "Radiação solar média da estação (MJ/m²/dia, ~5 a 25): "),
        "nivel_inicial_m": ler_float("Nível atual do rio (metros): "),
    }


# =====================================================================
# 8. PROGRAMA PRINCIPAL
# =====================================================================
def main():
    if "--validar" in sys.argv:
        validar()
        return

    cenario = coletar_cenario()
    n = ler_int("Número de rodadas Monte Carlo (ex.: 10000): ", minimo=100)
    print("\nSimulando...")
    res = monte_carlo(cenario, n)
    mostrar_resultado(res)


if __name__ == "__main__":
    main()