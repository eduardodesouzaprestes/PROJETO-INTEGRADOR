import time

def simular_enchente():
    print("=== SIMULADOR DE RISCO DE INUNDAÇÃO - ROLANTE ===")
    
    # Entradas do Usuário
    chuva_mm = float(input("Precipitação total prevista (mm): "))
    duracao_h = int(input("Duração estimada da chuva (horas): "))
    nivel_atual = float(input("Nível inicial do rio (metros): "))
    solo_saturado = input("Solo já está encharcado? (s/n): ").strip().lower() == 's'

    # Parâmetros calibrados
    cota_transbordo = 4.5  # metros
    c_runoff = 0.70 if solo_saturado else 0.40  # Coeficiente de escoamento
    taxa_subida_por_mm = 0.035 * c_runoff  # Fator simplificado de elevação
    vazão_esgotamento = 0.15  # mm/h equivalente que o rio consegue escoar
    
    chuva_por_hora = chuva_mm / duracao_h
    nivel_rio = nivel_atual
    pico_atingido = nivel_atual
    
    print("\n--- INICIANDO SIMULAÇÃO HORA A HORA ---")
    
    for hora in range(1, duracao_h + 1):
        # Entrada vs Saída
        efetivo_mm = max(0.0, chuva_por_hora - vazão_esgotamento)
        variacao = efetivo_mm * taxa_subida_por_mm
        nivel_rio += variacao
        
        if nivel_rio > pico_atingido:
            pico_atingido = nivel_rio

        # Barra visual ASCII
        barra_tamanho = int(nivel_rio * 4)
        grafico_barra = "█" * barra_tamanho
        
        status = "NORMAL"
        if nivel_rio >= cota_transbordo:
            status = "ALERTA DE INUNDAÇÃO!"
        elif nivel_rio >= cota_transbordo - 1.0:
            status = "ATENÇÃO"

        print(f"Hora {hora:02d}h | Rio: {nivel_rio:.2f}m | [{grafico_barra:<20}] | Status: {status}")
        time.sleep(0.3)

    print("\n=== RESULTADO DA SIMULAÇÃO ===")
    print(f"Nível Máximo Atingido: {pico_atingido:.2f} metros")
    
    if pico_atingido >= 5.5:
        print("Impacto: Inundação severa. Áreas centrais e periféricas atingidas.")
    elif pico_atingido >= cota_transbordo:
        print("Impacto: Transbordamento do leito. Ruas baixas próximas ao rio afetadas.")
    else:
        print("Impacto: Sem risco imediato de transbordamento urbano.")

if __name__ == "__main__":
    simular_enchente()