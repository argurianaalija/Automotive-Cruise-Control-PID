import pandas as pd
import matplotlib.pyplot as plt

# 1. Leggiamo il dataset generato dal nostro sistema in C
try:
    df = pd.read_csv('cruise_log.csv')
except FileNotFoundError:
    print("Errore: Il file 'cruise_log.csv' non è stato trovato. Assicurati di aver avviato prima il simulatore in C!")
    exit()

# 2. Creiamo una figura con due grafici impilati (condividono l'asse del tempo)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# --- GRAFICO SUPERIORE: La dinamica della Velocità ---
ax1.plot(df['Tempo_s'], df['Velocita_ms'], label='Velocità Reale', color='blue', linewidth=2.5)
ax1.plot(df['Tempo_s'], df['V_Riferimento_ms'], label='Riferimento (Rampa)', color='red', linestyle='--', linewidth=2)
ax1.set_title('Simulatore Cruise Control - Risposta del Sistema', fontsize=14, fontweight='bold')
ax1.set_ylabel('Velocità [m/s]', fontsize=12)
ax1.legend(loc='lower right')
ax1.grid(True, linestyle=':', alpha=0.7)

# --- GRAFICO INFERIORE: Lo sforzo dell'attuatore (Saturazione) ---
ax2.plot(df['Tempo_s'], df['Forza_N'], label='Azione di Controllo (Forza)', color='green', linewidth=2)
# Disegniamo la linea del limite fisico del motore
ax2.axhline(y=5000, color='darkorange', linestyle='-.', linewidth=1.5, label='Saturazione Motore (5000 N)')
ax2.fill_between(df['Tempo_s'], df['Forza_N'], 5000, where=(df['Forza_N'] >= 4999), color='orange', alpha=0.3, label='Windup Zone')

ax2.set_xlabel('Tempo [secondi]', fontsize=12)
ax2.set_ylabel('Forza Motore [N]', fontsize=12)
ax2.legend(loc='upper right')
ax2.grid(True, linestyle=':', alpha=0.7)

# 3. Rendiamo l'interfaccia elegante e la mostriamo a schermo
plt.tight_layout()
plt.show()