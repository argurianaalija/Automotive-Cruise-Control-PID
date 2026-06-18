/**
 * ============================================================
 *  CRUISE CONTROL SIMULATOR — Versione Avanzata
 *  Controllore PID con Anti-Windup, Rampa di Riferimento,
 *  Saturazione Simmetrica, Metriche di Performance e Log CSV
 * ============================================================
 *
 *  Modello del veicolo (Plant lineare del 1° ordine):
 *      m * dv/dt = F_motore - b * v
 *
 *  Controllore: PID discreto con Back-Calculation Anti-Windup
 *      u[k] = Kp*e[k] + Ki*I[k] + Kd*(e[k]-e[k-1])/dt
 *
 *  Integrazione numerica: metodo di Eulero in avanti
 *
 *  Output:
 *      - Tabella formattata su stdout
 *      - File CSV "cruise_log.csv" per analisi esterna
 *      - Report finale con metriche ISE, overshoot, tempo di assestamento
 * ============================================================
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* ──────────────────────────────────────────────────────────
 *  1. PARAMETRI DEL VEICOLO
 * ────────────────────────────────────────────────────────── */
#define MASSA_KG        1500.0   /* Massa del veicolo [kg]          */
#define ATTRITO_NS_M      50.0   /* Coefficiente di attrito viscoso  */

/* ──────────────────────────────────────────────────────────
 *  2. PARAMETRI DEL CONTROLLORE PID
 *     Sintonizzazione basata su metodo Ziegler-Nichols
 *     adattato a plant del 1° ordine con ritardo trascurabile.
 * ────────────────────────────────────────────────────────── */
#define KP              800.0   /* Guadagno proporzionale            */
#define KI               40.0   /* Guadagno integrale                */
#define KD               20.0   /* Guadagno derivativo               */

/* Anti-windup: back-calculation gain  Kaw = 1/sqrt(Ki)  (euristica) */
#define KAW               0.16  /* Gain anti-windup                  */

/* ──────────────────────────────────────────────────────────
 *  3. LIMITI FISICI DELL'ATTUATORE
 * ────────────────────────────────────────────────────────── */
#define FORZA_MAX_N    5000.0   /* Forza massima erogabile [N]       */
#define FORZA_MIN_N       0.0   /* Nessuna trazione negativa (no ABS)*/

/* ──────────────────────────────────────────────────────────
 *  4. PROFILO DI RIFERIMENTO — RAMPA + GRADINO
 *     La velocità target non cambia a gradino (urto meccanico)
 *     ma segue una rampa con slew-rate limitato.
 * ────────────────────────────────────────────────────────── */
#define V_TARGET_MS      10.0   /* Velocità di crociera [m/s]        */
#define SLEW_RATE_MS2     2.0   /* Rampa massima sul riferimento [m/s²] */

/* ──────────────────────────────────────────────────────────
 *  5. PARAMETRI DI SIMULAZIONE
 * ────────────────────────────────────────────────────────── */
#define DT               0.01   /* Passo di integrazione [s] (100 Hz) */
#define TEMPO_TOTALE    20.0    /* Durata simulazione [s]             */
#define PRINT_DECIMATION   10   /* Stampa ogni N passi (→ 10 Hz)     */

/* ──────────────────────────────────────────────────────────
 *  6. SOGLIA DI ASSESTAMENTO
 * ────────────────────────────────────────────────────────── */
#define BANDA_SETTLING   0.02   /* ±2 % di V_TARGET                  */

/* ──────────────────────────────────────────────────────────
 *  STRUTTURA: stato completo del controllore PID
 * ────────────────────────────────────────────────────────── */
typedef struct {
    double kp, ki, kd, kaw;
    double integrale;
    double errore_precedente;
    double u_sat_prec;          /* uscita saturata al passo k-1 (anti-windup) */
} PID;

/* ──────────────────────────────────────────────────────────
 *  STRUTTURA: metriche di performance
 * ────────────────────────────────────────────────────────── */
typedef struct {
    double ise;                  /* Integral Squared Error             */
    double iae;                  /* Integral Absolute Error            */
    double overshoot_pct;        /* Overshoot massimo [%]              */
    double tempo_salita;         /* Rise time: 10%→90% del target [s] */
    double tempo_assestamento;   /* Settling time [s]                  */
    double velocita_max;         /* Picco di velocità                  */
    int    salita_10_raggiunto;
    double t_10;
    int    settled;              /* Flag: assestato definitivamente     */
    int    fuori_banda;          /* Campioni consecutivi fuori banda    */
} Metriche;

/* ──────────────────────────────────────────────────────────
 *  FUNZIONI DI UTILITÀ
 * ────────────────────────────────────────────────────────── */

/** Saturazione simmetrica con ritorno del valore saturato */
static double satura(double val, double min_v, double max_v) {
    if (val > max_v) return max_v;
    if (val < min_v) return min_v;
    return val;
}

/** Inizializza il PID */
static void pid_init(PID *c,
                     double kp, double ki, double kd, double kaw) {
    c->kp  = kp;  c->ki  = ki;
    c->kd  = kd;  c->kaw = kaw;
    c->integrale          = 0.0;
    c->errore_precedente  = 0.0;
    c->u_sat_prec         = 0.0;
}

/**
 * Calcola l'uscita del PID discreto con back-calculation anti-windup.
 *
 *   Back-calculation:
 *     I[k] = I[k-1] + Ki*e[k]*dt + Kaw*(u_sat[k-1] - u_raw[k-1])*dt
 *
 *   In questo modo l'integrale "scarica" se l'attuatore è saturo.
 */
static double pid_step(PID *c, double errore, double dt) {
    /* Termine derivativo (differenza all'indietro, niente filtro per semplicità) */
    double derivata = (errore - c->errore_precedente) / dt;

    /* Uscita PID non saturata */
    double u_raw = c->kp * errore
                 + c->ki * c->integrale
                 + c->kd * derivata;

    /* Saturazione */
    double u_sat = satura(u_raw, FORZA_MIN_N, FORZA_MAX_N);

    /* Aggiornamento integrale con anti-windup (back-calculation) */
    c->integrale        += dt * (errore + c->kaw * (u_sat - c->u_sat_prec));
    c->errore_precedente = errore;
    c->u_sat_prec        = u_sat;

    return u_sat;
}

/** Aggiorna il riferimento con slew-rate limiter */
static double rampa_riferimento(double v_ref_attuale,
                                double v_target,
                                double slew, double dt) {
    double delta = v_target - v_ref_attuale;
    double passo_max = slew * dt;
    if (fabs(delta) <= passo_max) return v_target;
    return v_ref_attuale + (delta > 0.0 ? passo_max : -passo_max);
}

/** Aggiorna le metriche ad ogni passo */
static void aggiorna_metriche(Metriche *m,
                               double t,
                               double v,
                               double errore,
                               double v_target) {
    double e2 = errore * errore;
    m->ise += e2 * DT;
    m->iae += fabs(errore) * DT;

    if (v > m->velocita_max) m->velocita_max = v;

    /* Overshoot rispetto al target finale (non al riferimento rampato) */
    double os = (v - v_target) / v_target * 100.0;
    if (os > m->overshoot_pct) m->overshoot_pct = os;

    /* Rise time 10 % */
    if (!m->salita_10_raggiunto && v >= 0.10 * v_target) {
        m->salita_10_raggiunto = 1;
        m->t_10 = t;
    }
    /* Rise time 90 % */
    if (m->salita_10_raggiunto && m->tempo_salita < 0.0 &&
        v >= 0.90 * v_target) {
        m->tempo_salita = t - m->t_10;
    }

    /* Settling time: ultima volta che esce dalla banda */
    double banda = BANDA_SETTLING * v_target;
    if (fabs(v - v_target) > banda) {
        m->settled        = 0;
        m->fuori_banda   = 1;
        m->tempo_assestamento = t;   /* aggiornato continuamente */
    } else if (!m->settled) {
        m->settled = 1;
    }
}

/* ──────────────────────────────────────────────────────────
 *  MAIN
 * ────────────────────────────────────────────────────────── */
int main(void) {

    /* --- Stato del sistema --- */
    double velocita_attuale = 0.0;
    double v_riferimento    = 0.0;   /* riferimento rampato             */

    /* --- Controllore --- */
    PID ctrl;
    pid_init(&ctrl, KP, KI, KD, KAW);

    /* --- Metriche --- */
    Metriche met;
    memset(&met, 0, sizeof(met));
    met.tempo_salita        = -1.0;
    met.tempo_assestamento  =  0.0;

    /* --- File CSV --- */
    FILE *csv = fopen("cruise_log.csv", "w");
    if (!csv) {
        fprintf(stderr, "[ERRORE] Impossibile creare cruise_log.csv\n");
        return EXIT_FAILURE;
    }
    fprintf(csv, "Tempo_s,Velocita_ms,V_Riferimento_ms,"
                 "Errore_ms,Forza_N,Accelerazione_ms2\n");

    /* --- Header console --- */
    printf("\n");
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║         CRUISE CONTROL SIMULATOR — PID + Anti-Windup    ║\n");
    printf("╠══════════════════════════════════════════════════════════╣\n");
    printf("║  Target: %.1f m/s (%.1f km/h)   Durata: %.0f s   dt: %.0f ms  ║\n",
           V_TARGET_MS, V_TARGET_MS * 3.6, TEMPO_TOTALE, DT * 1000.0);
    printf("╠══════════════════════════════════════════════════════════╣\n");
    printf("║  Kp=%.0f  Ki=%.0f  Kd=%.0f  Kaw=%.2f                      ║\n",
           KP, KI, KD, KAW);
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    printf("  %-8s %-14s %-14s %-12s %-12s %-10s\n",
           "t [s]", "v [m/s]", "v_ref [m/s]",
           "errore [m/s]", "F [N]", "a [m/s²]");
    printf("  %s\n",
           "──────────────────────────────────────────────────────────────");

    /* ── CICLO PRINCIPALE ── */
    int passi_totali = (int)(TEMPO_TOTALE / DT) + 1;

    for (int k = 0; k < passi_totali; k++) {
        double t = k * DT;

        /* 1. Aggiorna il riferimento con slew-rate limiter */
        v_riferimento = rampa_riferimento(v_riferimento,
                                          V_TARGET_MS,
                                          SLEW_RATE_MS2, DT);

        /* 2. Calcola errore rispetto al riferimento rampato */
        double errore = v_riferimento - velocita_attuale;

        /* 3. PID → forza motore */
        double forza_motore = pid_step(&ctrl, errore, DT);

        /* 4. Modello fisico del veicolo (Eulero in avanti)
         *    m·dv/dt = F - b·v  ⟹  a = (F - b·v) / m           */
        double accelerazione = (forza_motore
                                - ATTRITO_NS_M * velocita_attuale)
                               / MASSA_KG;
        velocita_attuale += accelerazione * DT;

        /* 5. Metriche (errore rispetto a target finale, non rampa) */
        double errore_finale = V_TARGET_MS - velocita_attuale;
        aggiorna_metriche(&met, t, velocita_attuale,
                          errore_finale, V_TARGET_MS);

        /* 6. Log CSV (ogni passo) */
        fprintf(csv, "%.3f,%.4f,%.4f,%.4f,%.2f,%.4f\n",
                t, velocita_attuale, v_riferimento,
                errore, forza_motore, accelerazione);

        /* 7. Stampa console ogni PRINT_DECIMATION passi */
        if (k % PRINT_DECIMATION == 0) {
            /* Indicatore visivo della saturazione */
            char sat_flag = (forza_motore >= FORZA_MAX_N) ? '!' :
                            (forza_motore <= FORZA_MIN_N) ? 'v' : ' ';
            printf("  %7.2f  %7.3f m/s    %7.3f m/s   %+8.4f     "
                   "%7.1f N%c   %+6.3f\n",
                   t, velocita_attuale, v_riferimento,
                   errore, forza_motore, sat_flag, accelerazione);
        }
    }

    fclose(csv);

    /* ── REPORT FINALE ── */
    printf("\n");
    printf("╔══════════════════════════════════════════════════════════╗\n");
    printf("║                  REPORT DI PERFORMANCE                  ║\n");
    printf("╠══════════════════════════════════════════════════════════╣\n");
    printf("║  ISE  (Integral Squared Error)  : %10.4f             ║\n",
           met.ise);
    printf("║  IAE  (Integral Absolute Error) : %10.4f             ║\n",
           met.iae);
    printf("║  Overshoot massimo              : %+9.3f %%            ║\n",
           met.overshoot_pct);
    printf("║  Velocità picco                 : %10.4f m/s         ║\n",
           met.velocita_max);

    if (met.tempo_salita >= 0.0)
        printf("║  Rise time  (10%%→90%%)          : %10.3f s           ║\n",
               met.tempo_salita);
    else
        printf("║  Rise time  (10%%→90%%)          :    n/d               ║\n");

    printf("║  Settling time (±%.0f%%)          : %10.3f s           ║\n",
           BANDA_SETTLING * 100.0, met.tempo_assestamento);
    printf("╠══════════════════════════════════════════════════════════╣\n");
    printf("║  Log salvato in: cruise_log.csv                         ║\n");
    printf("╚══════════════════════════════════════════════════════════╝\n\n");

    return EXIT_SUCCESS;
}