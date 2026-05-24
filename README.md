# Monte Carlo Simulation (GARCH)

Simulerar framtida prisbanor för att utvärdera en swing-trade (lång eller kort) med stop-loss och trailing stop.

## Översikt

Skriptet kör 5 000 simuleringar över 20 handelsdagar för att besvara:
- Vad är det förväntade utfallet (EV) för denna trade?
- Hur ofta vinner/förlorar strategin?
- Är payoff ratio tillräcklig för att vara lönsam?

## Metodik

### GJR-GARCH(1,1,1) — Volatilitetsmodellering

Volatilitet är inte konstant. Efter en stor rörelse tenderar fler stora rörelser att följa ("volatility clustering"). GJR-GARCH-modellen fångar detta, inklusive att negativa chocker ofta ger *högre* volatilitet än positiva (leverage-effekten):

```
σ²(t) = ω + (α + γ·I(ε<0))·ε²(t-1) + β·σ²(t-1)
```

där σ² är variansen, ε är gårdagens avvikelse från medel, I(ε<0) är 1 vid negativ avkastning (leverage-term), och ω, α, γ, β skattas från historisk data. Residualerna modelleras med **Student-t-fördelning** för att bättre fånga feta svansar.

### Filtered Historical Simulation (Bootstrapping)

Istället för att anta normalfördelade avkastningar drar skriptet slumpmässigt från aktiens faktiska historiska "chocker" (standardiserade GARCH-residualer). Detta bevarar:
- **Feta svansar** — verkliga krascher/rusningar är vanligare än normalfördelningen antar
- **Skevhet** — asymmetri mellan upp- och nedrörelser

### Stop-mekanik

Varje simulerad bana avslutas om:
1. **Stop-loss** — priset passerar en fast nivå (default: 2×ATR(14) från start)
2. **Trailing stop** — priset passerar en rörlig nivå som följer extrempunkter (default: 2×ATR(14))

Stop-storlek skalar med tickerns volatilitet så att tight stops inte slår ut normala dagliga rörelser. Overrides: `NxATR` (t.ex. `1.5xatr`), procent (`0.03`), eller absolut pris (t.ex. `23500`). Detta matchar `summary.py` som också använder 2×ATR för stop-sizing.

Riktningen på stop-logiken beror på positionstyp:

| | Lång (default) | Kort (`--short`) |
|--|----------------|-------------------|
| Stop-loss | Priset faller under `start × (1 - pct)` | Priset stiger över `start × (1 + pct)` |
| Trailing stop | Trackar toppar, stoppar under `peak × (1 - pct)` | Trackar bottnar, stoppar över `trough × (1 + pct)` |
| Vinst | Priset stiger | Priset faller |
| Default target | `start × 1.07` (+7%) | `start × 0.93` (-7%) |

## Nyckeltal i rapporten

| Mått | Tolkning |
|------|----------|
| **EV (Expected Value)** | Genomsnittlig vinst/förlust per trade. Positivt = statistisk edge. |
| **Win Rate** | Andel lönsamma trades. Trendföljande strategier har ofta 35–45%. |
| **Payoff Ratio** | Medianvinst / Medianförlust (PnL) över alla simulerade paths. Median valt över medel för att inte låta en handfull fat-tail-utfall blåsa upp ration. Kompenserar låg win rate om Payoff > 1.5. Skiljs från "R/R" i `levels.py`, som baseras på faktisk target/stop från S/R-nivåer — Payoff Ratio kommer ur fördelningen, inte en pre-trade plan.|
| **Break-even win rate** | Minsta win rate för att gå ±0, givet payoff: `1 / (1 + Payoff)` |
| **Half Kelly** | Konservativ positionsstorlek = 50% av Kelly = `0.5 × (win_rate × payoff − loss_rate) / payoff`, där `payoff` = Payoff Ratio ovan. |
| **Median Exit (winners)** | Typisk exitnivå vid vinst — användbart som kursmål/target. |
| **Median Exit (losers)** | Typisk exitnivå vid förlust — visar var stopparna biter. |
| **Target Price** | Ingångsparameter, inte en prediktion. Default: +7% (lång) / -7% (kort) från start. |
| **P(≥ Target) / P(≤ Target)** | Sannolikhet att nå målpriset (riktning beror på lång/kort). |
| **Stopped paths** | Andel simuleringar som träffade stop-loss eller trailing stop. |

### Median Exit (winners) som target price

Av de tillgängliga måtten är **Median Exit (winners)** det bästa att använda som kursmål/target price:

- **Mean Final Price** dras upp av enstaka extrema simuleringar (feta svansar). Det ger ett mål som ser mer optimistiskt ut än vad majoriteten av utfallen faktiskt levererar.
- **Median Exit (winners)** filtrerar bort förlorande simuleringar och ger mittpunkten av de banor som faktiskt gick rätt väg — ett realistiskt mål för en lyckad trade.

Det viktiga: Median Exit (winners), det är där typiska vinnande paths toppar innan trailing stop slår till. Att sätta target vid denna nivå betyder att du plockar ut ~max av fördelningen — varje krona därutöver är extremt dyrköpt sannolikhetsmässigt.

### Är strategin lönsam?

En strategi är lönsam om `Win Rate > Break-even Win Rate`, där:

```
Break-even = 1 / (1 + Payoff Ratio)
```

Exempel: Payoff = 1.8 → `1 / (1 + 1.8) = 35.7%`. Om din win rate är 42% har du en edge på 6.3 procentenheter.

## Grafen

Grafen består av tre paneler:

### Huvudpanel (uppe till vänster) — prisbanor

- **Grönt fält** — 50% konfidensintervall (p25–p75)
- **Blått fält** — 90% konfidensintervall (p5–p95)
- **Lila linje** — Median (50:e percentilen)
- **Blå streckad** — Medelvärde (dras upp av stora vinster om positiv skevhet)
- **Grön horisontell** — Startpris
- **Orange streckad** — Målpris (target)
- **Röd prickad** — Fast stop-loss-nivå
- **Rosa streckad-prickad** — Trailing stop (median av aktiva banor)
- **Tunna grå linjer** — sampel av enskilda banor (upp till 15 stoppade + 15 överlevande). Stoppade banor klipps vid sin stop-dag så det syns var de träffar SL/trailing. Ger en känsla för spridningen som percentilbanden döljer.

Konfidensbanden beräknas endast på *aktiva* (ej stoppade) banor — annars skulle stoppade banor frysa nere vid stop-nivån och dra ned percentilerna artificiellt. När väldigt få banor återstår (< 0.5%) blankas banden ut.

### Stoppanel (nere till vänster) — kumulativ knockout

Visar hur stor andel av alla simuleringar som har stoppats ut vid varje given dag, uppdelat på stop-typ:

- **Röd area (Fixed SL)** — banor som stoppades utan att någonsin ha rört sig gynnsamt. Inkluderar både rena träffar på den fasta stop-loss-nivån och trailing-stoppar som utlöstes innan priset gick över start. Alltid förlust.
- **Orange area (Trailing)** — banor där priset först rörde sig gynnsamt (extrempunkten passerade start) och sedan retracerade in i trailing-stoppen. Kan vara vinst (priset hann nå en nivå klart över start innan vändning) eller förlust (toppen var bara marginellt över start, trailing-nivå hamnade under start).
- **Svart linje** — total andel stoppade

Areorna är *staplade*, så den övre kanten av orange = total stop-out. Läs t.ex. "dag 10: 35%" som "35% av alla simuleringar har stoppats ut senast dag 10". Brant lutning = stop biter hårt; platt = strategin överlever.

Trailing-andelen som fångas vid vinst respektive förlust syns i textrapporten och i statistikrutan i grafen (`30%V/15%F` = 30% vinst, 15% förlust av trailing-stoppade).

### Histogrampanel (höger) — slutprisfördelning

Horisontellt histogram över priserna på sista dagen:

- **Grön** — banor som slutade i vinst
- **Röd** — banor som slutade i förlust
- Stop-loss/start/target-linjerna sträcker sig över både huvudpanel och histogram för att visa var massan ligger relativt nivåerna

Skev fördelning åt höger (lång svans uppåt) = positiv skevhet; mean ligger då högre än median.

### Tolka median vs startpris (survivorship-bias)

**Viktigt:** Den lila medianlinjen visar medianen av *aktiva* banor — alla som stoppats ut är borttagna ur beräkningen. Det är alltså en **villkorad median**: "givet att traden fortfarande lever, var ligger 50:e percentilen?"

Konsekvenser att vara medveten om:

- Linjen tippar systematiskt uppåt över tid — inte för att marknaden förbättras, utan för att stop-loss plockar bort förlorarna ur urvalet (survivorship bias by design)
- "Medianen ligger över startpris dag 15" betyder **inte** att >50% av alla simuleringar har vinst — det betyder att av dem som överlevt till dag 15 ligger 50% över start
- För helhetsbilden måste medianen läsas *tillsammans* med stoppanelen: t.ex. "median +4% över start, men 40% kumulativt stoppade" → faktisk vinstandel är ungefär 60% × 50% = ~30%

Den exakta vinstandelen för hela populationen finns i textrapporten som **Win Rate**, och i statistikrutan i grafen.

### Median vs EV

En strategi kan ha **låg/negativ Win Rate men positivt EV** om de få vinsterna är tillräckligt stora — typiskt för trendföljande strategier med tight stop-loss.

**Exempel — 100 trades med 3% stop-loss:**

| Utfall | Antal | Resultat |
|--------|-------|----------|
| Förlust (stoppas ut) | 60 | -3% |
| Liten vinst | 15 | +4% |
| Stor vinst (trenden höll) | 25 | +15% |

- **Win Rate:** 40% (40 av 100 i vinst)
- **EV:** `(60 × -3%) + (15 × +4%) + (25 × +15%) = +2.55%` per trade

Strategin förlorar 60% av gångerna, men de 25 trades där trenden drar iväg kompenserar för alla små förluster. Stop-lossen begränsar nedsidan till ett känt belopp, medan uppsidan är obegränsad. Därför är EV det viktigaste måttet — det säger vad du tjänar *över tid*.

## Modelldiagnostik

Rapporten inkluderar tre diagnostiska tester som validerar GARCH-modellens kvalitet:

| Test | Vad det mäter | Godkänt |
|------|---------------|---------|
| **Persistence (α+γ/2+β)** | Volatilitetens uthållighet. ≥ 1.0 innebär icke-stationär variansprocess. | < 1.0 |
| **Ljung-Box(20)** | Kvarvarande autokorrelation i kvadrerade residualer. Lågt p-värde = modellen missar mönster. | p ≥ 0.05 |
| **ARCH-LM(10)** | Kvarvarande heteroskedasticitet. Lågt p-värde = modellen fångar inte all volatilitetsdynamik. | p ≥ 0.05 |

Om ett test misslyckas visas en varning (⚠) i rapporten. Simuleringsresultaten bör då tolkas med försiktighet.

## Begränsningar

1. **Gap-risk** — Simuleringen bevarar det simulerade gap-priset som exit-kurs när banan hoppar förbi stoppet, vilket approximerar slippage. Verkliga gap (t.ex. vid rapport) kan fortfarande vara större än de bootstrappade dagsavkastningarna.
2. **Courtage/spread** — Ingår ej. Dra av mentalt.
3. **Regimskifte** — Bootstrapping bygger på historik. Fundamentalt nya marknadsförhållanden fångas inte.

## Användning

```bash
# Standardkörning: 20 dagar, 2×ATR stop-loss, 2×ATR trailing (lång position)
python montecarlo.py EQT.ST

# Annan aktie med custom parametrar (procent eller ATR-multipel)
python montecarlo.py VOLV-B.ST --days 30 --stop-loss 0.05 --trailing-stop 0.04
python montecarlo.py VOLV-B.ST --stop-loss 1.5xatr --trailing-stop 3xatr

# Fast stop utan trailing
python montecarlo.py ASSA-B.ST --no-trailing

# Specifik ingångskurs och målkurs
python montecarlo.py SAND.ST --start-price 250 --target 280

# Kort position (short) — inverterad stop/target-logik
python montecarlo.py VOLV-B.ST --short

# Kort position med specifikt målpris
python montecarlo.py ERIC-B.ST --short --target 60 --stop-loss 0.04

# Reproducerbar körning med fast seed
python montecarlo.py EQT.ST --seed 42

# Begränsa GARCH-fit till senaste N handelsdagar (default: hela historiken)
python montecarlo.py SIVE.ST --lookback 250
```

### `--lookback` — när och varför

Default fittas GARCH på hela CSV-historiken. För aktier som nyligen bytt regim (t.ex. parabolisk rusning från låg-vol-period) blir den långa historiken missvisande: konstant-mean-skattningen drar μ mot noll och simuleringarna ser inte den nya driften.

`--lookback N` trimmar serien till de senaste N handelsdagarna innan fit. Riktlinjer:

- **N ≥ 200** — säker zon, GARCH-stationaritet brukar hålla
- **N = 100–200** — ofta OK men kontrollera att Persistence (α+β) < 1.0 i diagnostiken
- **N < 100** — riskabelt, GARCH blir lätt icke-stationär (α+β ≈ 1.0) och resultaten exploderar; varningen i diagnostiken fångar detta

Modelldiagnostiken i rapporten är ditt skyddsnät — om Persistence ≥ 1.0 eller Ljung-Box/ARCH-LM faller, öka lookback eller kör utan flaggan.

## Output

- `montecarlo_TICKER.txt` — Statistikrapport
- `montecarlo_TICKER.png` — Visualisering med konfidensband och histogram
