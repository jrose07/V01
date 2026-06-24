1.  Woher stammen kosmische Myonen?:
	 - Kosmische Strahlung (hochenergetische Protonen) wechselwirken mit Atmosphäre und dabei entstehen hautpsächlich Pionen und teilweise Kaonen.  \
		Diese zerfallen nach kurzer Zeit dann einfach in Myonen und myonen-neutrinos
	- Entstehen in ca. 10 km Höhe und zerfallen schwach in Elektronen [[Q1_ATP2.pdf#page=218&selection=55,0,55,81|Q1_ATP2, page 218]] []()
 2.  Lebensdauer $\tau$ ist $N(t) = N_0 \cdot exp(-t/\tau)$ , aus Halbwertszeit $T_{1/2} = \tau \ln(2)$ .
    Für Myonen ist $\tau \approx 2.197 \, \mathrm{µs}$  
    
 3. Myonen bewegen sich mit quasi Lichtgeschwindigkeit $v\approx c$.
		 1. Klassisch:  
		    $S = c * \tau = 2.197 µs \cdot 3e8 m/s \approx 659,1 m$ 
		 2. Relativistisch:
		    $\gamma = E_\mu / (m_\mu c^2) = \frac{10Gev}{105,658Mev} \approx 94,64$ 
		    $\implies S = \gamma \tau \cdot c \approx 62km$.

4.  Der Fluss der Myonen auf Meereshöhe liegt bei ca. $1 cm^-2 min^-1$  [[Q1_ATP2.pdf#page=217&selection=27,0,31,11|Q1_ATP2, page 217]]  
		Szintilatortank: $V = A \cdot h = \pi r^2 \cdot 2r = 2\pi r^3 \leftrightarrow r = \left(\frac{V}{2 \pi}\right)^{1/3} \implies r \approx  19.96 \mathrm{cm} \; \text{and} \; A =  1251 \mathrm{cm}^2$   
		Damit :
			$R = \Phi \cdot A = 1 \cdot 1251 \, \mathrm{min^-1} \approx 20 Hz$ 

5. Szintillator: Detektor ionisierender Strahlung der bei Energiedisposition der Strahlung Lichtblitze erzeugt.	Unterscheidungen und verschiedene Typen:
	1. Anorganische Szintillationskristalle:
	   Szinitillationsmechanismus hängt stark von Bandstruktur ab. 
	    - Vorteile:
		    - Deutlich höhere Energieausbeute und Energieauflösung -> Spektroskopisch relevant
		- Nachteile:
			- Langsamere Lichtsignale, länger andauernd (Schlechtere Zeitauflösung) und teurere Materialien.
	2. Organische Szinitllationsmaterialien:
	   Plastik, Flüssigkeiten, auf Kohlenstoff basiert
		- Vorteile:
			- Günstig, man kann also große Flächen damit versehen
			- Schnelle Lichtisgnale (ns-Bereich) -> Gute Zeitauflösung
		- Nachteile:
			- Geringere Lichtausbeute und damit geringere Energieauflösung

6. Maßnahmen:
	1. Diskriminator mit variabler Schwelle:
	   Die Signale der PMT's an Ein- und Ausgang werden mit einem Diskriminator in einer *Koinzidenzschaltung*  verbunden. An dem Diskriminator kann dann eingstellt werden, wie verzögert die Signale sein dürfen. 
	2. Suchzeit $T_s$ kann am Monoflop vorgegeben werden. Es wird die Zeit am TAC gemessen, wenn das Myon eintritt und wenn es zerfällt und maximal-Suchwert ist $T_s$ (maximale Zeit, damit nicht noch mehr Events gezählt werden als sollten.)
- Untergrundrate:
  $P(X>=1) = 1 - e^{-R\cdot T_s}$ mit $R$: Konizidenzrate, $T_s \approx 50 \mathrm{ns}$: Suchzeit.
	$R \approx R_\mu (aus 4.) \approx 20 \mathrm{Hz} \implies P = 1e-6$  (Sehr klein!)
7. Vielkanalanalysator:
   VKA misst Impulshöhe des Signals und ordnet diese je nach Kalibration einem Kanal zu
		$\implies$ Es misst also hier, da das Signal aus dem TAC kommt die Messzeiten (Zerfalllszeiten) der Myonen und stellt sie histogrammatisch auf.
	Wahrscheinliche Verteilung (Myon ist MIP):
	![[Pasted image 20260608092828.png]]

8. Myonen Stoppen in Material -> Zerfallen in Elektronen nach gewisser Verzögerungszeit -> Neues Signal 