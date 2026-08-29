# Projektanweisungen

## Lieferung von Artefakten

Wenn ein Auftrag Abbildungen oder andere Artefakte verlangt, ist er erst
abgeschlossen, nachdem die Dateien nachweislich existieren, nicht leer sind und
dem Benutzer tatsaechlich angezeigt wurden. Ein Pfad oder die Behauptung einer
Erzeugung genuegt nicht. Vor jeder Abschlussmeldung ist die vollstaendige
Lieferliste des Auftrags nochmals gegen die tatsaechlich vorhandenen Artefakte
zu pruefen.

Konkret vor jeder Abschlussmeldung:

1. Auftragstext erneut lesen und jede geforderte Lieferung einzeln auflisten.
2. Jede Datei auf dem Dateisystem pruefen: Existenz unter absolutem Pfad,
   Dateigroesse groesser null, gueltiger Dateityp, bei Bildern zusaetzlich die
   Abmessungen.
3. Erfolgreichen Abschluss (Exit-Code) des erzeugenden Skripts oder Programms
   pruefen.
4. Die Abbildungen in der Antwort als sichtbare Bildvorschau ausgeben. Eine
   Liste von Pfaden ist kein Ersatz.
5. Scheitert die Anzeige technisch, den exakten Werkzeugfehler melden und den
   Auftrag nicht als erfolgreich beenden.

Ein bestandener Testlauf ersetzt eine geforderte Lieferung nicht.

## Reproduzierbarkeit von Abbildungen

Abbildungen entstehen ausschliesslich aus versionierten Skripten in `python/`
und aus Daten, die ein Werkzeug aus `tools/` in ein Ergebnisverzeichnis unter
`results/` geschrieben hat. Keine improvisierte Einmalgrafik ohne
reproduzierbare Eingangsdaten.
