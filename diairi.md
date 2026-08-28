# Dies abans de comemçar el diari 

# Transportation Data

## Metro

### TMB
- Clean segment-by-segment track data:
  https://developer.tmb.cat/data/recs#trams_linia_metro

### Barcelona Open Data
- Includes metro stations, but the data is less structured and does not include trajectories:
  https://opendata-ajuntament.barcelona.cat/data/ca/dataset/transports

### OpenStreetMap (OSMnx)
- Contains station locations.
- Data quality is lower (e.g. some station names are missing).

---

## Bus

### TMB

**Stops**
- General stop dataset (does not include the lines serving each stop):
  https://developer.tmb.cat/data/parades
- Separate dataset with one line per stop:
  https://developer.tmb.cat/data/parades#parades_linia

**Routes**
- Complete route geometries (not split into segments):
  https://developer.tmb.cat/data/recs#recs

### AMB

Stops can be downloaded as an Excel file containing:
- Stop information
- Geographic coordinates
- Served lines

https://www.ambmobilitat.cat/Principales/busquedaparadas.aspx

### Notes

- AMB includes more stops because it also covers municipalities such as L'Hospitalet and Esplugues.
- TMB is missing some of these peripheral stops.
- Plan:
  - Use the AMB stop dataset.
  - Exclude stops that are only served by night buses (lines starting with `N`).
- Route geometries are not available for non-TMB AMB lines (e.g. line 87), but these lines have relatively low ridership, so this is not considered a major issue.

---

## FGC

### Stops

FGC Open Data provides stop information, but not the list of lines serving each stop:

https://dadesobertes.fgc.cat/explore/assets/gtfs_stops/

### Routes

Route geometries for each line (whole routes, not split into segments):

https://fgc.opendatasoft.com/explore/assets/gtfs_routes/view/

### Processing decisions

- Metro and FGC stations that share the same physical station will be merged into a single node.
- Many suburban FGC services share the same infrastructure within Barcelona. These will either be represented by a single line or merged together.

Equivalent lines:

```text
S3, S4, S8, R50, S9, R53, R5, R6, R60, R63, L8
```

Equivalent lines:

```text
S1, S2
```

Lines kept:

```python
['L6', 'L7', 'L8', 'S1','L12']
```


### Shared Metro–FGC stations

| Metro station | FGC line(s) |
|--------------|-------------|
| Plaça Catalunya | L6, L7 |
| Diagonal / Provença | L6, L7 |
| Av. Carrilet | L8 |
| Espanya | L8 |
| Europa \| Fira | L8 |

estan a l'overleaf

# 20 juliol

Crear nodes del metro. 
Afegir Linies de fgc a parades de tipus metro. 
Decidir que les parades de fgc seran nodes diferents, ecara que linies de fgc puguin estar en parades de metro

He descobert que la L6 no arriba a reina eliseda, pero a les dades de trajectòries de la L6 si que hi arriba -> corregir.
Per tant hem d'afegir la L12 encara que siguin només dos parades, ja que arriba a una destinació que no arriben altres linies
Dels nodes de fgc ens carreguem vallvidriera inferior i superir que 

Considerant si afegir o no les parades fora de bcn. Per exemple a l'hospitalet hi ha certes zones d'oficines. Segurmanet decidir un cop tingui els punts d'interès

# 21 Juliol

Rebut dades de AMB, encara no accessibles degut a servei no disponible
Construir nodes tram. Les dades de ide amb no son bones perque estan desactualitzades -> ocontruirem el map paparada - linia manual igual que el de fgc
Dubte: excloure les parades de Tram que tenen metro de la solució (segurmanet si)
Si que ho fem, pero no modifiquem la seva geometria com en el cas del fgc, ja que pot alterar a les rutes d'enllaç i la llargada de les trajectòries.

['Ernest Lluch','Zona Universitària','Palau Reial','Maria Cristina']
['Fòrum','El Maresme','Selva de Mar','Glòries','Marina','Ciutadella | Vila Olímpica']
['Gorg','Sant Roc','Besòs','Verdaguer']


Començat a buscar dades de punts d'interès de bcn
https://www.institutmetropoli.cat/en/estudi/mobility-survey-on-weekday-2024-emef/?utm_source=chatgpt.com
https://portalrecerca.uab.cat/es/publications/walking-short-distances-the-socioeconomic-drivers-for-the-use-of-/?utm_source=chatgpt.com
https://ce.atm.cat/en/web/observatori/w/repartiment-modal-tipus-despla-comarca-residencia?utm_source=chatgpt.com

https://arxiv.org/abs/2103.15638?utm_source=chatgpt.com
https://www.sciencedirect.com/science/article/pii/S0965856414002419?utm_source=chatgpt.com

lo millor trobat és això: https://www.amb.cat/ca/web/mobilitat/gestio-i-organitzacio/dades/dades-de-fluxos
dades de fluxe per districte, pero els districtes son massa generals


Potser és una bona opció només considerar parades de bus com a possibles parades, i considerar tram/metro nomñes si no hi ha alguna parada de bus a prop
Dades d'AMB són molt semblants a les de TMB pero arriben millor a l'hospitalet i sant adria. Hi ha bastantes més parades per poc més coverage potser renta no fer-ho servir
Començat a generar nodes per bus, pero falten les parades per linia de hos

## 22 de juliol

el paper -A city of cities: Measuring how 15-minutes urban accessibility shapes human mobility in Barcelona- és força bo, i té dades sorbe les zones de destí i les seves frequencues
Les dades no són públiques ja que són de vodafone, però poden ser compartides escrivint-li a aquest senyor https://ajuntament.barcelona.cat/es/organigrama-municipal/organ/Oficina%20Municipal%20de%20Datos

preguntar a ferran -> https://cit.upc.edu/en/portfolio-item/second_phase_vml/

OD data by 10 districts, not granular enough -> https://ropenspain.github.io/spanishoddata/

Idea hexagons + scores
faltaria afegir les geometries de tota la AMB


# 23 de juilol 

Hexagons -> només els creem a bcn i a les de les vores ja hi arribarem, així que les parades candidates seràn properes al centre.
Potser fins i tot filtrar la capa d'hexagonns més propera a la frontera
Score parades -> no contar algo tenint en compte les linies de la parada + distancia altres parades + les linies que tingui. Així una parada amb més linies sera beneficiada. No contar una línia dos cops si té dos parades dins el buffer. 

Aquest paper tenia bona pinta pero fa al reves -> busca ocupar les zones on no hi ha parades : https://onlinelibrary.wiley.com/doi/epdf/10.1155/2020/8853872

Aquest té bona pinta pero no el puc llegir : https://www.sciencedirect.com/org/science/article/abs/pii/S2324993522004699
Demano accés als autors

Aquest tb té pintaca: https://www.emerald.com/jtran/article-abstract/doi/10.1680/jtran.25.00025/1303212/Design-method-of-feeder-bus-network-based-on-belt?redirectedFrom=fulltext



# 24 de juliol

llegit aquest paper i buscar-ne més https://onlinelibrary.wiley.com/doi/full/10.1155/atr/9728885
es bastant útil, pero no considerar integració, només parades als punts de demanda

he inventat 3 destinacions
he creat els edges de les linies de metro i la primera linia de fgc (L6)
no estandaritzem els punts e fgc i les parades de metro, ja que putejen a la projecció de les linies de fgc a lhora de crear els segments

# 25 jul - 7 ago

Escriure correus per conseguir dades
Les de moventis no tenen dades de les destinacions finals, sinó de com arriben a les parades i a les destinacions
Les del paper, encara estic esperant reposta, un cop enviat el document que justifica lús de les dades
Idea -> mantenir stop id, per tal de crear els edges 'recursius després'


# 8 agost

Afegir la columna de stop_id, per aconseguir fer mes facils els temps d'espera del trasbord
A trambesos li falten parts dels recorreguts -> enviar correu a OPENDATA demanant-los
la trajectoria de la linia 7 de la font amb està malament, el que fà que no es porjectin bé les parades
fa falta mirar si les altres tb estan malament i potser es poden 'imputar' desde tmb algunes
les trajectories cicliques no es porjecten bé, ja que al fer-ho el bucle no sap quines aprades són d'anada i quines de tornada


# 9 agost

He provat de fer servir les dades de gtfs de amb, però tenen molt poques linies dins de barcelona, i sobretot tenen línies de la perfifèria
He eliminat algunes linies de bus problmatiques (+2 sentits) i que es trobaven a la periferia.
Altres que si que aportaven valor les he manitingut (88)
a les dades de amb falten parades com les de la linia 115 (eg 1771), parada 32 de la linia 7...
fer un merge amb tmb

el bus 127 te les trajectories d'anada i tornada malament

linies on he trobat probles:

PR3, M27, 65, 34 (només lúltim punt), 114, 155 (i falten parades), 121, 127 (disjoint), 133, 180 es raro crec q falten punts, 196, 197, 34 (nms lultim), 39,52, 59, 65
7, 86, 87, b14, b15, b18, b23, b25, b81, h10, h12, h14, lh1, lh2, m27, v11, v17, v25, v3, v31, v5 (nomes ultima), v7, v9, x1


el que fare sera mantenir la geometria de tmb si es troba en els dos datasets. En el proces tb afegim alguna linia més. No eliminar cap linia de moemnt

------

Per tractar amb les parades que faltaven simplement he fet servir les geografies de l'excel, enlloc de fer un merge excel tmb + amb


# 11 d'agost


hem arreglat el problema de r5py instalant : https://adoptium.net/es/download?link=https%3A%2F%2Fgithub.com%2Fadoptium%2Ftemurin21-binaries%2Freleases%2Fdownload%2Fjdk-21.0.12%252B8%2FOpenJDK21U-jdk_aarch64_mac_hotspot_21.0.12_8.pkg&vendor=Adoptium


Començat el tutorial de r5py. Hi ha una funció que dona temps de trajecte de O a D per diferents modalitats. Un dels inputs importants és l'elevació del terreny. Buscarem dades com ho va fer el ferran per el treball de DMT 

altres alternatives : https://visors.icgc.cat/appdownloads/?c=fmevtmet

r5py requereix un fitxer en format pbf, l'estic buscant a aqui https://download.bbbike.org/osm/bbbike/Barcelona i aqui : https://download.geofabrik.de/europe/spain/cataluna.html
Hem conractat a la gent del tram. Falta part del recorregut peruqe estan en obres. Eliminarem les parades que no estan operatives. EL T4 torna a estar operatiu el 14/09

# 12 Agost

Plantejat fer servir r5py per calcular ruutes mes curtes entre parades - destinacions
Fa falta gtfs, el que comporta inconvenients:
1- GTFS AMB no té les millors linies
2- Solució dependria de la hora

Farem sevir r5py per calcular els camins entre parades de bus properes i per bus interurbà.
Veig que fa linestrings entre vertexs de carrers, no punt a punt exacte, crec que no hi ha solucio

Considerem fer servir iscorones de 5 min, enllco de buffers de 300 m per considerar una parada com a propera. La funcií isochornes retorna un linestring enlloc de multipolygon, la qual cosa no mola. Hem de mirar com convertirho per fer un within.

# 13 d'agost


Hem trobat quer per walk podem fer el multilinestring -> polygon facil
Hem afegit l'elevation model i comprovat que el temps de A -> B != B -> A per pendents
La parada de metro de can serra tenia la ubicació mal posada, enmig d'un edifici, el que feia que no crees be lisochrone, l'he corregit i ara va bé. Quantes més n'hi haura així??
La parada de santa eulalia estvaa sobre una teulada, el que feia que no es creeisn els isochrones
Arreglat els punts que donaven problemes amb els isocrones
Els isocrones no semblen del tot bons ni realistes amb r5py, amb osmxn semblem millors.

També seran millors les rutes amb osmxn?

# 14 agost

provar de fer els iochrones amb osmxn
queden millor, però al fer les routes després superen els 5 minuts
amb         point_grid_resolution=35, millora la forma dels isochornes de r5py, pero tarden molt. SI redueixo molt peta el kernel
solució fer-los amb osmxn i dp filtar les que surtin amb + 5 min

falta executar tot


# 15 d'agost

el punt de la parada 110024 estava malament i feia que no pogues calcular ruta amb alguna parada que te a prop
Algunes parades de bus de amb i de metro tmb tenien els mateixos ids


veig que a vegades r5py calcula temps raros de trajecte i no son reals. He tornat a executar el mateix codi i ha donat temps realistes
he corregit la speed walk a 5 kmh 

He acabat de calcular tots els dedges, dema tocara fer els edges del bus IU i els de egress (almenys preparar)


# 17 d'agist

he creat els edges egress
next edges IU, wait time, afegir pesos a egress, exchange wait...


# 19 d'agost

Fem els edges del bus IU, primer de la city gate a tot i dp entre ells i bidreeccionals. No fem servir el temps que ens donen els outputs de r5py, sino que el calculem ambla distancia ja que podem adaptar el temps als parametres del paper

Tenim alguns problemes amb les parades candidates, al havernhi moltes (sobretot bus) no podem calcular tots els edges entre elles pq tardaria massa. Solucins -> cluster, considerar parades amb +1 linia, eliminar algines línies... 

He reduit les de metro, fgc, i tram peqrque estiguin només a barcelona, però el major probelma és el bus

em peto les linies de bus barri de les solucions, no dels edges

passem de 2235 parades a 1631, el temps esperat de calcular totes les rutes seguei sent 8h, l'hem reduit a la meitat, pero no es suficient. Em carrego les linies interurbanes també, comencen amb L,M,B


# 24 agost

Seguir eliminanr paradees
Alguns busos barri que ens vam deixar
He eliminat les parades de les linies de bus barri i les que tenien només una linia de bus
El temps esperat d'execució és de 2 hores


# 25 agost

he executat les rutes del bus IU
el kenrel peta al aplicar la geometria amb wkt.loads, chucnkejarho pot solucinar-ho
eliminem rutes massa curtes o massa llarges

# 28 d'agost

Aconseguit calcular edges IU bus, però segurament fara faltar reduir mes les conexions pq n'hi ha casi 1M. 
Idees :
 - trahectes > 5 min i < 20.
 - Petar mes parades de la zona de dlat i només pearades amb +2 linies enlloc de +1 

To - do aquesta tarda -> executar i guardae egress
I wait times després

Comencem amb els wait times, els fem per cada mode de tranport per spearat pq alguns ids coincideixen

Els gtfs del metro no se ben bé com estimar wait times per a les parades que tenen més d'una linia
Pensaba que serien els trip ids pero per exempl id 216 no te cap trip id que surti de la l5
Tb hi ha duplciats de direccio i temps, parada com si surtissin +1 a la vegada


Demà continuar amb wait tiemes, fer un comput global de temps i provar de passar ja a neo4j per calcular millor ruta Soltion -> Dest (no fa falta passar IU trajs)