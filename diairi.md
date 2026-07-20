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

- TMB includes more stops because it also covers municipalities such as L'Hospitalet and Esplugues.
- AMB is missing some of these peripheral stops.
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
['L6', 'L7', 'L8', 'S1']
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


