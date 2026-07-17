    PROPERTY GRAPH

    ---------------------------------------------------------------------------

    Nodes

    - id: Gate/Metro/Bus/Tram/Ferro/Destination + number 
    - name
    - type : gate/metro/bus/

    * Origin
        - id : O1
        - name: mataró
        - geo primera parada/sortida

    * Gate
        - id : G1
        - name : glòries (potser probar altres)
        - geo

    * Stops
        - id M/B/T/F + num
        - name
        - type : metro/bus/tram/ferro
        - is possible sol : T/F
        - geo

    * Destination
        - id: d + num
        - name
        - geo

    ---------------------------------------------------------------------------

    Edges

    * autopista
        - id a + num
        - origin : Origin
        - destination : Gate
        - time

    * intra-city - precalc with r5py osmxn...
        - id
        - origin : gate | node stop wehre is possible stop : T
        - destination :  node stop wehre is possible stop : T
        - speed: 18 km/h (maybe no need if graph if we obtain time)
        - distance
        - time = distance / speed

    * exchange (per trasbords dins d'un mateix metro assumim)
        - id
        - origin: node of type Stop
        - destination : same node | another node of type stop
        - time
        - exchange weight
        - wait time
        - wait weight
        - percieved time = time * exhange weight
        - percieved wait time = wait time * wait weight

    * metro travel
        - id 
        - origin : node stop of type metro
        - destination : node stop of type metro
        - line : L1/L2/L3....
        - distance
        - speed : 25 km/h
        - time : distance / speed
        - time weight
        - percieved travel time = time * time weight

    * bus travel
        - id 
        - origin : node stop of type bus
        - destination : node stop of type bus
        - line : H2/D20/7....
        - distance
        - speed : 12 km/h
        - time : distance / speed
        - time weight
        - percieved travel time = time * time weight

    *  tram travel
        - id 
        - origin : node stop of type tram
        - destination : node stop of type tram
        - line : T1/T2/T3...
        - distance
        - speed : 22 km/h
        - time : distance / speed
        - time weight
        - percieved travel time = time * time weight


    *  ferro travel
        - id 
        - origin : node stop of type ferro
        - destination : node stop of type ferro
        - line : L6/L7/L8....
        - distance
        - speed : 25 km/h (pot cambiar)
        - time : distance / speed
        - time weight
        - percieved travel time = time * time weight

    * walk travel (consider walk if Destination is in a buffer of ~500 m of a stop) - precalc with r5py osmxn...
        - id
        - origin: a node of type Stop
        - destination : node of type Destination
        - speed : 5km/h (maybe no need if graph if we obtain time)
        - distance 
        - time = distance /speed
        - egress weight
        - percieved time  = egress weight * time
