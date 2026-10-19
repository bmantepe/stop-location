CALL gds.graph.drop('transport-routing')
YIELD graphName
RETURN graphName;

CALL gds.graph.project(
  'transport-routing',
  ['StopSol', 'Stop', 'POI'],
  {
    EXCHANGE: {properties: 'cost',orientation: 'NATURAL'},
    TRAVEL_BUS: {properties: 'cost',orientation: 'NATURAL'},
    TRAVEL_RAIL: {properties: 'cost',orientation: 'UNDIRECTED'},
    EGRESS: {properties: 'cost',orientation: 'NATURAL'}

  }
);