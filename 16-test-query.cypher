MATCH (n)
RETURN labels(n), keys(n), count(*) AS count
LIMIT 20;

// LOAD CSV imports every field as text. Convert the existing relationship costs
// before projecting them into GDS.
MATCH ()-[r:EXCHANGE|TRAVEL_TO|EGRESS]->()
WHERE r.cost IS NOT NULL AND trim(toString(r.cost)) <> ''
SET r.cost = toFloat(r.cost);

CALL gds.graph.drop('transport-routing', false)
YIELD graphName
RETURN graphName;

CALL gds.graph.project(
  'transport-routing',
  ['StopSol', 'Stop', 'PoI'],
  {
    EXCHANGE: {properties: 'cost'},
    TRAVEL_TO: {properties: 'cost'},
    EGRESS: {properties: 'cost'}
  }
);