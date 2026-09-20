MATCH (source:StopSol),
      (destination:POI {poi_name: 'Zona Universitària'})

CALL gds.shortestPath.dijkstra.stream('transport-routing', {
  sourceNode: source,
  targetNode: destination,
  relationshipWeightProperty: 'cost'
})
YIELD totalCost, nodeIds, costs

WITH
  source AS route_source,
  destination AS route_destination,
  totalCost,
  nodeIds,
  costs

UNWIND range(0, size(nodeIds) - 2) AS i

WITH
  route_source,
  route_destination,
  i + 1 AS step,
  totalCost,
  gds.util.asNode(nodeIds[i]) AS origin,
  gds.util.asNode(nodeIds[i + 1]) AS destination,
  costs[i + 1] AS cumulativeCost

// Prefer the physical relationship in the path direction. Only fall back to
// the reverse relationship for relationships projected as UNDIRECTED.
OPTIONAL MATCH (origin)-[forward]->(destination)

WITH
  route_source,
  route_destination,
  step,
  totalCost,
  origin,
  destination,
  cumulativeCost,
  collect(forward)[0] AS forward_relationship

OPTIONAL MATCH (origin)<-[reverse]-(destination)

WITH
  route_source,
  route_destination,
  step,
  totalCost,
  origin,
  destination,
  cumulativeCost,
  forward_relationship,
  collect(reverse)[0] AS reverse_relationship

WITH
  route_source,
  route_destination,
  step,
  totalCost,
  origin,
  destination,
  cumulativeCost,
  CASE
    WHEN forward_relationship IS NOT NULL THEN forward_relationship
    ELSE reverse_relationship
  END AS relationship

RETURN
  // add route_id as source-destination string
  route_source.id as from,
  route_destination.id as to,
  coalesce(route_source.id, route_source.poi_name) + ' to ' + coalesce(route_destination.id, route_destination.poi_name) AS route_id,
  step,
  coalesce(origin.id, origin.poi_name) AS origin,
  coalesce(destination.id, destination.poi_name) AS destination,
  relationship.tram as tram,
  relationship.type AS edge_type,
  relationship.cost AS edge_cost,
  cumulativeCost,
  totalCost,
  relationship.geometry AS geometry

ORDER BY route_id, step;

