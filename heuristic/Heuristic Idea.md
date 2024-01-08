# Heuristic Idea

- _Idea 1_: Create disks for launching/returning points for each one of the graphs to be visited, parameterizing the disks by center and radius.
- _Idea 2_: Discretize the space of feasible centers on a regular grid:
  1. Solve a TSP using, for each graph, the closest point to the origin. This solution determines the initial sequence of visits.
  2. Generate a subgrid of lattice points for these vertices and choose randomly the launching / returning points within the lattice generated.
  3. Fixing the sequence, solve the restricted truck-drone problem for launching / returning points on the neighborhoods.
  4. Perform interchange heuristic on the neighborhoods of the sequence: interchange the positions, ciclic displacement, move centers of disks to points in the grid.

