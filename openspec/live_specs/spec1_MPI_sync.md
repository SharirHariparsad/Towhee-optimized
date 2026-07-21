## OBJECTIVE
Eliminate redundant MPI synchronization overhead and convert blocking communication patterns into non-blocking, asynchronous sequences to reduce parallel latency and prevent serial bottlenecks.

## TARGET FILES
- `control.c` (Focus: lines 70–110)
- `towhee.c` (Focus: lines 130–170)
- `rex.c`     (Focus: lines 100–130)

## STEP-BY-STEP IMPLEMENTATION PLAN

### 1.Eliminate Redundant Barriers in control.c (Tramonto Mode Path)
- Locate the two code blocks where `MPI_Bcast` handles the `logicflag` synchronization (approx. lines 77–79 and 101–105).
- Identify the triple-synchronization pattern: `MPI_Barrier` -> `MPI_Bcast` -> `MPI_Barrier`.
- **Action:** Delete both the preceding and succeeding `MPI_Barrier(MPI_COMM_WORLD);` calls for both blocks. Keep only the core `MPI_Bcast` statement, as it inherently acts as a collective synchronization point.

### 2.Remove Duplicate Teardown Barriers in towhee.c
- Locate the termination/dispatch functions leading into `MPI_Finalize` (approx. lines 136–138 and 161–169).
- Identify instances where two back-to-back `MPI_Barrier(MPI_COMM_WORLD);` calls are invoked.
- **Action:** Remove the second, duplicate `MPI_Barrier` call in all three dispatch paths, leaving exactly one barrier to guarantee that all ranks have completed their workload before finalization.

### 3. Convert Sequential Handshakes to Async Isend in rex.c (Replica Exchange Path)
- Locate the loop sending lambda keys sequentially from the server to the client nodes (approx. lines 115–119):
  ```c
  for (k = 1; k <= paramwrap->num_nodes; k++) {
      m = J[k];
      MPI_Send(&m, 1, MPI_INT, k, LAMBDA_KEY, MPI_COMM_WORLD);
  }
  
** Action:** Refactor this block to use non-blocking communication to prevent serialization:

-Declare an array of MPI requests: MPI_Request reqs[MAX_NODES]; (Ensure MAX_NODES or an equivalent upper bound is available, or dynamically allocate based on paramwrap->num_nodes).
-Replace MPI_Send with MPI_Isend, passing the tracking request: &reqs[k-1].
-Immediately following the loop, inject a blocking MPI_Waitall(paramwrap->num_nodes, reqs, MPI_STATUSES_IGNORE); to guarantee all data is safely out of the buffer before proceeding.

## CONSTRAINTS & GUARDRAILS
-Do NOT modify or alter the actual payload logic, variable assignments, or data values being passed (e.g., logicflag, m, J[k]).
-Do NOT touch any code paths outside of these specific MPI message passing layers.
-Maintain strict standard C compliance; ensure any dynamically allocated memory for tracking arrays (if utilized) is appropriately freed.
