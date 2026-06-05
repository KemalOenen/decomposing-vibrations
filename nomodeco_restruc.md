# Plan: nomodeco Memory Efficiency & Tracking

## Context

nomodeco computes vibrational energy distributions for molecules. The current codebase stores internal coordinates (bonds, angles, dihedrals) as Python lists of string tuples like `("H1", "O2")`, runs O(n³) triple-nested Python loops to fill 3D density matrices, and has no memory tracking. The goal is to (1) replace string-tuple ICs with integer-index tuples for compact representation, (2) add `tracemalloc`-based memory profiling behind a CLI flag, and (3) vectorize the hotspot loops as a bonus computational win.

---

## Phase 1: Integer-Index Internal Coordinates

**What changes**: IC tuples go from `("H1", "O2")` to `(0, 1)` — integers indexing into the `Molecule` atom list. This propagates through every file that generates or consumes IC tuples.

### 1a. `nomodeco/libraries/nomodeco_classes.py`

Add two helpers to `Molecule`:
- `symbol_to_index(self) -> dict[str, int]` — `{atom.symbol: i for i, a in enumerate(self)}`
- `index_to_symbol(self) -> list[str]` — `[atom.symbol for atom in self]`

Convert all generation methods to return integer tuples:
- `covalent_bonds()`, `intermolecular_h_bond()`, `intermolecular_acceptor_donor()` → `list[tuple[int, int]]`
- `generate_angles()` → `(list[tuple[int,int,int]], list[tuple[int,int,int]])`
- `generate_dihedrals()`, `generate_out_of_plane()`, `generate_oop_planar_subunits()` → `list[tuple[int,...]]`

Inside each method: build `idx = self.symbol_to_index()` once at the top, then replace every appended symbol string with `idx[symbol]`.

Add `bond_angle_by_index(self, i, j, k)` and `actual_length_by_index(self, i, j)` that use `self[i].coordinates` directly — call these from the generation methods instead of the symbol-based versions.

`InternalCoordinates.add_coord_diff` and `add_coord_diff_linear` require no structural changes; they use set operations and length checks on tuples, both of which work on integer tuples.

### 1b. `nomodeco/libraries/bmatrix.py`

Currently builds `atom_index = {a.symbol: i for i, a in enumerate(atoms)}` and then does `atom_index[a]` for each element of each IC tuple. With integer tuples, `a` is already the index:
- Remove the `atom_index` dict construction
- Replace `atom_index[a] * 3` with `a * 3`
- Replace `coordinates[atom_index[a]]` with `coordinates[a]`

Apply uniformly in both `b_matrix()` and `b_matrix2()` across bond, angle, linear_angle, outofplane, and dihedral loop bodies.

### 1c. `nomodeco/libraries/specifications.py`

- `is_string_in_tuples(string, list_of_tuples)` → rename to `is_index_in_tuples(idx: int, ...)` and update callers to pass `atoms.index_of(atom)` instead of `atom.symbol`
- `bound_to_atom(query_atom, bonds, central_atom)` → parameters become integer indices
- The `"multiplicity"` dict in `calculation_specification()` changes to `dict[tuple[int,int], int]`
- `"equivalent_atoms"` value changes to `list[list[int]]` (pymatgen already returns position-based sets)

### 1d. `nomodeco/libraries/icsel.py`

- `remove_enumeration()` / `remove_enumeration_tuple()` — no callers found; can be removed or left as no-ops
- `all_atoms_can_be_superimposed_*()` functions compare tuple elements with `==`; string→int comparison is faster, no structural change
- `number_terminal_bonds(mult_list)` — `mult_list` element type changes from `tuple[str, int]` to `tuple[int, int]`; the logic `atom_and_mult[1] == 1` is unchanged

### 1e. `nomodeco/libraries/topology.py` (minimal changes only)

Only 4 types of changes needed across this 4465-line file:

1. **Lines ~2648, ~2964, ~3642, ~4001** — `removed_bond[0].strip(string.digits) == "H"` pattern. Replace with `atoms_list[removed_bond[0]].symbol.strip(string.digits) == "H"` using the already-available `atoms_list`.

2. **Line ~1482** (`extract_atoms_of_submolecules`) — `atom_obj.symbol in component` where `component` is now a set of integers from `nx.connected_components`. Replace with `enumerate`-based loop checking `i in component`.

3. **`get_multiplicity()` callers** — now receives integer atom index; `if atom_name == atom_and_mult[0]` becomes int==int, which is correct.

4. All `nx.Graph.add_edges_from(bonds)` calls work identically with integer-pair tuples — no change needed.

### 1f. `nomodeco/nomodeco.py`

- **Cache repeated generate calls** (lines 927–968): `generate_angles(cov_bonds)` and `generate_dihedrals(cov_bonds)` are called multiple times on the same bond sets. Cache results before the `add_coordinate` / `add_coord_diff` calls.

- **Display conversion** (line ~1491): `", ".join(internal)` fails on integer tuples. Add:
  ```python
  symbol_list = atoms.index_to_symbol()
  def ic_to_label(ic): return "(" + ", ".join(symbol_list[i] for i in ic) + ")"
  all_internals_string = [ic_to_label(ic) for ic in all_internals]
  ```

- **Backward compat for `--graph_ic` / `--nomodeco_coords` file reading**: After parsing IC tuples from text files, check `isinstance(first_element, str)` and if so apply `atoms.symbol_to_index()` to convert. Old files continue to work transparently.

### 1g. `nomodeco/libraries/icset_opt.py` and `icset_opt_multproc.py`

No structural changes needed — these consume IC tuples only as opaque keys for counting and dict lookups, which work identically on integer tuples.

---

## Phase 2: Memory Tracking

### 2a. New file: `nomodeco/libraries/memprofile.py`

```python
import tracemalloc

class MemoryProfiler:
    def __init__(self, enabled: bool): ...
    def start(self): ...
    def checkpoint(self, label: str): ...  # logs current + peak MiB
    def stop(self): ...
    def report(self): ...  # prints summary table
```

Standard library only — no new dependencies.

### 2b. `nomodeco/libraries/arguments.py`

Add `--memory_profile` boolean flag (store_true).

### 2c. Checkpoint insertion in `nomodeco/nomodeco.py`

```
mem = MemoryProfiler(enabled=args.memory_profile)
mem.start()
```

Checkpoints at:
1. After molecule construction
2. After all IC generation (`Total_IC_dict` fully populated)
3. After B-matrix construction
4. After G-matrix construction
5. Just before P/T/E tensor allocation (shows pre-tensor baseline)
6. After P/T/E fill loops
7. After `ved_matrix` and `contribution_matrix` computed

`mem.report()` at the end of `main()`.

### 2d. `nomodeco/libraries/icset_opt.py`

Pass `mem` as optional parameter to `find_optimal_coordinate_set()`. Add checkpoints before/after P/T/E allocation per IC set candidate.

### 2e. `nomodeco/libraries/icset_opt_multproc.py`

Within each worker `process_ic_set()`, use a local `tracemalloc` start/stop and return peak memory as an extra return value. Aggregate in the orchestrator.

---

## Phase 3 (Bonus): Vectorize O(n³) Loops

Replace the triple-nested loops in `nomodeco.py` (lines 1406–1475), `icset_opt.py`, and `icset_opt_multproc.py` with numpy broadcasting. The same pattern applies in all three files:

```python
# current (O(n³) Python loop):
for i in range(R):
    for m in range(N):
        for n in range(N):
            k = i + num_rottra
            P[i,m,n] = D[m,k] * InternalF_Matrix[m,n] * D[n,k] / eigenvalues[k]

# vectorized:
k_idx = np.arange(num_rottra, num_rottra + R)
D_k = D[:, k_idx]                                  # shape (N, R)
D_outer = D_k.T[:, :, np.newaxis] * D_k.T[:, np.newaxis, :]  # (R, N, N)
P = D_outer * InternalF_Matrix[np.newaxis] / eigenvalues[k_idx, np.newaxis, np.newaxis]
T = D_outer * G_inv[np.newaxis]
E = 0.5 * (T + P)
```

Additional simplifications:
- `ved_matrix = P.sum(axis=2)` replaces the triple loop for VED (lines 1431–1435)
- `Diag_elements = P[:, :n_internals, :n_internals].diagonal(axis1=1, axis2=2)` replaces the `np.diag(P[i])` loop (lines 1448–1450)
- `nu = np.einsum('mi,mn,ni->n', D_k, InternalF_Matrix[:n_internals,:n_internals], D_k)` replaces the triple loop for intrinsic frequencies (lines 1469–1473)

After Phase 3 is verified, extract the P/T/E block into a shared `nomodeco/libraries/ped_calc.py` to eliminate code triplication across the three files.

---

## Verification

1. **Regression test**: Run on a known molecule (any file in `test_calculations/`). Compare `ved_matrix`, `contribution_matrix`, and `nu_final` values between old and new code. They must be numerically identical (within floating-point tolerance).

2. **Memory tracking**: Run with `--memory_profile` and confirm the profiling section prints MiB values at each checkpoint. Run without the flag and confirm zero overhead.

3. **Vectorization correctness**: The normalization sums (`sum_check_PED`, `sum_check_KED`, `sum_check_TED`) already exist in the code — confirm they still pass after Phase 3.

4. **Backward compat**: Run with an old-format `--graph_ic` file (string tuples) and confirm the auto-conversion shim kicks in without errors.

---

## Implementation Order

1. Phase 1 — data structure migration (highest impact on memory, touches the most files)
2. Phase 2 — memory tracking (independent of Phase 1, can be done in parallel)
3. Phase 3 — loop vectorization (depends on Phase 1 being stable)