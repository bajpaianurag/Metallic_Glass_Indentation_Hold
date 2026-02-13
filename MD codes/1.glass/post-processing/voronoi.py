# voronoi_ovito_radical.py
import os
import numpy as np
from collections import defaultdict, Counter

from ovito.io import import_file
from ovito.modifiers import VoronoiAnalysisModifier, PythonScriptModifier

from ase import Atoms
from ase.io import write

import logging

# Define a desired radius (Å) mapping per element. (Examples use metallic radii.)
DEFAULT_RADII = {
    "Cu": 1.28,
    "Al": 1.43,
    "Ti": 1.47,
    "Zr": 1.60,
    "Ag": 1.44,
}

# Configure logging to capture a reproducible record of the analysis.
logging.basicConfig(
    filename="voronoi_analysis.log",   # Output log file
    filemode="w",                      # 'a' appends; 'w' overwrites
    level=logging.INFO,                 # Minimum severity to record
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logging.info("Voronoi analysis started")

# Fallback radius used for species not present in the mapping.
DEFAULT_RADIUS_FALLBACK = 1.35


def add_radii_property_modifier(radii_map=None, default_radius=DEFAULT_RADIUS_FALLBACK):
    """
    Create a PythonScriptModifier that sets/updates per-particle 'Radius'.
    These radii are consumed by Radical Voronoi (use_radii=True).
    """
    if radii_map is None:
        radii_map = DEFAULT_RADII.copy()

    def _set_radii(frame, data):
        # Resolve each particle's chemical name from its 'Particle Type',
        # then assign the corresponding radius (with fallback).
        types = data.particles.particle_types
        type_ids = data.particles["Particle Type"].array
        radii = np.empty(len(type_ids), dtype=float)
        for i, tid in enumerate(type_ids):
            name = types.type_by_id(int(tid)).name
            radii[i] = radii_map.get(name, default_radius)
        # Create or overwrite the 'Radius' property used by Radical Voronoi.
        data.particles_.create_property("Radius", data=radii)

    return PythonScriptModifier(function=_set_radii)


def load_with_ovito(filename, radii_map=None):
    """
    Build and execute the OVITO pipeline:
      1) Inject the per-particle Radius property.
      2) Run VoronoiAnalysisModifier with use_radii=True (Radical Voronoi).
    Returns the evaluated DataCollection.
    """
    pipeline = import_file(filename)

    # Inject radii first so the Voronoi modifier can use them.
    pipeline.modifiers.append(add_radii_property_modifier(radii_map=radii_map))

    # Radical Voronoi (accounts for particle radii when tessellating).
    pipeline.modifiers.append(
        VoronoiAnalysisModifier(
            compute_indices=True,
            generate_bonds=True,
            use_radii=True,         # Key switch for Radical Voronoi
            # Optional cleanup if needed:
            # edge_threshold=0.0,
            # face_threshold=0.0,
        )
    )
    data = pipeline.compute()
    return data


def get_species_names_per_particle(data):
    """Return a list of element names for each particle."""
    types = data.particles.particle_types
    type_ids = data.particles["Particle Type"].array
    return [types.type_by_id(int(tid)).name for tid in type_ids]


def voronoi_index_string(voro_vec):
    """
    Convert OVITO's Voronoi index array into a compact "<n3,n4,n5,...>" string,
    trimming trailing zeros for readability.
    """
    arr = np.asarray(voro_vec, dtype=int)
    arr = arr[2:]  # Drop coordination number and number of faces if present at front
    if arr.size == 0:
        return "<>"
    nz = np.nonzero(arr)[0]
    trimmed = [] if nz.size == 0 else arr[: nz[-1] + 1].tolist()
    return "<" + ",".join(str(x) for x in trimmed) + ">"


def bonds_topology_neighbors_of(data, i):
    """Return unique bonded neighbor indices of particle i from the bonds topology."""
    topo = data.particles.bonds["Topology"].array
    neigh = []
    mask0 = topo[:, 0] == i
    mask1 = topo[:, 1] == i
    if np.any(mask0):
        neigh.extend(topo[mask0, 1].tolist())
    if np.any(mask1):
        neigh.extend(topo[mask1, 0].tolist())
    return sorted(set(int(x) for x in neigh))


def is_well_inside_by_frac(data, i, frac_margin=0.1):
    """
    Heuristic: Check if particle i is well inside the simulation cell by requiring
    its fractional coordinates to be at least `frac_margin` away from each face.
    Helps avoid exporting clusters that wrap across periodic boundaries.
    """
    pos = data.particles["Position"][i]
    cell = data.cell
    a, b, c = cell[:, 0], cell[:, 1], cell[:, 2]
    origin = cell[:, 3]
    M = np.column_stack((a, b, c))
    frac = np.linalg.solve(M, pos - origin)
    return np.all((frac_margin < frac) & (frac < (1.0 - frac_margin)))


def export_local_cluster_xyz(data, center_id, neighbor_ids, species_names, filename):
    """Export a local cluster (center + bonded neighbors) to an .xyz file using ASE."""
    pos_all = data.particles["Position"].array
    cell = data.cell
    cell3 = np.column_stack((cell[:, 0], cell[:, 1], cell[:, 2]) )
    ids = [center_id] + list(neighbor_ids)
    symbols = [species_names[j] for j in ids]
    positions = [pos_all[j] for j in ids]
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell3, pbc=True)
    write(filename, atoms)
    logging.info(f"✔ Saved: {filename}")


def analyze_voronoi_by_species_ovito(data, target_species):
    """
    For a chosen element symbol (e.g., 'Cu'), compute the distribution of
    Voronoi index strings among atoms of that species and collect example atoms
    and their bonded neighbors for each index.
    Returns (Counter of indices, mapping index->list of (center_id, neighbor_ids)).
    """
    species_names = get_species_names_per_particle(data)
    voro_prop = data.particles["Voronoi Index"].array
    index_map = defaultdict(list)
    index_counter = Counter()
    N = data.particles.count
    for i in range(N):
        if species_names[i] != target_species:
            continue
        idx_str = voronoi_index_string(voro_prop[i])
        neighbor_ids = bonds_topology_neighbors_of(data, i)
        index_map[idx_str].append((i, neighbor_ids))
        index_counter[idx_str] += 1
    return index_counter, index_map


# --- Main ---
if __name__ == "__main__":
    # Input structure file readable by OVITO (e.g., LAMMPS .cfg).
    filename = "final.cfg"

    # Optional: customize radii per species for Radical Voronoi.
    custom_radii = {
        "Cu": 1.28,
        "Al": 1.43,
        "Ti": 1.47,
        "Zr": 1.60,
        "Ag": 1.44,
    }
    data = load_with_ovito(filename, radii_map=custom_radii)

    # Prepare output directory for local cluster exports.
    os.makedirs("voronoi_clusters", exist_ok=True)
    species_names = get_species_names_per_particle(data)

    # Analyze selected center species and log distributions and example clusters.
    for element in ["Cu", "Zr"]:
        logging.info(f"\n===== {element}-centered Voronoi index distribution =====")
        index_counter, index_map = analyze_voronoi_by_species_ovito(data, element)
        total = sum(index_counter.values()) or 1
        for index_str, count in index_counter.most_common():
            percent = 100.0 * count / total
            logging.info(f"{index_str:15s}: {count:5d} atoms ({percent:.2f}%)")

        logging.info(f"\n--- Saving clusters for top 10 indices ({element}) ---")
        top10 = index_counter.most_common(10)
        for _, (index_str, _) in enumerate(top10, start=1):
            examples = index_map[index_str]
            for atom_id, neighbors in examples:
                if is_well_inside_by_frac(data, atom_id):
                    outname = (
                        f"voronoi_clusters/{element}_voro_"
                        f"{index_str.replace('<','').replace('>','').replace(',','-')}"
                        ".xyz"
                    )
                    export_local_cluster_xyz(data, atom_id, neighbors, species_names, outname)
                    break  # Save a representative example and move to next index

