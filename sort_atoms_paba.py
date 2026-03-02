import numpy as np

def get_dist(p1, p2):
    return np.linalg.norm(p1 - p2)

def sort_paba_custom(input_xyz, output_xyz, cutoff=1.7):
    with open(input_xyz, 'r') as f:
        lines = f.readlines()

    if not lines: return

    num_atoms_str = lines[0].strip()
    lattice_line = lines[1].strip()
    atom_data = lines[2:]

    clusters = {}
    for line in atom_data:
        parts = line.split()
        if len(parts) < 6: continue
        c_id = int(parts[0])
        species = parts[2]
        pos = np.array([float(x) for x in parts[3:6]])
        clusters.setdefault(c_id, []).append({'species': species, 'pos': pos})

    final_output = []
    global_id = 1

    for c_id in sorted(clusters.keys()):
        atoms = clusters[c_id]
        for i, a in enumerate(atoms): a['local_id'] = i

        # 1. Identify Key Atoms
        oxygens = [a for a in atoms if a['species'] == 'O']
        hydrogens = [a for a in atoms if a['species'] == 'H']
        carbons = [a for a in atoms if a['species'] == 'C']
        n_amine = next(a for a in atoms if a['species'] == 'N')

        o_oh = next(o for o in oxygens if any(get_dist(o['pos'], h['pos']) < 1.2 for h in hydrogens))
        o_co = next(o for o in oxygens if o['local_id'] != o_oh['local_id'])
        h_hydroxyl = next(h for h in hydrogens if get_dist(h['pos'], o_oh['pos']) < 1.2)

        c_carboxyl = next(a for a in carbons if sum(1 for o in oxygens if get_dist(a['pos'], o['pos']) < cutoff) == 2)
        ring_carbons = [a for a in carbons if a['local_id'] != c_carboxyl['local_id']]
        ring_center = np.mean([a['pos'] for a in ring_carbons], axis=0)

        # 2. Define the Viewing Plane (Up = N, Bottom = COOH, Right = H_hydroxyl side)
        # Vertical Axis (Y)
        v_y = n_amine['pos'] - c_carboxyl['pos']
        v_y /= np.linalg.norm(v_y)

        # Temp vector to define the "right" side
        v_to_h = h_hydroxyl['pos'] - c_carboxyl['pos']

        # Normal to the plane (Z)
        v_z = np.cross(v_y, v_to_h)
        v_z /= np.linalg.norm(v_z)

        # Horizontal Axis (X) - pointing right toward the H side
        v_x = np.cross(v_z, v_y)

        # 3. Find the "Starting Atom" for the sweep
        # Requirement: Atom in ring/amine closest to the Hydroxyl Hydrogen
        candidates = ring_carbons + [n_amine] + [h for h in hydrogens if h['local_id'] != h_hydroxyl['local_id']]
        start_atom = min(candidates, key=lambda a: get_dist(a['pos'], h_hydroxyl['pos']))

        # 4. Sorting Function: Counter-Clockwise in the XY plane
        def get_ccw_angle(pos):
            rel = pos - ring_center
            # Project 3D vector onto our 2D local X and Y axes
            proj_x = np.dot(rel, v_x)
            proj_y = np.dot(rel, v_y)
            # Calculate angle in the 2D plane
            angle = np.arctan2(proj_y, proj_x)

            # Adjust so start_atom is the 0-point for the sweep
            start_rel = start_atom['pos'] - ring_center
            offset = np.arctan2(np.dot(start_rel, v_y), np.dot(start_rel, v_x))

            final_angle = angle - offset
            if final_angle < 0: final_angle += 2 * np.pi
            return final_angle

        # 5. Apply the Sorting to Subgroups
        h_ring_amine = [h for h in hydrogens if h['local_id'] != h_hydroxyl['local_id']]
        h_ring_amine.sort(key=lambda h: get_ccw_angle(h['pos']))

        ring_carbons.sort(key=lambda c: get_ccw_angle(c['pos']))

        # 6. ASSEMBLE FINAL SEQUENCE
        ordered_list = [o_oh, o_co, n_amine, h_hydroxyl] + h_ring_amine + [c_carboxyl] + ring_carbons

        for a in ordered_list:
            final_output.append(f"{c_id} {global_id} {a['species']} {a['pos'][0]:.8f} {a['pos'][1]:.8f} {a['pos'][2]:.8f}\n")
            global_id += 1

    with open(output_xyz, 'w') as f:
        f.write(f"{num_atoms_str}\n{lattice_line}\n")
        f.writelines(final_output)
    print(f"Sorted {len(clusters)} molecules with CCW sweep (H-hydroxyl on right).")

if __name__ == "__main__":
    sort_paba_custom("gamma_16mol_frame36.xyz", "topology_sorted.xyz")
