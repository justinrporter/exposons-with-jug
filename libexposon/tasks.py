import jug

from . import util


def map_sasa_core(trajectoryfile, topologyfile, probe_radius):
    import mdtraj as md

    print(
        "computing using threading ", probe_radius, "nm sasa for",
        trajectoryfile, "using topology", topologyfile)

    if topologyfile:
        trj = md.load(trajectoryfile, top=topologyfile)
    else:
        trj = md.load(trajectoryfile)
    print("loaded", trajectoryfile, 'with topology', topologyfile)
    sasas = md.shrake_rupley(trj, probe_radius=probe_radius)

    return sasas


@jug.TaskGenerator
def map_sasa_sparse(trajectoryfile, topologyfile, probe_radius, out):

    import os
    from scipy import sparse

    if not os.path.exists(os.path.dirname(out)):
        os.makedirs(os.path.dirname(out), exist_ok=True)
    # assert os.path.exists(os.path.dirname(out)), \
    #     "Directory doesn't exist for output %s" % out

    sasas = map_sasa_core(trajectoryfile, topologyfile, probe_radius)

    sparse.save_npz(
        file=out,
        matrix=sparse.csr_matrix(sasas)
    )

    return out


def condense_sidechain_sasas_core(sasas, top):
    from tqdm import tqdm
    import time
    import numpy as np

    assert top.n_atoms == sasas.shape[1], '%s != %s' % (top.n_atoms,
                                                        sasas.shape[1])

    SELECTION = ('not (name N or name C or name CA or name O or '
                 'name HA or name H or name H1 or name H2 or name '
                 'H3 or name OXT)')

    sc_ids = [top.select('resid %s and ( %s )' % (i, SELECTION))
              for i in range(top.n_residues)]

    rsd_sasas = np.zeros((sasas.shape[0], len(sc_ids)), dtype='float32')

    for i in tqdm(range(len(sc_ids))):
        try:
            # squeeze is necessary because scipy matrices return (n, 1)
            rsd_sasas[:, i] = sasas[:, sc_ids[i]].sum(axis=1).squeeze()
        except:
            print('condensing residue', i, 'of', top.n_residues)
            print(sc_ids[i])
            print(sasas.shape)
            print(rsd_sasas.shape)
            raise

    return rsd_sasas


@jug.TaskGenerator
def condense_sparse_sidechain_sasas(sasas_file, topology):
    import time
    import mdtraj as md
    from scipy import sparse

    print(topology)
    top = md.load(topology).top

    t0 = time.perf_counter()
    sasas = sparse.load_npz(sasas_file)
    print("loading sparse array took", time.perf_counter() - t0)

    return condense_sidechain_sasas_core(sasas, top)


@jug.TaskGenerator
def assemble_sasa_h5(sasas, filename):

    import os
    import tables
    from tqdm import tqdm

    if not os.path.isdir(os.path.dirname(filename)):
        os.mkdir(os.path.dirname(filename))

    if os.path.isfile(filename):
        raise FileExistsError(f"File '{filename}' already exists.")

    compression = tables.Filters(complevel=9, complib='zlib', shuffle=True)
    n_zeros = len(str(len(sasas))) + 1

    print(filename)
    with tables.open_file(filename, 'a') as handle:
        shape = None

        for i, sasa in enumerate(tqdm(sasas)):

            data = jug.value(sasa.t)

            atom = tables.Atom.from_dtype(data.dtype)
            tag = 'sasas_' + str(i).zfill(n_zeros)

            if tag in handle.root:
                logger.warn(
                    'Tag %s already existed in %s. Overwriting.',
                    tag, filename)
                handle.remove_node('/', name=tag)

            if shape is None:
                shape = data.shape
            elif len(shape) > 1:
                assert shape[1] == data.shape[1], "We had %s residues, but then loaded trajectory %s and it had %s." % (shape[1], i, data.shape[1])

            node = handle.create_carray(
                where='/', name=tag, atom=atom,
                shape=data.shape, filters=compression)
            node[:] = data

            sasa.t.unload()

    return filename


@jug.TaskGenerator
def cluster_features(features, assignments, distances, center_features,
                     center_indices, cluster_radius,
                     cluster_distance='euclidean',
                     algorithm='khybrid',
                     cluster_iterations=5):
    import os
    import enspara
    import subprocess

    center_features = util.set_ext(center_features, '.npy')
    cluster_executable = util.enspara_path('apps/cluster.py')

    for f in [assignments, center_indices, distances, center_features]:
        os.makedirs(os.path.dirname(f), exist_ok=True)

    args = [
        'mpiexec', 'python', '-u',
        cluster_executable,
        '--algorithm', algorithm,
        '--features', features,
        '--assignments', assignments,
        '--distances', distances,
        '--center-features', center_features,
        '--center-indices', center_indices,
        '--cluster-distance', cluster_distance,
        '--cluster-radius', cluster_radius,
    ]

    if algorithm == 'khybrid':
        args += ['--cluster-iterations', str(cluster_iterations)]

    print(" ".join(args))
    subprocess.run(args, check=True)

    return assignments, distances, center_features, center_indices


@jug.TaskGenerator
def write_struct_ctrs(trajectoryfiles, topology, ctr_inds_file,
                      ctr_structs_file, run_after=None):

    import os
    import pickle
    import mdtraj as md
    from enspara.util.load import load_as_concatenated

    print("Loading center indices at", ctr_inds_file)

    if os.path.isfile(ctr_structs_file):
        print("Overwriting", ctr_structs_file)
    #     print("Refusing to overwrite", ctr_structs_file)
    #     return ctr_structs_file

    try:
        with open(ctr_inds_file, 'rb') as f:
            ctr_inds = pickle.load(f)
    except pickle.UnpicklingError:
        import numpy as np
        ctr_inds = np.load(ctr_inds_file)

    print("Loaded %s center indices." % len(ctr_inds))
    top = md.load(topology).top

    try:
        lengths, xyz = load_as_concatenated(
            filenames=[trajectoryfiles[tr] for tr, fr in ctr_inds],
            args=[{'frame': fr, 'top': top} for tr, fr in ctr_inds],
            processes=4
        )
    except IndexError:
        print(len(trajectoryfiles), len(ctr_inds),
              max([tr for tr, fr in ctr_inds]))
        raise

    ctr_structs = md.Trajectory(xyz=xyz, topology=top)
    ctr_structs.save(ctr_structs_file)

    return ctr_structs_file


@jug.TaskGenerator
def implied_timescales(assignments, plot_path, method='prior_counts',
                       lag_times='5:200:2'):
    import os
    import subprocess
    timescales_app_location = util.enspara_path('apps/implied_timescales.py')

    os.makedirs(os.path.dirname(plot_path), exist_ok=True)

    subprocess.run(
        ['python', '-u',
         timescales_app_location,
         '--assignments', assignments,
         '--lag-times', lag_times,
         '--symmetrization', method,
         '--plot', plot_path,
         '--logscale'],
        check=True)

    return plot_path


def msm_core(assignments, lag_time, method='prior_counts', **kwargs):

    from enspara import msm
    from functools import partial

    max_n_states = kwargs.get('max_n_states', assignments.max())

    if method == 'prior_counts':
        method = partial(
            msm.builders.normalize,
            prior_counts=1 / max_n_states)
    else:
        method = getattr(msm.builders, method)

    mkv = msm.MSM(lag_time=lag_time, method=method)

    assert hasattr(assignments, 'shape')
    mkv.fit(assignments)

    return mkv


@jug.TaskGenerator
def msm2file(filename, assignments, lag_time, **kwargs):

    import os
    import shutil

    os.makedirs(os.path.dirname(filename), exist_ok=True)

    assignments = util.load_assignments(assignments)

    msm = msm_core(assignments, lag_time, **kwargs)

    if os.path.isdir(filename) or os.path.isfile(filename):
        print(f"Path {filename} already existed. Moving "
              f"it to {filename}.bk before building an MSM there.")
        shutil.move(filename, filename + ".bk")

    msm.save(filename)

    return filename
