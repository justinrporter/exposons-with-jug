import jug
from jug.io import NoLoad

from . import tasks
from . import util


def featurize(
        tag, trajectories, topology, probe_radius=2.8):
    """Run the whole pockets pipeline
    """

    probe_radius = '%.2f' % probe_radius

    DATA_STEM = util.data_stem(tag, probe_radius)

    sasas = [
        tasks.map_sasa_sparse(
            trj, topology, float(probe_radius) / 10,
            out=DATA_STEM + 'atomic-%06d.sparse.npz' % i)
        for i, trj in enumerate(trajectories)]

    sidechain_sasas = [tasks.condense_sparse_sidechain_sasas(sasa, topology)
                       for sasa in sasas]

    SC_SASA_FILE = DATA_STEM + 'sasas-sidechains.h5'
    sasa_sidechain_h5 = tasks.assemble_sasa_h5(
        sasas=[NoLoad(r) for r in sidechain_sasas],
        filename=SC_SASA_FILE)

    return sasa_sidechain_h5


def cluster(
        tag, trajectories, topology, sasa_sidechain_h5, cluster_radii,
        cluster_algorithm='khybrid', cluster_distance='euclidean',
        kmedoids_updates=5):

    import os
    from collections import namedtuple

    sc_sasa_filename = jug.bvalue(sasa_sidechain_h5)

    def make_cluster_name(sc_sasa_filename, suffix):
        f = sc_sasa_filename \
            .rstrip('sasas-sidechains.h5') \
            .replace('/features', '/cluster/') + \
            "-".join(['sidechains', cluster_algorithm, radius,
                      cluster_distance] + (
                          ['kmedoids' + str(kmedoids_updates)]
                          if cluster_algorithm == 'khybrid'
                          else []))
        return f

    ClusterFiles = namedtuple(
        'ClusterFiles',
        ['assignments', 'distances', 'center_indices',
         'structure_centers', 'feature_centers'])

    cluster_results = []
    for radius in cluster_radii:

        CLUSTER_STEM = sc_sasa_filename \
            .replace('sasas-sidechains.h5', '') \
            .replace('/features/', '/cluster/') + \
            "-".join(['sidechains', cluster_algorithm, radius,
                      cluster_distance] + (
                          ['kmedoids' + str(kmedoids_updates)]
                          if cluster_algorithm == 'khybrid'
                          else []))

        ASSIGNMENTS_FILE = CLUSTER_STEM + '-assignments.h5'
        DISTANCES_FILE = CLUSTER_STEM + '-distances.h5'
        FEATURE_CENTERS_FILE = CLUSTER_STEM + '-feature-centers.npy'
        CENTER_INDS_FILE = CLUSTER_STEM + '-center-inds.npy'

        assignments, distances, center_features, center_indices = \
            jug.iteratetask(tasks.cluster_features(
                sasa_sidechain_h5,
                ASSIGNMENTS_FILE,
                DISTANCES_FILE,
                FEATURE_CENTERS_FILE,
                CENTER_INDS_FILE,
                radius,
                'euclidean',
                cluster_algorithm,
                cluster_iterations=kmedoids_updates
            ), n=4)

        ctr_structs = tasks.write_struct_ctrs(
            trajectories, topology, center_indices,
            CLUSTER_STEM + "-structure-centers.h5")

        result = ClusterFiles(
            assignments=assignments,
            distances=distances,
            center_indices=center_indices,
            feature_centers=center_features,
            structure_centers=ctr_structs)
        cluster_results.append(result)

        PLOT_PATH = os.path.join(
            'figures', 'implied',
            os.path.basename(CLUSTER_STEM) + 'implied-timescales.png')

        tasks.implied_timescales(
            assignments=assignments,
            plot_path=PLOT_PATH
        )

    return cluster_results
