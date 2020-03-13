import os
import json

import jug

from lib.pipeline import featurize, cluster
from lib.tasks import msm2file


with open(os.path.join(os.path.dirname(__file__),
                       'proteins.json'), 'r') as f:
    CONFIGS = json.load(f)

prior_count = 1

cluster_results = {}

for protein, cfg in CONFIGS.items():
    assert cfg['trajectories'], "No trjs for %s@%s" % protein

    sc_sasa_h5 = featurize(
        tag=protein,
        trajectories=cfg['trajectories'],
        topology=cfg['topology'],
    )

    cluster_results[protein] = cluster(
        tag=protein,
        trajectories=cfg['trajectories'],
        topology=cfg['topology'],
        sasa_sidechain_h5=sc_sasa_h5,
        cluster_distance=cfg['cluster_distance_metric'],
        cluster_radii=cfg['cluster_radii'],
        kmedoids_updates=10
    )


selected_cluster_results = {}
for protein, cluster_result_list in cluster_results.items():
    if "lag_time" not in CONFIGS[protein]:
        continue

    for radius, cluster_result in zip(cfg['cluster_radii'],
                                      cluster_result_list):
        if CONFIGS[protein]["model_cluster_radius"] == radius:
            selected_cluster_results[protein] = cluster_result
            lag_time = CONFIGS[protein]["lag_time"]

            if cluster_result.assignments.can_load():
                dirname, assigs_file = os.path.split(
                    jug.bvalue(cluster_result.assignments))

                dirname = os.path.join(os.path.split(dirname)[0], 'models')
                msm_filename = assigs_file.replace('-assignments.h5',
                                                   '-%02dprior-%slt-msm' %
                                                   (prior_count, lag_time))

                assignments = cluster_result.assignments
                msm2file(os.path.join(dirname, msm_filename),
                         assignments,
                         lag_time=lag_time)
