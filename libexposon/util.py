import os


def enspara_path(path):
    import enspara
    return os.path.join(enspara.__path__[0],
                        path)


def set_ext(p, ext):
    import os
    path, ext = os.path.splitext(p)
    return path + ext if ext[0] == '.' else '.' + ext


def data_stem(tag, probe_radius):
    return f'data/{tag}/features/{tag}-{probe_radius}A-sub1-'


def subset_trj(xtc, top, selection):
    import mdtraj as md
    trj = md.load(xtc, top)
    return trj.atom_slice(trj.top.select(selection))


def load_assignments(assignments):

    from enspara.util import array as ra
    from tables import NoSuchNodeError

    if not hasattr(assignments, 'shape'):
        print('loading msm assignments from', assignments)
        try:
            assignments = ra.load(assignments, keys=None)
        except NoSuchNodeError:
            assignments = ra.load(assignments, keys=...)

    return assignments
