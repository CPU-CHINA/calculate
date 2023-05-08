"""
此文件是计算 Lennard-Jones potential 的纯 Python 实现
"""

import numpy as np


def get_distance(a, b):  # a,b should be a tuple or list with froms (x_axe,y_axe,z_axe)
    a, b = map(lambda x: np.array(x), [a, b])
    return np.linalg.norm(a - b)


def ca_need_atom(tree_dict: dict, atom_set: set, atom: int):
    for angle in tree_dict["[angles]"]:
        '''
        angle=
        [6, 15]
        [5, 6, 15]
        []
        '''
        if atom in angle:
            angle.remove(atom)
            atom_set -= set(angle)
        atom_set -= {atom}
    return atom_set


def get_mols_LJ(mols: list, forest: dict, mole_loc: dict, LJ_cutoff: bool, CC_lib: dict):
    e_LJ = 0.0
    for mol in mols:
        atom_set = set(map(int, forest[mol]["[atoms]"].keys()))
        # atom_set={1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20}
        for atom1 in forest[mol]["[atoms]"]:
            # atom1 = 1,2,3,4,5,6,7... 是key
            atom_set = ca_need_atom(forest[mol], atom_set, atom1)
            for atom2 in atom_set:
                dis = get_distance(mole_loc[mol][atom1 - 1][-1], mole_loc[mol][atom2 - 1][-1])
                # print(f"A={mole_loc[mol][atom1 - 1][-1]}")
                # print(f"B={mole_loc[mol][atom2 - 1][-1]}")
                if LJ_cutoff and dis >= 1.4:
                    continue
                else:
                    atomtype1 = str(forest[mol]["[atoms]"][atom1][1])
                    atomtype2 = str(forest[mol]["[atoms]"][atom2][1])
                    c12 = CC_lib[f'{atomtype1}_{atomtype2}'][0]
                    c6 = CC_lib[f'{atomtype1}_{atomtype2}'][1]
                    e_LJ += c12 / (dis ** 12) - c6 / (dis ** 6)
    return e_LJ


if __name__ == "__main__":
    # 读取
    with open("c.mols.txt", "r") as f:
        mols = eval(f.read().strip("dict_keys(").strip(")"))

    with open("c.forest.txt", "r") as f:
        forest = eval(f.read())

    with open("c.mole_loc.txt", "r") as f:
        mole_loc = eval(f.read())

    LJ_cutoff = True

    import pickle

    with open("CC_lib.pkl", "rb") as f:
        CC_lib = pickle.load(f)

    import time

    start = time.time()
    print(get_mols_LJ(mols, forest, mole_loc, LJ_cutoff, CC_lib))
    end = time.time()
    print(f'Time spend in Python: {end - start}s')
