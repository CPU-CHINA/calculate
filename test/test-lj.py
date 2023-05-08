"""
此文件是计算 Lennard-Jones potential 的 Rust 实现调用示例
"""

import calculate

if __name__ == "__main__":
    # 读取
    with open("c.mols.txt", "r") as f:
        mols = eval(f.read().strip("dict_keys(").strip(")"))

    with open("c.forest.txt", "r") as f:
        forest = eval(f.read())
        # No need in LJ calculate
        # forest["complex"] = forest.pop("complex_rpp_dop")

    with open("c.mole_loc.txt", "r") as f:
        mole_loc = eval(f.read())

    LJ_cutoff = True

    import pickle

    with open("CC_lib.pkl", "rb") as f:
        CC_lib = pickle.load(f)

    CalculateLjPotential = calculate.CalculateLjPotential(CC_lib)
    import time

    start = time.time()
    print(CalculateLjPotential.calculate(mols, forest, mole_loc, LJ_cutoff))
    end = time.time()
    print(f'Time spend in Rust'
          f': {end - start}s')
