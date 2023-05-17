import calculate
import numpy as np

if __name__ == "__main__":
    with open("rpp_xyz.txt", "r") as f:
        xyz = eval(f.read())

    with open("rpp_zmat.txt", "r") as f:
        zmat = eval(f.read())

    print("xyz=")
    # for i in xyz:
    #     print(xyz[i])
    # print(xyz)
    for i in range(1, max(xyz.keys())+1):
        print(xyz[i])

    print("zmat=")
    # print(zmat)
    for i in zmat:
        print(zmat[i])

    print()
    print("--------------------")

    out = calculate.zmat2xyz(zmat)
    # print(out)
    print("转换结果：")
    for i in range(1, max(out.keys())+1):
        print(out[i])




