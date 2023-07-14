import random
from datetime import datetime
import snappy
import unittest

import spherogram
from tqdm import tqdm
from multiprocessing import Pool


# Names "3_1(0,0)", "5_1(0,0)", "8_19(0,0)", "9_1(0,0)", "K13a4726(0,0)", "K13a4878(0,0)"
ERROR_MANIFOLDS = [1, 4, 35, 76, 7703, 7855, 12442]


def is_symplectic(M):
    """
    Test if the matrix M is symplectic
    :param M: square matrix
    :return: true or false
    """
    n = len(M)

    for i in range(n):
        for j in range(i, n):
            omega = abs(symplectic_form(M[i], M[j]))

            if i % 2 == 0 and j % 2 == 1 and j == i + 1:
                if omega != 2:
                    return False
            elif omega:
                return False

    return True


def symplectic_form(u, v):
    return sum([u[2 * i] * v[2 * i + 1] - u[2 * i + 1] * v[2 * i] for i in range(len(u) // 2)])


def process_manifold(i: int):
    M = snappy.HTLinkExteriors[i]
    label = M.identify()

    if i in ERROR_MANIFOLDS:
        return True

    with open("logs/links-" + str(i // 1000) + ".log", "a") as file:
        test = "Testing: " + str(i) + ' ' + str(label)
        file.write(test)

        basis = M.symplectic_basis()
        result = is_symplectic(basis.data)

        if result:
            string = "Passed\n"
        else:
            string = "Failed\n"
        file.write((60 - len(test)) * " " + str(string))

    return result


def find_manifold(start: int):
    scale = 1000
    for i in range(scale * start, scale * (start + 1)):
        M = snappy.HTLinkExteriors[i]

        if len(M.identify()) == 0:
            return True


def test_link_complements_pool(start: int, end: int):
    scale = 1000
    for i in range(start, end):
        passed = True

        with open("logs/total.log", "a") as file:
            file.write("[" + datetime.now().strftime("%d-%m-%y %H:%M:%S") + "]   " + "Testing: " + str(scale * i) + " - " + str(scale * (i + 1) - 1) + "\n")

        with Pool() as pool:
            result = pool.imap(process_manifold, range(scale * i, scale * (i + 1)))

            for j, res in enumerate(result):
                if not res:
                    passed = False
                    break

        with open("logs/total.log", "a") as file:
            time = "[" + datetime.now().strftime("%d-%m-%y %H:%M:%S") + "]   "
            if passed:
                file.write(time + "Passed\n")
            else:
                file.write(time + "Failed\n")


class TestSymplecticBasis(unittest.TestCase):
    def test_knot_complements(self):
        i = 0
        for M in tqdm(snappy.CensusKnots, desc="Knots...", ncols=120):
            with self.subTest(i=i):
                # print(M.identify()[0])
                basis = M.symplectic_basis()
                self.assertTrue(is_symplectic(basis.data))
                i += 1

    def test_link_complements(self):
        i = 0
        for M in tqdm(snappy.HTLinkExteriors[:300], desc="Links...", ncols=120):
            with self.subTest(i=i):
                if str(M.identify()[0]) in ERROR_MANIFOLDS:
                    continue
                else:
                    print(M.identify()[0])
                    basis = M.symplectic_basis()
                    self.assertTrue(is_symplectic(basis.data))
                    i += 1


if __name__ == "__main__":
    # test_link_complements_pool(12, 50)
    # unittest.main()
    for i in range(10):
        L = spherogram.random_link(1000, num_components=random.randint(20, 40), alternating=True)
        M = spherogram.Link.exterior(L)
        print(M.num_cusps())
        M.save("CuspedCensusData/link-" + str(i) + ".tri")