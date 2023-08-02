import random
from datetime import datetime, timedelta
import snappy
import unittest

import spherogram
from tqdm import tqdm
from multiprocessing import Pool


def is_symplectic(M):
    """
    Test if the matrix M is symplectic
    :param M: square matrix
    :return: true or false
    """
    n = M.dimensions()

    for i in range(n[0]):
        for j in range(i, n[1]):
            omega = abs(symplectic_form(M.data[i], M.data[j]))

            if i % 2 == 0 and j % 2 == 1 and j == i + 1:
                if omega != 2:
                    return False
            elif omega:
                return False

    return True


def symplectic_form(u, v):
    return sum([u[2 * i] * v[2 * i + 1] - u[2 * i + 1] * v[2 * i] for i in range(len(u) // 2)])


def process_manifold(i: int):
    # index = random.randint(1, 200000)
    M = snappy.HTLinkExteriors[i]

    if len(M.identify()) > 0:
        label = M.identify()[0]
    else:
        label = ""

    if i == 0:
        return True

    basis = M.symplectic_basis()
    result = is_symplectic(basis)

    if result:
        string = "Passed"
    else:
        string = "Failed"

    with open("logs/links-" + str(i // 1000) + ".log", "a") as file:
        file.write(f"Testing: {str(index)} {(10 - len(str(index))) * ' '} {str(label)} {(30 - len(str(label))) * ' '} {string} \n")

    return result


def find_manifold(start: int):
    scale = 1000
    for i in range(scale * start, scale * (start + 1)):
        M = snappy.HTLinkExteriors[i]

        if len(M.identify()) == 0:
            return True


def random_link_exteriors(n: int, n_tet: int, n_cusps: int):
    for i in range(n):
        L = spherogram.random_link(n_tet, n_cusps, alternating=True)
        M = spherogram.Link.exterior(L)
        print(M.num_cusps())
        M.save(f"CuspedCensusData/link-{n_tet}-{n_cusps}-{i}.tri")


def test_link_complements_pool(start: int, end: int):
    scale = 1000
    for i in range(start, end):
        with open("logs/total.log", "a") as file:
            file.write(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Testing: {str(scale * i)} - {str(scale * (i + 1) - 1)}\n")

        # for i in range(scale * i, scale * (i + 1)):
        #     process_manifold(i)

        with Pool() as pool:
            result = pool.imap(process_manifold, range(scale * i, scale * (i + 1)))

            with open("logs/total.log", "a") as file:
                file.write(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Passed: {sum(result)} / {len(result)}\n")


class TestSymplecticBasis(unittest.TestCase):
    def test_knot_complements(self):
        i = 0
        for M in tqdm(snappy.CensusKnots, desc="Knots...", ncols=120):
            with self.subTest(i=i):
                # print(M.identify()[0])
                basis = M.symplectic_basis()
                self.assertTrue(is_symplectic(basis), str(M.identify()[0]))
                i += 1

    @unittest.skip
    def test_link_complements(self):
        i = 0
        for M in tqdm(snappy.HTLinkExteriors[1:5000], desc="Links...", ncols=120):
            with self.subTest(i=i):
                # print(M.identify()[0])
                basis = M.symplectic_basis()
                self.assertTrue(is_symplectic(basis))
                i += 1

    @unittest.skip
    def test_random_links(self):
        iterations = 10

        for i in tqdm(range(iterations), desc="Random Links...", ncols=120):
            with self.subTest(i=i):
                L = spherogram.random_link(100, num_components=random.randint(3, 10), alternating=True)
                M = spherogram.Link.exterior(L)
                basis = M.symplectic_basis()
                self.assertTrue(is_symplectic(basis))


if __name__ == "__main__":
    test_link_complements_pool(0, 5)
    # unittest.main()
    # M = snappy.HTLinkExteriors[159285]
    # M.symplectic_basis()
