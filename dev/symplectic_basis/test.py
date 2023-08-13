import random
from datetime import datetime
import snappy
import unittest

import spherogram
from tqdm import tqdm
from multiprocessing import Pool
import itertools


if len(snappy.HTLinkExteriors(crossings=15)) == 0:
    file_name = "links-linux"
else:
    file_name = "links"

with open(file_name, "r") as file:
    lst = file.readline()
    print(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Building test set")
    manifolds = list(set([int(x) for x in lst.split(',')[:-1]]))
    manifolds_tri = [snappy.HTLinkExteriors[i] for i in manifolds]
    manifolds_labels = [M.identify()[0] for M in manifolds_tri if len(M.identify()) > 0]


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
    return sum([u[2 * i] * v[2 * i + 1] - u[2 * i + 1] * v[2 * i]
                for i in range(len(u) // 2)])


def save_manifold(index: int):
    M = snappy.HTLinkExteriors[index]
    M.save(f"CuspedCensusData/link-{index}.tri")


def process_manifold(i: int, output: bool = True):
    if output is False:
        M = snappy.HTLinkExteriors[i]
        index = i
        label = M.identify()[0] if len(M.identify()) > 0 else ""
    else:
        M = manifolds_tri[i]
        index = manifolds[i]
        label = manifolds_labels[i]

    if index == 0:
        return True

    basis = M.symplectic_basis()
    result = is_symplectic(basis)

    if result:
        string = "Passed"
    else:
        string = "Failed"

    if output:
        with open("logs/links-0.log", "a") as file:
            file.write(f"Testing: {str(index)} {(20 - len(str(index))) * ' '} {str(label)} {(40 - len(str(label))) * ' '} {string}\n")

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


def generate_tests(output_file: str):
    print(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Generating symplectic basis tests")
    for _ in range(2000):
        index = random.randint(1, len(snappy.HTLinkExteriors) - 1)
        process_manifold(index, output=False)
        print(index)

        with open(output_file, "a") as file:
            file.write(f"{index},")


def test_link_complements(start: int, end: int, manifolds_list: bool):
    scale = 1000
    with open("logs/total.log", "a") as file:
        file.write(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Testing: {str(scale * start)} - {str(scale * end - 1)}\n")

    if manifolds_list:
        result = [process_manifold(i) for i in range(len(manifolds))]
    else:
        result = [process_manifold(i) for i in range(scale * start, scale * end)]

    with open("logs/total.log", "a") as file:
        file.write(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Passed: {sum(result)} / {len(result)}\n")


def test_link_complements_pool(start: int, end: int, manifolds_list: bool):
    scale = 1000
    with open("logs/total.log", "a") as file:
        file.write(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Testing: {str(scale * start)} - {str(scale * end - 1)}\n")

    with Pool(maxtasksperchild=25) as pool:
        if manifolds_list:
            print(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Testing from Manifold List")
            result = pool.imap(process_manifold, range(len(manifolds)))
        else:
            print(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Testing from {scale * start} to {end * scale}")
            result = pool.imap(process_manifold, range(start * scale, end * scale))

        for _ in range(start, end):
            lst = list(itertools.islice(result, scale))

            #with open("logs/total.log", "a") as file:
            print(f"[{datetime.now().strftime('%d-%m-%y %H:%M:%S')}]  Passed: {sum(lst)} / {len(lst)}")


class TestSymplecticBasis(unittest.TestCase):
    def test_knot_complements(self):
        i = 0
        for M in tqdm(snappy.CensusKnots, desc="Knots...", ncols=120):
            with self.subTest(i=i):
                # print(M.identify()[0])
                basis = M.symplectic_basis()
                self.assertTrue(is_symplectic(basis), str(M.identify()[0]))
                i += 1

    # @unittest.skip
    def test_link_complements(self):
        i = 0
        for M in tqdm(snappy.HTLinkExteriors[1:1000], desc="Links...", ncols=120):
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
    test_link_complements_pool(0, 1, True)
    # test_link_complements(0, 1, True)
    # generate_tests("links-linux")
    # unittest.main()
