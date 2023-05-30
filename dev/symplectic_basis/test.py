import snappy
import unittest
from tqdm import tqdm
import sys


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


class TestSymplecticBasis(unittest.TestCase):
    def test_knot_complements(self):
        i = 0
        for M in tqdm(snappy.CensusKnots, desc="Knots...", ncols=120):
            with self.subTest(i=i):
                # print(M.identify()[0])
                M = snappy.CensusKnots[i]
                basis = M.symplectic_basis()
                self.assertTrue(is_symplectic(basis.data))
                i += 1

    def test_link_complements(self):
        with open('dev/symplectic_basis/test.log', 'w') as file:
            i = 0
            initial_pos = file.tell()
            for M in tqdm(snappy.HTLinkExteriors[:1000], desc="Knots...", ncols=120, file=file):
                with self.subTest(i=i):
                    # M = snappy.HTLinkExteriors[i]
                    if str(M.identify()[0]) in ["3_1(0,0)", "5_1(0,0)", "8_19(0,0)", "9_1(0,0)"]:
                        continue
                    else:
                        # print(M.identify()[0])
                        basis = M.symplectic_basis()
                        self.assertTrue(is_symplectic(basis.data))

                    file.seek(initial_pos)
                    i += 1


if __name__ == "__main__":
    unittest.main()
