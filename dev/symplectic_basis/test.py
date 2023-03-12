import snappy
import unittest


class TestSymplecticBasis(unittest.TestCase):
    MANIFOLDS = [
        (
            # Figure 8 knot
            snappy.manifold('4_1'),
            [
                [0, -1, 1, 0],
                [-2, 0, 2, 0],
                [0, 0, -2, 0],
                [1, -1, 1, -1]
            ]
        )
    ]

    def test_symplectic_basis(self):
        for m, sb in self.MANIFOLDS:
            basis = m.symplectic_basis()
            self.assertEqual(basis, sb)


if __name__ == "__main__":
    unittest.main()
