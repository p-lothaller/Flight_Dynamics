import unittest
import Cit_par
import numpy as np
import warnings

warnings.filterwarnings('ignore', category=PendingDeprecationWarning)

class TestStateSpace(unittest.TestCase):

    def test_system_eigenvalues(self):

        poles_a, A_a, poles_s, A_s = Cit_par.main(0)
        print(A_a)
        assert all(poles_a) == all(np.linalg.eig(A_a)[0])
        assert all(poles_s) == all(np.linalg.eig(A_s)[0])


if __name__ == '__main__':
    unittest.main()
