import unittest
import warnings


if __name__ == '__main__':
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        suite = unittest.TestLoader().discover('.')
        unittest.TextTestRunner(verbosity=2).run(suite)
