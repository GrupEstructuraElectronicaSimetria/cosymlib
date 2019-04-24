import tests


if __name__ == '__main__':
    testsuite = tests.TestLoader().discover('.')
    tests.TextTestRunner(verbosity=2).run(testsuite)
