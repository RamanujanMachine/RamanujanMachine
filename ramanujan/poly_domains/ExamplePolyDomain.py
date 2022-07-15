from .CartesianProductPolyDomain import CartesianProductPolyDomain


class ExampleDomain(CartesianProductPolyDomain):
    """
    This class is an example to how one may create a new structure for an and bn polynomials.

    We'll define the series to take the following forms:
    a(n) = c_0 * ((n+1)^2 + n^2) + c_1
    b(n) = c_3 * n^4

    c_0, c_1, c_2, c_3 are three coefficient that can get any range of integers.

    This is a descendant of CartesianProductPolyDomain since each coefficient is independent from the others.
    """

    def __init__(
        self, a_coefs_ranges, b_coef_range, a_deg, a_coef_range, b_deg, *args, **kwargs
    ):
        """
        Under this function you'll need to create 2 objects:
        1. a_coef_range
        2. b_coef_range
        both are lists that contains the range for each coefficient. so in our case, an should look like this
        [[c_0_min, c_0_max], [c_1_min, c_1_max]], and bn like [[c_3_min, c_3_max]]
        :param a_coefs_ranges:  ranges for an coefficients, as described above
        :param b_coef_range:  bn has only one coefficient, so this is the range for it - [c_3_min, c_3_max]
        """
        super().__init__(a_deg, a_coef_range, b_deg, b_coef_range, *args, **kwargs)
        self.a_coef_range = a_coefs_ranges
        # bn's coefficient range is converted to the expected form.
        self.b_coef_range = [b_coef_range]

        # make sure to call this function at the end of __init__. It will calculate an_size and bn_size, and
        # other data relevant to the algorithm
        self._setup_metadata()

    def get_calculation_method(self):
        """
        Here you must specify exactly how to use the coefficients generated. e.g. give a formula for an and bn
        :return: two iterators that generate series based on the given coefficients
        """

        def an_iterator(coefs, max_runs, start_n=0):
            for i in range(start_n, max_runs):
                yield coefs[0] * (i**2 + (i + 1) ** 2) + coefs[1]

        def bn_iterator(coefs, max_runs, start_n=0):
            for i in range(start_n, max_runs):
                yield coefs[0] * (i**4)

        return an_iterator, bn_iterator
