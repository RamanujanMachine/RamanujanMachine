from .CartesianProductPolyDomain import CartesianProductPolyDomain


class ExplicitCartesianProductPolyDomain(CartesianProductPolyDomain):
    """
    This domain is the equivalent to CartesianProductPolyDomain, but the user may
    state the required range for each coef individually.
    It's used mostly in splitting a cartesian domain to BOINC clients.
    """

    def __init__(self, a_coefs=((0, 0),), b_coefs=((0, 0),), *args, **kwargs):
        super().__init__(a_deg=len(a_coefs) - 1, b_deg=len(b_coefs) - 1, a_coef_range=[0, 0], b_coef_range=[0, 0],
                         *args, **kwargs)

        self._setup_metadata()

    def _setup_metadata(self):
        super()._setup_metadata()

        # This domain's degree is not constant and needs to be calculated 
        self.a_deg = len(self.a_coef_range) - 1
        self.b_deg = len(self.b_coef_range) - 1
