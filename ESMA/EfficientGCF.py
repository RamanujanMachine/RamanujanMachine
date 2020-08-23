import ramanujan.utils.mobius


class EfficientGCF(ramanujan.utils.mobius.EfficientGCF):
    """
    In MITM we refer to the first item in bn as b1, while in ESMA it is b0.
    This class allows ESMA to use EfficientGCF without enforcing the use of b1 as the first item.

    TODO - change ESMA to use b1 as the first item, and use ramanujan.utils.mobius.EfficientGCF without This patch
    """
    def __init__(self, a_, b_):
        self.prev_A = 0
        self.A = 1
        self.prev_B = 1
        self.B = a_[0]

        for i in range(1, len(a_)):
            tmp_a = self.A
            tmp_b = self.B
            self.A = a_[i] * self.A + b_[i - 1] * self.prev_A
            self.B = a_[i] * self.B + b_[i - 1] * self.prev_B
            self.prev_A = tmp_a
            self.prev_B = tmp_b
