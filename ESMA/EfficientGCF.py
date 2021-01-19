import ramanujan.utils.mobius


class EfficientGCF(ramanujan.utils.mobius.EfficientGCF):
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
