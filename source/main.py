import mobius
import massey
import mpmath
import time

start = time.time()
with mpmath.workdps(10000):
    an = mobius.create_simple_continued_fraction(mpmath.phi, 50)
    print("\tsimple continued fraction of phi:\n{}".format(an))
    shift_reg = massey.slow_massey(an, 199)
    print("\tmassey shift register:\n{}\nlen: {}".format(shift_reg, len(shift_reg)))

    an = mobius.create_simple_continued_fraction(mpmath.e, 1000)
    print("\tsimple continued fraction of e:\n{}".format(an))
    shift_reg = massey.slow_massey(an, 5657)
    print("\tmassey shift register:\n{}\nlen: {}".format(shift_reg, len(shift_reg)))

    an = mobius.create_simple_continued_fraction(mpmath.pi, 2000)
    print("\tsimple continued fraction of pi:\n{}".format(an))
    shift_reg = massey.slow_massey(an, 5657)
    print("\tmassey shift register:\n{}\nlen: {}".format(shift_reg, len(shift_reg)))

end = time.time()

print("This took {0:.2f} seconds".format(end-start))
