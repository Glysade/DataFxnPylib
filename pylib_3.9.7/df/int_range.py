# Custom argparse type representing a bounded int
# from https://pyquestions.com/in-python-using-argparse-allow-only-positive-integers

import argparse


class IntRange:

    def __init__(self, imin=None, imax=None):
        self.imin = imin
        self.imax = imax

    def __call__(self, arg):
        try:
            value = int(arg)
        except ValueError:
            raise self.exception()
        if (self.imin is not None and value < self.imin) \
                or (self.imax is not None and value > self.imax):
            raise self.exception()
        return value

    def exception(self):
        if self.imin is not None and self.imax is not None:
            return argparse.ArgumentTypeError(f"Must be an integer in the"
                                              f" range [{self.imin},"
                                              f" {self.imax}]")
        elif self.imin is not None:
            return argparse.ArgumentTypeError(f"Must be an integer >="
                                              f" {self.imin}")
        elif self.imax is not None:
            return argparse.ArgumentTypeError(f"Must be an integer <="
                                              f" {self.imax}")
        else:
            return argparse.ArgumentTypeError("Must be an integer")
