# Custom argparse type representing a bounded float
# from https://pyquestions.com/in-python-using-argparse-allow-only-positive-integers

import argparse


class FloatRange:

    def __init__(self, fmin=None, fmax=None):
        self.fmin = fmin
        self.fmax = fmax

    def __call__(self, arg):
        try:
            value = float(arg)
        except ValueError:
            raise self.exception()
        if (self.fmin is not None and value < self.fmin) \
                or (self.fmax is not None and value > self.fmax):
            raise self.exception()
        return value

    def exception(self):
        if self.fmin is not None and self.fmax is not None:
            return argparse.ArgumentTypeError(f"Must be a float in the"
                                              f" range [{self.fmin},"
                                              f" {self.fmax}]")
        elif self.fmin is not None:
            return argparse.ArgumentTypeError(f"Must be an float >="
                                              f" {self.fmin}")
        elif self.fmax is not None:
            return argparse.ArgumentTypeError(f"Must be an float <="
                                              f" {self.fmax}")
        else:
            return argparse.ArgumentTypeError("Must be a float")
