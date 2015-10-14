"""Abstraction of mathematical operators to be between aospy objects."""


class Operator(object):
    def __init__(self, operator, objects):
        self.operator = operator
        self.objects = objects

    def __str__(self):
        return ("Operator object with operator '%s' and objects '%s'"
                % (self.operator, self.objects))

    __repr__ = __str__
