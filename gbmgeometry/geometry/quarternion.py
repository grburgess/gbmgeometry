import numpy as np


class Quaternion(object):
    def __init__(self, q):

        assert len(q) == 4

        # self._q = q / np.dot(q,q)
        self._q = q
        # convert q to (n,a) --- not used
        # if n is None and a is None:
        #    self.a = 2*np.arccos(self._q[0])
        #    self.n = self._q[1:] / np.sin(self.a/2)

    @property
    def w(self):
        return self._q[0]

    @property
    def x(self):
        return self._q[1]

    @property
    def y(self):
        return self._q[2]

    @property
    def z(self):
        return self._q[3]

    @property
    def q(self):
        return self._q

    @property
    def vector(self):
        """ get the vector part """
        return self._q[1:]

    @property
    def scalar(self):
        """ get the scalar part """
        return self._q[0]

    @property
    def magnitude(self):
        """ get the magnitude """
        return np.sqrt(np.dot(self._q, self._q))

    def is_unit(self, tol=1e-10):
        """ check for unit quaternion """
        return np.abs(1.0 - self.magnitude) < tol

    @property
    def normalized(self):
        """ normalize a quaternion """
        if not self.is_unit():
            n = self.magnitude
            if n > 0:
                # self._q /= n # this will change the main object
                return self.__class__(q=self._q / n)

        return self._q

    @property
    def conjugate(self):
        """ 
        get conjugate of a quaternion
        """
        return self.__class__(
            q=np.concatenate((np.array([self.scalar]), -self.vector), axis=0)
        )

    @property
    def inverse(self):
        """ get inverse of a quaternion """
        ss = np.dot(self._q, self._q)
        if ss > 0:
            d = self.conjugate
            return self.__class__(q=d.q / ss)
        else:
            raise ZeroDivisionError("a zero quaternion cannot be inverted")

    @property
    def log(self):
        """ get log of a quaternion """
        v = self.vector
        s = self.scalar
        z = (v / np.sqrt(np.dot(v, v))) * np.arccos(s / self.magnitude)
        r = np.concatenate(([np.log(self.magnitude)], z), axis=0)
        return self.__class__(q=r)

    @property
    def exp(self):
        """ get exp of a quaternion """
        v = self.vector
        vn = np.sqrt(np.dot(v, v))
        s = self.scalar
        r = np.exp(s) * np.concatenate(([np.cos(vn)], (np.sin(vn) / vn) * v), axis=0)
        return self.__class__(q=r)

    def axisangle(self):
        """ quaternion to axis-angle """
        self.a = 2 * np.arccos(self._q[0])
        self.n = self._q[1:] / np.sin(self.a / 2)
        return self.n, self.a

    #### OPERATIONS ####
    def add(self, other):
        return self.__class__(q=self._q + other.q)

    def sub(self, other):
        return self.__class__(q=self._q - other.q)

    def mul(self, other):
        q1 = self._q
        q2 = other.q
        w = q1[0] * q2[0] - np.dot(q1[1:], q2[1:])
        v = q2[0] * q1[1:] + q1[0] * q2[1:] + np.cross(q1[1:], q2[1:])
        m = np.concatenate((np.array([w]), v), axis=0)
        return self.__class__(q=m)

    def div(self, other):
        q2i = other.inverse
        return self.mul(q2i)

    def rotatePoint(self, p):
        """
        rotate a point using P = q P q^-1
        """
        # convert the point to a quaternion format
        P = np.concatenate((np.array([0.0]), p), axis=0)
        P = self.__class__(q=P)
        Pn = self.mul(P)
        Pr = Pn.mul(self.inverse)
        return Pr.vector

    def __str__(self):
        """ the format when we print the quaternion """
        return "[{:.3f} {:.3f} {:.3f} {:.3f}]".format(
            self._q[0], self._q[1], self._q[2], self._q[3]
        )

    def __repr__(self):
        """ the format in the command line, using repr() """
        return "Quaternion({}, {}, {}, {})".format(
            repr(self._q[0]), repr(self._q[1]), repr(self._q[2]), repr(self._q[3])
        )
