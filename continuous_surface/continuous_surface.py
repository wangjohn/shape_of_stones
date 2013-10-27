from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

from pylab import *
from numpy import fft
from numpy import linalg
from scipy import integrate
from scipy import interpolate
from numpy.polynomial import chebyshev


# -----------------------------------------------------------------------------


PLOT_RESOLUTION = 1E3
T = 0.1
N_STEPS = T/0.01
N_PTS = 20
H_MIN = 0.2
METHOD = ["remove_points", "remove_lowest"][0]
SHAPE = ["circle", "ellipse", "blob"][2]

def remove_points(shape):
    return arctan((shape.curvature()-1)*5) + pi/2

def remove_lowest(shape):
    shape.extrema()

def circle(s):
    x = cos(s)
    y = sin(s)
    return SpectralShape(vstack([x, y]).T, s)

def ellipse(s):
    x = 2*cos(s)
    y = sin(s)
    return SpectralShape(vstack([x, y]).T, s)

def blob(s):
    x = cos(s)*chebyshev.chebval((s-pi)/pi, [2,0,1])
    y = sin(s)*chebyshev.chebval((s-pi)/pi, [2,0,0,0,1])
    return SpectralShape(vstack([x, y]).T, s)


# -----------------------------------------------------------------------------

DIGITS = 53;
TOL = sqrt(0.5**DIGITS)
def ellipse_circumference(a, b):
    """
    Compute the circumference of an ellipse with semi-axes a and b.
    Require a >= 0 and b >= 0.  Relative accuracy is about 0.5^53.
    """
    x = max(a, b)
    y = min(a, b)
    
    if DIGITS*y < TOL*x: 
        return 4 * x

    s = 0
    m = 1
    while x-y > tol*y:
        x = 0.5 * (x + y) 
        y = sqrt(x * y)
        m *= 2 
        s += m * (x - y)**2
    return pi * ((a + b)**2 - s) / (x + y)


def frobenius_norm(x):
    return sqrt(sum(x**2, axis=-1))

def apply_to_cols(f):
    def func(mat, *args, **kwargs):
        mat = mat.copy().reshape(len(mat), -1)

        for i in arange(mat.shape[1]):
            try:
                out = vstack([out, f(mat[:,i], *args, **kwargs)])
            except NameError:
                out = f(mat[:,i], *args, **kwargs)

        return out.T

    return func

@apply_to_cols
def get_spectral(x, s):
    n_pts = len(s)
    k = arange(n_pts/2+1)

    mat = exp(-1j*s*k[:,newaxis])

    return mat.dot(x)

@apply_to_cols
def get_values(x_hat, s):
    x_hat_full = hstack([x_hat, conjugate(x_hat[-2:0:-1])])
    n_pts = int((len(x_hat)-1)*2)

    k = hstack([arange(n_pts/2+1), arange(-n_pts/2+1, 0)])
    mat = (1/n_pts)*exp(1j*s[:,newaxis]*k)

    return real(mat.dot(x_hat_full))

@apply_to_cols
def spectral_interpolate(x, s, s_new):
    return get_values(get_spectral(x, s), s_new)

@apply_to_cols
def spectral_derivative(x, s, n=1):
    n_pts = len(s)

    k = arange(int(n_pts/2+1))
    x_hat = get_spectral(x, s) 
    w_hat = x_hat * (1j*k)**n
    if n % 2 == 1:
        w_hat[-1] = 0

    return get_values(w_hat, s)

@apply_to_cols
def spectral_integral(x, s, n=1, endpoint=False):
    n_pts = len(s)

    k = arange(int(n_pts/2+1))
    x_hat = get_spectral(x, s)
    w_hat = x_hat * hstack([0, (1/(1j*k[1:])**n)])
    if n % 2 == 1:
        w_hat[-1] = 0
    
    if endpoint:
        s = linspace(0, 2*pi, n_pts+1)
    else:
        s = linspace(0, 2*pi, n_pts, endpoint=False)

    const = (1/n_pts)*s*real(x_hat[0])
    
    return get_values(w_hat, s) + const

@apply_to_cols
def improve_resolution(x, s, min_pts):
    s_new = linspace(0, 2*pi, max([min_pts, len(x)]), endpoint=False)
    return spectral_interpolate(x, s, s_new)

def plot_spectral(s, x):
    s_fine = linspace(0, 2*pi, PLOT_RESOLUTION)
    x_fine = spectral_interpolate(x, s, s_fine)
    plot(s_fine, x_fine)


class SpectralShape(object):
    def __init__(self, x, s, h_min=H_MIN):
        self.h_min = h_min
        self.s = s
        self.x = x

    def __len__(self):
        return len(self.s)

    def index_to_full(self, index):
        return nonzero(self._valid_points)[0][index]

    def ignore_point(self, index):
        self._valid_points[self.index_to_full(index)] = False

    def convert_to_partial(self, mat):
        return mat[self._valid_points]

    def convert_to_full(self, mat):
        out = self._x.copy()
        out[self._valid_points] = mat
        return out

    @property
    def s(self):
        return self.convert_to_partial(self._s)

    @s.setter
    def s(self, value):
        self._s = value.copy()
        if not hasattr(self, '_valid_points'):
            self._valid_points = ones(len(self._s), dtype=bool)

    @property
    def x(self):
        return self.convert_to_partial(self._x)

    @x.setter
    def x(self, value):
        self._x = value.copy().reshape(-1, 2)
        self.redistribute()

    def surface_normal(self):
        x_dot = spectral_derivative(self.x, self.s, n=1)
        x_dot_n = x_dot[:,(1,0)]
        x_dot_n *= [-1,1]
        x_dot_n /= frobenius_norm(x_dot_n)[:,newaxis]
        return x_dot_n

    def curvature(self):
        x_dot = spectral_derivative(self.x, self.s, n=1)
        x_ddot = spectral_derivative(self.x, self.s, n=2)
            
        k = cross(x_dot, x_ddot) / frobenius_norm(x_dot)**3
        return k

    def redistribute(self):
        v = frobenius_norm(spectral_derivative(self.x, self.s, n=1))
        l = spectral_integral(v, self.s)
        dl = diff(l)

        while True:
            small = nonzero(dl < self.h_min)[0]
            if len(small) == 0:
                break

            index = small[0]
            dl[index+1] += dl[index]

            dl = delete(dl, index)
            self.ignore_point(index+1)

    def dxdt(self, method):
        return method(self)[:,newaxis] * self.surface_normal()

    def plot(self):
        s_fine = linspace(0, 2*pi, PLOT_RESOLUTION, endpoint=False)
        x_fine = spectral_interpolate(self.x, self.s, s_fine)

        plot(x_fine[:,0], x_fine[:,1])
        plot(self.x[:,0], self.x[:,1], 'x')
        axis('equal')


# -----------------------------------------------------------------------------


def run_simulation(shape, t_steps, method):

    def func(x, t):
        shape.x = x
        return shape.convert_to_full(shape.dxdt(method)).flatten()

    x_simulation = integrate.odeint(func, shape.x.flatten(), t_steps).reshape(len(t_steps), -1, 2)

    print(shape._valid_points)
    for i in arange(N_STEPS, step=int(N_STEPS/4)):
        shape.x = x_simulation[i]
        shape.plot()

    show()


method = {"remove_points": remove_points, "remove_lowest": remove_lowest}[METHOD]
s = linspace(0, 2*pi, N_PTS, endpoint=False)
shape = {"circle":circle, "ellipse": ellipse, "blob": blob}[SHAPE](s)
t = linspace(0, T, N_STEPS)

# run_simulation(shape, t, method)



plot_spectral(linspace(0, 2*pi, N_PTS-3, endpoint=False), sin(delete(s, [5, 8, 10])))
show()