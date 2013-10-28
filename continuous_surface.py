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

from matplotlib import rc
rc("text", usetex=True)

from mpltools import style
style.use('ggplot')

figure(figsize=(5, 4))
plot_name = "remove_lowest_blob_bad_points_zoom"
filename = "/home/jonathan/Dropbox/Documents/Coursework/Project Lab in Mathematics/shape_of_stones/continuous_surface/figures/{}.pdf".format(plot_name)

# -----------------------------------------------------------------------------

SHOW_X = False
SHOW_LINE = True
SHOW_POINTS = False
RESOLUTION_FACTOR = 2**(0)
PLOT_NPTS_FACTOR = 1E1
T = 0.2
N_STEPS = T/1E-2
N_PTS = 20
H_MIN = 0.2
METHOD = ["remove_points", "remove_lowest"][1]
SHAPE = ["circle", "ellipse", "blob"][2]

def remove_points(shape):
    return arctan((shape.curvature()-1)*5) + pi/2

def remove_lowest(shape):
    y_min = amin(get_values(increase_spectral_points(shape.x_hat, 1E1))[:,1], axis=0)
    g = 1/(10*absolute(shape.x[:,1] - y_min)+1E-1)
    g -= 1E-1
    g[g<0] = 0
    return g

def circle(s):
    x = cos(s)
    y = sin(s)
    return SpectralShape(vstack([x, y]).T)

def ellipse(s):
    x = 2*cos(s)
    y = sin(s)
    return SpectralShape(vstack([x, y]).T)

def blob(s):
    x = cos(s)*chebyshev.chebval((s-pi)/pi, [2,0,1])
    y = sin(s)*chebyshev.chebval((s-pi)/pi, [2,0,0,0,1])
    return SpectralShape(vstack([x, y]).T)


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
def get_spectral(x):
    return fft.rfft(x)

@apply_to_cols
def get_values(x_hat):
    return real(fft.irfft(x_hat))

def spectral_npts(x_hat):
    return (len(x_hat)-1)*2

@apply_to_cols
def spectral_derivative(x_hat, n=1):
    k = arange(len(x_hat))
    w_hat = x_hat * (1j*k)**n
    if n % 2 == 1:
        w_hat[-1] = 0

    return w_hat

@apply_to_cols
def spectral_integral(x_hat, n=1):
    k = arange(len(x_hat))
    w_hat = x_hat * hstack([0, (1/(1j*k[1:])**n)])
    if n % 2 == 1:
        w_hat[-1] = 0
    w_hat[0] = linspace(0, 2*pi, spectral_npts(x_hat), endpoint=False)
    
    return w_hat

@apply_to_cols
def increase_spectral_points(x_hat, factor):
    n_pts = spectral_npts(x_hat)
    if factor > 1:
        x_hat = append(x_hat, zeros((n_pts*(factor - 1))/2))

    elif factor < 1:
        x_hat = x_hat[:int((n_pts*factor)/2+1)]

    x_hat *= spectral_npts(x_hat) / n_pts
    return x_hat

def plot_spectral(x_hat):
    n_pts = spectral_npts(x_hat) * PLOT_NPTS_FACTOR
    s_fine = linspace(0, 2*pi, n_pts)
    x_fine = get_values(increase_spectral_points(x_hat, PLOT_NPTS_FACTOR))
    plot(s_fine, x_fine)


class SpectralShape(object):
    def __init__(self, x):
        self.x = x
        self.resolution_factor = RESOLUTION_FACTOR

    def __len__(self):
        return int(spectral_npts(self._x_hat) * self.resolution_factor)

    @property
    def x(self):
        return get_values(self.x_hat)

    @x.setter
    def x(self, value):
        self._x_hat = get_spectral(value)

    @property
    def x_hat(self):
        return increase_spectral_points(self._x_hat, self.resolution_factor)

    @x_hat.setter
    def x_hat(self, value):
        self._x_hat = value.copy()
        # self._x_hat[4:,0] = 0

    def surface_normal(self):
        x_dot = get_values(spectral_derivative(self.x_hat, n=1))
        x_dot_n = x_dot[:,(1,0)]
        x_dot_n *= [-1,1]
        x_dot_n /= frobenius_norm(x_dot_n)[:,newaxis]
        return x_dot_n

    def curvature(self):
        x_dot = get_values(spectral_derivative(self.x_hat, n=1))
        x_ddot = get_values(spectral_derivative(self.x_hat, n=2))
            
        k = cross(x_dot, x_ddot) / frobenius_norm(x_dot)**3
        return k

    def dxdt(self, method):
        dx_hatdt = get_spectral(method(self)[:,newaxis] * self.surface_normal())
        
        # x_dot = get_values(spectral_derivative(self.x_hat, n=1))
        # v = frobenius_norm(x_dot)
        # x_ddot = get_values(spectral_derivative(self.x_hat, n=2))
        
        # a_t = x_dot / v[:,newaxis]
        # a_t *= sum(x_ddot * x_dot, axis=1)[:,newaxis] / 5

        # dx_hatdt += get_spectral(a_t)*5

        return dx_hatdt[:len(self._x_hat)]

    def plot(self, label=None):
        n_pts = len(self) * PLOT_NPTS_FACTOR
        s_fine = linspace(0, 2*pi, n_pts, endpoint=False)
        
        x_fine = get_values(increase_spectral_points(self._x_hat, PLOT_NPTS_FACTOR))

        if SHOW_X:
            plot_spectral(self._x_hat[:,0])

        if SHOW_LINE:
            plot(x_fine[:,0], x_fine[:,1], label=label)
        
        if SHOW_POINTS:
            plot(self.x[:,0], self.x[:,1], 'x')

        axis('equal')

# -----------------------------------------------------------------------------

def complex_to_real(mat):
    out = hstack([real(mat), imag(mat)])
    return out

def real_to_complex(mat):
    out = mat[:,:2] + 1j*mat[:,2:]
    return out

def run_simulation(shape, t_steps, method):
    def func(x_hat_real, t):
        shape.x_hat = real_to_complex(x_hat_real.reshape(-1,4))
        return complex_to_real(shape.dxdt(method)).flatten()

    x_hat_real = complex_to_real(shape._x_hat).flatten()
    x_hat_simulation = integrate.odeint(func, x_hat_real, t_steps)
    x_hat_simulation = x_hat_simulation.reshape(len(t_steps), -1, 4)

    for i in arange(N_STEPS, step=int(N_STEPS/4)):

        shape.x_hat = real_to_complex(x_hat_simulation[i])
        shape.plot("t = {:.2f}".format(t_steps[i]))
    legend()
    savefig(filename)
    show()


method = {"remove_points": remove_points, "remove_lowest": remove_lowest}[METHOD]
s = linspace(0, 2*pi, N_PTS, endpoint=False)
shape = {"circle":circle, "ellipse": ellipse, "blob": blob}[SHAPE](s)
t = linspace(0, T, N_STEPS)

run_simulation(shape, t, method)