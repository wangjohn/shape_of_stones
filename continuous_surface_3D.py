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

import os

from matplotlib import rc
rc("text", usetex=True)

from mpltools import style
style.use('ggplot')

figure(figsize=(5, 4))
ax = subplot()

# -----------------------------------------------------------------------------
# Plotting
# --------

SHOW_PLOT = True
SAVE_PLOT = False
PLOT_NAME = ""

PLOT_LINE = True
PLOT_POINTS = False

filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), "paper/figures", PLOT_NAME) + ".pdf"

# -----------------------------------------------------------------------------
# Simulations
# -----------

T = 0.5
N_STEPS = T/1E-2
N_PTS = 100
METHOD = ["remove_lowest"][0]
SHAPE = ["sphere"][0]


# -----------------------------------------------------------------------------
# Methods
# -------

def remove_lowest(shape):
    y_min = amin(irdft(increase_spectral_points(shape.x_hat, 10))[...,1])
    g = 1/(5*absolute(shape.x[...,1] - y_min) + 0.5) - 0.5
    g[g<0] = 0
    return g

method = {"remove_lowest": remove_lowest}[METHOD]

def sphere(s):
    x = cos(s)
    y = sin(s)
    z = (s - pi)/pi
    return SpectralShape(vstack([x, y, z]).T)

shape_func = {"sphere": sphere}[SHAPE]

# -----------------------------------------------------------------------------

def insert_dims(mat, index, n_dims):
    for i in arange(n_dims):
        mat = expand_dims(mat, index)
    
    return mat

def vec_to_mat(vec, dim, n_dim):
    vec = insert_dims(vec, 0, dim - 1)
    vec = insert_dims(vec, vec.ndims, n_dim - vec.ndim)
    return vec

def vnorm(x):
    return sqrt(sum(x**2, axis=-1))

def vdot(a, b):
    return sum(a * b, axis=-1)

def vcross(a, b):
    return cross(a, b)

def change_n(x_hat, n):
    n_old = x_hat.shape[:-1]
    for i, n_i in enumerate(n):
        if n_i > n_old[i]:
            shape = x_hat.shape
            shape[i] *= n_i - n_old[i]
            x_hat = insert(x_hat, int(n_old[i]/2), zeros(shape), axis=i)
        else:
            x_hat = take(x_hat, indices=fft_k(n_i), axis=i)

    return (prod(n) / prod(n_old)) * x_hat

def fft_k(n):
    return hstack([arange(n/2+1), arange(-n/2+1, 0)])

def fft_theta(n):
    return linspace(0, 2*pi, n, endpoint=False)

def fft_axes(x):
    return arange(x.ndim-1)

def spectral_derivative(x_hat, p):
    n = x_hat.shape[:,-1]
    for i, p_i  in enumerate(p):
        k = fft_k(n[i])
        w = (1j*k)**p[i]
        if p[i] % 2 == 1:
            w[n[i]/2] = 0

        x_hat *= vec_to_mat(w, i, len(p))

    return x_hat

# def spectral_integral(x_hat, n=1):
#     k = arange(len(x_hat))
#     w_hat = x_hat * hstack([0, (1/(1j*k[1:])**n)])
#     if n % 2 == 1:
#         w_hat[-1] = 0

#     w_hat[0] = linspace(0, 2*pi, rdft_n_shape(x_hat), endpoint=False)
    
#     return w_hat

PLOT_N = 1E1
def plot_spectral(x_hat):
    s_fine = fft_theta(len(x_hat))
    x_fine = real(fft.ifftn(change_n(x_hat, PLOT_N), axes=fft_axes(x_hat)))
    plot(s_fine, x_fine)


class SpectralShape(object):
    def __init__(self, x):
        self.x = x

    def __len__(self):
        return self.x_hat.shape

    @property
    def x(self):
        return fft.ifftn(self.x_hat, axes=fft_axes(self.x_hat))

    @x.setter
    def x(self, value):
        self.x_hat = fft.fftn(value, axes=fft_axes(value))

    def x_dot():

    def x_ddot():

    def surface_normal(self):
        x_dot = real(fft.ifftn(spectral_derivative(self.x_hat), axes=fft_axes(self.x_hat)))
        x_dot_n = x_dot[:,(1,0)] * [-1,1]
        x_dot_n /= vnorm(x_dot_n)[:,newaxis]
        return x_dot_n

    def surface_tangent(self):
        x_dot = real(fft.ifftn(spectral_derivative(self.x_hat, n=1)))
        x_dot /= vnorm(x_dot)[:,newaxis]
        return x_dot

    def centroid(self):
        x_dot = irdft(spectral_derivative(self.x_hat, n=1))
        area_hat = rdft(-sign(x_dot[:,0])*self.x[:,1])
        xy_hat = rdft(-sign(x_dot[:,0])*self.x[:,0]*self.x[:,1])
        yx_hat = rdft(sign(x_dot[:,1])*self.x[:,0]*self.x[:,1])



    def curvature(self):
        x_dot = irdft(spectral_derivative(self.x_hat, n=1))
        x_ddot = irdft(spectral_derivative(self.x_hat, n=2))
            
        k = vcross(x_dot, x_ddot) / vnorm(x_dot)**3
        return k

    def dxdt(self, method):
        g = method(self)
        dx_hatdt = rdft(g[:,newaxis] * self.surface_normal())
        
        x_ddot = irdft(spectral_derivative(self.x_hat, n=2))
        a_t = vdot(x_ddot, self.surface_tangent())
        a_t *= norm(g) / norm(a_t)
        dx_hatdt += rdft(a_t[:,newaxis] * self.surface_tangent())

        return dx_hatdt

    def plot(self, label=None):
        n_pts = len(self) * PLOT_N
        s_fine = linspace(0, 2*pi, n_pts, endpoint=False)
        
        x_fine = irdft(increase_spectral_points(self.x_hat, PLOT_N))

        color = ax._get_lines.color_cycle.next()
        if PLOT_LINE:
            plot(x_fine[:,0], x_fine[:,1], color, label=label)
        
        if PLOT_POINTS:
            plot(self.x[:,0], self.x[:,1], "x", color="{}".format(color))

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
        shape.x_hat = real_to_complex(x_hat_real.reshape(-1,4)).copy()
        return complex_to_real(shape.dxdt(method)).flatten()

    x_hat_real = complex_to_real(shape.x_hat).flatten()
    x_hat_simulation = integrate.odeint(func, x_hat_real, t_steps)
    x_hat_simulation = x_hat_simulation.reshape(len(t_steps), -1, 4)

    for i in arange(N_STEPS, step=int(N_STEPS/4)):
        shape.x_hat = real_to_complex(x_hat_simulation[i])
        shape.plot(label="t = {:.2f}".format(t_steps[i]))

    legend()
    savefig(filename)
    show()


s = linspace(0, 2*pi, N_PTS, endpoint=False)
# shape = shape_func(s)
# t = linspace(0, T, N_STEPS)

# run_simulation(shape, t, method)
