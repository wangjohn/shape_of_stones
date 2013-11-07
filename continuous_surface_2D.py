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
PLOT_POINTS = True

filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), "paper/figures", PLOT_NAME) + ".pdf"

# -----------------------------------------------------------------------------
# Simulations
# -----------

T = 0.1
STEPS = T/1E-2
N = 20
METHOD = ["remove_points", "remove_lowest"][1]
SHAPE = ["circle", "ellipse", "blob"][2]


# -----------------------------------------------------------------------------
# Methods
# -------

def remove_points(shape):
    return arctan((shape.curvature()-1)*5) + pi/2

def remove_lowest(shape):
    n_fine = 1E2
    y_min = amin(real(fft.ifft(change_n(shape.x_hat, n_fine), axis=0))[:,1])
    g = 1/(5*absolute(shape.x[:,1] - y_min) + 0.5) - 0.5
    g[g<0] = 0
    return g

method = {"remove_points": remove_points, "remove_lowest": remove_lowest}[METHOD]

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

shape_func = {"circle":circle, "ellipse": ellipse, "blob": blob}[SHAPE]


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

def vnorm(x):
    return sqrt(sum(x**2, axis=-1))

def vdot(a, b):
    return sum(a * b, axis=-1)

def vcross(a, b):
    return vcross(a, b)

def change_n(x_hat, n):
    n_old = len(x_hat)
    if n > n_old:
        x_hat = insert(x_hat, int(n_old/2), zeros([n-n_old, 2]), axis=0)

    else:
        x_hat = take(x_hat, indices=fft_k(n), axis=0)

    return (n / n_old) * x_hat

def fft_k(n):
    return hstack(arange(n))

def fft_theta(n):
    return linspace(0, 2*pi, n, endpoint=False)

def spectral_derivative(x_hat, p=1):
    n = len(x_hat)
    k = fft_k(n)[:,newaxis]
    w_hat = x_hat * (1j*k)**p
    if p % 2 == 1:
        w_hat[n/2] = 0

    return w_hat

# @apply_to_cols
# def spectral_integral(x_hat, n=1):
#     k = arange(len(x_hat))
#     w_hat = x_hat * hstack([0, (1/(1j*k[1:])**n)])
#     if n % 2 == 1:
#         w_hat[-1] = 0

#     w_hat[0] = fft_theta(len(x_hat))
    
#     return w_hat

PLOT_N = 1E2
def plot_spectral(x_hat):
    s_fine = fft_theta(len(x_hat))
    x_fine = real(fft.ifft(change_n(x_hat, PLOT_N), axis=0))
    plot(s_fine, x_fine)


class SpectralShape(object):
    def __init__(self, x):
        self.x = x

    def __len__(self):
        return len(self.x_hat)

    @property
    def x(self):
        return real(fft.ifft(self.x_hat, axis=0))

    @x.setter
    def x(self, value):
        self.x_hat = fft.fft(value, axis=0)

    def surface_normal(self):
        x_dot = real(fft.ifft(spectral_derivative(self.x_hat), axis=0))
        x_dot_n = x_dot[:,(1,0)] * [-1,1]
        x_dot_n /= vnorm(x_dot_n)[:,newaxis]
        return x_dot_n

    def surface_tangent(self):
        x_dot = real(fft.ifft(spectral_derivative(self.x_hat), axis=0))
        x_dot /= vnorm(x_dot)[:,newaxis]
        return x_dot

    def curvature(self):
        x_dot = real(fft.ifft(spectral_derivative(self.x_hat), axis=0))
        x_ddot = real(fft.ifft(spectral_derivative(self.x_hat, p=2), axis=0))
            
        kappa = vcross(x_dot, x_ddot) / vnorm(x_dot)**3
        return kappa

    def dxdt(self, method):
        g = method(self)
        dx_hatdt = g[:,newaxis] * self.surface_normal()
        
        x_ddot = real(fft.ifft(spectral_derivative(self.x_hat, p=2), axis=0))
        a_t = vdot(x_ddot, self.surface_tangent())
        a_t *= norm(g) / norm(a_t)
        dx_hatdt += a_t[:,newaxis] * self.surface_tangent()

        return dx_hatdt

    def plot(self, label=None):
        x_fine = real(fft.ifft(change_n(self.x_hat, PLOT_N), axis=0))

        color = ax._get_lines.color_cycle.next()
        if PLOT_LINE:
            ax.plot(x_fine[:,0], x_fine[:,1], color, label=label)
        
        if PLOT_POINTS:
            ax.plot(self.x[:,0], self.x[:,1], "x", color="{}".format(color))

        axis('equal')

# -----------------------------------------------------------------------------

def run_simulation(shape, t_steps, method):
    def func(x, t):
        shape.x = x.reshape(-1,2)
        return shape.dxdt(method).flatten()

    x_simulation = integrate.odeint(func, shape.x.flatten(), t_steps)
    x_simulation = x_simulation.reshape(len(t_steps), -1, 2)

    for i in arange(STEPS, step=int(STEPS/2)):
        shape.x = x_simulation[i]
        shape.plot(label="t = {:.2f}".format(t_steps[i]))

    legend()
    savefig(filename)
    show()


s = fft_theta(N)
shape = shape_func(s)
t = linspace(0, T, STEPS)

run_simulation(shape, t, method)