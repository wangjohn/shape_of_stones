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

from mpl_toolkits.mplot3d import Axes3D

import os

from matplotlib import rc
rc("text", usetex=True)

from mpltools import style
style.use('ggplot')

fig = figure(figsize=(5, 4))
ax = fig.add_subplot(111, projection='3d')

# -----------------------------------------------------------------------------
# Plotting
# --------

SHOW_PLT = True
SAVE_PLT = False
PLT_NAME = ""

PLOT_NPTS = 1E2

filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                        "paper/figures", PLT_NAME) + ".pdf"

# -----------------------------------------------------------------------------
# Simulations
# -----------

T = 0.5
SIM_NSTEPS = T/1E-2
SIM_NPTS = 10
METHOD = ["remove_lowest"][0]
SHAPE = ["sphere"][0]

FFT_NDIM = 2
FFT_AXES = (0, 1)

# -----------------------------------------------------------------------------
# Methods
# -------

def remove_lowest(shape, dim=1):
    x_min = amin(real(fft.ifft2(change_n(shape.x_hat, 1E2 * ones(FFT_NDIM)), 
                 axes=FFT_AXES))[...,dim])

    g = 1/(5*absolute(shape.x[...,dim] - x_min) + 0.5) - 0.5
    g[g<0] = 0
    return g

method = {"remove_lowest": remove_lowest}[METHOD]

def sphere(s):
    x = ( cos(s[1]) )[...,newaxis]
    y = ( sin(s[0]) * sin(s[1]) )[...,newaxis]
    z = ( cos(s[0]) )[...,newaxis]

    return SpectralShape(concatenate([x, y, z], axis=-1))

shape_func = {"sphere": sphere}[SHAPE]

# -----------------------------------------------------------------------------

def insert_dims(mat, index, n):
    for _ in arange(n):
        mat = expand_dims(mat, index)
    return mat

def vec_to_mgrid(vec, dim, n):
    vec = insert_dims(vec, 0, dim-1)
    return insert_dims(vec, -1, n - dim)

def vnorm(a):
    return sqrt(sum(a**2, axis=-1))

def vdot(a, b):
    return sum(a * b, axis=-1)

def vcross(a, b):
    return cross(a, b)

def e_vec(i, npts):
    vec = zeros(npts)
    vec[i] = 1
    return vec

def change_npts(x_hat, npts):
    npts_old = x_hat.shape[:-1]
    for i in arange(FFT_NDIM):
        if npts[i] > npts_old[i]:
            x_hat = insert(x_hat, int(npts_old[i]/2) * ones(npts[i] - npts_old[i]), 0, axis=i)
        else:
            x_hat = take(x_hat, indices=fft_k(npts[i]), axis=i)

    return (prod(npts) / prod(npts_old)) * x_hat

def fft_k(npts):
    if size(npts) == 1:
        return hstack([arange(npts/2+1, dtype=int), arange(-npts/2+1, 0, dtype=int)])
    else:
        k = [fft_k(npts_i) for npts_i in npts]
        return meshgrid(*k)

def fft_s(npts):
    if size(npts) == 1:
        s = linspace(0, 2*pi, npts, endpoint=False)
        return s + (s[1] - s[0])/2

    else:
        s = [fft_s(npts_i) for npts_i in npts]
        return meshgrid(*s)

def spectral_derivative(x_hat, p):
    npts = x_hat.shape[:-1]
    for i in arange(FFT_NDIM):
        k = fft_k(npts[i])
        w = (1j*k)**p[i]
        if p[i] % 2 == 1:
            w[npts[i]/2] = 0

        x_hat *= vec_to_mgrid(w, i, FFT_NDIM)

    return x_hat

def plot_spectral(x_hat):
    x_fine = real(fft.ifft2(change_npts(x_hat, PLOT_NPTS * ones(FFT_NDIM)), axes=FFT_AXES))
    s_fine = fft_s(x_fine.shape)

    contour(s_fine[...,0], s_fine[...,1], x_hat)


class SpectralShape(object):
    def __init__(self, x):
        self.x = x

    def shape(self):
        return self.x_hat.shape

    @property
    def x(self):
        return real(fft.ifft2(self.x_hat, axes=FFT_AXES))

    @x.setter
    def x(self, value):
        self.x_hat = fft.fft2(value, axes=FFT_AXES)

    def x_dot(self, p=1):
        return concatenate([real(fft.ifft2(spectral_derivative(self.x_hat, 
                    p=e_vec(i, FFT_NDIM)), axes=FFT_AXES))[...,newaxis] for i in arange(FFT_NDIM)], axis=-1)

    def surface_normal(self):
        x_dot = self.x_dot(p=1)
        x_dot_n = vcross(x_dot[...,0], x_dot[...,1])
        print(vnorm(x_dot_n))
        exit()
        x_dot_n /= vnorm(x_dot_n)[...,newaxis]
        return x_dot_n

    # def surface_tangent(self):
    #     x_dot = real(fft.ifft2(spectral_derivative(self.x_hat, n=1)))
    #     x_dot /= vnorm(x_dot)[:,newaxis]
    #     return x_dot

    # def centroid(self):
    #     x_dot = irdft(spectral_derivative(self.x_hat, n=1))
    #     area_hat = rdft(-sign(x_dot[:,0])*self.x[:,1])
    #     xy_hat = rdft(-sign(x_dot[:,0])*self.x[:,0]*self.x[:,1])
    #     yx_hat = rdft(sign(x_dot[:,1])*self.x[:,0]*self.x[:,1])



    # def curvature(self):
    #     x_dot = irdft(spectral_derivative(self.x_hat, n=1))
    #     x_ddot = irdft(spectral_derivative(self.x_hat, n=2))
            
    #     k = vcross(x_dot, x_ddot) / vnorm(x_dot)**3
    #     return k

    # def dxdt(self, method):
    #     g = method(self)
    #     dx_hatdt = rdft(g[:,newaxis] * self.surface_normal())
        
    #     x_ddot = irdft(spectral_derivative(self.x_hat, n=2))
    #     a_t = vdot(x_ddot, self.surface_tangent())
    #     a_t *= norm(g) / norm(a_t)
    #     dx_hatdt += rdft(a_t[:,newaxis] * self.surface_tangent())

        dxdt = self.surface_normal()

        return dxdt

    def plot(self, label=None):
        x_fine = real(fft.ifft2(change_npts(self.x_hat, PLOT_NPTS * ones(FFT_NDIM)), axes=FFT_AXES))

        ax.plot_surface(x_fine[...,0], x_fine[...,1], x_fine[...,2], rstride=int(PLOT_NPTS/10), cstride=int(PLOT_NPTS/10))
        # ax.scatter(x_fine[...,0], x_fine[...,1], x_fine[...,2], c='k')
        axis('equal')

# -----------------------------------------------------------------------------

def run_simulation(shape, t_steps, method):
    def func(x, t):
        shape.x = x.reshape(-1,2)
        return shape.dxdt(method).flatten()

    x_simulation = integrate.odeint(func, shape.x.flatten(), t_steps)
    x_simulation = x_simulation.reshape(size(t_steps), -1, 2)

    for i in arange(SIM_NSTEPS, step=int(SIM_NSTEPS/4), dtype=int):
        shape.x = x_simulation[i]
        shape.plot(label="t = {:.2f}".format(t_steps[i]))

    legend()
    savefig(filename)
    show()


s = fft_s(SIM_NPTS * ones(FFT_NDIM))
shape = shape_func(s)
# shape.surface_normal()

shape.plot()
show()

# t = linspace(0, T, SIM_NSTEPS)
# run_simulation(shape, t, method)
