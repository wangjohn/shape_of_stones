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
METHOD = ["remove_points", "remove_lowest"][0]
SHAPE = ["ellipse", "blob"][1]

def remove_points(shape):
    return arctan((shape.curvature()-1)*5) + pi/2

def remove_lowest(shape):
    shape.extrema()

def ellipse(theta):
    x = 2*cos(theta)
    y = sin(theta)
    return SpectralShape(vstack([x, y]).T)

def blob(theta):
    x = cos(theta)*chebyshev.chebval((theta-pi)/pi, [2,0,1])
    y = sin(theta)*chebyshev.chebval((theta-pi)/pi, [2,0,0,0,1])
    return SpectralShape(vstack([x, y]).T)


# -----------------------------------------------------------------------------


def frobenius_norm(x):
    return sqrt(sum(x**2, axis=-1))

def get_spectral(x):
    return fft.rfft(x)

def get_values(x_hat, theta):
    k = arange(len(x_hat))
    mat = (1/len(x_hat))*exp(complex(0,1)*theta[:,newaxis]*k)
    return mat.dot(x_hat)

def spectral_derivative(x, n=1):
    n_pts = x.shape[0]
    x_col = x.reshape(n_pts, -1)
    w = zeros(x_col.shape)

    for i in arange(x.shape[-1]):
        w_hat = fft.rfft(x_col[:,i]) * (complex(0,1)*arange(int(n_pts/2+1)))**n
        if n % 2 == 1:
            w_hat[-1] = 0
        w[:,i] = fft.irfft(w_hat)

    return w.reshape(x.shape)

def spectral_integral(x, n=1):
    n_pts = x.shape[0]
    x_col = x.reshape(n_pts, -1)
    w = zeros(x_col.shape)

    for i in arange(x_col.shape[1]):
        w_hat = fft.rfft(x_col[:,i]) * (1/(complex(0,1)*arange(int(n_pts/2+1)))**n)
        if n % 2 == 1:
            w_hat[-1] = 0
        w[:,i] = fft.irfft(w_hat)
    return w.reshape(x.shape)

def spectral_interpolation(x, min_n_pts):
    x_hat = fft.rfft(x)
    if len(x) < min_n_pts:
        x_hat = append(x_hat, zeros((min_n_pts - len(x))/2))
    
    x_fine = fft.irfft(x_hat)
    x_fine *= len(x_fine)/len(x)
    return x_fine

def plot_spectral(x):
        resolution = maximum(len(self.x), self.PLOT_RESOLUTION)
        s_fine = linspace(0, 2*pi, resolution, endpoint=False)
        x_fine = spectral_interpolation(self.x, resolution)
        plot(s_fine, x_fine)


class SpectralShape:
    def __init__(self, x):
        self.set_x(x)

    def set_x(self, x):
        self.x = x.reshape(-1, 2)

    def surface_normal(self):
        x_dot = spectral_derivative(self.x, 1)
        x_dot_n = x_dot[:,(1,0)]
        x_dot_n *= [-1,1]
        x_dot_n /= frobenius_norm(x_dot_n)[:,newaxis]
        return x_dot_n

    def curvature(self):
        x_dot = spectral_derivative(self.x, 1)
        x_ddot = spectral_derivative(self.x, 2)
            
        k = cross(x_dot, x_ddot) / frobenius_norm(x_dot)**3
        return k

    def redistribute(self, n_pts=N_PTS):
        x_dot = spectral_derivative(self.x, 1)
        v = frobenius_norm(x_dot)
        L = spectral_integral(v)

        v *= (2*pi)/L # Rescale path length to 2*pi
        theta_hat = get_spectral(spectral_integral(1/v))
        theta = get_values(theta_hat, linspace(0, 2*pi, n_pts, endpoint=False))

        self.set_x(get_values(get_spectral(self.x), theta))

    def dxdt(self, method):
        return method(self)[:,newaxis] * self.surface_normal()

    def plot(self):
        resolution = maximum(len(self.x), PLOT_RESOLUTION)
        x_fine = zeros([resolution, 2])
        for i in arange(2):
            x_fine[:,i] = spectral_interpolation(self.x[:,i], resolution)
        
        x_fine = vstack([x_fine, x_fine[0]])
        plot(x_fine[:,0], x_fine[:,1])
        plot(self.x[:,0], self.x[:,1], 'x')
        axis('equal')


# -----------------------------------------------------------------------------


def run_simulation(shape, t, method):
    def func(x, t):
        shape.set_x(x)
        return shape.dxdt(method).flatten()

    x_simulation = integrate.odeint(func, shape.x.flatten(), t).reshape(len(t), -1, 2)
    

    for i in arange(N_STEPS, step=int(N_STEPS/4)):
        shape.set_x(x_simulation[i])
        shape.plot()

    show()


method = {"remove_points": remove_points, "remove_lowest": remove_lowest}[METHOD]
theta = linspace(0, 2*pi, N_PTS, endpoint=False)
shape = {"ellipse": ellipse, "blob": blob}[SHAPE](theta)
t = linspace(0, T, N_STEPS)

#run_simulation(shape, t, method)

shape.redistribute()
shape.plot()
