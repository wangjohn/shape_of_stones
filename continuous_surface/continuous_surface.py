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


def frobenius_norm(x):
    return sqrt(sum(x**2, axis=-1))

def spectral_derivative(x, n=1):
    N = len(x)
    w = zeros(x.shape)

    for i in arange(2):
        w_hat = fft.rfft(x[:,i]) * (complex(0,1)*arange(int(N/2+1)))**n
        if n % 2 == 1:
            w_hat[-1] = 0
        w[:,i] = fft.irfft(w_hat)

    return w

def spectral_integral(x, n=1):
    N = len(x)
    w = zeros(x.shape)

    for i in arange(2):
        w_hat = fft.rfft(x[:,i]) * (1/(complex(0,1)*arange(int(N/2+1)))**n)
        if n % 2 == 1:
            w_hat[-1] = 0
        w[:,i] = fft.irfft(w_hat)

    return w

def curve_normal(x):
    out = x[:,(1,0)]
    out *= [-1,1]
    out /= frobenius_norm(out)[:,newaxis]
    return out

PLOT_RESOLUTION = 100
def plot_spectral(x):
    x_hat = fft.rfft(x)
    x_fine = fft.irfft(append(x_hat, zeros(PLOT_RESOLUTION/2)))
    
    x_fine *= len(x_fine)/len(x)
    plot(x_fine)

def plot_curve(x):
    x_fine = zeros([len(x)+PLOT_RESOLUTION, 2])
    for i in arange(2):
        x_hat = fft.rfft(x[:,i])
        x_fine[:,i] = fft.irfft(append(x_hat, zeros(PLOT_RESOLUTION/2)))
    
    x_fine *= len(x_fine)/len(x)
    plot(x_fine[:,0], x_fine[:,1])
    axis('equal')

def curvature(x):
    x_dot = spectral_derivative(x, 1)
    x_ddot = spectral_derivative(x, 2)
        
    k = cross(x_dot, x_ddot) / frobenius_norm(x_dot)**3
    return k

def integrate_surface(x_init, t, g):
    def func(x, t):
        x_col = x.reshape(-1, 2)
        dxdt = g(curvature(x_col))[:,newaxis] * curve_normal(spectral_derivative(x_col))
        return dxdt.flatten()
    
    return integrate.odeint(func, x_init.flatten(), t).reshape(len(t), -1, 2)


def g(k):
    return arctan((k-1)*5) + pi/2

def main():
    T = 0.1
    N = T/0.01
    s = linspace(0, 2*pi, 20, endpoint=False)
    x_init = vstack([2*cos(s), sin(s)]).T
    x_init = vstack([cos(s)*(chebyshev.chebval((s-pi)/pi, [0,0,1])+2), sin(s)*(chebyshev.chebval((s-pi)/pi, [0,0,0,0,1])+2)]).T

    t = linspace(0, T, N)
    x = integrate_surface(x_init, t, g)

    for i in arange(N, step=int(N/4)):
        plot_curve(x[i])
        plot(x[i,:,0],x[i,:,1], 'x')

    show()

main()
