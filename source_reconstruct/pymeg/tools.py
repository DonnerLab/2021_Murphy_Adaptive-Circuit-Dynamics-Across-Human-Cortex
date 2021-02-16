from numpy import asarray, iscomplexobj, zeros, newaxis, zeros_like, real
import multiprocessing
import os
import pickle
import socket
import warnings
import logging

nthread = multiprocessing.cpu_count()

try:
    import pyfftw
    pyfftw.interfaces.cache.enable()
except:
    pyfftw = None
# Load wisdom from previous plans if it exists

#hostname = socket.gethostname()
#cache = os.path.join('%s_fftw.wisdom.pickle' % hostname)
#try:
#    wisdom = pickle.load(open(cache))
#    pyfftw.import_wisdom(wisdom)
#except IOError:
#    print('Did not load wisdom cache')
#    pass


def fft(X, shape=None):
    X = asarray(X)
    if shape is None:
        shape = X.shape
    inp = X.astype('complex128')
    out = pyfftw.empty_aligned(shape, dtype='complex128')
    fft_object = pyfftw.FFTW(out, out, flags=['FFTW_MEASURE'], threads=nthread,
                             planning_timelimit=1.)
    return fft_object(X)


def ifft(X, shape=None):
    if shape is None:
        shape = X.shape
    inp = pyfftw.empty_aligned(shape, dtype='complex128')
    out = pyfftw.empty_aligned(shape, dtype='complex128')
    ifft_object = pyfftw.FFTW(inp, out, flags=['FFTW_MEASURE'], threads=nthread,
                              planning_timelimit=1., direction='FFTW_BACKWARD')
    return ifft_object(X)


def fftconvolve(inp, kernel):
    shape = (len(inp) + len(kernel) - 1, )
    t = zeros(shape)
    t[:len(inp)] = inp
    fin = fft(t)
    t = 0 * t
    t[:len(kernel)] = kernel
    fK = fft(t)
    return real(ifft(fin * fK))[(len(kernel) / 2) - 1:], t


def hilbert(x, N=None):
    """
    Compute the analytic signal, using the Hilbert transform.
    The transformation is done along the last axis by default.
    -> From scipy but replace fft with fftw

    Note: Inpute must be (N_channels x Time)!
    Note: Input array is destroyed!
    """
    axis = 1
    x = asarray(x)
    if iscomplexobj(x):
        raise ValueError("x must be real.")
    if N is None:
        N = x.shape[axis]
    if N <= 0:
        raise ValueError("N must be positive.")

    ax = [0] * len(x.shape)
    xshape, xdim = x.shape, x.ndim
    Xf = fft(x)
    del x
    h = zeros(N)
    if N % 2 == 0:
        h[0] = h[N // 2] = 1
        h[1:N // 2] = 2
    else:
        h[0] = 1
        h[1:(N + 1) // 2] = 2
    if len(xshape) > 1:
        ind = [newaxis] * xdim
        ind[axis] = slice(None)
        h = h[ind]

    x = ifft(Xf * h)
    del Xf
    try:
        cPickle.dump(pyfftw.export_wisdom(), open(cache, 'w'))
    except:
        print('Did not save wisdom cache')
    return x


def deprecated(func):
    '''This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.'''
    def new_func(*args, **kwargs):
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func


import collections
import functools


class memoized(object):
    '''Decorator. Caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned
    (not reevaluated).
    '''

    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if not isinstance(args, collections.Hashable):
            # uncacheable. a list, for instance.
            # better to not cache than blow up.
            return self.func(*args)
        if args in self.cache:
            return self.cache[args]
        else:
            value = self.func(*args)
            self.cache[args] = value
            return value

    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__

    def __get__(self, obj, objtype):
        '''Support instance methods.'''
        return functools.partial(self.__call__, obj)
