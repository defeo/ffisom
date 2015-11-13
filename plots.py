#!/usr/bin/env sage

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from sage.all import is_prime, prime_powers

def primes_arith_prog(bound=10**4, out='arith_prog.pdf', d=3):
    '''
    Run through all prime powers r up to `bound`, and select for each
    the smallest u such that ur+1 is prime.
    
    Then do degree `d` linear regression, and plot everything into 
    `out`.
    '''
    # Find primes
    x, y = np.array([(p, next(u for u in xrange(p) if is_prime(u*p+1)))
                     for p in prime_powers(3,bound)]).T
    # linear regression
    lstsq = np.polyfit(np.log10(x), y, d)

    # Figure settings
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif='Computer Modern')
    plt.figure(figsize=(8,4))

    # Plot a 2D histogram of the data, with abscissa in log scale
    plt.hist2d(np.log10(x), y, bins=(1.5*max(y),max(y)/2),
               weights=np.log10(x)/x,
               cmap=cm.jet, norm=LogNorm())
    plt.colorbar()
    # Plot the linear regression
    x = np.linspace(*plt.xlim(), num=100)
    plt.plot(x, np.polyval(lstsq, x), c='black', label='degree %d polynomial regression' % d)
    plt.legend(loc='upper left')

    # More figure settings
    locs, _ = plt.xticks()
    plt.xticks(locs, map(lambda l: '$10^%d$' % l, locs))
    plt.tight_layout()
    # Write to file
    plt.savefig(out, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    # this takes ~2min to complete
    primes_arith_prog(10**8)
