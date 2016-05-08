#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import re

def prime_powers(bound):
    '''
    Very naive sieve outputting the lists of primes and prime powers up to bound.
    '''
    
    primes = list(range(2, bound + 1))
    ppowers = list(range(2, bound + 1))
    for p in primes:
        for n in range(2*p, bound+1, p):
            try:
                primes.remove(n)
            except:
                try:
                    ppowers.remove(n)
                except:
                    pass
    return primes, ppowers

def parse_timings(datasets=('01/01', '01/02',
                                '02/01', '02/02', '02/03', '02/04', '02/05', '02/06',
                                '03/01', '03/02', '03/03', '03/04', '03/05', '03/06')):
    '''
    Parse timing data from files. Returns a pandas DataFrame.
    '''
    def read_file(f):
        d = pd.read_table(f, sep=r'[,\)]? \(?', comment='#', engine='python',
                            names=['prime', 'degree', 'rains_aux', 'allombert_aux',
                                       't_rains', 't_ellrains', 't_pari', 't_javad'],
                            na_values='0')
        d['dataset'] = f
        return d
    d = pd.concat([read_file('benchdata/%s.dat' % run) for run in datasets])
    return d

def plot_rains(d, size=(10,10)):
    '''
    Plot cyclotomic Rains' vs elliptic Rains' algorithm.
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.loglog(basex=2, basey=2)
    ax.set_xlabel('degree $r$')
    ax.set_ylabel('seconds')

    # Plot cyclotomic Rains', one line per auxiliary extension degree up to 8
    df = d.groupby([d.rains_aux, d.degree]).mean().loc[:9]
    ax.set_xlim(df.index.min()[1], df.index.max()[1])
    ax.set_ylim(df.t_rains.min(), df.t_rains.max())
    for key, g in df.groupby(level=0):
        ax.plot(g.index.get_level_values('degree'), g.t_rains, label=key, color='blue')
    # Plot elliptic Rains', scatter style
    d1 = d.sort_values('degree')
    ax.plot(d1.degree, d1.t_ellrains, 'r.', ms=2)

    # Legend
    cyclo = mlines.Line2D([], [], color='blue', label="Cyclotomic Rains'")
    ell = mlines.Line2D([], [], color='red', label="Elliptic Rains'")
    ax.legend(handles=[cyclo, ell], loc=2)

    return fig

def plot_allombert(d, size=(10,10)):
    '''
    Plot Pari/GP implementation of Allombert's algorithm vs ours.
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.loglog(basex=2, basey=2)
    ax.set_xlabel('degree $r$')
    ax.set_ylabel('seconds')

    # Plot Pari vs Javad's implementation of Allombert's algorithm,
    # one line per auxiliary extension degree up to 6
    df = d.groupby([d.allombert_aux, d.degree]).mean().loc[:7]
    ax.set_xlim(df.index.min()[1], df.index.max()[1])
    ax.set_ylim(df.t_javad.min(), df.t_javad.max())
    for key, g in df.groupby(level=0):
        ax.plot(g.index.get_level_values('degree'), g.t_pari, 'r-', label=key)
        ax.plot(g.index.get_level_values('degree'), g.t_javad, '--', label=key)

    # Legend
    pari = mlines.Line2D([], [], color='black', label="Pari/GP")
    javad = mlines.Line2D([], [], color='black', ls='--', label="Our implementation")
    ax.legend(handles=[pari, javad], loc=2)

    return fig

def plot_allombert_vs_ell(d, size=(10,10)):
    '''
    Scatter plot of running time ratios Allombert / elliptic Rains.
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'order of $q$ mod $r$')
    ax.set_ylabel("ratio")

    # Compute average running times for elliptic Rains'
    reference = d.groupby('degree')[['t_ellrains']].mean()
    df = d[['degree', 'allombert_aux', 't_pari', 't_javad']].merge(reference, left_on='degree', right_index=True)
    # Plot running time ratios
    df = df.sort_values('allombert_aux')
    ax.plot(df.allombert_aux, df.t_pari / df.t_ellrains, 'b.', label="Pari", ms=2)
    ax.plot(df.allombert_aux, df.t_javad / df.t_ellrains, 'rx', label="Our implementation", ms=2)
    ax.plot(df.allombert_aux, np.ones_like(df.allombert_aux), 'k-', label="Elliptic Rains'")

    # Legend
    ax.legend(loc=2)
    
    return fig

def plot_all(d, size=(10,10)):
    '''
    Comparison of all algorithms by degree.
    '''
    
    # Figure settings
    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    ax.loglog(basex=2, basey=2)
    ax.set_xlabel('degree $r$')
    ax.set_ylabel('seconds')

    # Plot cyclotomic Rains', auxiliary degree 1 and 2
    df = d.groupby([d.rains_aux, d.degree]).mean().loc[:2]
    for key, g in df.groupby(level=0):
        ax.plot(g.index.get_level_values('degree'), g.t_rains, label="Cyclotmic Rains' $s=%d$" % key)
    # Plot Pari, auxiliary degree 1, 2 and 10
    df = d.groupby([d.allombert_aux, d.degree]).mean().loc[[1,2,10],:]
    for key, g in df.groupby(level=0):
        ax.plot(g.index.get_level_values('degree'), g.t_pari, '--', label="Pari/GP $s=%d$" % key)
    # Plot elliptic Rains', scatter style
    df = d.sort_values('degree')
    ax.plot(df.degree, df.t_ellrains, '.', ms=2, label="Elliptic Rains'")
    
    # Legend
    ax.legend(loc=2)
    
    return fig

def pdf_plots(prefix='bench-'):
    '''
    Plot benchmarks.
    '''
    d = parse_timings()
    
    # Take out non-prime-powers to limit noise
    primes, ppowers = prime_powers(d.degree.max())
    d = d[d.degree.isin(ppowers) & (d.prime > 2**10) & (d.prime < 2**13)]

    # Global figure configurations
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', serif='Computer Modern')
    plt.rc('legend', fontsize=10)

    # Save pdf plots
    plot = plot_rains(d, (4,3))
    plot.savefig(prefix + 'rains.pdf', bbox_inches='tight')

    plot = plot_allombert(d, (4,3))
    plot.savefig(prefix + 'allombert.pdf', bbox_inches='tight')
    
    plot = plot_allombert_vs_ell(d, (4,3))
    plot.savefig(prefix + 'allombert-vs-ell.pdf', bbox_inches='tight')

    allombert = plot_all(d, (5,4))
    allombert.savefig(prefix + 'all.pdf', bbox_inches='tight')

    
if __name__ == '__main__':
    pdf_plots()
