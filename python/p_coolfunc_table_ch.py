import sys

import time

import numpy as np

from scipy.optimize import brentq

import mpmath as mp

import h5py

from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD

mp.mp.dps = 100

m_e       = 9.1093826e-28
m_i       = 1.67262171e-24 # proton mass --- really should be mean nuclear weight
sig_T     = 6.65e-25
chi       = 1.20
ln_Lambda = 20.0

def f(th_e, A):
    tmp = A - chi*(m_e/m_i)*th_e
    if (tmp <= 0.0):
        tmp = 1.0e-20

    th_e = mp.mpf(th_e)
    th_i = mp.mpf(tmp)

    tmp = (1/(mp.besselk(2, 1/th_e) * mp.besselk(2, 1/th_i))) * (((2*(th_e + th_i)**2 + 1)/(th_e + th_i)) * mp.besselk(1, (th_e + th_i)/(th_e*th_i)) + 2*mp.besselk(0, (th_e + th_i)/(th_e*th_i)))

#   print repr(tmp)

    return float(tmp)

def eqn(th_e, A, B, C):
    return ((3.0/8.0)*(m_e/m_i)*ln_Lambda)*(A - (m_e/m_i)*(1.0 + chi)*th_e)*f(th_e, A) - B*th_e*(1.0 + 4.0*th_e) + (1.0/4.0)*B*C

# ----------+---------- start: program execution ----------+----------

rank = COMM_WORLD.Get_rank()
size = COMM_WORLD.Get_size()

print 'hello from ' + repr(rank)
sys.stdout.flush()
sys.stderr.flush()

A_grid = np.logspace(-15, 15, 301, endpoint = True)
B_grid = np.logspace(-15, 15, 301, endpoint = True)
C_grid = np.logspace(-6,  1,  71,  endpoint = True)
C_grid = np.insert(C_grid, 0, 0.0)

# ----------+---------- root process ----------+----------
if (rank == 0):
    # start the clock
    start_time = time.time()

    th_e = np.zeros((len(A_grid), len(B_grid), len(C_grid)))

    todo = []
    for i in range(0, len(A_grid)):
        for j in range(0, len(B_grid)):
            for k in range(0, len(C_grid)):
                todo.append([i, j, k])

    job_index = 0
    for i in range(1, size):
        if (job_index < len(todo)):
            COMM_WORLD.send([job_index, todo[job_index]], dest = i, tag = 0)
            job_index += 1

    finished = 0
    while (finished != len(todo)):
        result    = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)
        finished += 1
        print repr(finished) + '/' + repr(len(todo))
        # record result
        th_e[result[0], result[1], result[2]] = result[3]
        # send next job to this processor
        if (job_index < len(todo)):
            COMM_WORLD.send([job_index, todo[job_index]], dest = result[4], tag = 0)
            job_index += 1

    # kill *all* processors
    for i in range(1, size):
        COMM_WORLD.send([-1, -1], dest = i, tag = 0)

    f = h5py.File('coolfunc_table_ch.h5', 'w')
    f.create_dataset('A_grid', (len(A_grid),), dtype = 'f')
    f['/A_grid'][...] = A_grid
    f.create_dataset('B_grid', (len(B_grid),), dtype = 'f')
    f['/B_grid'][...] = B_grid
    f.create_dataset('C_grid', (len(C_grid),), dtype = 'f')
    f['/C_grid'][...] = C_grid
    f.create_dataset('th_e', (len(A_grid), len(B_grid), len(C_grid)), dtype = 'f')
    f['/th_e'][...] = th_e
    f.close()

    # stop the clock, record time elapsed
    end_time = time.time()
    elapsed  = end_time - start_time
    print 'Time elapsed: ' + '{0:.2f}'.format(elapsed/3600) + 'h.'

# ----------+---------- worker process ----------+----------
else:
    while True:
        # receive job
        tmp = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)

        # break out if sent done signal
        if (tmp[0] == -1):
            break

        i = tmp[1][0]
        j = tmp[1][1]
        k = tmp[1][2]
    
        try:
            th_e = brentq(eqn, 1.0e-20, (1.0/chi)*(m_i/m_e)*A_grid[i], args = (A_grid[i], B_grid[j], C_grid[k]))
        except ValueError:
            th_e = (1.0/chi)*(m_i/m_e)*A_grid[i]
        else:
            th_i = A_grid[i] - (m_e/m_i)*chi*th_e
            if (th_i <= 0.0):
                th_e = (1.0/chi)*(m_i/m_e)*A_grid[i]
#       print 'th_e = ' + repr(th_e)
#       print 'th_i = ' + repr(A_grid[i] - (m_e/m_i)*chi*th_e)

        reply = (i, j, k, th_e, rank)

        COMM_WORLD.send(reply, dest = 0, tag = 0)

