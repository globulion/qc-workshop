#!/usr/bin/python3
"""
 Test the 4-index transformation approaches: sequential and full
"""
import numpy
from tutor.psithon import util
from time import time

# number of basis functions
n = 40
# number of occupied orbitals
o = 5

# generate some random LCAO-MO matrix and ERI's in AO basis
C = numpy.random.random((n, o))
eri_ao = numpy.random.random((n,n,n,n))

# compare
t1 = time()
eri_mo_1 = util.four_index_transform(eri_ao, C, C, C, C)
t2 = time()
eri_mo_2 = util.four_index_transform_full(eri_ao, C, C, C, C)
t3 = time()

# compute times
T1 = t2-t1
T2 = t3-t2

# check if we got the same result
similar = numpy.allclose(eri_mo_1, eri_mo_2, rtol=1e-09, atol=1e-09)
assert(similar is True), "There is an error in the implementation!"

print(" 4-index transformation: sequential  t = %10.4f" % T1)
print(" 4-index transformation: full        t = %10.4f" % T2)
print(" Sequential method is %4.1f times faster than full method" % (T2/T1))
