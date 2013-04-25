#
# Copyright 2013 Thomas Trogdon
# This software is distributed under the terms of the GNU General Public License
#

import sys
import string
from libhill import *
import pickle

# Number of Fourier Modes
N = int(sys.argv[1])
# Period Length
L = eval(sys.argv[2])
# Number of Periods
P = int(sys.argv[3])
# Vector Size
V = int(sys.argv[4])
# Number of mu values
Z = int(sys.argv[5])
# String containing coefficients, then formatted
U = reformat(sys.argv[6],V)
# Boolean for Eigenvectors
vecbool = sys.argv[7]
# Boolean for using exixting Fourier Coefficients
coefbool = sys.argv[8]
# Initial variable for Fourier Coefficients
C = None

# Look to see if non-local coefficients are present, and format it
try:
    R = reformat(sys.argv[9],V)
    print 'Non-Local'
except:
    R = None
    print 'Local'

# Convert vecbool to a true boolean
if vecbool == 'True':
    vecflag = True
else:
    vecflag = False

if mod(Z,2) == 0:
    Z += -1

MU = linspace(-pi/(L*P),pi/(L*P),Z)



# If Fourier Coefficients are already present in a file.
if coefbool == 'True':
    pyDat = open("pyDat.dat")
    C = pickle.load(pyDat)
    pyDat.close()
    # Compute Spectrum Depending on local/non-local
    # returns  A = evals,lhs-fcoeff,evects(,rhs-fcoeff)
    if len(C) == 2:
        A = hill_comp(MU,L,N,P,U,R,vecflag,C[0],C[1])
    else:
        A = hill_comp(MU,L,N,P,U,R,vecflag,C[0])
else:
    # returns  A = evals,lhs-fcoeff,evects(,rhs-fcoeff)
    A = hill_comp(MU,L,N,P,U,R,evecs=vecflag)


# Open all storage files
EigenValuesReal = open("EigenValuesReal.dat","w")
EigenValuesIm = open("EigenValuesIm.dat","w")
EigenVectorsReal = open("EigenVectorsReal.dat","w")
EigenVectorsIm = open("EigenVectorsIm.dat","w")
mu = open("mu.dat","w")

print 'Writing Data'
for i in range(len(MU)):
    for j in range(len(A[0][i])):
        mu.write(str(MU[i])+'\n')
        EigenValuesReal.write(str(real(A[0][i][j]))+'\n')
        EigenValuesIm.write(str(imag(A[0][i][j]))+'\n')        
        #Efficient way to turn vector into delimited string
        EigenVectorsReal.write(string.join([`num` for num in real(A[2][i][:,j])],'\t')+'\n')
        EigenVectorsIm.write(string.join([`num` for num in imag(A[2][i][:,j])],'\t')+'\n')

# If Fourier coefficients were not given
if C is None:
    print 'Writing Coefficients'
    pyDat = open("pyDat.dat","w")
    if R is not None:
        pickle.dump([A[1],A[3]],pyDat)
    else:
        pickle.dump([A[1]],pyDat)
        
# All files automagically closed.
print "xyzzy"
#print MU*P*L/(2*pi)

