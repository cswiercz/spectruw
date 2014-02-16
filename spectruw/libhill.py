#
# Copyright 2013 Thomas Trogdon
# This software is distributed under the terms of the GNU General Public License
#


from numpy import *
from scipy.integrate import *
from scipy.linalg import norm, lu_factor, lu_solve, solve, inv
from scipy.lib.lapack.flapack import zggev, zgeev 
from scipy.special import *
import sys
import cPickle
import logging

###############################################
###     Functions for testing purposes      ###
###############################################


def printmat(A,string):
	f = open(string,'w')
	S = len(A)
	for i in range(S):
		for j in range(S):
			f.write(str(A[i,j])+' ')
		f.write('\n')
	f.close()

def savemat(A,string):
	f = file(string,'w')
	cPickle.dump(A,f,1)
	f.close()

##############################################
###Syntax checking and reformatting strings###
##############################################
def chk(a):
	x = 1
	k = 0
	try:
		eval(a)
		return
	except:
		b = sys.exc_info()[0]
		return str(b)

def syntax_chk(a):
	left_paren,right_paren = 0,0
	for item in a[:]:
		if item is '^':
			ab,sep,aa = a.partition('^')
			a = ab + '**' + aa
		if item is 'i':
			ab,sep,aa = a.partition('*i')
			if len(sep) > 0:
				a = ab + '*1j' + aa
			ab,sep,aa = a.partition('i*')
			if len(sep) > 0:
				a = ab + '1j*' + aa	
			if a is 'i':
				a = '1j'
		if item is '(':
			left_paren += 1
		if item is ')':
			right_paren += 1
	
	if right_paren is not left_paren:
		return "exceptions. Parentheses mismatch." 	
	
	b = chk(a)
	if b == None:
		return a
	else:
		return b

def reformat(coef_strlist,V):
	'''Accepts a list of coefficient functions separated by ";" for each 
	scalar problem and ":" for each coefficient function within a particular 
	scalar problem. V represents the dimension of the vector problem.
	
	Function returns a nested list such that first level of nesting represents 
	the ith equation in system, second level the jth function and finally the
	kth coefficient. All numbering done in Python syntax, i.e., starting 
	from 0. The returned list is in proper Python syntax.'''
	coef_strlist = coef_strlist.split(';')
	vec_coef_list = []
	l = len(coef_strlist)
	coef_list = []
	i = 0.0
	for scalar_strlist in coef_strlist:

		scalar_strlist = scalar_strlist.split(':')
		coef = []
		
		for coef_str in scalar_strlist:
			coef_str = syntax_chk(coef_str)
			if 'exceptions' in coef_str:
				return coef_str
			else:
				coef.append(coef_str)
	
		coef_list.append(coef)
		i = i + 1.

		if i%V == 0:
			vec_coef_list.append(coef_list)
			coef_list = []
		
	
	return vec_coef_list		

def insert_parameter(coef_strlist,param,num):
	'''Accepts a list of string coefficient functions and replaces every instance
	of param string with num.

	Function returns the modified list of string coefficient functions.'''
	
	return coef_strlist.replace(param, str(num))

###############################################################
###Matrix creation and E-value/E-vector estimation functions###
###############################################################
def construct_ef(v,L,N,P,mu,x):
	q = len(v)/(2*N+1)
	v = v/norm(v,ord=inf)
	Ev = zeros([q,len(x)],complex)
	for j in range(2*N+1):
		for n in range(q):
			Ev[n,:] += exp(1j*mu*x)*v[n*(2*N+1)+j]*exp(1j*(j-N)*pi*x/(L*P))
	print 'Efunc Real = ',norm(real(Ev),ord=inf)
	print 'Efunc Imag = ',norm(imag(Ev),ord=inf)
	return Ev

def invpower(e,A,B=None):
	if B is None:
		B = eye(len(A))
	K = A - B*e
	G = lu_factor(K)
	x = ones([len(A),1],complex)/len(A)
	iter = 0
	error = 10
	#for i in range(8):
	while error > 1e-8 and iter < 20:
		try:
			x = lu_solve(G,x)
		    
		except:
			print 'LU Exception'
			x = solve(K,x)
		    
		x = x/norm(x,ord=inf)
		error = norm(dot(A,x)-e*dot(B,x))
		iter = iter +1
		print 'invpower error = ',error
	x = x*conj(x[0])
	print 'Eval = ',e
	print 'Evect Real = ',norm(real(x))
	print 'Evect Imag = ',norm(imag(x))
	return x



# Find fourier coefficents
# The function F is passed in as a string
# The FFT version should be used instead
def fcoef_old(L,N,F):
	print 'Computing Fourier Coef'
	C = zeros(2*N+1,complex)
	L = .5*L
       	for i in range(2*N+1):
		    #Let F = U + iV
		U = lambda x: real(eval(F))
		V = lambda x: imag(eval(F))
		Ur = quad(U,-L,L,epsabs=1e-14,full_output=1,weight='cos',wvar=pi*(i-N)/L)[0]
		Vr = quad(V,-L,L,epsabs=1e-14,full_output=1,weight='sin',wvar=pi*(i-N)/L)[0]
		Vi = quad(V,-L,L,epsabs=1e-14,full_output=1,weight='cos',wvar=pi*(i-N)/L)[0]
		Ui = quad(U,-L,L,epsabs=1e-14,full_output=1,weight='sin',wvar=pi*(i-N)/L)[0]
		C[i] = 1/(2.0*L)*(Ur+Vr+1j*(Vi-Ui))
	#print F,sum(C)

	#Uncomment to output the same matricies as SpectrUW 2.0
	#for i in xrange(len(C)):
	#    if (i < N/2.) or (i > 3./2.*N):
	#	    C[i] = 0


	return C


def fcoef(L,N,F):

    logging.info('Computing Fourier Coef using FFT')
    C = zeros(2*N+1,complex)
    p = 2.0
    while p<N:
        p = p*2

    x = r_[0:2*p]*L*.5/p-L*.5
    F = eval(F)
    if (type(F) == type(1)) or (type(F)==type(1j)) or (type(F)==type(1.0)):
        F = ones_like(x)*F
    #f = F[where(x==0)]
    F = fft.fft(fft.fftshift(F))
    F = F/(2.0*p)
    C[:] = r_[F[-N:],F[:N+1]]

    #print 'sums =',sum(C),f
    
    #Uncomment to output the same matricies as SpectrUW 2.0
    #for i in xrange(len(C)):
    #    if (i < N/2.) or (i > 3./2.*N):
    #	    C[i] = 0

    return C

# Pass in f-coefficents, order, mu, L, N
# improved matrix construction: uses 1/4 the cpu_time
def hillmat(C,k,mu,L,N,P):
	J = xrange(2*N+1)
	for i in J:
		if mod(i-2*N,P) == 0:
			try:
				M += C[2*N-i/P]*eye(2*N+1,2*N+1,i-2*N)
			except:
				M = C[2*N-i/P]*eye(2*N+1,2*N+1,i-2*N)
			if i-2*N is not 0:
				try:
					M += C[i/P]*eye(2*N+1,2*N+1,2*N-i)
				except:
					M = C[i/P]*eye(2*N+1,2*N+1,2*N-i)
	for i in J:
		M[:,i] = (1.j*(mu+(2*pi*(i-N)/(P*L))))**k*M[:,i]
	return M

# If supplied with fourier coef then H is computed and returned
# Otherwise H and C are computed and returned
def scalarmat(U,mu,L,N,P,C='0'):
	a = 0
	if C is '0':
		a = 1
		C = zeros([len(U),2*N+1],complex)
	k = 0
	H = zeros([2*N+1,2*N+1],complex)
	for U_val in U:
		if U_val is not '0':
			if a is 1:
				C[k,:] = fcoef(L,N,U_val)

			H += hillmat(C[k,:],k,mu,L,N,P)		
		k = k + 1
		
	if a is 1:
		return [H,C]
	else:
		return H


def vecmat(U,mu,L,N,P,*FD):
	'''Main entry point for the Hill Method. Constructs the Hill Matrix as 
	explained in Computing spectra of linear operators using Hill's method by 
	Bernard Deconinck and J. Nathan Kutz (J. Comp. Physics 219, 296-321, 2006).

	The function returns the calculated eigen-values using LAPACK routines for 
	a given set of Floquet parameters. Input variables consist of
	
	U   Nested list of coefficient functions formed using reformat() 
		or otherwise
	mu  List of mu values (Floquet resolution)
	FD  Fourier coefficients may be given as a nested dictionary for 
		a given mu value. This will help speed up calculation for repeated 
		evaluations of eigenvalues.'''
	s = len(U)
	V = 2*N+1
	H = zeros([s*V,s*V],complex)
	#E = zeros([len(mu),s*V],complex)
	if FD:
		fd_exist=True
		C = FD[0]
	else:
		fd_exist=False
		C = {}
	
	#for mu_val in mu:
	logging.info('Constructing Matrix')
	for i in range(s):
		if fd_exist is False:
			C[i]={}
		for j in range(s):
			if fd_exist is False:
				[H[i*V:(i+1)*V,j*V:(j+1)*V],C[i][j]] = scalarmat(U[i][j],mu,L,N,P)
			else:
				H[i*V:(i+1)*V,j*V:(j+1)*V] = scalarmat(U[i][j],mu,L,N,P,C[i][j])
		
	
	#	if fd_exist is False:
	#		print 'Computing Evals'
	#		E[p,:] = eigvals(H)
	#	p = p + 1
	#print 'RNorm ' + str(norm(real(H)))
	#print 'INorm ' + str(norm(imag(H)))
	if fd_exist:
		return H
	else:
		return [H,C]


def hill_comp(MU,L,N,P,U,R=None,evecs=False,*fcoefs):
	''' This one doesn't suck'''
	# MU is an array of mu values
	s = len(U)
	n = 2*N+1
	E = zeros([len(MU),s*n],complex)
	V = zeros([len(MU),s*n,s*n],complex)
	# If fcoefs exist then get them
	if fcoefs:
		fcoefs_exist = True
		C = fcoefs[0]
		try:
			K = fcoefs[1]
			print 'Non-Local'
		except:
			print 'Local'
	else:
		fcoefs_exist = False

	# loop for each mu value
	# For the first iteration, when the fcoefs aren't provided
	# we need to compute fcoefs, else don't compute them
	
	p = 0
	print 'Computing Eigenvalues'
	for mu_val in MU:
		if fcoefs_exist:
			H = vecmat(U,mu_val,L,N,P,C)
			if R is not None:
				M = vecmat(R,mu_val,L,N,P,K)	
		else:
			H,C = vecmat(U,mu_val,L,N,P)
			if R is not None:
				M,K = vecmat(R,mu_val,L,N,P)
			fcoefs_exist = True
		if evecs:
			logging.info('Computing Eigenvectors')
			if R is not None:
				print 'Generalized Problem'
				a,b,vl,V[p,:,:],info = zggev(H,M,0,1)
				if info > 0:
					print 'Trying to Invert'
					E[p,:],vl,V[p,:,:],info = zgeev(dot(inv(M),H),0,1)
					print info,'Eigenvalues Did Not Converge'
				else:
					E[p,:] = a/b
			else:
				E[p,:],vl,V[p,:,:],info = zgeev(H,0,1)
				if info > 0:
					print info,'Eigenvalues Did Not Converge'
		else:
			if R is not None:
				print 'Generalized Problem'
				a,b,vl,vr,info = zggev(H,M,0,0)
				# print info,'Eigenvalues Did Not Converge'
				if info > 0:
					print 'Trying to Invert'
					E[p,:],vl,vr,info = zgeev(dot(inv(M),H),0,0)
					print info,'Eigenvalues Did Not Converge'
				else:
					E[p,:] = a/b
			else:
				E[p,:],vl,vr,info = zgeev(H,0,0)
				if info > 0:
					print info,'Eigenvalues Did Not Converge'
		p = p + 1
	if R is not None:
		return E,C,V,K,H
	else:
		return E,C,V,H
	
		       
#####################################################
###Additional functions not defined in Scipy/Numpy###
#####################################################
def sech(x):
	return 1/cosh(x)
	
def coth(x):
	return 1/tanh(x)
	
def csch(x):
	return 1/sinh(x)

def sn(x,k):
	return ellipj(x,k)[0]

def cn(x,k):
	return ellipj(x,k)[1]

def dn(x,k):
	return ellipj(x,k)[2]

def ph(x,k):
	return ellipj(x,k)[3]
	
