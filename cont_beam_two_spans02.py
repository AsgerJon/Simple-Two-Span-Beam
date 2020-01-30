import numpy as np
# from numba import njit
import matplotlib.pyplot as plt
from beam_drawing import *

# @njit
def load(a, b, L, q, R1, N=101):
	N += (N+1)%2
	Rq = q * (L - a - b)
	Rxa = a + Rq/2/q
	Rxb = b + Rq/2/q
	R2 = 2*Rq/L*(Rxb)-2*R1
	R3 = Rq - R2 - R1
	dx = (N+2)*L/N**2
	p = np.zeros((N))
	x = np.linspace(-L/N,L*(N+1)/N,N)
	for i in range(N):
		if x[i] > a and x[i] < L - b:
			p[i] = -q
	return p*L/N
		
# @njit
def shear(a, b, L, q, R1, N=100):
	Rq = q * (L - a - b)
	Rxa = a + Rq/2/q
	Rxb = b + Rq/2/q
	R2 = 2*Rq/L*(Rxb)-2*R1
	R3 = Rq - R2 - R1
	V = np.zeros((N))
	x = np.linspace(0,L,N)
	X = np.array([a, L/2, L-b, L])
	Vs = np.array([R1, R1, R1 - q*(L/2-a) + R2, R1+R2-Rq])
	for i in range(N):
		if x[i] < X[0]:
			V[i] = Vs[0]
		elif x[i] < X[1]:
			V[i] = Vs[1] - q*(x[i] - X[0])
		elif x[i] < X[2]:
			V[i] = Vs[2] - q*(x[i] - X[1])
		elif i < N - 1:
			V[i] = Vs[3]
		else:
			V[i] = Vs[3]
	return V
	
# @njit
def moment(a, b, L, q, R1, N=101):
	Rq = q * (L - a - b)
	Rxa = a + Rq/2/q
	Rxb = b + Rq/2/q
	R2 = 2*Rq/L*(Rxb)-2*R1
	R3 = Rq - R2 - R1
	M = np.zeros((N))
	V = shear(a, b, L, q, R1, N)
	x = np.linspace(0, L, N)
	for i in range(N):
		if x[i] < a:
			M[i] = R1 * x[i]
		elif x[i] < L/2:
			M0 = R1 * a
			h = R1
			dx = x[i] - a
			M[i] = M0 + dx/2*(2*h - q*dx)
		elif x[i] < L - b:
			h = R1
			dx = L/2 - a
			M0 = R1 * a + dx/2*(2*h-q*dx)
			h = R1 + R2 - q*(L/2-a)
			dx = x[i] - L/2
			M[i] = M0 + dx/2*(2*h - q*dx)
		else:
			h = R1
			dx = L/2 - a
			M0 = R1 * a + dx/2*(2*h-q*dx)
			h = R1 + R2 - q*(L/2-a)
			dx = L - b - L/2
			M0 += dx/2*(2*h-q*dx)
			dx = x[i] - (L - b)
			M[i] = M0 - R3 * dx
	return M

# @njit
def energi(a, b, L, q, R1, N=101):
	N += (N+1)%2
	M = moment(a, b, L, q, R1, N)
	u = M[0]**2 + M[-1]**2
	for i in range(1,N-1):
		u += (2+(i%2)*2)*M[i]**2
	return u

# @njit
def get_R1(a, b, L, q, N=101):
	r1 = 0
	r1max = q*L
	r1inc = q*L/20
	r1min = r1
	minu = energi(a, b, L, q, r1, N)
	while r1inc > 1e-06:
		while r1 < r1max:
			utest = energi(a, b, L, q, r1, N)
			if utest < minu:
				r1min = r1
				minu = utest
			r1 += r1inc
		r1 = r1min - r1inc
		r1max = r1min + r1inc
		r1inc /= 10
	return r1min
	
def get_R2(a, b, L, q, R1):
	Rq = q * (L - a - b)
	Rxa = a + Rq/2/q
	Rxb = b + Rq/2/q
	R2 = 2*Rq/L*(Rxb)-2*R1
	return R2
	
def get_R3(a, b, L, q, R1):
	Rq = q * (L - a - b)
	Rxa = a + Rq/2/q
	Rxb = b + Rq/2/q
	R2 = 2*Rq/L*(Rxb)-2*R1
	R3 = Rq - R2 - R1
	return R3
	
def get_x_F0(L, F):
	N = F.shape[0]
	x = np.linspace(0,L,N)
	xF0 = []
	for i in range(1,N-1):
		if F[i-1]*F[i] < 0:
			if abs(F[i-1]) < abs(F[i]):
				xF0.append(x[i-1])
			else:
				xF0.append(x[i])
	return xF0
	
def draw_beam(a, b, L, q, N=1001):
	R1 = get_R1(a, b, L, q, N)
	R2 = get_R2(a, b, L, q, R1) *1e-03
	R3 = get_R3(a, b, L, q, R1) *1e-03
	x = np.linspace(0, L, N)
	V = shear(a, b, L, q, R1, N) *1e-03
	M = moment(a, b, L, q, R1, N) *1e-06
	R1 *= 1e-03
	xM0 = get_x_F0(L, M)
	xV0 = get_x_F0(L, V)
	xV0 = [i for i in xV0 if abs((i-L/2))/L > 0.1 ]
	xt = sorted([0, L ]+xV0+xM0)
	s = L*15e-03
	plt.figure(dpi=100, figsize=(6,9), )
	ax1 = plt.subplot(411)
	basic_beam(L)
	add_pinned_support(0, s)
	add_roller_support(L/2, s)
	add_roller_support(L, s)
	add_distributed_load(a, L-b, s)
	plt.yticks([])
	plt.subplots_adjust(hspace = 0)
	ax1.set_xticks([])
	ax2 = plt.subplot(412, sharex = ax1)
	plt.plot(x, V, 'k')
	plt.plot(x, np.zeros((N)), '--k', lw = 0.5)
	for i in xV0:
		plt.plot((i,i), (np.amin(V), np.amax(V)), '--k', lw = 0.75)
	for i in xM0:
		plt.plot((i,i), (np.amin(V), np.amax(V)), '--k', lw = 0.75)
	plt.ylabel("V(x)\n[kN]", fontsize='x-large')
	ax3 = plt.subplot(413, sharex = ax1)
	plt.plot(x, M, 'k')
	plt.plot(x, np.zeros((N)), '--k', lw = 0.5)
	for i in xV0:
		plt.plot((i,i), (np.amin(M), np.amax(M)), '--k', lw = 0.75)
	for i in xM0:
		plt.plot((i,i), (np.amin(M), np.amax(M)), '--k', lw = 0.75)
	plt.gca().invert_yaxis()
	plt.ylabel("M(x)\n[kNm]", fontsize='x-large')
	ax3.set_xticks(xt)
	plt.setp(ax3.get_xticklabels(), rotation=90, horizontalalignment='right')
	ax4 = plt.subplot(414)
	ax4.axis('off')
	topmargin = 0.4
	linemargin = 0.125
	leftmargin = 0.02
	textlines = ["Venstre Reaktion: "]
	textlines.append("Midter Reaktion: ")
	textlines.append("Højre Reaktion: ")
	textlines.append("Maksimalt Positivt Moment: ")
	textlines.append("Maksimalt Negativt Moment: ")
	textvals = [R1, R2, R3, np.amax(M), np.amin(M)]
	textlength = max([len(i) for i in textlines])
	textlength += int(textlength/5)
	for (i,j) in enumerate(textlines):
		if textvals[i] < 0:
			c = 1
		else: 
			c = 0
		cc = 0
		while abs(textvals[i]) > 10**(cc):
			cc += 1
		cc -= 1
		if abs(textvals[i]) < 1:
			cc = 0
		while len(textlines[i]) < textlength - c - cc:
			textlines[i] += "."
		textlines[i] += " "
	for (i,j) in enumerate(textlines):
		textlines[i] += "%.3f" %(textvals[i])
	textunits = [" kN", " kN", " kN", " kNm", " kNm"]
	for (i,j) in enumerate(textlines):
		textlines[i] += textunits[i]
	for (i,j) in enumerate(textlines):
		plt.text(leftmargin, 1 - topmargin - i*linemargin, j,
			family = 'monospace')