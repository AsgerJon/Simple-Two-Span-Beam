import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets

def basic_beam(L):
	plt.axis('equal')
	plt.plot((0, L), (0, 0), 'k', lw = 2)

def add_roller_support(x0, s=25):
	# Upper Triangle
	Dx = [x0, x0+2*s, x0-2*s, x0]
	Dy = [0, -3*s, -3*s, 0]
	plt.plot(Dx,Dy,'k')
	# Circular Rollers
	N = 24
	t = np.linspace(0, 2*np.pi, N)
	r = 0.9 * s
	xc1 = x0-s
	xc2 = x0+s
	yc = -4*s
	Dx1 = [xc1 + r*np.cos(t[i]) for i in range(N)]
	Dx2 = [xc2 + r*np.cos(t[i]) for i in range(N)]
	Dy = [yc + r*np.sin(t[i]) for i in range(N)]
	plt.plot(Dx1,Dy,'k')
	plt.plot(Dx2,Dy,'k')
	# Horizontal Line under Rollers
	plt.plot([x0-2*s, x0+2*s],[-5*s, -5*s],'k')
	# Lower Inclined Lines
	for i in range(5):
		plt.plot([x0+(2-i)*s, x0+(1-i)*s],[-5*s, -6*s],'k')
		
def add_pinned_support(x0, s=25):
	# Upper Triangle
	plt.axis('equal')
	Dx = [x0, x0+2*s, x0-2*s, x0]
	Dy = [0, -3*s, -3*s, 0]
	plt.plot(Dx,Dy,'k')
	# Lower Inclined Lines
	for i in range(5):
		plt.plot([x0+(2-i)*s, x0+(1-i)*s],[-3*s, -4*s],'k')
		
def add_distributed_load(x0, x1, s=25):
	# Draw Load Box
	Dx = [x0, x1, x1, x0, x0]
	Dy = [s/2, s/2, 7*s/2, 7*s/2, s/2]
	plt.plot(Dx, Dy, 'k')
	narrows = int((x1 - x0) // (2*s))
	arrow_dist = (x1 - x0) / (narrows + 1)
	for i in range(narrows):
		xa = x0 + arrow_dist*(i+1)
		plt.arrow(xa, 7*s/2, 0, -3*s, 
			width = 1,
			head_width = s,
			head_length = s,
			overhang = 0,
			head_starts_at_zero = False,
			length_includes_head = True,
			color = 'k',
			)
def add_concentrated_load(x0, s=25):
	plt.arrow(x0, 11*s/2, 0, -5*s, 
		width = 2,
		head_width = s,
		head_length = s,
		overhang = 0,
		head_starts_at_zero = False,
		length_includes_head = True,
		color = 'k',
		)
	