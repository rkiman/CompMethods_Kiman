import numpy as np
import matplotlib.pyplot as plt

def Verlet(star1,star2,star3,t,h):
	'''
	This functions runs the Verlet algorithm and change the positions of 
	the stars according to the gravitational force
	Needs:
		Three object of the class star
		t time
		h step size
	'''
	
	for star in [star1,star2,star3]:
		star.r = star.r + h*star.v05
	
	k1 = h*f(star1,star2,star3)
	k2 = h*f(star2,star1,star3)
	k3 = h*f(star3,star2,star1)
	
	for star,k in zip([star1,star2,star3],[k1,k2,k3]):
		star.v = star.v05 + 0.5*k
		star.v05 = star.v05 + k
	return 0
	
def f(star1,star2,star3):
	G = 1
	r1 = star1.r
	m1 = star1.m
	r_int = [r1 - star2.r,r1 - star3.r]  #vector to calculate force in r1
	m_int = [star2.m, star3.m]
	fx,fy,fz = [0,0,0]
	for r,m in zip(r_int,m_int):
		x0 = r[0]
		y0 = r[1]
		z0 = r[2]
		r_norm = np.linalg.norm(r)
		fx = fx - G*m*x0/(r_norm**3)
		fy = fy - G*m*y0/(r_norm**3)
		fz = fz - G*m*z0/(r_norm**3)
	return np.array([fx,fy,fz],float)

class star():
	def __init__(self,r0,v0,v05,m0):
		self.r = r0
		self.v = v0
		self.v05 = v05
		self.m = m0
	
def calc_energy(star1,star2,star3):
	G = 1
	T = 0.5 * star1.m * np.linalg.norm(star1.v)**2 + 0.5 * star2.m * np.linalg.norm(star2.v)**2 + 0.5 * star3.m * np.linalg.norm(star3.v)**2
	V = 0
	for stari in [star1,star2,star3]:
		for starj in [star1,star2,star3]:
			if(stari.r[0] != starj.r[0]):
				V = V - 0.5*G*stari.m*starj.m/np.linalg.norm(stari.r - starj.r) 
	Ltot = np.linalg.norm(star1.m*np.cross(star1.r,star1.v) + star2.m*np.cross(star2.r,star2.v) + star3.m*np.cross(star3.r,star3.v))
	Ptot = np.linalg.norm(star1.m * star1.v + star2.m * star2.v + star3.m * star3.v)
	return T + V,T,V,Ltot,Ptot	

def define_step_size(r10,v10,r20,v20,r30,v30,m1,m2,m3):
	'''
	This functions uses the initial conditions to repeat the verlet integration
	making the step size smaller until the orbit converges to the solution.
	'''
	
	tfin = 0
	delta = 0
	while(tfin == 0):
		try:
			tfin = float(input('Define final time to check best step size: (1000.0 recomended)'))
		except:
			a = 0
	while(delta == 0):
		try:
			delta = float(input('Define delta to improve best step size: (10.0 recomended)'))
		except:
			a = 0

	star1 = star(r10,v10,v10,m1)
	star2 = star(r20,v20,v20,m2)
	star3 = star(r30,v30,v30,m3)

	h = 0.1
	t = 0
	dist1_points,dist2_points,t_points = [],[],[]

	while(t < tfin):
		Verlet(star1,star2,star3,t,h)
		t = t + h
		t_points.append(t)
		Xcm = (star1.m*star1.r + star2.m*star2.r)/(star1.m + star2.m)
		dist1,dist2 = [np.linalg.norm(star1.r - star2.r),np.linalg.norm(Xcm - star3.r)]
		dist1_points.append(dist1)
		dist2_points.append(dist2)

	for h in [0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005,0.00001]:
		
		star1 = star(r10,v10,v10,m1)
		star2 = star(r20,v20,v20,m2)
		star3 = star(r30,v30,v30,m3)
		
		t = 0
		delta1,delta2 = [0,0]
		delta1_points,delta2_points = [],[]
		dist1_replace,dist2_replace = [],[]
		i_compare = 0
		t_compare = t_points[i_compare]
		
		while(t < (tfin -0.1) ):
			Verlet(star1,star2,star3,t,h)
			t = t + h
			if(round(t,4) == round(t_compare,4)):
				Xcm = (star1.m*star1.r + star2.m*star2.r)/(star1.m + star2.m)
				dist1,dist2 = [np.linalg.norm(star1.r - star2.r),np.linalg.norm(Xcm - star3.r)]
				delta1,delta2 = [abs(dist1 - dist1_points[i_compare]),abs(dist2 - dist2_points[i_compare])]
				delta1_points.append(delta1)
				delta2_points.append(delta2)
				dist1_replace.append(dist1)
				dist2_replace.append(dist2)
				i_compare += 1
				t_compare = t_points[i_compare]
		dist1_points = dist1_replace
		dist2_points = dist2_replace
		if(all(x < delta for x in delta1_points) and all(y < delta for y in delta2_points)):
			print('Best step size: {}'.format(h))
			break
	return h

def set_Xcm_zero(star1,star2,star3):
	Xcm = (star1.m*star1.r + star2.m*star2.r + star3.m*star3.r)/(star1.m + star2.m + star3.m)
	star1.r -= Xcm
	star2.r -= Xcm
	star3.r -= Xcm

def plot(j):
	size = 2
	name = 'Results/Stability_Distance_'+ str(j) + '.txt'
	data = np.loadtxt(name)
	t = data[:,0]
	dist1 = data[:,1]
	dist2 = data[:,2]
	E = data[:,3]
	T = data[:,4]
	V = data[:,5]
	Ltot = data[:,6]
	Ptot = data[:,7]
	name2 = 'Results/'+str(j)+'_Distance.jpg'
	name3 = 'Results/'+str(j)+'_Energy.jpg'
	name4 = 'Results/'+str(j)+'_Momentum.jpg'
	
	fig1 = plt.figure()
	plt.plot(t,dist1,'r.',markersize=size,label = 'First binary')
	plt.plot(t,dist2,'g.',markersize=size,label = 'Second binary')
	plt.legend()
	plt.xlabel('t')
	plt.ylabel('Distance')
	fig1.savefig(name2)

	fig2 = plt.figure()
	plt.plot(t,E,'r.',markersize=size,label = 'Energy')
	plt.plot(t,V,'g.',markersize=size,label = 'Potential')
	plt.plot(t,T,'b.',markersize=size,label = 'Kinetic energy')
	plt.legend()
	plt.xlabel('t')
	fig2.savefig(name3)

	fig3 = plt.figure()
	plt.plot(t,Ptot,'r.',markersize=size,label = 'Total Linear Momentum')
	plt.plot(t,Ptot,'g.',markersize=size,label = 'Total Angular Momentum')
	plt.legend()
	plt.xlabel('t')
	fig3.savefig(name4)

def plot_trajectories():
	from mpl_toolkits.mplot3d import Axes3D

	data = np.loadtxt('Coordinates.txt')
	r1x = data[:,1]
	r1y = data[:,2]
	r1z = data[:,3]
	r2x = data[:,4]
	r2y = data[:,5]
	r2z = data[:,6]
	r3x = data[:,7]
	r3y = data[:,8]
	r3z = data[:,9]

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	ax.plot(r1x,r1y,r1z)
	ax.plot(r2x,r2y,r2z)
	ax.plot(r3x,r3y,r3z)
	plt.show()
