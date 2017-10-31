import numpy as np
import matplotlib.pyplot as plt
import time

def Verlet(r,v,t,h):
	'''
	needs r in t, v in t+h/2
	returns r in t+h, v in t + h, v in t+3/2*h
	'''
	r_new = r + h*v
	k = h*f(r_new)
	v_new = v + 0.5*k
	v_new2 = v + k
	return r_new,v_new,v_new2 

def f(r):
	x0 = r[0]
	y0 = r[1]
	r_norm = np.linalg.norm(r)
	x = -G*M*x0/(r_norm**2*np.sqrt(r_norm**2+L**2/4))
	y = -G*M*y0/(r_norm**2*np.sqrt(r_norm**2+L**2/4))
	return np.array([x,y],float)


global G
global M
global L
G = 1
M = 10  
L = 2

x = 1.0 
y = 0.0
vx = 0.0
vy = 1.0

h = 0.01

print('Space garbage\n')
cont = 0
while(cont != 20):
	while(cont != 1):
		print('Choose action:\n \t 1) Plot result \n \t 2) Movie \n \t 3) Exit')
		choise = input('')
		try:
			choise_num = int(choise)
			if(choise_num==1 or choise_num==2 or choise_num==3):
				cont = 1
		except:
			cont = 0
	if(choise_num == 2):  
		r = np.array([x,y],float)
		v = np.array([vx,vy],float)

		x_points = []
		y_points = []
		t_points = []
		t = 0.0
		tmax = 10.0

		N = int((tmax - t)/h)
		x_points.append(x)
		y_points.append(y)
		t_points.append(t)          
		plt.ion() # set plot to animated
		line, = plt.plot(x_points,y_points)
		plt.ylim(-1.0,1.0)
		plt.xlim(-1.0,1.0)
		plt.xlabel('x position')
		plt.ylabel('y position')
		plt.gca().set_aspect('equal', adjustable='box')
		for i in range(0,N):
			r_new,v_new,v_new2 = Verlet(r,v,t,h)
			r = r_new
			v = v_new2
			x_points.append(r_new[0])
			y_points.append(r_new[1])
			t_points.append(t+h)
			t = t + h
			line.set_xdata(x_points)
			line.set_ydata(y_points)  # update the data
			plt.draw() # update the plot
			plt.pause(0.000001)
		plt.ioff()
		plt.show()
		cont = 0
	elif(choise_num == 1):
		r = np.array([x,y],float)
		v = np.array([vx,vy],float)

		x_points = []
		y_points = []
		t_points = []
		t = 0.0
		tmax = 10.0

		N = int((tmax - t)/h)
		x_points.append(x)
		y_points.append(y)
		t_points.append(t)
		for i in range(0,N):
			r_new,v_new,v_new2 = Verlet(r,v,t,h)
			r = r_new
			v = v_new2
			x_points.append(r_new[0])
			y_points.append(r_new[1])
			t_points.append(t+h)
			t = t + h
		plt.plot(x_points,y_points)  
		plt.gca().set_aspect('equal', adjustable='box')  
		plt.xlabel('x position')
		plt.ylabel('y position')    
		plt.show()
		cont = 0
	elif(choise_num == 3):
		cont = 20
