import numpy as np
import matplotlib.pyplot as plt

global G
global R 
global M 
global m 
global w2 

G = 6.674*10**(-11) #m**3Kg**(-1)s**(-2)
R = 3.844*10**8 #m
M = 5.974*10**(24)
m = 7.348*10**(22) #kg
w2 = (2.66*10**(-6))**2*R**3/G #s**(-1)

def f(x):
	func = M/(x**2) - m/((1-x)**2) - w2*x
	return func
def df(x):
	func = -2*M/x**3 - 2*m/(1-x)**3 - w2
	return func
	
def Newton(x):	
	func = x - float(f(x))/df(x)
	return func
	

print('The Lagrange Point \n')

print('Equation we want to solve: G*M/r**2 - G*m/(R-r)**2 - w**2*r = 0')
print('Find r')

stop = str(1)
while(stop == str(1)):
	print('Choose one option:\n \t 1 - plot f(r) = G*M/r**2 - G*m/(R-r)**2 - w**2*r \n \t 2 - Calculate r using Newton method')
	choose = input('')

	if(choose == str(1)):
		x = np.linspace(-0.1,1.1,10**6)
		y = f(x)

		plt.plot(x,y)
		plt.xlabel('r/R')
		plt.ylabel('f(r)')
		plt.grid()
		#plt.ylim(-10**(16),10**30)
		plt.show()
		print('1 - Menu 2 - Exit')
		stop = input('')
		
	if(choose == str(2)):
		
		r = input('Choose starting r:')
		try:
			r = float(r)
		except:
			print('Not valid r')
			continue
		N = input('Choose N, number of steps:')
		try:
			N = int(N)
		except:
			print('Not valid N')
			continue
		for i in range(0,N):
			r_new = Newton(r)
			r = r_new
			delta = abs(f(r_new) - 0)
			if(delta < 0.00001):
				i = N

		
		print('r = {} Km'.format(int(r*R/1000)))
		print('Tabulated value = {} Km'.format(326390))

		x = np.linspace(r-0.1,r+0.1,10**6)
		y = f(x)

		plt.plot(x,y)
		plt.plot(r,f(r),'.',color = 'r')
		plt.xlabel('r/R')
		plt.ylabel('f(r/R)')
		#plt.ylim(-10**16,10**16)
		plt.grid()
		plt.show()

		print('1 - Menu 2 - Exit')
		stop = input('')

