#intergral = 6.49394
import numpy as np
import matplotlib.pyplot as plt

def f(x):
	g = x**3/(np.exp(x)-1)
	return g


print('The Stefan-Boltzmann constant \n')

stop = 0

while(stop != 1):
	print("What do you want to do? \n \t 1 - Calculate Stefan-Boltzmann constant \n \t 2 - Exit")
	integration = float(input(""))
	while(integration != 1 and integration !=2):
		integration = float(input("Select 1 or 2: "))	
	while(integration == 1):
		kb = 1.38064852 * 10**(-23) #m^2 kg/(s^2 K)
		c = 299792458  #m/s
		h = 6.62607004 * 10**(-34)/(2*np.pi) # m^2 kg / s
		a = 0.00000000001
		b = 200
	
		print("Choose an integration method: \n \t 1 - Simpson's rule \n \t 2 - Gaussian quadrature \n \t 3 - Return to menu")
		select = float(input(""))
		while(select != 1 and select != 2 and select != 3):
			select = float(input("Select 1, 2 or 3: "))	
		if(select == 3):
			continue
		elif(select == 1):
			simpson = 1
			while(simpson == 1):
				#Simpson's Rule
				N = float(input("Choose number of slices: ")) #10000 #Has to be even
				while(N % 2 != 0 or N - int(N) != 0 or N <= 0 ):
					if(N % 2 != 0 and N - int(N) == 0):
						N = float(input("Please choose an even number: "))	
					if(N - int(N) != 0):
						N = float(input("Please choose an integer: "))
					if(N <= 0):
						N = float(input("Please choose a positive number: "))		
				delta = (b-a)/N
				F = delta/3*(f(a) + f(b) + 4*sum([f(a + (2*k - 1)*delta) for k in range(1,int(N/2)+1)]) \
				+ 2*sum([f(a + 2*k*delta) for k in range(1,int(N/2))]))

				sb = kb**4/(4*np.pi**2*c**2*h**3)*F
				
				delta = (b-a)/(N/2)
				F1 = delta/3*(f(a) + f(b) + 4*sum([f(a + (2*k - 1)*delta) for k in range(1,int(N/2)+1)]) \
				+ 2*sum([f(a + 2*k*delta) for k in range(1,int(N/2))]))
				error = 1/15*abs(F - F1)
				print('Stefan-Boltzmann constant: {:.3e} +/- {:.3e} kg s^(-3) K^(-4)'.format(sb,error))
				print('Tabulated value: {:.3e}  kg s^(-3) K^(-4)'.format(5.670*10**(-8)))
				
				print("What do you want to do now?: \n \t 1 - Repeat method \n \t 2 - Choose another method \n \t 3 - Return to menu")
				select = float(input(""))
				while(select != 1 and select != 2 and select != 3):
					select = float(input("Select 1, 2 or 3: "))	
				if(select == 1): 
					continue
				elif(select == 2):
					simpson = 0
					continue
				elif(select == 3):
					simpson = 0
					integration = 0 
					continue
		elif(select == 2):
			#Gaussian quadrature
			gaussian = 1
			while(gaussian == 1):
				from gaussxw_func import *
				order = float(input("Choose order for gaussian integration: "))
				while(order - int(order) != 0 or order <= 0):
					order = float(input("Please choose a positive integer: "))
				order = int(order)
				x, w = gaussxwab(order,a,b)

				F2 = sum([f(x0)*w0 for x0,w0 in zip(x,w)])

				sb = kb**4/(4*np.pi**2*c**2*h**3)*F2
			
				x, w = gaussxwab(order*2,a,b)
				F3 = sum([f(x0)*w0 for x0,w0 in zip(x,w)])
				error = abs(F2 - F3)
				print('Stefan-Boltzmann constant: {:.3e} +/- {:.3e} kg s^(-3) K^(-4)'.format(sb,error))
				print('Tabulated value: {:.3e} kg s^(-3) K^(-4)'.format(5.670 *10**(-8)))
				
				print("What do you want to do now?: \n \t 1 - Repeat method \n \t 2 - Choose another method \n \t 3 - Return to menu")
				select = float(input(""))
				while(select != 1 and select != 2 and select != 3):
					select = float(input("Select 1, 2 or 3: "))	
				if(select == 1): 
					continue
				elif(select == 2):
					gaussian = 0
					continue
				elif(select == 3):
					gaussian = 0
					integration = 0 
					continue
	if(integration == 2):
		stop = 1
			

			
			
