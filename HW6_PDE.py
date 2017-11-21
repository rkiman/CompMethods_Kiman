import numpy as np
import matplotlib.pyplot as plt
import sys

M = 100

#Boundary conditions
V_wall = 0.0
V_capacitor_l = 1.0
V_capacitor_r = -1.0

print('Calculate potential of a capacitor in an isolated box with a final precision of 10^6 volts, using relaxation method\ndelta:')
target = 1e-6 # Target accuracy

phi = np.zeros([M+1,M+1],float)
phi[20:80,20] = V_capacitor_l
phi[20:80,80] = V_capacitor_r
phiprime = np.empty([M+1,M+1],float)
# Main loop
delta = 1.0
while delta>target:
	for i in range(M+1):
		for j in range(M+1):
			if i==0 or i==M or j==0 or j==M:
				phiprime[i,j] = phi[i,j]
			elif(20 < i and i < 80 and j == 20):
				phiprime[i,j] = phi[i,j]
			elif(20 < i and i < 80 and j == 80):
				phiprime[i,j] = phi[i,j]
			else:
				phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])/4
	# Calculate maximum difference from old values
	delta = np.max(abs(phi-phiprime))
	# Swap the two arrays around
	phi,phiprime = phiprime,phi
	sys.stdout.write("\r%f" % delta)
	sys.stdout.flush()

sys.stdout.write("\r%f\n" % delta)

plt.imshow(phi)
plt.colorbar(label='Potential [Volts]')
plt.title('Density plot for the potential')
plt.xlabel('$x$ [mm]')
plt.ylabel('$y$ [mm]')
plt.show()

