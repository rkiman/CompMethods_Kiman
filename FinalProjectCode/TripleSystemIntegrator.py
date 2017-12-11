from VerletCode import *
import numpy as np

#Load initial conditions
grid = np.loadtxt('grid.txt')

r01_points,v01_points,m1_points = [np.array(grid[:,0]), np.array(grid[:,1]),np.array(grid[:,2])]
r02_points,v02_points,m2_points = [np.array(grid[:,3]),np.array(grid[:,4]),np.array(grid[:,5])] 
r03_points,v03_points,m3_points = [np.array(grid[:,6]),np.array(grid[:,7]),np.array(grid[:,8])] 
stability = np.array(grid[:,9])

mask = (stability == 0)

#Main Loop
i=0
for r01i,v01i,m1i,r02i,v02i,m2i,r03i,v03i,m3i in zip(r01_points[mask],v01_points[mask],m1_points[mask],r02_points[mask],v02_points[mask],m2_points[mask],r03_points[mask],v03_points[mask],m3_points[mask]):
	r10,v10,m1 = [np.array([r01i,0,0]),np.array([0,v01i,0]),m1i]
	r20,v20,m2 = [np.array([r02i,0,0]),np.array([0,v02i,0]),m2i]
	r30,v30,m3 = [np.array([r03i,0,0]),np.array([v03i,0,0]),m3i] 
	
	tfin = 2000.0
	
	#h = define_step_size(r10,v10,r20,v20,r30,v30,m1,m2,m3)
	h = 0.001 #If you already run the define_step_size you can put it here.

	#Define the three stars of the system with the class star	
	star1 = star(r10,v10,v10,m1)
	star2 = star(r20,v20,v20,m2)
	star3 = star(r30,v30,v30,m3)

	set_Xcm_zero(star1,star2,star3)
	
	N = int(float(tfin)/h)
	t = 0
	
	name = 'Results/Stability_Distance_'+ str(i) + '.txt'
	file_coord = open(name,'a')
	file_coord.write('#t\tDist1\tDist2\tE\tT\tV\tL\tPtot\n' )

	for step in range(0,N):
		Verlet(star1,star2,star3,t,h)
		t = t + h
		#Calculate the new coordinates useful to analyze the system
		Xcm = (star1.m*star1.r + star2.m*star2.r)/(star1.m + star2.m)
		dist1,dist2 = [np.linalg.norm(star1.r - star2.r),np.linalg.norm(Xcm - star3.r)]
		E,T,V,Ltot,Ptot = calc_energy(star1,star2,star3)
		file_coord.write(str(t) + '\t' + str(dist1) + '\t' + str(dist2) + '\t' + str(E) + '\t' + str(T) + '\t' + str(V) + '\t' + str(Ltot) + '\t' + str(Ptot) + '\n')
		
	file_coord.close()
	
	plot(i)
	
	i += 1
