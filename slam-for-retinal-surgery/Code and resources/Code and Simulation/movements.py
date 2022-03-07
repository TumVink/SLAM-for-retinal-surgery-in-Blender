import numpy as np
import math
import matplotlib.pyplot as plt
def create_shit_movement(points=100): #100
    global p_insert

    n=points
    
    #define the spiral orbit
    theta = np.linspace(start=0,stop=10*math.pi,num=100)
    r = 0+0.1*theta  #0.01 is the distance of neighbor spiral orbits
    
    x = np.multiply(r,np.cos(theta))
    y = np.multiply(r,np.sin(theta))
    z = np.linspace(start = 2, stop = 7,num = 100)
    
        
    p0=p_insert
    positions=list()
    for i in range(n):
        p1=[x[i],y[i],z[i]]
        diff=[p0[0]-p1[0],p0[1]-p1[1],p0[2]-p1[2]]
        if diff[0]==0:
            angle_z=90
        else:
            angle_z=math.atan(diff[1]/diff[0])
            angle_z=angle_z*180.0/math.pi
        if diff[0]<0:
            angle_z+=180

        angle_1=math.atan(math.sqrt((diff[0]*diff[0]+diff[1]*diff[1]))/diff[2])
        angle_1=angle_1*180.0/math.pi
        positions.append([x[i],y[i],z[i],angle_1,0,angle_z])

    #pickle.dump(positions,open("/home/jingsong/blender_file/Master Thesis-EyeSpotlight/square.pk","wb")) 
    return positions


p_insert=[6,0,22]
positions = create_shit_movement()

print(positions)