import bpy
import cv2
import numpy as np
import math
import pickle
import os

def create_helix_movement(points=100):
    global p_insert
    angle = np.linspace(0,7*np.pi,points)
    z=angle/4+2
    x=np.sin(angle)*2
    y=np.cos(angle)*2

    p0=p_insert
    positions=list()
    for i in range(points):
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

    pickle.dump(positions,open("C:/positions/helix.pk","wb")) 
    return positions

def create_square_movement(points=100):
    global p_insert
    positions=list()
    x=np.zeros(points)
    y=np.zeros(points)
    z=np.zeros(points)
    p_0=[-2.5,-2.5,3]
    p_1=[2.5,-2.5,4]
    p_2=[2.5,2.5,5]
    p_3=[-2.5,2.5,6]
    n=points
    for i in range(n//4):
        x[i]=p_0[0]+(p_1[0]-p_0[0])/(n//4)*i
        y[i]=p_0[1]+(p_1[1]-p_0[1])/(n//4)*i
        z[i]=p_0[2]+(p_1[2]-p_0[2])/(n//4)*i

    for i in range(n//4):
        x[i+n//4]=p_1[0]+(p_2[0]-p_1[0])/(n//4)*i
        y[i+n//4]=p_1[1]+(p_2[1]-p_1[1])/(n//4)*i
        z[i+n//4]=p_1[2]+(p_2[2]-p_1[2])/(n//4)*i

    for i in range(n//4):
        x[i+2*n//4]=p_2[0]+(p_3[0]-p_2[0])/(n//4)*i
        y[i+2*n//4]=p_2[1]+(p_3[1]-p_2[1])/(n//4)*i
        z[i+2*n//4]=p_2[2]+(p_3[2]-p_2[2])/(n//4)*i

    for i in range(n//4):
        x[i+2*n//4+n//4]=p_3[0]+(p_0[0]-p_3[0])/(n//4)*i
        y[i+2*n//4+n//4]=p_3[1]+(p_0[1]-p_3[1])/(n//4)*i
        z[i+2*n//4+n//4]=p_3[2]+(p_0[2]-p_3[2])/(n//4)*i

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

    pickle.dump(positions,open("C:/positions/square.pk","wb")) 
    return positions




def write_txt_output(predicted_positions,real_positions,folder="C:/tmp"):
    #Real positions
    if os.path.exists(folder+"/real_positions.txt"):
      os.remove(folder+"/real_positions.txt")


    f = open(folder+"/real_positions.txt", "a")
    for real_position in real_positions:
        f.writelines(str(real_position[0])+","+str(real_position[1])+","+str(real_position[2])+","+str(real_position[3])+","+str(real_position[4])+","+str(real_position[5])+"\n")
    f.close()

    #Predicted positions
    if os.path.exists(folder+"/prediction.txt"):
      os.remove(folder+"/prediction.txt")


    f = open(folder+"/prediction.txt", "a")
    for prediction in predicted_positions:
        f.writelines(str(prediction[0])+","+str(prediction[1])+","+str(prediction[2])+","+str(prediction[3])+","+str(prediction[4])+","+str(prediction[5])+"\n")
    f.close()

    #Difference
    if os.path.exists(folder+"/error.txt"):
      os.remove(folder+"/error.txt")

    f = open(folder+"/error.txt", "a")
    for i in range(len(predicted_positions)):
        f.writelines(str(abs(real_positions[i][0]-predicted_positions[i][0]))+","+str(abs(real_positions[i][1]-predicted_positions[i][1]))+","+str(abs(real_positions[i][2]-predicted_positions[i][2]))+","+str(abs(real_positions[i][3]-predicted_positions[i][3]))+","+str(abs(real_positions[i][4]-predicted_positions[i][4]))+","+str(abs(real_positions[i][5]-predicted_positions[i][5]))+"\n")
    f.close()


def get_vertex(elp,beta):

    a=(elp[1][1]) #major axis
    b=(elp[1][0]) #minor axis
    
    rotation=elp[2] #angle between the major axis and the x-axis

    #angle=alpha as introduced in the paper
    angle=math.asin(math.sqrt(1-((b/2)*(b/2))/((a/2)*(a/2)))*math.cos(math.radians(beta/2)))
    angle=angle*180.0/math.pi
    

    x=(a/2)*math.sin(math.radians(2*angle))/(math.sin(math.radians(beta))) #(s1=x)
    h=(a/2)*(math.cos(math.radians(2*angle))/math.sin(math.radians(beta))+1/math.tan(math.radians(beta))) #(s2=h)
    
    #cordinates:
    x_center=elp[0][0]
    y_center=elp[0][1]
    #two possible vertex positions (symmetry of the ellipse)
    x_vertex=x_center+math.sin(math.radians(elp[2]))*x
    y_vertex=y_center+math.cos(math.radians(elp[2]))*x
    x_vertex2=x_center-math.sin(math.radians(elp[2]))*x
    y_vertex2=y_center-math.cos(math.radians(elp[2]))*x
    
    return h,x_vertex,y_vertex,x_vertex2,y_vertex2


def heigth_based_on_distance_to_center(y1):
    global f,d0,r

    if y1==0:
        real_x1=0
    else:
        a=1+(f*f)/(y1*y1)
        b=(2*d0*(-1*f)/y1)-2*r*(-1*f)/y1
        c=d0*d0-2*d0*r
        real_x1=(-b+math.sqrt(b*b-4*a*c))/(2*a)
    return real_x1


def contour_to_position(contour):
    global f,d0,r,pixels

    num_points=len(contour)
    x_e=np.zeros(num_points)
    y_e=np.zeros(num_points)
    z_e=np.zeros(num_points)

    i=0
    #reconstruct shape based on the surface
    for coord in contour:
        point_x=coord[0][0]
        point_y=coord[0][1]
        x_from_center=point_x-(pixels-1)/2
        y_from_center=point_y-(pixels-1)/2
        distance_to_center_in_pixels=math.sqrt(x_from_center*x_from_center+y_from_center*y_from_center)
        distance_to_center_in_mm=(sensor_size/pixels)*distance_to_center_in_pixels
        real_x1=heigth_based_on_distance_to_center(distance_to_center_in_mm)

        real_x1_on_sensor_in_mm=real_x1/d0*f
        real_x1_on_sensor_in_pixel=real_x1_on_sensor_in_mm/sensor_size*pixels
        
        if x_from_center==0:
            x_new=0
        else:
            x_new=real_x1_on_sensor_in_pixel/math.sqrt(1+(y_from_center*y_from_center)/(x_from_center*x_from_center))
        
        if x_from_center<0:
            x_new=0-x_new
        
        if x_from_center==0:
            if y_from_center>=0:
                y_new=real_x1_on_sensor_in_pixel
            else:
                y_new=-real_x1_on_sensor_in_pixel
        else:
            y_new=y_from_center/x_from_center*x_new
        

        coord[0][0]=(pixels-1)/2+x_new
        coord[0][1]=(pixels-1)/2+y_new

        height_in_mm=0.0125-math.sqrt(0.0125*0.0125-real_x1*real_x1)
        heigth_in_px=((height_in_mm/d0*f)/sensor_size)*pixels

        x_e[i]=coord[0][0]
        y_e[i]=coord[0][1]
        z_e[i]=heigth_in_px

        i=i+1


    #get x and y center of the ellipse
    min_x=np.amin(x_e)
    max_x=np.amax(x_e)
    min_y=np.amin(y_e)
    max_y=np.amax(y_e)
    min_z=np.amin(z_e)
    max_z=np.amax(z_e)

    z_offset=(max_z-min_z)/2+min_z
    z_offset=z_offset/pixels*sensor_size/f*d0
    
    ellipse_center_x=(max_x+min_x)/2
    ellipse_center_y=(max_y+min_y)/2
    
    x_from_center=ellipse_center_x-(pixels-1)/2
    y_from_center=ellipse_center_y-(pixels-1)/2

    distance_to_center_in_pixels=math.sqrt(x_from_center*x_from_center+y_from_center*y_from_center)
    distance_to_center_in_mm=(sensor_size/pixels)*distance_to_center_in_pixels	
    real_x1=heigth_based_on_distance_to_center(distance_to_center_in_mm)
    
    #angle of a tangent plane at the position of the ellipse center
    angle_sphere=math.acos(math.sqrt(real_x1*real_x1)/0.0125)
    angle_sphere=90-angle_sphere*180.0/math.pi
    
    #move ellipse into origin:
    x_c=np.zeros(num_points)
    y_c=np.zeros(num_points)
    z_c=np.zeros(num_points)

    for i in range(num_points):
        x_c[i]=x_e[i]-ellipse_center_x
        y_c[i]=y_e[i]-ellipse_center_y
        z_c[i]=z_e[i]-min_z-(max_z-min_z)/2

    #rotate ellipse around Z-axis:
    if x_from_center==0:
        if y_from_center<0:
            angle_around_z=-90
        else:
            angle_around_z=90
    else:
        angle_around_z=math.atan(y_from_center/x_from_center)
        angle_around_z=-angle_around_z*180.0/math.pi
    
    rotation_matrix_z=[math.cos(math.radians(angle_around_z)),-math.sin(math.radians(angle_around_z)),0,math.sin(math.radians(angle_around_z)),math.cos(math.radians(angle_around_z)),0,0,0,1]
    x_cr=np.zeros(num_points)
    y_cr=np.zeros(num_points)
    z_cr=np.zeros(num_points)

    for i in range(num_points):
        x_cr[i]=rotation_matrix_z[0]*x_c[i]+rotation_matrix_z[1]*y_c[i]+rotation_matrix_z[2]*z_c[i]
        y_cr[i]=rotation_matrix_z[3]*x_c[i]+rotation_matrix_z[4]*y_c[i]+rotation_matrix_z[5]*z_c[i]
        z_cr[i]=rotation_matrix_z[6]*x_c[i]+rotation_matrix_z[7]*y_c[i]+rotation_matrix_z[8]*z_c[i]
    
    #rotate ellipse around Y-axis:
    if ellipse_center_x<(pixels-1)/2:
        angle_y=-angle_sphere
    else:
        angle_y=angle_sphere

    rotation_matrix=[math.cos(math.radians(angle_y)),0,math.sin(math.radians(angle_y)),0,1,0,-math.sin(math.radians(angle_y)),0,math.cos(math.radians(angle_y))]
    x_cr2=np.zeros(num_points)
    y_cr2=np.zeros(num_points)
    z_cr2=np.zeros(num_points)

    for i in range(num_points):
        x_cr2[i]=rotation_matrix[0]*x_cr[i]+rotation_matrix[1]*y_cr[i]+rotation_matrix[2]*z_cr[i]
        y_cr2[i]=rotation_matrix[3]*x_cr[i]+rotation_matrix[4]*y_cr[i]+rotation_matrix[5]*z_cr[i]
        z_cr2[i]=rotation_matrix[6]*x_cr[i]+rotation_matrix[7]*y_cr[i]+rotation_matrix[8]*z_cr[i]

        contour[i][0][0]=x_cr2[i]+ellipse_center_x
        contour[i][0][1]=y_cr2[i]+ellipse_center_y
    

    #fit ellipse into reuslting shape
    elps = cv2.fitEllipse(contour)
    h,x_vertex,y_vertex,x_vertex2,y_vertex2=get_vertex(elps,spot_light_angle)	

    #move vertex into origin(x,y)	
    x_vertex_2=x_vertex2-ellipse_center_x
    y_vertex_2=y_vertex-ellipse_center_y		
    x_vertex_4=x_vertex-ellipse_center_x
    y_vertex_4=y_vertex2-ellipse_center_y 

    z_vertex=h

    #perform inverse rotation around the y-axis
    angle_y=-angle_y

    rotation_matrix=[math.cos(math.radians(angle_y)),0,math.sin(math.radians(angle_y)),0,1,0,-math.sin(math.radians(angle_y)),0,math.cos(math.radians(angle_y))]

    x_vertex_r2=rotation_matrix[0]*x_vertex_2+rotation_matrix[1]*y_vertex_2+rotation_matrix[2]*z_vertex
    y_vertex_r2=rotation_matrix[3]*x_vertex_2+rotation_matrix[4]*y_vertex_2+rotation_matrix[5]*z_vertex
    z_vertex_r2=rotation_matrix[6]*x_vertex_2+rotation_matrix[7]*y_vertex_2+rotation_matrix[8]*z_vertex

    x_vertex_r4=rotation_matrix[0]*x_vertex_4+rotation_matrix[1]*y_vertex_4+rotation_matrix[2]*z_vertex
    y_vertex_r4=rotation_matrix[3]*x_vertex_4+rotation_matrix[4]*y_vertex_4+rotation_matrix[5]*z_vertex
    z_vertex_r4=rotation_matrix[6]*x_vertex_4+rotation_matrix[7]*y_vertex_4+rotation_matrix[8]*z_vertex

    #perform inverse rotation around the z-axis
    angle_around_z=-angle_around_z
    rotation_matrix_z=[math.cos(math.radians(angle_around_z)),-math.sin(math.radians(angle_around_z)),0,math.sin(math.radians(angle_around_z)),math.cos(math.radians(angle_around_z)),0,0,0,1]
    
    x_vertex_r22=rotation_matrix_z[0]*x_vertex_r2+rotation_matrix_z[1]*y_vertex_r2+rotation_matrix_z[2]*z_vertex_r2
    y_vertex_r22=rotation_matrix_z[3]*x_vertex_r2+rotation_matrix_z[4]*y_vertex_r2+rotation_matrix_z[5]*z_vertex_r2
    z_vertex_r22=rotation_matrix_z[6]*x_vertex_r2+rotation_matrix_z[7]*y_vertex_r2+rotation_matrix_z[8]*z_vertex_r2

    x_vertex_r24=rotation_matrix_z[0]*x_vertex_r4+rotation_matrix_z[1]*y_vertex_r4+rotation_matrix_z[2]*z_vertex_r4
    y_vertex_r24=rotation_matrix_z[3]*x_vertex_r4+rotation_matrix_z[4]*y_vertex_r4+rotation_matrix_z[5]*z_vertex_r4
    z_vertex_r24=rotation_matrix_z[6]*x_vertex_r4+rotation_matrix_z[7]*y_vertex_r4+rotation_matrix_z[8]*z_vertex_r4

    #undo translations
    reconstructed_x2=ellipse_center_x+x_vertex_r22-(pixels-1)/2
    reconstructed_y2=ellipse_center_y+y_vertex_r22-(pixels-1)/2
    reconstructed_z2=z_vertex_r22
    
    reconstructed_x4=ellipse_center_x+x_vertex_r24-(pixels-1)/2
    reconstructed_y4=ellipse_center_y+y_vertex_r24-(pixels-1)/2
    reconstructed_z4=z_vertex_r24 

   
    print("reconstructed_position1:",1000*reconstructed_x2/pixels*sensor_size/f*d0,-1000*reconstructed_y2/pixels*sensor_size/f*d0,1000*(reconstructed_z2/pixels*sensor_size/f*d0) +1000*z_offset)
    print("reconstructed_position2:",1000*reconstructed_x4/pixels*sensor_size/f*d0,-1000*reconstructed_y4/pixels*sensor_size/f*d0,1000*(reconstructed_z4/pixels*sensor_size/f*d0) +1000*z_offset)
    
    #calcluate positions of the two vertex positions [mm] 
    x_1=reconstructed_x2/pixels*sensor_size/f*d0
    y_1=-reconstructed_y2/pixels*sensor_size/f*d0
    z_1=(reconstructed_z2/pixels*sensor_size/f*d0+z_offset)
    x_2=reconstructed_x4/pixels*sensor_size/f*d0
    y_2=-reconstructed_y4/pixels*sensor_size/f*d0
    z_2=(reconstructed_z4/pixels*sensor_size/f*d0+z_offset)
    
    return x_1,y_1,z_1,x_2,y_2,z_2

def get_angles(p0,p1):
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
    return angle_1,0,angle_z

def test_scenario(spotlight1,spotlight2,spotlight3,x=0,y=0,z=5,angle_1=0,angle_2=0,angle_3=0,img=0,single=True):
    global p_insert
    #output location
    file="C:/tmp/temp"+str(img)+".png"

    #adjust position and rotation of the spotlight(s)
    spotlight1.location=(x,y,z)
    
    if single:
        spotlight1.rotation_euler=(math.radians(angle_1),math.radians(angle_2),math.radians(angle_3))
        spotlight2.location=(x+10,y,z)
        spotlight3.location=(x+10,y,z)
    else:
        spotlight1.rotation_euler=(math.radians(angle_1-10),math.radians(angle_2-6),math.radians(angle_3))
        spotlight2.location=(x,y,z)
        spotlight3.location=(x,y,z)
    
    spotlight2.rotation_euler=(math.radians(angle_1),math.radians(angle_2+11),math.radians(angle_3))
    spotlight3.rotation_euler=(math.radians(angle_1+9),math.radians(angle_2-6),math.radians(angle_3))
    
    #render current scene
    bpy.data.scenes["Scene"].render.filepath=file
    bpy.ops.render.render(write_still=True)
    
    z_offset=0
        
    #load rendered image
    img = cv2.imread(file,0)
    img = cv2.medianBlur(img,9)
    #img = cv2.GaussianBlur(img,(9,9),0)
    #ret,img = cv2.threshold(img,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    #ret, img = cv2.threshold(img, 170, 255, cv2.THRESH_BINARY)
    ret, img = cv2.threshold(img, 50, 255, cv2.THRESH_BINARY)

    #detect contours in the image
    contours, hierarchy = cv2.findContours(img, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    
    #retrieve vertex positions for each complex contour (length > 5) 
    positions=list()
    for contour in contours:
        if len(contour)>5:
            x1,y1,z1,x2,y2,z2=contour_to_position(contour)

            #choose point closest to the real position
            dist_point1=math.sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))
            dist_point2=math.sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2))

            if dist_point1<dist_point2:
                positions.append((x1,y1,z1))
            else:
                positions.append((x2,y2,z2))
                
                
                
    
    #calcluate average position based on all output positions
    avg_x=0
    avg_y=0
    avg_z=0
    for position in positions:
        avg_x+=position[0]
        avg_y+=position[1]
        avg_z+=position[2]
    avg_x=avg_x/ len(positions)
    avg_y=avg_y/ len(positions)
    avg_z=avg_z/ len(positions)
    dist=math.sqrt((x-avg_x) *(x-avg_x) + (y-avg_y) *(y-avg_y) + (z-avg_z) *(z-avg_z))

    #return the single position
    if single:
        return avg_x,avg_y,avg_z,dist*1000,[]
    
    #version with 3 spotlights
    #calcluate distance to the average positon for each resulting position
    dist_0=math.sqrt((positions[0][0]-avg_x)*(positions[0][0]-avg_x)+(positions[0][1]-avg_y)*(positions[0][1]-avg_y)+(positions[0][2]-avg_z)*(positions[0][2]-avg_z))
    dist_1=math.sqrt((positions[1][0]-avg_x)*(positions[1][0]-avg_x)+(positions[1][1]-avg_y)*(positions[1][1]-avg_y)+(positions[1][2]-avg_z)*(positions[1][2]-avg_z))
    dist_2=math.sqrt((positions[2][0]-avg_x)*(positions[2][0]-avg_x)+(positions[2][1]-avg_y)*(positions[2][1]-avg_y)+(positions[2][2]-avg_z)*(positions[2][2]-avg_z))
    
    multiple_positions=[positions[0][0],positions[0][1],positions[0][2],positions[1][0],positions[1][1],positions[1][2],positions[2][0],positions[2][1],positions[2][2]]
    
    #return position closest to the average position 
    if dist_0 <= dist_1 and dist_0 <= dist_2:
        return positions[0][0],positions[0][1],positions[0][2],1000*math.sqrt((positions[0][0]-x)*(positions[0][0]-x)+(positions[0][1]-y)*(positions[0][1]-y)+(positions[0][2]-z)*(positions[0][2]-z)), multiple_positions
    if dist_1 <= dist_0 and dist_1 <= dist_2:
        return positions[1][0],positions[1][1],positions[1][2],1000*math.sqrt((positions[1][0]-x)*(positions[1][0]-x)+(positions[1][1]-y)*(positions[1][1]-y)+(positions[1][2]-z)*(positions[1][2]-z)),multiple_positions
    else:
        return positions[2][0],positions[2][1],positions[2][2],1000*math.sqrt((positions[2][0]-x)*(positions[2][0]-x)+(positions[2][1]-y)*(positions[2][1]-y)+(positions[2][2]-z)*(positions[2][2]-z)),multiple_positions
    

def run_test(spotlight1,spotlight2,spotlight3,real_positions,single,output_txt=True,output_pickle=True):
    global p_insert
    prediction=list()
    errors=list()
    multiple_predictions=list()

    i=0
    print(len(real_positions))
    for position in real_positions:
        x=position[0]
        y=position[1]
        z=position[2]
        angle1=(position[4])
        angle2=(position[3])
        angle3=(position[5])

        x_p,y_p,z_p,error,multiple_positions=test_scenario(spotlight1,spotlight2,spotlight3,x=x/1000,y=y/1000,z=z/1000,angle_1=angle1,angle_2=angle2,angle_3=angle3,img=i,single=single)
        p0=p_insert
        p1=[1000*x_p,1000*y_p,1000*z_p]
        angle_y,angle_x,angle_z=get_angles(p0,p1)
        prediction.append((1000*x_p,1000*y_p,1000*z_p,angle_y,angle_x,angle_z))
        errors.append(error)
        i+=1
    
    if output_pickle:
        pickle.dump(prediction,open("C:/tmp/predicted_positions.pk","wb"))  
        pickle.dump(errors,open("C:/tmp/errors.pk","wb"))  
        if not single:
            pickle.dump(multiple_predictions,open("C:/tmp/multiple_predictions.pk","wb"))  

    if output_txt:
        write_txt_output(prediction,real_positions)

    print(sum(errors)/len(errors))
    print(max(errors))
    
      

#Start Simulation
spot_light_angle=18 #degree
spot_light_blur=1  #0=none 1=max
pixels=5000
f=300 #focal length in mm
d0=0.11 #distance of focal point to surface in cm
sensor_size=36 #in mm
p_insert=[6,0,22]
r=0.0125
bpy.context.scene.camera = bpy.data.objects["Camera"]
bpy.context.scene.render.resolution_x=pixels
bpy.context.scene.render.resolution_y=pixels
bpy.context.scene.camera.data.lens=f
bpy.context.scene.camera.location[0]=0
bpy.context.scene.camera.location[1]=0
bpy.context.scene.camera.location[2]=d0
bpy.context.scene.camera.data.sensor_width=sensor_size

spotlight1 = bpy.data.objects["Spot1"]
spotlight1.data.spot_size=math.radians(spot_light_angle)
spotlight1.data.spot_blend=spot_light_blur

spotlight2 = bpy.data.objects["Spot2"]
spotlight2.data.spot_size=math.radians(spot_light_angle)
spotlight2.data.spot_blend=spot_light_blur


spotlight3 = bpy.data.objects["Spot3"]
spotlight3.data.spot_size=math.radians(spot_light_angle)
spotlight3.data.spot_blend=spot_light_blur

#real_positions = create_helix_movement()
real_positions = create_square_movement()
#real_positions = pickle.load(open("C:/positions/helix.pk","rb"))
#real_positions = pickle.load(open("C:/positions/square.pk","rb"))

#single: [True=single spotlight,False=three spotlights]
run_test(spotlight1,spotlight2,spotlight3,real_positions,single=False,output_pickle=False,output_txt=True)

