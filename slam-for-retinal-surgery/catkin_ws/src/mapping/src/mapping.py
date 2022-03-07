#!/usr/bin/env python3

import rospy
#from std_msgs.msg import String, Int16
from visualization_msgs.msg import Marker
import cv2
import struct
from sensor_msgs import point_cloud2
from sensor_msgs.msg import PointCloud2, PointField
from std_msgs.msg import Header
from numpy import loadtxt
import numpy as np


class Mapping():
    def __init__(self):
        rospy.init_node('mapping', anonymous=True)
        rate = rospy.Rate(20)

        # Publish msg to rviz
        self.pub_pc2 = rospy.Publisher('/pcl', PointCloud2, queue_size=2)
        self.pub_marker = rospy.Publisher('/marker', Marker, queue_size=5)

        rospy.sleep(.3)

        self.marker = Marker()
        self.sensor_size = 36 # in mm
        self.pixels = 5000
        self.mode = 1  #0 for center, 1 for center_dist
        self.scale = 100 # for better visualization in rviz
        self.radius = 12.5 # in mm


        self.create_marker()
        self.pc2 = self.create_pc2()
        #self.out_errors()
        #self.read_text()

    def create_marker(self):
        #self.marker = Marker()
        self.marker.header.frame_id = "my_frame"
        self.marker.header.stamp = rospy.get_rostime()
        self.marker.ns = "window"
        self.marker.action = self.marker.ADD
        self.marker.pose.orientation.w = 1.0
        self.marker.id = 1
        self.marker.type = self.marker.SPHERE
        self.marker.scale.x = 0.025*self.scale
        self.marker.scale.y = 0.025*self.scale
        self.marker.scale.z = 0.025*self.scale
        self.marker.pose.position.x = 0
        self.marker.pose.position.y = 0
        self.marker.pose.position.z = 0.0125*self.scale
        self.marker.color.a = 0.5
        self.marker.color.r = 1.0
        self.marker.color.g = 1.0
        self.marker.color.b = 1.0

    def create_pc2(self):

        x,y,z = self.read_text()
        self.out_errors(x, y, z)
        x = x/1000*self.scale # in ros default system is m
        y = y/1000*self.scale
        z = z/1000*self.scale

        r,g,b = self.extract_rgb(x, y)
        #print(g)
        points = []
        #print(color)

        for i in range(len(x)):

            a = 255
            rgb = struct.unpack('I', struct.pack('BBBB', b[i], g[i], r[i], a))[0]
            #print(hex(rgb))
            pt = [x[i], y[i], z[i], rgb]
            points.append(pt)

        fields = [PointField('x', 0, PointField.FLOAT32, 1),
                  PointField('y', 4, PointField.FLOAT32, 1),
                  PointField('z', 8, PointField.FLOAT32, 1),
                  # PointField('rgb', 12, PointField.UINT32, 1),
                  PointField('rgba', 12, PointField.UINT32, 1),
                  ]

        header = Header()
        header.frame_id = "my_frame"
        pc2 = point_cloud2.create_cloud(header, fields, points)
        pc2.header.stamp = rospy.Time.now()

        return pc2
    ### read the points' coordinates from text ###
    ### input:
    ### output: arrays x, y, z
    def read_text(self):

        file = ""
        if self.mode == 1:
            file = "/home/jingsong/Documents/SA/dataset/results_shit/center_dist.txt"
        if self.mode == 0:
            file = "/home/jingsong/Documents/SA/dataset/results_shit/center.txt"
        lines = loadtxt(file)
        #print(lines.shape)
        x = lines[:,0]
        y = lines[:,1]
        z = lines[:,2]
        if self.mode == 0:
            z = z
            x = x - ((5000.0 - 1) / 2)
            y = y - ((5000.0 - 1) / 2)

            z = 1000 * z / 5000 * 36 / 300 * 0.11
            x = 1000 * x / 5000 * 36 / 300 * 0.11
            y = 1000 * y / 5000 * 36 / 300 * 0.11

        return x, y, z

    ### extract the rgb values of the points in img ###
    ### input: x, y, z Type:array
    ### output: r, g, b Type: array
    def extract_rgb(self, x, y):
        #according to line 388 in my_codes.py
        #test: (0,0) -> (2500,2500)
        #test: (18,18) -> (5000,5000)
        u = (x * self.pixels/self.sensor_size + (self.pixels-1)/2).astype(int)
        v = (y * self.pixels/self.sensor_size + (self.pixels-1)/2).astype(int)

        image = cv2.imread("/home/jingsong/Documents/SA/Code and resources/images/5000.png")
        color = (image[u,v]).astype(int) # color in order bgr
        #print(color.shape)

        return color[:,2], color[:,1], color[:,0]

    def out_errors(self,x,y,z):
        sum = x**2+y**2+(self.radius-z)**2
        sum = np.sqrt(sum)
        ave_error = np.average(sum-12.5)
        print("ave_error:" + str(ave_error))

if __name__ == '__main__':
    mapping = Mapping()
    rospy.sleep(0.2)
    mapping.pub_marker.publish(mapping.marker)
    mapping.pub_pc2.publish(mapping.pc2)
    rospy.spin()

