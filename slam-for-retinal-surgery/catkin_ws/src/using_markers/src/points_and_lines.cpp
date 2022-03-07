#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <iostream>
#include <fstream>
#include <cmath>

int main( int argc, char** argv )
{
  ros::init(argc, argv, "points_and_lines");
  ros::NodeHandle n;
  ros::Publisher marker_pub = n.advertise<visualization_msgs::Marker>("visualization_marker", 10);

  ros::Rate r(30);
 

  float f = 0.0;
  while (ros::ok())
  {

    visualization_msgs::Marker points,marker;
    points.header.frame_id =marker.header.frame_id= "/some_tf_frame";
    points.header.stamp =marker.header.stamp = ros::Time::now();
    points.ns  = marker.ns ="points_and_lines";
    points.action  =marker.action = visualization_msgs::Marker::ADD;
    points.pose.orientation.w =marker.pose.orientation.w = 1.0;



    points.id = 0;
    marker.id = 1;
    //line_list.id = 2;



    points.type = visualization_msgs::Marker::POINTS;
	marker.type = visualization_msgs::Marker::SPHERE;
    //line_strip.type = visualization_msgs::Marker::LINE_STRIP;
    //line_list.type = visualization_msgs::Marker::LINE_LIST;



    // POINTS markers use x and y scale for width/height respectively
    points.scale.x = 0.05;
    points.scale.y = 0.05;
	
	marker.scale.x = 25;
	marker.scale.y = 25;
	marker.scale.z = 25;
    // LINE_STRIP/LINE_LIST markers use only the x component of scale, for the line width
    //line_strip.scale.x = 0.1;
    //line_list.scale.x = 0.1;
	marker.pose.position.x = 0;
	marker.pose.position.y = 0;
	marker.pose.position.z = 12.5;


    // Points are green
    points.color.g = 1.0f;
    points.color.a = 1.0;

	marker.color.a = 0.2; // Don't forget to set the alpha!
	marker.color.r = 1.0;
	marker.color.g = 0.0;
	marker.color.b = 0.0;
	//marker.mesh_resource = "https://www.sim.informatik.tu-darmstadt.de/~kohlbrecher/ros_trac/sphere_mesh/sphere_32_16.dae";

    // Line strip is blue
    //line_strip.color.b = 1.0;
    //line_strip.color.a = 1.0;

    // Line list is red
    //line_list.color.r = 1.0;
    //line_list.color.a = 1.0;



    // Create the vertices for the points and lines

    //marker_pub.publish(points);
    marker_pub.publish(marker);
    //marker_pub.publish(line_list);

    r.sleep();

    //f += 0.04;
  }
}
