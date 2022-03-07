#include <ros/ros.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_conversions/pcl_conversions.h>
//#include <opencv2>

#include <visualization_msgs/Marker.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

typedef pcl::PointCloud<pcl::PointXYZRGB> PointCloud;

int main(int argc, char** argv)
{


  ros::init (argc, argv, "pub_pcl");
  ros::NodeHandle nh;
  ros::Publisher pub = nh.advertise<PointCloud> ("points2", 100);
  ros::Publisher marker_pub = nh.advertise<visualization_msgs::Marker>("visualization_marker", 40);


  PointCloud::Ptr msg (new PointCloud);
  msg->header.frame_id = "my_frame";
  msg->height = msg->width = 1;

  visualization_msgs::Marker marker;
  marker.header.frame_id= "my_frame";
  marker.header.stamp = ros::Time::now();
  marker.ns ="points_and_lines";
  marker.action = visualization_msgs::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.id = 1;
  marker.type = visualization_msgs::Marker::SPHERE;
  marker.scale.x = 0.025;
  marker.scale.y = 0.025;
  marker.scale.z = 0.025;
  marker.pose.position.x = 0;
  marker.pose.position.y = 0;
  marker.pose.position.z = 0.0125;
  marker.color.a = 0.8; // Don't forget to set the alpha!
  marker.color.r = 1.0;
  marker.color.g = 0.0;
  marker.color.b = 0.0;

// read the data from .txt
  int num_points = 100;

  int mode =0; //1 for center_dist; 0 for center, since center need additional process.
  
  std::string file_name;
  float x,y,z;
  pcl::PointXYZRGB p= pcl::PointXYZRGB(x,y,z);

  if(mode==0)
  {
	file_name = "/home/jingsong/Documents/SA_documents/dataset/results_shit/center.txt";
  } 
  if(mode==1)
  {
	file_name = "/home/jingsong/Documents/SA_documents/dataset/results_shit/center_dist.txt";
  }
  std::ifstream file(file_name.c_str());

  if(file.is_open())
  {
	for(int i=0;i<num_points;i++){

		file >> x>>y>>z;

		//data process funcs
		if(mode==0)
		{
			z = z;	
			x = x-((5000.0-1)/2);
			y = y-((5000.0-1)/2);
			
			z = 1000*z/5000*36/300*0.11;
			x = 1000*x/5000*36/300*0.11;
			y = 1000*y/5000*36/300*0.11;
		}
	// adpat for ROS system: meter as default
		x = x/1000;
		y = y/1000;
		z = z/1000;

		std::cout<<x<<" "<<y<<" "<<z;
		std::cout<<"\n";
	//assign xyz coordinates
		p.x=x;
		p.y=y;
		p.z=z;
	//assign rgb coordiantes(working)
		p.r = 0;
		p.g = 255;
		p.b = 0;
		//push back to pointcloud
  		msg->push_back (p);		
	}		
  }


  //ROS_INFO("%d",msg.points.size());

  ros::Rate loop_rate(1);
  while (nh.ok())
  {
    pcl_conversions::toPCL(ros::Time::now(), msg->header.stamp);
    pub.publish (msg);
    marker_pub.publish(marker);
    ros::spinOnce ();
    loop_rate.sleep ();
  }
}
