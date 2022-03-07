

#include <iostream>
#include <fstream>
#include <math.h> 
#include <opencv2>

int main(){
	float Points[100][3];
	float sum = 0;
	float ave_error=0.0;
	int mode = 1; //mode=0 for center, mode=1 for center_dist
	std::ifstream file("/home/jingsong/Documents/SA_documents/dataset/results_shit/center_dist.txt");
	if(file.is_open())
	{
		for(int i=0;i<100;i++){
			for(int j=0;j<3;j++)
			{
				file >>Points[i][j];
				//std::cout<<Points[i][j]<<" ";
			}
			//printf("\n");		
			if(mode==0){
			Points[i][2] = Points[i][2];	
			Points[i][0] = Points[i][0]-((5000.0-1)/2);
			Points[i][1] = Points[i][1]-((5000.0-1)/2);		
			}
		}		
	}
	//printf("%f\n",pow(12.5,2));

	for(int i=0;i<100;i++)
	{	
		if(mode==0){
		for(int j=0;j<3;j++)
		{
			
			Points[i][j] = 1000*Points[i][j]/5000*36/300*0.11;
					
		}
		}			

		sum=pow(Points[i][0],2)+pow(Points[i][1],2)+pow((12.5-Points[i][2]),2);
		sum = sqrt(sum);
		ave_error = ave_error + (12.5-sum);
		//printf("\n");
		sum=0;		
	}
	ave_error = ave_error / 100;
	std::cout<<ave_error;  //-0.02 for center
						   //-0.07 for center_dist
}
