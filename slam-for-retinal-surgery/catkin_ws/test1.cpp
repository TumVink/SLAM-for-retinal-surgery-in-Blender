

#include <iostream>
#include <fstream>
#include <math.h> 

int main(){
	float Points[100][3];
	float sum = 0;
	std::ifstream file("/home/jingsong/Documents/SA_documents/dataset/center_square_100.txt");
	if(file.is_open())
	{
		for(int i=0;i<100;i++){
			for(int j=0;j<3;j++)
			{
				file >>Points[i][j];
				//std::cout<<Points[i][j]<<" ";
				Points[i][j] = Points[i][j]/1000;
			}
			//printf("\n");	
			Points[i][2] = 12.5-Points[i][2];	
	
		}		
	}
	//printf("%f\n",pow(12.5,2));

	for(int i=0;i<100;i++)
	{	
		for(int j=0;j<3;j++)
		{
		sum += pow(Points[i][j],2);
		}			
		printf("%f",pow(12.5,2)-sum);
		printf("\n");
		sum=0;		
	}
}
