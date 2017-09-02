#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iterator>
using namespace std;

class Point
{
public:
  double x;  
  double y;
   
  friend istream& operator>>(istream& input, Point& p);
 
  double getSqX(void);
  
  double getSqY(void);
  
 double LengthSquared(void);
};

  double Point::getSqX(void){
	  return pow(x,2);
	  }

  double Point::getSqY(void){
	  return pow(y,2);
	  }
	  
 double Point::LengthSquared(){ return getSqX() + getSqY(); }

istream& operator>>(istream& input, Point& p)
{
  char c;
  input >> c; // Read open parenthesis
  input >> p.x;
  input >> c; // Read comma
  input >> p.y;
  input >> c; // Read closing parenthesis
  
  return input;
};



vector<vector<Point> > LoadFFT(string path){
	string Filename;   
    vector<vector<Point> > matrix;	// placeholder for the file
    Filename.append(path);                                                              // path of the directory

    Filename.append("FFTChannel1.txt"); 
	ifstream fileFFT(Filename.c_str());
	string raw_text;
	while(getline(fileFFT, raw_text)){
		
		vector<Point> row;
		istringstream strm(raw_text);
        Point p;
		while( strm >> p ){
			row.push_back(p);
			
		}
		
	matrix.push_back(row);
	}
		
return(matrix);
}

vector<double> LoadDistance(string path){
	//" This function loads the distance data"
	string Filename;
	double temp;
	vector<double> TrainingDistance;
	Filename.append(path);
	Filename.append("distance.txt");

	ifstream distOb(Filename.c_str());

	while(distOb.good()){
		distOb >> temp;
		TrainingDistance.push_back(temp);
	}
	distOb.close();

return TrainingDistance;
}

struct ReturnObject{
  double value1;
  double value2;
};

ReturnObject LinearRegressor(const vector<double> &TrainingX, const vector<double> &TrainingY){
  /* Stochastic Gradient descent algorithm is used to estimate the regression parameters*/
  float alpha=pow(10,-5),p,err,a;
  
  int m=TrainingX.size();                                       // Find the size of the vector
  float w=0;                                                   // initialize slope
  float b=0;                                                   // initialize bias term

  for(int i = 0 ; i < (TrainingX.size() * pow(10,6)); i++){                // The loop passes through each training data to run the update, it has 1 million epochs 
        int idx= i % m;                                       // This enables the code to go through each data file in the training set and update the weights based on data from  individual data observation

        p= b + w * TrainingX[idx];
        err = p - TrainingY[idx];
        b = b- alpha *err;
        w = w- alpha * err * TrainingX[idx] ;
        if(i % 10000 ==0){
         cout<<"W: "<<w<<">>>>>>>>>>>"<<"b: "<<b<<">>>>>>>>>>>"<<">>>>>>>>>>>"<<"error: "<<err<<endl;
       }
	   }

  ReturnObject result={w,b};

  return(result);
}


int main(){
	double temp,PowerDBtemp;
	vector<Point> tempVector;
	vector<vector<double> > PowerTotal;
	vector<double> power,PowerDB;

	vector<vector<Point> > FFTfile=LoadFFT("C:\\Users\\Spandan Mishra\\Documents\\GitHub\\GopherCppcodes\\");
	vector<double> dist=LoadDistance("C:\\Users\\Spandan Mishra\\Documents\\GitHub\\GopherCppcodes\\");
	
	for (int i = 0; i < FFTfile.size(); i++){
		tempVector= FFTfile[i];
		for (int j = 0; j < tempVector.size(); j++){
			 temp=tempVector[j].LengthSquared();													// caluclate power for each signal
			 power.push_back(temp);
		}
        PowerDBtemp = *max_element(power.begin(),power.end());	
        PowerDB.push_back(10*  log10(PowerDBtemp));	
		
		PowerTotal.push_back(power);																// saving all the power vector  in powerTotal
		power.clear();		
	}
	
	ReturnObject prediction= LinearRegressor(PowerDB, dist);										// Calls the gradient descent for regression
	/*
	for(vector<double>::iterator itr=PowerDB.begin(); itr!= PowerDB.end(); ++itr){
		cout<< *itr<<" ";
	}
	*/
	return(0);
}