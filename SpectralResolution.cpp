#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <complex>
#include <valarray>
#include <vector>
#include <algorithm>
#include <dirent.h>
#include <string>
#include <typeinfo>
#define PI 3.14159265359
using namespace std;
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
// Declare the variables outside the main function
float dataObject1[1000];
float dataObject2[1000];
float dataObject3[1000];
float dataObject4[1000];
float dataObject5[1000];
float dataObject6[1000];
float dataObject7[1000];
float dataObject8[1000];
// bandpass filtered obtained from Jeff
int hamming(float *hamming_coeff, int order)
{
  int i,nhlf;

  if(order%2 == 1)
  {
    nhlf = (order+1)/2;
  }
  nhlf = order/2;
  for(i = 0; i<nhlf; i++)
  {
    hamming_coeff[i] = (float)(0.54-0.46*cos((float)(i+1)*2*PI/order));
    if((order-i)>nhlf)
    {
      hamming_coeff[order-1-i] = hamming_coeff[i];
    }
  }
  return 0;
}

int bandpass_symmetric_filter(float *data, int n_data, int samprate, int center_freq, int order)
{
  float lfreq, ufreq, fl, fh, c1;
  float *filter_coeff, *hamming_wind, *temp_data,tmp;
  int nhlf,i,j,k;

  lfreq = (center_freq-center_freq*0.4)/(0.5*samprate);
  ufreq = (center_freq+center_freq*0.4)/(0.5*samprate);


  if(ufreq < lfreq || ufreq > 1 || lfreq < 0)
  {
    printf("error filter parameters are not valid\n");
    return -1;
  }

  //order should be odd to create a symmetric filter
  if(order%2 == 1)order = order+1;

  fl = lfreq/2;
  fh = ufreq/2;
  nhlf = (order+1)/2;
  c1 = fh-fl;

  filter_coeff= new float[order];
  hamming_wind= new float[order];
  temp_data   = new float[n_data];

  if(filter_coeff ==  NULL || hamming_wind == NULL || temp_data == NULL)
  {
    printf("error: memory allocation failed\n");
    return -1;
  }


  filter_coeff[nhlf-1] = c1;

  //create the coefficients
  for(i = 1; i<nhlf; i++)
  {
    filter_coeff[nhlf-1-i] = (float)(2*sin(PI*i*c1)*cos(PI*i*(fh+fl))/(PI*i));
    filter_coeff[nhlf+i-1] = filter_coeff[nhlf-1-i];
  }

  //create a hamming window
  hamming(hamming_wind,order);

  for(k=0; k<n_data; k++)
  {
    temp_data[k] = data[k];
  }

  for(k = 0; k < n_data; k++)
  {
    tmp = 0;
    for(i = 0; i<order; i++)
    {
      j = k-(nhlf-1)+i;
      if(j>=0 && j<n_data)
      {
        tmp+=filter_coeff[i]*hamming_wind[i]*temp_data[j];
      }
    }
    data[k] = tmp;
  }

  free(filter_coeff);
  free(hamming_wind);
  free(temp_data);
  //cout<< "Data Filtered"<< endl;
  return 0;

}
void LoadData(string path, string  filename)
{

  string Filename;                                                                    // placeholder for the file
  Filename.append(path);                                                              // path of the directory
  //cout<<Filename<< endl;
  Filename.append(filename);                                                           //name of the file
 //  cout<< Filename<<endl;

  ifstream file(Filename.c_str());                                              //c_str returns a const char* that points to a null-terminated string (i.e. a C-style string).
                                                                                //   It is useful when you want to pass the "contents"¹ of an std::string to a function that expects
                                                                                // to work with a C-style string.
 
  if(file.is_open())
  {
   

 // cout<<" File opened sucessfully!! now reading the data"<<endl;
    
  for(int i=0;i<8000;i++)          //8000 is the toal number of data present
    {                              // reads all data that is not end of file
      file >> dataObject1[i];      // column 1
      file >> dataObject2[i];      // column 2
      file >> dataObject3[i];      // column 3
      file >> dataObject4[i];      // column 4
      file >> dataObject5[i];      // column 5
      file >> dataObject6[i];
      file >> dataObject7[i];
      file >> dataObject8[i];

      i++;
    }
 }

  else{

   cout<< "Error !!! File couldn't be opened"<<endl;    
 }

file.close();

}



// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive

void fft(CArray& x)
{
    const size_t N = x.size();
    const double PI1 = 3.141592653589793238460;
    if (N <= 1) return;
 
    // divide
    CArray even = x[slice(0, N/2, 2)];                        //This class represents a valarray slice selector. It does not contain nor refers to any element - 
                                                              //it only describes a selection of elements to be used as an index in valarray::operator[].
    CArray  odd = x[slice(1, N/2, 2)];
 
    // conquer
    fft(even);
    fft(odd);
    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = polar(1.0, -2 * PI1 * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

vector<float> PowerCalculation(CArray &data){
  // This function calculates the power from FFT components

  vector<float> power;
  float temp1;
                for(int j=0;j < data.size(); j++)
               {
                  temp1=pow(real(data[j]),2) + pow(imag(data[j]),2);
                  power.push_back(temp1);
              }
              return(power);
}

// Lets open the names of all the files

vector<string> open1(string path){
  DIR *dir;
  dirent *pdir;
  vector<string> files;
  dir =opendir(path.c_str());
  while(pdir=readdir(dir))
{
 files.push_back(pdir->d_name);
}
return files;
}

double predict(double w, double b, float x){

    return(w * (x) +b);
 
}

// This is the return type of the object
struct ReturnObject{
  double value1;
  double value2;
};

ReturnObject LinearRegressor(const vector<float> &TrainingX, const vector<float> &TrainingY){
  /* Stochastic Gradient descent algorithm is used to estimate the regression parameters*/
  float alpha=pow(10,-5),p,err,a;
  //vector<float> yhat;
  // lets initialize the values of slope and bias term
  
  int m=TrainingX.size();                                       // Find the size of the vector
  float w=0;                                                   // initialize slope
  float b=0;                                                   // initialize bias term

  for(int i = 0 ; i < (63 * pow(10,6)); i++){                // The loop passes through each training data to run the update, it has 1 million epochs 
        int idx= i % m;                                       // This enables the code to go through each data file in the training set and update the weights based on data from  individual data observation

        p= b + w * TrainingX[idx];

        err = p - TrainingY[idx];

          //cout<<err;
          //cin>>a;

        b = b- alpha *err;

        w = w- alpha * err * TrainingX[idx] ;

        if(i % 10000 ==0){

         cout<<"W: "<<w<<">>>>>>>>>>>"<<"b: "<<b<<">>>>>>>>>>>"<<">>>>>>>>>>>"<<"error: "<<err<<endl;

       }

  }

  ReturnObject result={w,b};
  return(result);
}

float* diff(float *power_data,int n){

  // delete NewArray in client fucntion
  float* NewArray=new float[n-1];
  for(size_t i=0;i< n-1 ;i++){
    *(NewArray + i)= power_data[i+1] - power_data[i];     
  }
  return(NewArray);
 }
///////////////////////////
int main(){
 int n_data = 1000, samprate = 2048, center_freq = 90, order = 512;
                                                 
//float PowerDBtemp;
int dist;   
//float *dt1 , *dt3;                                                       // temporary folders
//vector<float> PowerDB, distance,predictedDistance;                                                           // storing the values of all the power
vector<string> Filesdistance,DirectoryNames;
string distance1File="50m\\";              // path to folder 1
DirectoryNames.push_back(distance1File);
string distance2File="100m\\";             // path to foldeer 2
DirectoryNames.push_back(distance2File);
string distance3File="150m\\";             // path to folder 3
DirectoryNames.push_back(distance3File);
vector<string> ListString;         
vector<float> power1;  

string token;
string pathName;                                                                                      // Will store the path to directory
string delimiter= "\\";                                                               // this will be used to stored parsed strings

for(vector<string>::iterator tempitr=DirectoryNames.begin(); tempitr !=DirectoryNames.end();++tempitr){
	
    string tempdirectory= *(tempitr);	
	size_t pos=0;
	while((pos=tempdirectory.find(delimiter)) != string::npos){
		token = tempdirectory.substr(0, pos);                                                         // components that make the directory filename
        ListString.push_back(token);                                                                   // we will store all theese components in a directory
        tempdirectory.erase(0, pos + delimiter.length());
		}
		 
	string Flag= ListString.back();                                                         // the last element  of the vector

    if(Flag=="50m"){
		pathName="50m\\"; 
        Filesdistance= open1(pathName);                                              // collect the names of files at distance 50 m
        dist=50;
		}      
    else if(Flag=="100m"){
       pathName="100m\\"; 
        Filesdistance= open1(pathName);                                              // collect the names of files at distance 100 m
        dist=100;
		}
    else if(Flag=="150m"){
		pathName="150m\\"; 
        Filesdistance= open1(pathName);                                              // collect the names of files at distance 150 m
        dist=150;
		}
	else{
		cout<<"Error loading file";
		break;
		}

    for(vector<string>::iterator itr= Filesdistance.begin()+2; itr != Filesdistance.end(); ++itr){  // we will iterate through each filename to esimate power of the the data
        //cout<<*(itr)<<endl;
		LoadData(pathName,*(itr));                                                                  // load the data
                
        bandpass_symmetric_filter(dataObject1, n_data, samprate, center_freq, order);    			// band pass filter   on channel 1 data
        bandpass_symmetric_filter(dataObject3, n_data, samprate, center_freq, order);    			// band passs filter on channel 3 data
		
		float* diff_dataObject1=diff(dataObject1,1000); 
		float* diff_dataObject3=diff(dataObject3,1000);

                
        Complex Data1[999];                                                            			//placeholder to convert real data type double into complex; needs #include <complex>
        Complex Data3[999];                                                            			//placeholder to convert real data type double into complex

		for(int i=0; i< 999;i++){
			Data1[i]= diff_dataObject1[i];                                             				    // convert a double[] array into Complex<double> datatype
			Data3[i]= diff_dataObject3[i];
		}              
        
		CArray data1 (Data1,999);                                                       			// convert data1 into Valarray<complex> datatype; needs valarray datatype
	    CArray data3 (Data3,999);                                                       			// convert data3 into valarray<complex> datatype
		
        fft(data1);                                                                        			// FFt of channel 1
        fft(data3);                                                                        			// FFT of channel 3
		
		// write the distance to a seperate file
		ofstream distObject("distance.txt",ios::out |ios::app);
		distObject << dist<<endl;
		distObject.close();
		
		//distance.push_back(dist); 

                                                                                   
        // save the FFT data in a text file
        ofstream FFTData1("FFTChannel1.txt",ios::out |ios::app);
        ofstream FFTData3("FFTChannel3.txt",ios::out |ios::app);
        for(int k=0;k<999;k++){
			FFTData1<<data1[k]<<" "; //Outputs array to txtFile
            FFTData3<<data3[k]<<" ";
			}
        FFTData1<<endl;
        FFTData3<<endl;

        FFTData1.close();
        FFTData3.close();   
		delete[] diff_dataObject1,diff_dataObject3;
		}
	}

return(0);


}
        


