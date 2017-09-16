1)  SpectralResolution.exe function estimates the FFT of the Channel 1 and Channel 3 of the audio signal. 
    It generates output files FFTChannel1.txt, FFTChannel3.txt and distance.txt.
	{ input: folders:50m, 100m and 150m
	 Ouputs: FFTChannel1.txt, FFTChannel3.txt and distance.txt}
	
2)  SpectralPower.exe computes the power of the signals using FFTChannel1.txt and FFTChannel3.txt as input files.
    It uses Stochastic gradient descent regression to obtain the coeffient of linear relationship
	between Power and distance. The output is saved in Predicted_distance.txt files
    {Input:FFTChannel1.txt, FFTChannel3.txt and distance.txt
	 Ouput: Predicted_distance.txt}