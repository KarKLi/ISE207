#include "FT.h"
int main()
{
	Complex x[5];
	x[0].real = 147;
	x[0].Imaginary = 0;
	x[1].real = 213;
	x[1].Imaginary = 0;
	x[2].real = 125;
	x[2].Imaginary = 0;
	x[3].real = 58;
	x[3].Imaginary = 0;
	x[4].real = 12;
	x[4].Imaginary = 0;
	cout << "Original data:" << endl;
	for (int i = 0; i < 5; i++)
	{
		cout << "x" << i;
		cout.setf(ios::left);
		cout.width(5);
		cout << ": " << x[i] << endl;
	}
	DFTOneDim dft(5);
	dft.LoadData(x);
	dft.DFT();
	int error = dft.GetLastErrorCode();
	if (error == 0)
		dft.WriteData(x);
	cout << "After processing data:" << endl;
	for (int i = 0; i < 5; i++)
	{
		cout << "x" << i;
		cout.setf(ios::left);
		cout.width(5);
		cout<< ": " << x[i] << endl;
	}
	dft.InverseDFT();
	error = dft.GetLastErrorCode();
	if (error == 0)
		dft.WriteData(x);
	cout << "After inverse process data:" << endl;
	for (int i = 0; i < 5; i++)
	{
		cout << "x" << i;
		cout.setf(ios::left);
		cout.width(5);
		cout << ": " << x[i] << endl;
	}
	dft.DFT();
	double Array[5];
	dft.GetAmplitudeArray(Array, 5);
	double Array2[5];
	dft.GetPhaseArray(Array2, 5);
	long double Array3[5];
	dft.GetPowerArray(Array3, 5);
	cout << "And its Amplitude data:" << endl;
	for (int i = 0; i < 5; i++)
	{
		cout << "x" << i;
		cout.setf(ios::left);
		cout.width(5);
		cout << ": " << Array[i] << endl;
	}
	cout << "And its Phase data:" << endl;
	for (int i = 0; i < 5; i++)
	{
		cout << "x" << i;
		cout.setf(ios::left);
		cout.width(5);
		cout << ": " << Array2[i] << endl;
	}cout << "And its Power data:" << endl;
	for (int i = 0; i < 5; i++)
	{
		cout << "x" << i;
		cout.setf(ios::left);
		cout.width(5);
		cout << ": " << Array3[i] << endl;
	}
  	return 0;
}