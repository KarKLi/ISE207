#ifndef _FT
#define _FT
/*
	FT.h, written by Kark Li, 2018. All rights reserved.
	This header file is provided for Fourier Transform.
	Because computer's element was discrete, so this header
	is used for DFT.
*/
#include <vector>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>
#define _ERRNO_SUCEESS 0
#define _ERRNO_FAIL 1
#define PI 3.1415926
using namespace std;
class Complex
{
public:
	double real;
	double Imaginary;
	Complex()
	{
		this->real = 0;
		this->Imaginary = 0;
	}
	Complex(double x,double y)
	{
		this->real = x;
		this->Imaginary = y;
	}
	Complex(const Complex &obj)
	{
		this->real = obj.real;
		this->Imaginary = obj.Imaginary;
	}
	Complex &operator=(double x)
	{
		real = x;
		Imaginary = 0.0;
		return *this;
	}
	friend Complex &operator+(Complex &obj1, Complex &obj2)
	{
		Complex tmp;
		tmp.real = obj1.real + obj2.real;
		tmp.Imaginary = obj1.Imaginary + obj2.Imaginary;
		return tmp;
	}
	friend Complex &operator-(Complex &obj1, Complex &obj2)
	{
		Complex tmp;
		tmp.real = obj1.real - obj2.real;
		tmp.Imaginary = obj1.Imaginary - obj2.Imaginary;
		return tmp;
	}
	friend Complex &operator*(Complex &obj1, Complex &obj2)
	{
		Complex tmp;
		tmp.real = obj1.real*obj2.real - obj1.Imaginary*obj1.Imaginary;
		tmp.Imaginary = obj1.real*obj2.Imaginary + obj1.Imaginary*obj2.real;
		return tmp;
	}
	friend Complex &operator/(Complex &obj1, double num)
	{
		_STL_VERIFY((num == 0), "Division cannot be zero!");
		Complex tmp;
		tmp.real = obj1.real / num;
		tmp.Imaginary = obj1.Imaginary / num;
		return tmp;
	}
	friend Complex &operator/(Complex &obj1, Complex &obj2)
	{
		_STL_VERIFY((obj2.real == obj2.Imaginary), "Division cannot be zero!");
		Complex tmp;
		tmp = (obj1 * obj2) / (pow(obj2.real, 2) + pow(obj2.real, 2));
		return tmp;
	}
	friend ostream &operator<<(ostream &os, const Complex &obj)
	{
		double epsilon = 1e-4;
		if ((obj.real - 0.0) < epsilon)
		{
			if (abs(obj.Imaginary - -1.0) < epsilon)
				os << "-i" << endl;
			else if (abs(obj.Imaginary - 1.0) < epsilon)
				os << "i" << endl;
			else if (abs(obj.Imaginary - 0.0) < epsilon)
				os << fixed << setprecision(3) << obj.real << endl;
			else
				os << fixed << setprecision(3) << obj.Imaginary << "i" << endl;
		}
		else
		{
			if (abs(obj.Imaginary - -1.0) < epsilon)
				os << fixed << setprecision(3) << obj.real << " - i" << endl;
			else if (abs(obj.Imaginary - 1.0) < epsilon)
				os << fixed << setprecision(3) << obj.real << " + i" << endl;
			else if (abs(obj.Imaginary - 0.0) < epsilon)
				os << fixed << setprecision(3) << obj.real << endl;
			else if (obj.Imaginary < 0.0)
				os << fixed << setprecision(3) << obj.real << " - " << setprecision(3) << abs(obj.Imaginary) << "i" << endl;
			else
				os << fixed << setprecision(3) << obj.real << " + " << setprecision(3) << obj.Imaginary << "i" << endl;
		}
		return os;
	}
	friend istream &operator>>(istream &is, Complex &obj)
	{
		is >> obj.real >> obj.Imaginary;
		return is;
	}
};
double sqrt(Complex obj)
{
	return std::sqrt(obj.real*obj.real + obj.Imaginary*obj.Imaginary);
}
class DFTOneDim
{
public:
	DFTOneDim(int length)
	{
		if (length <= 0)/*Length must bigger than 0*/
		{
			ErrorCode = 1;
			cerr << "In File: " << __FILE__ << ", line: " << __LINE__ << " find an error, ErrorCode: " << hex << ErrorCode;
			exit(ErrorCode);
		}
		else
		{
			this->length = length;
			Data.resize(length);
		}
	}
	DFTOneDim(const DFTOneDim& obj)
	{
		this->index = obj.index;
		this->ErrorCode = obj.ErrorCode;
		this->length = length;
		this->Data.resize(this->length);
		for (int i = 0; i < length; i++)
			this->Data[i] = obj.Data[i];
	}
	~DFTOneDim() {}
	int GetLastErrorCode()const { return ErrorCode; }
	void LoadData(const Complex *Array, int len=0)
	{
		if (len == 0)
			len = this->length;
		for (int i = 0; i < len; i++)
		{
			try {
				Data[i]=Array[i];
			}
			catch (exception &e)
			{
				ErrorCode = 2;
				return;
			}
		}
	}
	void DFT()
	{
		vector<Complex>DFTData(length);
		if (Data.empty())
		{
			ErrorCode = 3;
			return;
		}
		else
		{
			int N = length;
			for (int k = 0; k <= N-1; k++)
			{
				Complex temp = { 0,0 };
				for (int n = 0; n <= N - 1; n++)
				{
					double cosRes = cos(2 * PI*k*n / N);
					double sinRes = sin(2 * PI*k*n / N);
					temp.real += Data[n].real*cosRes+Data[n].Imaginary*sinRes;
					temp.Imaginary += Data[n].Imaginary*cosRes-Data[n].real*sinRes;
				}
				DFTData[k] = temp;
			}
		}
		hasDFTed = 1;
		Data = DFTData;
	}
	void InverseDFT()
	{
		vector<Complex>InverseData(length);
		if (hasDFTed == 0)
		{
			ErrorCode = 5;
			return;
		}
		if (Data.empty())
		{
			ErrorCode = 3;
			return;
		}
		else
		{
			int N = length;
			for (int n = 0; n <= N-1; n++)
			{
				Complex temp = { 0,0 };
				for (int k = 0; k <= N - 1; k++)
				{
					
					double cosRes = cos(2 * PI*k*n / N);
					double sinRes = sin(2 * PI*k*n / N);
					temp.real += (1.0/N)*(Data[k].real*cosRes - Data[k].Imaginary*sinRes);
					temp.Imaginary += (1.0/N)*(Data[k].Imaginary*cosRes + Data[k].real*sinRes);
				}
				InverseData[n] = temp;
			}
		}
		hasDFTed = 0;
		Data = InverseData;
	}
	void GetAmplitudeArray(double *Array, int length)
	{
		if (length > this->length)
		{
			ErrorCode = 4;
			return;
		}
		else
		{
			for (int i = 0; i < length; i++)
				Array[i] = sqrt(Data[i]);
		}
	}
	void GetPhaseArray(double *Array,int length)
	{
		if (length > this->length)
		{
			ErrorCode = 4;
			return;
		}
		else
		{
			for (int i = 0; i < length; i++)
				Array[i] = atan2(Data[i].Imaginary,Data[i].real);
		}
	}
	void GetFrequencyArray(double *Array, int length)
	{
		GetAmplitudeArray(Array, length);
	}
	void GetPowerArray(long double *Array, int length)
	{
		if (length > this->length)
		{
			ErrorCode = 4;
			return;
		}
		else
		{
			for (int i = 0; i < length; i++)
				Array[i] = Data[i].real*Data[i].real+Data[i].Imaginary*Data[i].Imaginary;
		}
	}
	int WriteData(Complex *Array, int length=-1)
	{
		if (length == -1)
			length = this->length;
		else if (length == 0)
			return length;
		if (length > this->length)
		{
			int SuccLen = 0;
			for (int i = 0; i < this->length; i++)
			{
				try
				{
					Array[i] = Data[i];
				}
				catch (exception &e)
				{
					ErrorCode = 4;
					return SuccLen;
				}
				SuccLen++;
			}
			return SuccLen;
		}
		else
		{
			int SuccLen = 0;
			for (int i = 0; i < length; i++)
			{
				Array[i] = Data[i];
				SuccLen++;
			}
			return length;
		}
	}
	int WriteData(vector<Complex>Array, int Beginindex = 0, int length = -1)
	{
		if (length == -1)
			length = this->length;
		else if (length == 0)
			return length;
		int SuccLen = 0;
		for (int i = 0; i < length; i++)
		{
			Array[i + Beginindex] = Data[i];
			SuccLen++;
		}
		return SuccLen;
	}
private:
	int index = 0;
	/*
		ErrorCode:
		0x00000001 Length Error
		0x00000002 Add Data Error
		0x00000003 Vector Empty
		0x00000004 Array Error
		0x00000005 IDFT before DFT
	*/
	int length;
	int hasDFTed = 0;
	vector<Complex>Data;
protected:
	int ErrorCode = 0;
};
class DFTTwoDim
{
public:
	DFTTwoDim(int x, int y)
	{
		if (x <= 0 || y <= 0)
		{
			this->ErrorCode = 1;
		}
		else
		{
			this->xLength = x;
			this->yLength = y;
			this->ErrorCode = 0;
		}
		Data.resize(x);
		for (int i = 0; i < x; i++)
			Data[i].resize(y);
	}
	DFTTwoDim(const DFTTwoDim & obj)
	{
		this->xLength = obj.xLength;
		this->yLength = obj.yLength;
		this->ErrorCode = obj.ErrorCode;
		this->Data = obj.Data;
	}
	~DFTTwoDim()
	{
		//do nothing
	}
	void DFT()
	{
		if (hasDFTed == 0)
		{
			ErrorCode = 5;
			return;
		}
		if (Data.empty())
		{
			ErrorCode = 3;
			return;
		}
		vector<vector<Complex>>DFTData;
		DFTData.resize(xLength);
		for (int i = 0; i < xLength; i++)
			DFTData[i].resize(yLength);
		int M = xLength;
		int N = yLength;
		/*
			利用可分离性，先在y方向上求傅里叶变换，将结果存入数组，将该数组在x方向上进行傅里叶变换
		*/
		//y方向
		for (int u = 0; u < xLength; u++)
		{
			int x = u;
			for (int v = 0; v < yLength; v++)
				for (int y = 0; y <= N - 1; y++)
				{
					//Fourier base
					Complex base = { cos(2 * PI*u*x + 2 * PI*v*y),-sin(2 * PI*u*x + 2 * PI*v*y) };
					DFTData[u][v] = (DFTData[u][v] + Data[x][y] * base);
				}
		}
		//x方向
		for (int v = 0; v < yLength; v++)
		{
			int y = v;
			for (int u = 0; u < yLength; u++)
				for (int x = 0; x <= M - 1; x++)
				{
					//Fourier base
					Complex base = { cos(2 * PI*u*x + 2 * PI*v*y),-sin(2 * PI*u*x + 2 * PI*v*y) };
					DFTData[u][v] = (DFTData[u][v] + Data[x][y] * base);
				}
		}
		Data = DFTData;
		hasDFTed = 1;
	}
	void InverseDFT()
	{
		vector<vector<Complex>>InverseData;
		InverseData.resize(yLength);
		for (int i = 0; i < yLength; i++)
		{
			InverseData[i].resize(xLength);
		}
		if (hasDFTed == 0)
		{
			ErrorCode = 5;
			return;
		}
		if (Data.empty())
		{
			ErrorCode = 3;
			return;
		}
		else
		{
			int M = xLength;
			int N = yLength;
			/*
				利用可分离性，先在y方向上求傅里叶变换，将结果存入数组，将该数组在x方向上进行傅里叶变换
			*/
			//y方向
			for (int x = 0; x < xLength; x++)
			{
				int u = x;
				for (int y = 0; y < yLength; y++)
					for (int v = 0; v <= N - 1; v++)
					{
						//Fourier base
						Complex base = { cos(2 * PI*u*x + 2 * PI*v*y),sin(2 * PI*u*x + 2 * PI*v*y) };
						base = base / (M*N);
						InverseData[x][y] = (InverseData[x][y] + Data[u][v] * base);
					}
			}
			//x方向
			for (int y = 0; y < yLength; y++)
			{
				int v = y;
				for (int x = 0; x < xLength; x++)
					for (int u = 0; u <= M - 1; u++)
					{
						//Fourier base
						Complex base = { cos(2 * PI*u*x + 2 * PI*v*y),sin(2 * PI*u*x + 2 * PI*v*y) };
						base = base / (M*N);
						InverseData[x][y] = (InverseData[x][y] + Data[u][v] * base);
					}
			}
		}
		hasDFTed = 0;
		Data = InverseData;
	}
	int WriteData(vector<vector<Complex>>Array, int xlength = -1, int ylength = -1, int Beginxindex = 0, int Beginyindex = 0)
	{
		if (xlength == -1)
			xlength = this->xLength;
		else if (xlength == 0)
			return xlength;
		if (ylength == -1)
			ylength = this->yLength;
		else if (ylength == 0)
			return ylength;
		int SuccLen = 0;
		for (int i = 0; i < yLength; i++)
		{
			for (int j = 0; j < xLength; j++)
			{
				Array[i + Beginyindex][j+Beginxindex] = Data[i][j];
				SuccLen++;
			}
		}
		return SuccLen;
	}
	void GetAmplitudeArray(vector<vector<Complex>>Array, int xlength,int ylength)
	{
		if (xlength > this->xLength)
		{
			ErrorCode = 4;
			return;
		}
		else if (ylength > this->yLength)
		{
			ErrorCode = 4;
			return;
		}
		else
		{
			for (int i = 0; i < ylength; i++)
				for (int j = 0; j < xlength; j++)
				{
					Array[i][j] = sqrt(Data[i][j]);
				}
		}
	}
	void GetPhaseArray(vector<vector<Complex>>Array, int xlength,int ylength)
	{
		if (xlength > this->xLength)
		{
			ErrorCode = 4;
			return;
		}
		else if (ylength > this->yLength)
		{
			ErrorCode = 4;
			return;
		}
		else
		{
			for (int i = 0; i < ylength; i++)
				for (int j = 0; j < xlength; j++)
				{
					Array[i][j] = atan2(Data[i][j].Imaginary, Data[i][j].real);
				}
		}
	}
	void GetFrequencyArray(vector<vector<Complex>>Array, int xlength,int ylength)
	{
		GetAmplitudeArray(Array, xlength,ylength);
	}
	void GetPowerArray(vector<vector<Complex>>Array, int xlength,int ylength)
	{
		if (xlength > this->xLength)
		{
			ErrorCode = 4;
			return;
		}
		else if (ylength > this->yLength)
		{
			ErrorCode = 4;
			return;
		}
		else
		{
			for (int i = 0; i < yLength; i++)
				for (int j = 0; j < xLength; j++)
				{
					Array[i][j] = Data[i][j].real*Data[i][j].real + Data[i][j].Imaginary*Data[i][j].Imaginary;
				}
		}
	}
private:
	int xLength;
	int yLength;
	int ErrorCode;
	vector<vector<Complex>>Data;
	int hasDFTed = 0;
};
#endif