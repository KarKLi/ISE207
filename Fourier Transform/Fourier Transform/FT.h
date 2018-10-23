#ifndef _FT
#define _FT
/*
	FT.h, written by Kark Li, 2018. All rights reserved.
	This header file is provided for Fourier Transform.
	Because computer's element was discrete, so this header
	is used for DFT.
	Please remember, when you centralize a DFT image, you cannot use the IDFT.
	So If you use DFT(true), be sure not use IDFT, or you will get strange Image.
	Why?
	In DFT, assuming the picture's size is N*N, the time complexity is given by this formula:O(N^3)
	If we I-centralize a picture to make sure it can be used to IDFT, we need to do the work with the time complexity of O(N^3),
	which time cost is very large.
*/
/*
			ErrorCode:
			0x00000001 Length Error
			0x00000002 Add Data Error
			0x00000003 Vector Empty
			0x00000004 Array Error
			0x00000005 IDFT before DFT
			0x00000006 Size not match
*/
#include <vector>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>
#define _ERRNO_SUCEESS 0
#define _ERRNO_FAIL 1
#define PI 3.14159265358979323846
typedef unsigned int uint;
typedef unsigned char uchar;
using namespace std;
namespace FT
{
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
		Complex(double x, double y)
		{
			this->real = x;
			this->Imaginary = y;
		}
		Complex(const Complex &obj)
		{
			this->real = obj.real;
			this->Imaginary = obj.Imaginary;
		}
		Complex operator=(double x)
		{
			real = x;
			Imaginary = 0.0;
			return *this;
		}
		friend Complex operator+(Complex obj1, Complex obj2)
		{
			Complex tmp;
			tmp.real = obj1.real + obj2.real;
			tmp.Imaginary = obj1.Imaginary + obj2.Imaginary;
			return Complex(tmp);
		}
		friend Complex operator-(Complex obj1, Complex obj2)
		{
			Complex tmp;
			tmp.real = obj1.real - obj2.real;
			tmp.Imaginary = obj1.Imaginary - obj2.Imaginary;
			return Complex(tmp);
		}
		friend Complex operator*(Complex obj1, double num)
		{
			Complex tmp;
			tmp.real = obj1.real*num;
			tmp.Imaginary = obj1.Imaginary*num;
			return Complex(tmp);
		}
		friend Complex operator*(Complex obj1, Complex obj2)
		{
			Complex tmp;
			tmp.real = obj1.real*obj2.real - obj1.Imaginary*obj2.Imaginary;
			tmp.Imaginary = obj1.real*obj2.Imaginary + obj1.Imaginary*obj2.real;
			return Complex(tmp);
		}
		friend Complex operator/(Complex obj1, double num)
		{
			_STL_VERIFY((num != 0), "Division cannot be zero!");
			Complex tmp;
			tmp.real = obj1.real / num;
			tmp.Imaginary = obj1.Imaginary / num;
			return Complex(tmp);
		}
		friend Complex operator/(Complex obj1, Complex obj2)
		{
			_STL_VERIFY((obj2.real != obj2.Imaginary), "Division cannot be zero!");
			Complex tmp;
			Complex tmp2 = obj2;
			tmp2.Imaginary = -obj2.Imaginary;
			tmp = (obj1 * tmp2) / (obj2.real*obj2.real + obj2.Imaginary + obj2.Imaginary);
			return Complex(tmp);
		}
		friend ostream &operator<<(ostream &os, const Complex obj)
		{
			double epsilon = 1e-4;
			if (abs(obj.real - 0.0) < epsilon)
			{
				if (abs(obj.Imaginary - -1.0) < epsilon)
					os << "-i";
				else if (abs(obj.Imaginary - 1.0) < epsilon)
					os << "i";
				else if (abs(obj.Imaginary - 0.0) < epsilon)
					os << fixed << setprecision(3) << abs(obj.real);
				else
					os << fixed << setprecision(3) << obj.Imaginary << "i";
			}
			else
			{
				if (abs(obj.Imaginary - -1.0) < epsilon)
					os << fixed << setprecision(3) << obj.real << " - i";
				else if (abs(obj.Imaginary - 1.0) < epsilon)
					os << fixed << setprecision(3) << obj.real << " + i";
				else if (abs(obj.Imaginary - 0.0) < epsilon)
					os << fixed << setprecision(3) << obj.real;
				else if (obj.Imaginary < -epsilon)
					os << fixed << setprecision(3) << obj.real << " - " << setprecision(3) << abs(obj.Imaginary) << "i";
				else
					os << fixed << setprecision(3) << obj.real << " + " << setprecision(3) << obj.Imaginary << "i";
			}
			return os;
		}
		friend istream &operator>>(istream &is, Complex obj)
		{
			is >> obj.real >> obj.Imaginary;
			return is;
		}
	};

	double sqrt(Complex obj)
	{
		return std::sqrt(obj.real*obj.real + obj.Imaginary*obj.Imaginary);
	}

	inline Complex W(int N, int k)
	{
		Complex tmp;
		tmp.real = cos(2 * PI*k / N);
		tmp.Imaginary = -sin(2 * PI*k / N);
		return tmp;
	}

	inline void Conjugation(Complex &obj)
	{
		obj.Imaginary = -obj.Imaginary;;
	}
	//typedef _complex Complex;

	void FFT(vector<Complex> a, vector<Complex> &dst, int len, bool IFFT = false)
	{
		int counts = 0;
		int num = len - 1;
		while (num > 0)
		{
			num = num >> 1;
			counts++;
		}
		int max = pow(2, counts);
		while (len < max)
			len++;
		a.resize(len);
		dst.resize(len);
		int *pos = new int[len];
		Complex *arr1 = new Complex[len];//���a[pos[i]]�ĸ���Ҷ�任
		Complex *arr2 = new Complex[len];//��Ż���DFT[a[pos[i]]�ĵ�ǰ�㸵��Ҷ�任
		string str;
		if (pos == NULL || arr1 == NULL || arr2 == NULL)
			return;
		for (int i = 0; i < len; i++)
		{
			int count = i;
			int result = 0;
			if (count == 0)
				str.append("0");
			while (count > 0)
			{
				int tmp = count % 2;
				if (tmp == 0)
					str.append("0");
				else
					str.append("1");
				count = count >> 1;
			}
			while (str.length() < counts)
				str.append("0");
			reverse(str.begin(), str.end());
			for (int i = 0; i < str.length(); i++)
				result += (str[i] - 48)*pow(2, i);
			pos[i] = result;
			str.clear();
		}
		//����a[pos[i]]�ĸ���Ҷ�任
		for (int i = 0; i < len; i++)
		{
			arr1[i] = a[pos[i]];
		}
		/*
		һ��Ҷ����Ϊ2^n�������������������Ϊn
		��a[pos[i]]�ĸ���Ҷ�任�У��������ö������ĵ�n�㣨��������Ϊ��1�㣩
		��ˣ��һ���Ҫѭ��for 1 to n-1����������һ�����ڵ�
		�����ٸ���Ҷ�任��ԭ�����ڣ���ֻ��Ҫ����log2(n)���������������õ�����ȫ�����ڵ㣬�Ϳ��Եõ��������x(n)�ĸ���Ҷ�任
		*/
		//���л���Ļ���arr1[i]�ĸ���Ҷ�任���õ���һ��
		int Vector = IFFT ? -1 : 1;
		for (int layer = log2(len); layer >= 1; layer--)
		{
			//ȷ��һ���
			int divide = len >> layer;
			//ȷ����ǰ���N
			int N = divide << 1;//N=divide*2
			//����������ÿ����������ӦΪgroup=pow(2,layer-1),��������㣬��Ҫ����������
			int group = pow(2, layer - 1);
			//����ÿ���������еĳ��ȣ���������㣬ÿ���������еĳ���Ϊ2����N
			int listlen = N;
			int count1 = 0;
			int count2 = 0;
			for (int k = 0; k < divide; k++)
			{
				for (int j = 0; j < pow(2, layer - 1); j++)
				{
					arr2[count1] = arr1[count2] + W(N, Vector*k)*arr1[count2 + 1];
					arr2[count1 + len / 2] = arr1[count2] - W(N, Vector*k)*arr1[count2 + 1];
					count1++;
					count2 += 2;
				}

			}
			for (int i = 0; i < len; i++)
			{
				arr1[i] = arr2[i];
			}
			memset(arr2, 0, 2 * sizeof(double)*len);
		}
		if (IFFT == false)
			Vector = 1;
		else
			Vector = len;
		for (int i = 0; i < len; i++)
			dst[i] = arr1[i] / Vector;
		delete arr1;
		delete arr2;
		arr1 = arr2 = nullptr;
		return;
	}

	void FFT2(vector<vector<Complex>> a, vector<vector<Complex>>&dst, int x, int y, bool IFFT = false)
	{
		vector<Complex>tempx(x);
		vector<Complex>tempy(y);
		if (!IFFT)
		{
			//��ʼ���ж���
			for (int i = 0; i < y; i++)
			{
				for (int j = 0; j < x; j++)
				{
					tempx[j] = a[i][j];
				}
				//����һ�У�����
				FFT(tempx, tempx, x, IFFT);
				for (int j = 0; j < x; j++)
				{
					a[i][j] = tempx[j];
				}
			}
			//��ʼ���ж���
			for (int i = 0; i < x; i++)
			{
				for (int j = 0; j < y; j++)
				{
					tempy[j] = a[j][i];
				}
				//����һ�У�����
				FFT(tempy, tempy, y, IFFT);
				for (int j = 0; j < y; j++)
				{
					a[j][i] = tempy[j];
				}
			}
		}
		else
		{
			for (int i = 0; i < x; i++)
			{
				for (int j = 0; j < y; j++)
				{
					tempy[j] = a[j][i];
				}
				//����һ�У�����
				FFT(tempy, tempy, y, IFFT);
				for (int j = 0; j < y; j++)
				{
					a[j][i] = tempy[j];
				}
			}
			for (int i = 0; i < y; i++)
			{
				for (int j = 0; j < x; j++)
				{
					tempx[j] = a[i][j];
				}
				//����һ�У�����
				FFT(tempx, tempx, x, IFFT);
				for (int j = 0; j < x; j++)
				{
					a[i][j] = tempx[j];
				}
			}
		}
		dst = a;
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
		void LoadData(const Complex *Array, int len = 0)
		{
			if (len == 0)
				len = this->length;
			for (int i = 0; i < len; i++)
			{
				try {
					Data[i] = Array[i];
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
				for (int k = 0; k <= N - 1; k++)
				{
					Complex temp = { 0,0 };
					for (int n = 0; n <= N - 1; n++)
					{
						double cosRes = cos(2 * PI*k*n / N);
						double sinRes = sin(2 * PI*k*n / N);
						temp.real += Data[n].real*cosRes + Data[n].Imaginary*sinRes;
						temp.Imaginary += Data[n].Imaginary*cosRes - Data[n].real*sinRes;
					}
					DFTData[k] = temp / N;
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
				for (int n = 0; n <= N - 1; n++)
				{
					Complex temp = { 0,0 };
					for (int k = 0; k <= N - 1; k++)
					{

						double cosRes = cos(2 * PI*k*n / N);
						double sinRes = sin(2 * PI*k*n / N);
						temp.real += (1.0 / N)*(Data[k].real*cosRes - Data[k].Imaginary*sinRes);
						temp.Imaginary += (1.0 / N)*(Data[k].Imaginary*cosRes + Data[k].real*sinRes);
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
		void GetPhaseArray(double *Array, int length)
		{
			if (length > this->length)
			{
				ErrorCode = 4;
				return;
			}
			else
			{
				for (int i = 0; i < length; i++)
					Array[i] = atan2(Data[i].Imaginary, Data[i].real);
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
					Array[i] = Data[i].real*Data[i].real + Data[i].Imaginary*Data[i].Imaginary;
			}
		}
		int WriteData(Complex *Array, int length = -1)
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
		int WriteData(vector<Complex>&Array, int Beginindex = 0, int length = -1)
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

		int GetLastErrorCode()const { return ErrorCode; }

		void LoadData(const vector<vector<Complex>> Array)
		{
			if (Array.size() != xLength)
			{
				ErrorCode = 2;
				return;
			}
			else
			{
				for (int i = 0; i < xLength; i++)
				{
					if (Array[i].size() != yLength)
					{
						ErrorCode = 2;
						return;
					}
				}
			}
			Data = Array;
		}
		
		void SwapMatrix(vector<vector<Complex>> &Array1,vector<vector<Complex>> &Array2, const uint x,const uint y)
		{
			if (
				(Array1.size() != Array2.size()) ||
				(Array1[0].size() != Array2[0].size()) ||
				(Array1.size() != y) ||
				(Array1[0].size() != x)
				)
			{
				ErrorCode = 6;
				return;
			}
			else
			{
				std::swap(Array1, Array2);
			}
		}

		void DFT(bool Central = false)
		{
			if (xLength != yLength || (xLength % 2 !=0) || (yLength % 2!=0))
				Central = false;
			if (hasDFTed == 1)
				return;
			if (Data.empty())
			{
				ErrorCode = 3;
				return;
			}
			vector<vector<Complex>>DFTData = this->Data;
			int M = xLength;
			int N = yLength;
			/*
				���ÿɷ����ԣ�����y����������Ҷ�任��������������飬����������x�����Ͻ��и���Ҷ�任
			*/
			if (Central == true)
			{
				for (int x = 0; x <= M - 1; x++)
				{
					for (int y = 0; y <= N - 1; y++)
					{
						Data[x][y] = Data[x][y] * pow(-1, x + y);
					}
				}
			}
			Complex Base;
			for (int u = 0; u < xLength; u++)
				for (int v = 0; v < yLength; v++)
				{
					Complex sum;
					for (int x = 0; x <= M - 1; x++)
					{
						for (int y = 0; y <= N - 1; y++)
						{
<<<<<<< HEAD
=======

>>>>>>> ac8337cffc1b28815bc76664cd8dc49c62a1b627
							Base.real = cos(2 * PI*u*x / M + 2 * PI*v*y / N);
							Base.Imaginary = -sin(2 * PI*u*x / M + 2 * PI*v*y / N);
							sum = sum + Data[x][y] * Base;
						}
					}
					sum = sum / (M*N);
					DFTData[u][v] = sum;
				}
			Data = DFTData;
			hasDFTed = 1;
			
		}

		void InverseDFT()
		{
			vector<vector<Complex>>InverseData;
			InverseData.resize(xLength);
			for (int i = 0; i < xLength; i++)
			{
				InverseData[i].resize(yLength);
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
					���ÿɷ����ԣ�����y����������Ҷ�任��������������飬����������x�����Ͻ��и���Ҷ�任
				*/
				Complex Base;
				for (int u = 0; u < xLength; u++)
					for (int v = 0; v < yLength; v++)
					{
						Complex sum;
						for (int x = 0; x <= M - 1; x++)
						{
							for (int y = 0; y <= N - 1; y++)
							{
								Base.real = cos(2 * PI*u*x / M + 2 * PI*v*y / N);
								Base.Imaginary = sin(2 * PI*u*x / M + 2 * PI*v*y / N);
								sum = sum + Data[x][y] * Base;
							}
						}
						InverseData[u][v] = sum;
					}
				/*if (hasCentraled)
				{
					for (int x = 0; x <= M - 1; x++)
					{
						for (int y = 0; y <= N - 1; y++)
						{
							InverseData[x][y] = InverseData[x][y] * pow(-1, x + y);
						}
					}
					hasCentraled = 0;
				}*/
				Data = InverseData;
				hasDFTed = 0;
			}
		}

		int WriteData(vector<vector<Complex>>&Array, int xlength = -1, int ylength = -1, int Beginxindex = 0, int Beginyindex = 0)
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
			for (int i = 0; i < xLength; i++)
			{
				for (int j = 0; j < yLength; j++)
				{
					Array[i + Beginxindex][j + Beginyindex] = Data[i][j];
					SuccLen++;
				}
			}
			return SuccLen;
		}

		void GetAmplitudeArray(vector<vector<double>>&Array, int xlength, int ylength)
		{
			double max = -1.0;
			double min = 32767.0;
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
				for (int i = 0; i < xlength; i++)
					for (int j = 0; j < ylength; j++)
					{
						Array[i][j] = sqrt(Data[i][j]);
						if (Array[i][j] > max)
							max = Array[i][j];
						if (Array[i][j] < min)
							min = Array[i][j];
					}
				for (int i = 0; i < xlength; i++)
					for (int j = 0; j < ylength; j++)
					{
						Array[i][j] = 255.0*((Array[i][j] - min) / (max - min));
					}
			}
		}

		void GetPhaseArray(vector<vector<double>>&Array, int xlength, int ylength)
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

		void GetFrequencyArray(vector<vector<double>>&Array, int xlength, int ylength)
		{
			GetAmplitudeArray(Array, xlength, ylength);
		}

		void GetPowerArray(vector<vector<double>>&Array, int xlength, int ylength)
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
				for (int i = 0; i < xLength; i++)
					for (int j = 0; j < yLength; j++)
					{
						Array[i][j] = Data[i][j].real*Data[i][j].real + Data[i][j].Imaginary*Data[i][j].Imaginary;
					}
			}
		}

		void Centralized()
		{
			hasCentraled = 1;
			for (uint i = 0; i < yLength / 2; i++)
				for (uint j = 0; j < xLength / 2; j++)
				{
					swap(Data[i][j], Data[i + yLength / 2][j + xLength / 2]);
				}
			for (uint i = 0; i < yLength / 2; i++)
			{
				for (uint j = xLength / 2; j < xLength; j++)
				{
					swap(Data[i][j], Data[i + yLength / 2][j - xLength / 2]);
				}
			}
		}

		inline void Uncentralized()
		{
			Centralized();
			hasCentraled = 0;
		}

		void _FFT2(bool Central=false)
		{
			if (xLength != yLength || (xLength % 2 != 0) || (yLength % 2 != 0))
				Central = false;
			if (hasDFTed == 1)
				return;
			if (Data.empty())
			{
				ErrorCode = 3;
				return;
			}
			_STL_VERIFY(xLength == yLength, "Two Dimension matrix must be square!");
			//vector<vector<Complex>>DFTData = this->Data;
			int M = xLength;
			int N = yLength;
			FFT2(Data, Data, M, N, false);
			if (Central == true)
			{
				hasCentraled = 1;
				/*
					After show or write your picture,please call the Uncentralized to restore the regular image.
				*/
				Centralized();
			}
			hasDFTed = 1;
		}

		void _IFFT2()
		{
			_STL_VERIFY(xLength == yLength, "Two Dimension matrix must be square!");
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
				if (hasCentraled == 1)
				{
					hasCentraled = 0;
					for (uint i = 0; i < yLength / 2; i++)
						for (uint j = 0; j < xLength / 2; j++)
						{
							swap(Data[i][j], Data[i + yLength / 2][j + xLength / 2]);
						}
					for (uint i = 0; i < yLength / 2; i++)
					{
						for (uint j = xLength / 2; j < xLength; j++)
						{
							swap(Data[i][j], Data[i + yLength / 2][j - xLength / 2]);
						}
					}
				}
				FFT2(Data, Data, M, N, true);
				hasDFTed = 0;
			}
		}

	private:
		int xLength;
		int yLength;
		int ErrorCode;
		vector<vector<Complex>>Data;
		int hasDFTed = 0;
		int hasCentraled = 0;
	};
}
#endif