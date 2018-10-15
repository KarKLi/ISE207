/*
kernel.h written by Kark.Li
Copyrights 2018, all rights reserved.
This header file is used for image processes' kernel,
which is a even square matrix used for convolution.
In this header, I use std::vector<vector<unsigned char>> to storage this matrix.
Because the size of matrix is constant when you initial it.
So the element of the std::vector is easy to locate.
The Image is storaged by two dimensions std::vector, and the pixel index can be found by this formula:
Pixel[Height][Width]
Use operator[] to access element.
In this header, I provide some methods to make convolution such as
Convolution();
which accept a enum element defined by kernel.h, to make the blur convolution,
sobel convolution, lapalacian convolution, etc.
Enjoy it!
This code will push on my GitHub page:Https://github.com/KarKLi/ISE207
Follow the agreement to download my source code!
*/
#ifndef _KERNEL_
#define _KERNEL_
enum ConvolutionOperation {Average,WeightedAverage,Medium,FirstDerivative,FirstDerivativex,FirstDerivativey,SecondDerivative,SecondDerivativex,SecondDerivativey,Lapacian};
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cassert>
#include <list>
#include <exception>
#include <algorithm>
#define TRUE 1
#define FALSE 0
#define ERRNO_FAIL 0
#define ERRNO_SUCCESS 1
typedef unsigned int UINT;
using namespace std;
inline bool comp(int i, int j) { return (i < j); }
class Kernel
{
public:
	Kernel(unsigned int ImageHeight,unsigned int ImageWidth, unsigned int KernelSize) : KernelSize(KernelSize), kerneldata(KernelSize,vector<double>(KernelSize)), ImageData(ImageHeight) ,
		                                                                                ImageWidth(ImageWidth),ImageHeight(ImageHeight),ImageInt(ImageHeight)
	{
		if (KernelSize % 2 == 0)
			ErrorCode = 1;// Kernel size must be even
		else
			ErrorCode = 0;
		for (UINT i = 0; i < ImageData.size(); i++)
			ImageData[i].resize(ImageWidth);
		for (UINT i = 0; i < ImageInt.size(); i++)
			ImageInt[i].resize(ImageWidth);
	};
	int checkErrorStat()const { return ErrorCode; }
	const char *Loadimage(unsigned char *Image);
	errno_t Convolution(ConvolutionOperation operators);
	errno_t Writeimage(unsigned char *dst);
private:
	unsigned int KernelSize;
	vector<vector<double>>kerneldata;
	vector<vector<unsigned char>>ImageData;
	vector<vector<int>>ImageInt;
	unsigned int ImageWidth;
	unsigned int ImageHeight;
	int ErrorCode = 0;
};
const char *Kernel::Loadimage(unsigned char *Image)
{
	if (ErrorCode == 1)
		return "Some matters occurred.";
	unsigned char *p = &Image[0];
	for (unsigned int i = 0; i < ImageHeight; i++)
	{
		for (unsigned int j = 0; j < ImageWidth; j++)
		{
			try {
				ImageData[i][j]=*p;
				ImageInt[i][j] = (int)*p;
				p++;
			}
			catch (exception &e){
				return e.what();
			}
		}
	}
}
errno_t Kernel::Convolution(ConvolutionOperation operators)
{
	if (ErrorCode == 1)
		return FALSE;
	//Create convolution kernel
	switch (operators)
	{
	case ConvolutionOperation::Average:
		for (unsigned int i = 0; i < KernelSize; i++)
		{
			for (unsigned int j = 0; j < KernelSize; j++)
			{
				kerneldata[i][j] = (1.0 / (double)(KernelSize*KernelSize));
			}
		}
		break;
	case ConvolutionOperation::WeightedAverage:
	{
		double sum = 0;
		for (unsigned int i = 0; i < KernelSize; i++)
		{
			int temp = 1;
			for (unsigned int j = 0; j < (KernelSize / 2) + 1; j++)
			{
				kerneldata[i][j] = temp;
				temp *= 2;
				sum += temp;
			}
			for (unsigned int j = (KernelSize / 2) + 2; j < KernelSize; j++)
			{
				kerneldata[i][j] = temp;
				temp /= 2;
				sum += temp;
			}
		}
		sum = 1.0 / sum;
		for (unsigned int i = 0; i < KernelSize; i++)
		{
			for (unsigned int j = 0; j < KernelSize; j++)
			{
				kerneldata[i][j] = sum * kerneldata[i][j];//Normalized
			}
		}
	}
	break;
	case ConvolutionOperation::Medium:
		//No kernel
		break;
	case ConvolutionOperation::FirstDerivativex:
		kerneldata[KernelSize / 2][KernelSize / 2] = 1;
		kerneldata[KernelSize / 2][KernelSize / 2 + 1] = 0;
		kerneldata[KernelSize / 2 + 1][KernelSize / 2] = -1;
		kerneldata[KernelSize / 2 + 1][KernelSize / 2 + 1] = 0;
		break;
	case ConvolutionOperation::FirstDerivativey:
		kerneldata[KernelSize / 2][KernelSize / 2] = 1;
		kerneldata[KernelSize / 2][KernelSize / 2 + 1] = -1;
		kerneldata[KernelSize / 2 + 1][KernelSize / 2] = 0;
		kerneldata[KernelSize / 2 + 1][KernelSize / 2 + 1] = 0;
		break;
	case ConvolutionOperation::FirstDerivative:
		/*Set Kernel.*/
		/*
		    first,this kernel given by next is the Gradient of X-axis:
			 1 0
			-1 0
			and the kernel give by next is the Gradient of Y-axis:
			1 -1
			0  0
			To reduce the calculation, we choose this formula:
			|G|≈|Gx|+|Gy|
			instead of
			|G|=(|Gx|+|Gy|)^1/2
		*/
		kerneldata[KernelSize/2][KernelSize/2] = 2;
		kerneldata[KernelSize/2][KernelSize/2 + 1] = -1;
		kerneldata[KernelSize/2 + 1][KernelSize/2] = -1;
		kerneldata[KernelSize/2 + 1][KernelSize/2 + 1] = 0;
		break;
	case ConvolutionOperation::SecondDerivative:
		//set kernel
		kerneldata[KernelSize / 2 - 1][KernelSize / 2] = 1;
		kerneldata[KernelSize / 2][KernelSize / 2 - 1] = 1;
		kerneldata[KernelSize / 2][KernelSize / 2] = -4;
		kerneldata[KernelSize / 2][KernelSize / 2 + 1] = 1;
		kerneldata[KernelSize / 2 + 1][KernelSize / 2] = 1;
		break;
	case ConvolutionOperation::SecondDerivativex:
		kerneldata[KernelSize / 2][KernelSize / 2] = -2;
		kerneldata[KernelSize / 2][KernelSize / 2 + 1] = 1;
		kerneldata[KernelSize / 2][KernelSize / 2-1] = 1;
		break;
	case ConvolutionOperation::SecondDerivativey:
		kerneldata[KernelSize / 2][KernelSize / 2] = -2;
		kerneldata[KernelSize / 2-1][KernelSize / 2] = 1;
		kerneldata[KernelSize / 2 + 1][KernelSize / 2] = 1;
		break;
	case ConvolutionOperation::Lapacian:
		//set Kernel
		kerneldata[KernelSize/2-1][KernelSize/2-1] = 1;
		kerneldata[KernelSize/2-1][KernelSize/2] = 1;
		kerneldata[KernelSize/2-1][KernelSize/2+1] = 1;
		kerneldata[KernelSize/2][KernelSize/2-1] = 1;
		kerneldata[KernelSize/2][KernelSize/2] = -8;
		kerneldata[KernelSize/2][KernelSize/2+1] = 1;
		kerneldata[KernelSize/2+1][KernelSize/2-1] = 1;
		kerneldata[KernelSize/2+1][KernelSize/2] = 1;
		kerneldata[KernelSize/2+1][KernelSize/2+1] = 1;
		break;
	default:
		break;
	}
	switch (operators)
	{
	case Average:
	case WeightedAverage:
	case Lapacian:
	case FirstDerivative:
	case FirstDerivativex:
	case FirstDerivativey:
	case SecondDerivative:
	case SecondDerivativex:
	case SecondDerivativey:
	{
		vector<vector<unsigned char>>Imagetemp=ImageData;
		Imagetemp = ImageData;
		int min = 32767;
		int max = 0;
		for (UINT i = KernelSize / 2; i < ImageHeight - 2 * (KernelSize / 2); i++)
		{
			for (UINT j = KernelSize / 2; j < ImageWidth - 2 * (KernelSize / 2); j++)
			{
				double cResult = 0;
				int ImgbxIndex = i - KernelSize / 2;//当前像素点相对于核的左上角
				int ImgbyIndex = j - KernelSize / 2;
				int KerbIndex = 0;
				while (ImgbxIndex <= i + KernelSize / 2)//左下角
				{
					//逐行扫描
					for (UINT k = 0; k < KernelSize; k++) {
						cResult += kerneldata[KerbIndex][k] * ImageData[ImgbxIndex][ImgbyIndex];
						ImgbyIndex++;
					}
					ImgbxIndex++;//扫描完一行，向下一行扫描
					KerbIndex++;
					ImgbyIndex = j-KernelSize/2;
				}
				cResult = round(cResult);
				if (cResult < 0)
					cResult=0;
				if (cResult > 255)
					cResult=255;
				Imagetemp[i][j] = (unsigned char)cResult;//卷积完的结果赋值回原来的像素点
			}
		}
		if (operators == Lapacian)
		{
			for(UINT i=0;i<ImageData.size();i++)
				for (UINT j = 0; j < ImageData[i].size(); j++)
				{
					int tmp1 = (int)ImageData[i][j] + (int)Imagetemp[i][j];
					if (tmp1 > 255)
						tmp1 = 255;
					else if (tmp1 < 0)
						tmp1 = 0;
					ImageData[i][j] = (unsigned char)tmp1;
		       }
			ImageData = Imagetemp;
		}
		else
			ImageData = Imagetemp;//拷贝
	}
	break;
	case ConvolutionOperation::Medium:
	{
		vector<vector<unsigned char>>Imagetemp=ImageData;
		int *temp = new int[KernelSize*KernelSize];
		memset(temp, 0, KernelSize*KernelSize);
		for (UINT i = KernelSize / 2; i < ImageHeight - 2 * (KernelSize / 2); i++)
		{
			for (UINT j = KernelSize / 2; j < ImageWidth - 2 * (KernelSize / 2); j++)
			{
				int ImgbxIndex = i - KernelSize / 2;//当前像素点相对于核的左上角
				int ImgbyIndex = j - KernelSize / 2;
				int KerbIndex = 0;
				int tmpIndex = 0;
				while (ImgbxIndex <= i + KernelSize / 2)//左下角
				{
					//逐行扫描
					for (UINT k = 0; k < KernelSize; k++) {
						temp[tmpIndex] = ImageData[ImgbxIndex][ImgbyIndex];
						tmpIndex++;
						ImgbyIndex++;
					}
					ImgbxIndex++;//扫描完一行，向下一行扫描
					KerbIndex++;
					ImgbyIndex = j - KernelSize / 2;
				}
				std::sort(temp, temp + KernelSize * KernelSize, comp);
				unsigned char result = temp[KernelSize*KernelSize / 2];
				Imagetemp[i][j] = result;
			}
		}
		delete temp;
		temp = nullptr;
		ImageData = Imagetemp;
	}
	break;
	default:
	break;
	}
}
errno_t Kernel::Writeimage(unsigned char *dst)
{
	if (ErrorCode == 1)
		return FALSE;
	auto p = &dst[0];
	for (UINT i = 0; i < ImageHeight; i++)
	{
		for (UINT j = 0; j < ImageWidth; j++)
		{
			*p = ImageData[i][j];
			p++;
		}
	}
}
class IMGINFO
{
public:
	IMGINFO() {/*do nothing*/}
	IMGINFO(UINT Height, UINT Width, bool _8_BITS, bool _24_BITS) :Height(Height),Width(Width),_8_BITS(_8_BITS),_24_BITS(_24_BITS) {}
	unsigned int Height;
	unsigned int Width;
	bool _8_BITS;
	bool _24_BITS;
	IMGINFO(const IMGINFO& imf)
	{
		this->Height = imf.Height;
		this->Width = imf.Width;
		this->_24_BITS = imf._24_BITS;
		this->_8_BITS = imf._8_BITS;
	}
};

errno_t ImageProcess(unsigned char *src1, unsigned char *src2, unsigned char *dst, IMGINFO inf1, IMGINFO inf2, bool add = false, bool minus = false)
{
	//accept 24 bits image
	if ((inf1._8_BITS&inf2._24_BITS) != 0)
		return ERRNO_FAIL;
	if ((inf1.Height != inf2.Height) || (inf1.Width != inf2.Width))
		return ERRNO_FAIL;
	if ((inf1._24_BITS != true) || (inf1._24_BITS != true))
		return ERRNO_FAIL;
	int index = 0;
	if (add)
	{
		for (UINT i = 0; i < inf1.Height*inf1.Width; i++)
		{
			int temp = src1[i] + src2[i];
			if (temp > 255)
				temp = 255;
			if (temp < 0)
				temp = 0;
			dst[i] = temp;
		}
		return ERRNO_SUCCESS;
	}
	else if (minus)
	{
		for(UINT i = 0; i < inf1.Height*inf1.Width; i++)
		{
			int temp = src1[i] - src2[i];
			if (temp > 255)
				temp=255;
			if (temp < 0)
				temp=0;
			dst[i] = temp;
		}
		return ERRNO_SUCCESS;
	}
	else return ERRNO_FAIL;
	return ERRNO_FAIL;//unexpected situation
}
#endif