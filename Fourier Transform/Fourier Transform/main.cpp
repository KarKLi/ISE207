<<<<<<< HEAD
#include "FT.h"
#include <opencv2\opencv.hpp>
#include <opencv2\highgui.hpp>
using namespace cv;
int main()
{
	Mat src = imread("C:\\Users\\KarK\\Desktop\\image\\peppers.bmp", 1);
	if (src.empty())
		exit(1);
	if (src.type() != CV_8U)
		cvtColor(src, src, COLOR_BGR2GRAY);
	vector<vector<FT::Complex>>Image;
	Image.resize(src.rows);
	for (uint i = 0; i < src.rows; i++)
		Image[i].resize(src.cols);
	for (uint i = 0; i < src.rows; ++i)
		for (uint j = 0; j < src.cols; ++j)
			Image[i][j] = (double)src.at<uchar>(i, j);
	FT::DFTTwoDim dfts(src.rows, src.cols);
	Mat dst;
	src.copyTo(dst);
	vector<vector<double>>Dst(src.rows, vector<double>(src.cols));
	vector<vector<FT::Complex>>Dstcomp(src.rows, vector<FT::Complex>(src.cols));
	dfts.LoadData(Image);
	dfts._FFT2(true);
	dfts.GetAmplitudeArray(Dst, src.rows, src.cols);
	for (int i = 0; i < src.rows; i++)
	{
		uchar *p = dst.ptr<uchar>(i);
		for (int j = 0; j < src.cols; j++)
		{
			int tmp = round(150.0*(255.0*log2(1.0 + Dst[i][j]/255)));
			if (tmp > 255)
				tmp = 255;
			p[j] = (uchar)tmp;
		}
	}
	/*namedWindow("Amplitude");
	resizeWindow("Ampulitude", 1024, 1024);
	imshow("Amplitude", dst);
	waitKey(0);
	destroyAllWindows();*/
	imwrite("C:\\Users\\KarK\\Desktop\\image\\DFTCentralAmplitude.bmp", dst);
	dfts.Uncentralized();
	dfts.GetAmplitudeArray(Dst, src.rows, src.cols);
	for (int i = 0; i < src.rows; i++)
	{
		uchar *p = dst.ptr<uchar>(i);
		for (int j = 0; j < src.cols; j++)
		{
			int tmp = round(150.0*(255.0*log2(1.0 + Dst[i][j] / 255)));
			if (tmp > 255)
				tmp = 255;
			p[j] = (uchar)tmp;
		}
	}
	imwrite("C:\\Users\\KarK\\Desktop\\image\\DFTAmplitude.bmp", dst);
	dfts.GetPhaseArray(Dst, src.rows, src.cols);
	for (int i = 0; i < src.rows; i++)
	{
		uchar *p = dst.ptr<uchar>(i);
		for (int j = 0; j < src.cols; j++)
		{
			p[j] = (uchar)round(Dst[i][j]);
		}
	}
	imwrite("C:\\Users\\KarK\\Desktop\\image\\DFTPhase.bmp", dst);
	dfts._IFFT2();
	dfts.WriteData(Dstcomp);
	for (int i = 0; i < src.rows; i++)
	{
		uchar *p = dst.ptr<uchar>(i);
		for (int j = 0; j < src.cols; j++)
		{
			p[j] = (uchar)round(Dstcomp[i][j].real);
		}
	}
	imwrite("C:\\Users\\KarK\\Desktop\\image\\IDFT.bmp", dst);
	return 0;
=======
#include "FT.h"
#include <opencv2\opencv.hpp>
#include <opencv2\highgui.hpp>
using namespace cv;
int main()
{
	Mat src = imread("C:\\Users\\Admin\\Desktop\\lena.bmp", 0);
	if (src.empty())
		exit(1);
	vector<vector<FT::Complex>>Image;
	Image.resize(src.rows);
	for (uint i = 0; i < src.rows; i++)
		Image[i].resize(src.cols);
	for (uint i = 0; i < src.rows; ++i)
		for (uint j = 0; j < src.cols; ++j)
			Image[i][j] = (double)src.at<uchar>(i, j);
	FT::DFTTwoDim dfts(src.rows, src.cols);
	dfts.LoadData(Image);
	dfts.DFT(false);
	Mat dst;
	src.copyTo(dst);
	vector<vector<double>>Dst(src.rows, vector<double>(src.cols));
	vector<vector<FT::Complex>>Dstcomp(src.rows, vector<FT::Complex>(src.cols));
	dfts.GetAmplitudeArray(Dst, src.rows, src.cols);
	for (int i = 0; i < src.rows; i++)
	{
		uchar *p = dst.ptr<uchar>(i);
		for (int j = 0; j < src.cols; j++)
		{
			int tmp=round(150.0*(255.0*log2(1.0 + Dst[i][j] / 255)));
			if (tmp > 255)
				tmp = 255;
			p[j] = (uchar)tmp;
		}
	}
	imwrite("C:\\Users\\Admin\\Desktop\\DFTAmplitude.bmp",dst);
	dfts.InverseDFT();
	dfts.WriteData(Dstcomp);
	for (int i = 0; i < src.rows; i++)
	{
		uchar *p = dst.ptr<uchar>(i);
		for (int j = 0; j < src.cols; j++)
			p[j] = (uchar)Dstcomp[i][j].real;
	}
	imwrite("C:\\Users\\Admin\\Desktop\\IDFT.bmp", dst);
  	return 0;
>>>>>>> ac8337cffc1b28815bc76664cd8dc49c62a1b627
}