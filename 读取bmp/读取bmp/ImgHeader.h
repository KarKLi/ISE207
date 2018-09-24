/*该头文件没有提供16位RGB555或RGB565的读取方式*/
#pragma once
#pragma warning(disable:4996)
#include <iostream>
#include <fstream>
#include <Windows.h>
#include <malloc.h>
#include <opencv2\highgui\highgui.hpp>
#include <math.h>
#define _CHANNEL_OF_COLOR 3
using namespace std;
using namespace cv;
class Image
{
public:
	Image(const char *filename)
	{
		ImgFilename = filename;
		memset(&bmpheader, 0, sizeof(BITMAPFILEHEADER));//清零
		memset(&bmpinfo, 0, sizeof(BITMAPINFOHEADER));//清零
	};
	Image(string filename)
	{
		ImgFilename = filename;
		memset(&bmpheader, 0, sizeof(BITMAPFILEHEADER));//清零
		memset(&bmpinfo, 0, sizeof(BITMAPINFOHEADER));//清零
	};
	bool BMPRead();
	void BMPCVShow();
	//bool BMPWrite(const char *filename);
	bool WriteInfoToText();
private:
	string ImgFilename;
	BITMAPFILEHEADER bmpheader;
	BITMAPINFOHEADER bmpinfo;
	//bool BMPNone24Write(const char *filename);
};
bool Image::BMPRead()
{
	/*bmpfp为读取文件指针*/
	FILE *bmpfp;
	bmpfp = fopen(ImgFilename.c_str(), "rb");
	if (!bmpfp)
		return false;
	if (!(fread(&bmpheader, sizeof(BITMAPFILEHEADER), 1, bmpfp) && fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, bmpfp)))//读取文件头，信息头
	{
		cerr << "Can't read image file, press any key to exit." << endl;
		cin.ignore();
		getchar();
		fclose(bmpfp);
		return false;
	}
	else if (bmpheader.bfType == 0x4d42)//前两个字头，4d=='B',42=='M'
	{
		cout << "The type of image is: BMP" << endl;
		cout << "其大小为：" << bmpheader.bfSize / 1024 << "KB" << endl;
		cout << "需要偏移的偏移量为：" << bmpheader.bfOffBits <<"字节"<< endl;
		cout << "其宽度为" << bmpinfo.biWidth << "像素" << endl;
		cout << "其高度为" << bmpinfo.biHeight << "像素" << endl;
		cout << "其位数为：" << bmpinfo.biBitCount << endl;
		switch (bmpinfo.biCompression)//压缩方式，用于16位BMP图像
		{
		case 0:
			cout << "没有压缩" << endl;
			break;
		case 1:
			cout << "采用8比特游程编码" << endl;
			break;
		case 2:
			cout << "采用4比特游程编码" << endl;
			break;
		default:
			cout << "采用其他压缩编码方式" << endl;
			break;
		}
		cout << "调色板中的颜色索引数为：" << bmpinfo.biClrUsed<<endl;
		fclose(bmpfp);
		cout << "Press any key and enter to continue..." << endl;
		cin.ignore();
		getchar();
		return true;
	}
	else
	{
		cerr << "This is not a bmp file! Press any key and enter to exit." << endl;
		cin.ignore();
		getchar();
		fclose(bmpfp);
		return false;
	}
}
//bool Image::BMPWrite(const char *filename)//用于真彩色图像的写入
//{
//	if (bmpinfo.biBitCount != 24)//非24位真彩色图像
//	{
//		return BMPNone24Write(filename);
//	}
//	FILE *file = fopen(filename, "wb");
//	if (!file)
//		return false;
//	//initalize the FILE pointer
//	FILE *bmpfp;
//	bmpfp = fopen(ImgFilename.c_str(), "rb");
//	if (!bmpfp)
//		return false;
//	//Get width and height
//	DWORD w, h;
//	w = bmpinfo.biWidth;
//	h = bmpinfo.biHeight;
//	//allocate the buffer memory
//	unsigned char *buff = new unsigned char[bmpinfo.biHeight*bmpinfo.biWidth*_CHANNEL_OF_COLOR];//相当于malloc
//	memset(buff, 0, sizeof(buff));
//	auto p = buff;
//	if (!(fread(&bmpheader, sizeof(BITMAPFILEHEADER), 1, bmpfp) && fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, bmpfp)))
//	{
//		cerr << "Can't read image file, press any key to exit." << endl;
//		cin.ignore();
//		getchar();
//		free(buff);
//		return false;
//	}
//	fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, file);
//	fwrite(&bmpinfo, sizeof(BITMAPINFO), 1, file);
// 	fseek(file,bmpheader.bfOffBits, SEEK_SET);//从0开始，漂移bfOffBits字节，类里的bfOffBits成员代表从0开始漂移到达像素点的字节量
//	fread(buff,w*h*_CHANNEL_OF_COLOR,1, bmpfp);//读取长*宽*颜色数（BGR）到缓冲区buff中
//	cout << "Writing to file, please wait..." << endl;
//	for (unsigned int j = 0; j<h; j++)
//	{
//		for (unsigned int i = 0; i<w * 3; i++)
//		{
//			fwrite(p++, 1, 1, file);//通过buff逐字节写入
//		}
//	}
//	delete buff;//相当于free
//	fclose(bmpfp);
//	fclose(file);
//	cout << "Done!" << endl;
//	cout << "The file name is: " << filename << endl;
//	cout << "Press enter to exit." << endl;
//	system("pause");
//	return true;
//}
void Image::BMPCVShow()
{
	//This function can only be executed in Release mode. I don't know why.
#ifndef _DEBUG
	Mat Img = imread(ImgFilename, CV_LOAD_IMAGE_UNCHANGED);//读入图像到Mat（矩阵）中
	if (Img.data == NULL)
		return;
	cout << "Showing picture..." << endl;
	namedWindow("MyWindow", CV_WINDOW_AUTOSIZE);
	//在MyWindow的窗中中显示存储在img中的图片
	startWindowThread();
	imshow("MyWindow", Img);
	//等待直到有键按下
	waitKey(0);
	destroyAllWindows();
#endif
}
//bool Image::BMPNone24Write(const char *filename)
//{
//	FILE *file = fopen(filename, "wb");
//	if (!file)
//		return false;
//	FILE *fp = fopen(ImgFilename.c_str(), "rb");//二进制读方式打开指定的图像文件
//	if (fp == 0)
//		return 0;
//	//fread(&bmpheader, sizeof(BITMAPFILEHEADER), 1, fp);
//	////定义位图信息头结构变量，读取位图信息头进内存，存放在变量head中
//	//fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, fp); //获取图像宽、高、每像素所占位数等信息
//	fseek(fp, 54, SEEK_SET);//漂移54字节,54为BITMAPFILEHEADER(14字节)+BITMAPINFOHEADER(40字节)的大小
//	LONG bmpWidth = bmpinfo.biWidth;
//	LONG bmpHeight = bmpinfo.biHeight;
//	WORD biBitCount = bmpinfo.biBitCount;//定义变量，计算图像每行像素所占的字节数（必须是4的倍数）
//	RGBQUAD pColorTable[256];
//	//int lineByte = (bmpWidth * biBitCount / 8 + 3) / 4 * 4;//每行字节归四化（内存对齐）
//	//申请调色板所需要的空间，读颜色表进内存
//	fread(&pColorTable, sizeof(RGBQUAD), pow(2,bmpinfo.biBitCount), fp);
//	//申请位图数据所需要的空间，读位图数据进内存
//	unsigned char *pBmpBuf = new unsigned char[lineByte * bmpHeight];
//	memset(pBmpBuf, 0, sizeof(pBmpBuf));
//	auto p = pBmpBuf;
//	//写入文件头和信息头
//	fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, file);
//	fwrite(&bmpinfo, sizeof(BITMAPINFO), 1, file);
//	fseek(file,54, SEEK_SET);//不要用sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER),会莫名其妙多出四个字节，应使用24位下的bmpheader.biOffset量来漂移，经测得，为54
//	//写入调色板
//	for (unsigned int i = 0; i < pow(2, bmpinfo.biBitCount); i++)
//	{
//		fwrite(&pColorTable[i].rgbBlue, sizeof(BYTE), 1, file);
//		fwrite(&pColorTable[i].rgbGreen, sizeof(BYTE), 1, file);
//		fwrite(&pColorTable[i].rgbRed, sizeof(BYTE), 1, file);
//		fwrite(&pColorTable[i].rgbReserved, sizeof(BYTE), 1, file);
//	}
//	fseek(file, bmpheader.bfOffBits,SEEK_SET);//bfOffBits应为1078个字节，因为1078-54=1024=2^10,2^10/4=2^8=256，对应pow(2,bmpinfo.biBitCount)
//	fread(pBmpBuf,lineByte * bmpHeight,1, fp);//行字节数*行数=总像素字节数
//	cout << "Writing to file, please wait..." << endl;
//	for (unsigned int j = 0; j<bmpHeight; j++)
//	{
//		for (unsigned int i = 0; i<lineByte; i++)
//		{
//			fwrite(p++, 1, 1, file);//这里p++的作用和上面是一样的。
//		}
//	}
//	delete pBmpBuf;
//	fclose(fp);//关闭文件
//	fclose(file);//关闭文件
//	cout << "Done!" << endl;
//	cout << "The file name is: " << filename << endl;
//	cout << "Press enter to exit." << endl;
//	system("pause");
//	return true;//写入文件成功
//}
bool Image::WriteInfoToText()
{
	cout << "Enter your filename that you want to write: ";
	char filename[100] = { 0 };
	cin >> filename;
	FILE *file = fopen(filename, "w");
	if (!file)
		return false;
	FILE *fp = fopen(ImgFilename.c_str(), "rb");//二进制读方式打开指定的图像文件
	if (fp == 0)
		return 0;
	fprintf(file,"The type of image is: BMP\n");
	fprintf(file, "其大小为：%d KB\n",bmpheader.bfSize / 1024);
	fprintf(file,"需要偏移的偏移量为：%d 字节\n",bmpheader.bfOffBits);
	fprintf(file,"其宽度为 %d 像素\n",bmpinfo.biWidth);
	fprintf(file,"其高度为 %d 像素\n",bmpinfo.biHeight);
	fprintf(file,"其位数为：%d\n",bmpinfo.biBitCount);
	switch (bmpinfo.biCompression)//压缩方式，用于16位BMP图像
	{
	case 0:
		fprintf(file,"没有压缩\n");
		break;
	case 1:
		fprintf(file, "采用8比特游程编码\n");
		break;
	case 2:
		fprintf(file, "采用4比特游程编码\n");
		break;
	default:
		fprintf(file,"采用其他压缩编码方式\n");
		break;
	}
	fprintf(file,"调色板中的颜色索引数为： %d\n",bmpinfo.biClrUsed);
	fprintf(file, "像素数据为：\n");
	//fread(&bmpheader, sizeof(BITMAPFILEHEADER), 1, fp);
//定义位图信息头结构变量，读取位图信息头进内存，存放在变量head中
//fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, fp); //获取图像宽、高、每像素所占位数等信息
	fseek(fp, 54, SEEK_SET);//漂移54字节
	//int counts = 0;
	LONG bmpWidth = bmpinfo.biWidth;
	LONG bmpHeight = bmpinfo.biHeight;
	WORD biBitCount = bmpinfo.biBitCount;//定义变量，计算图像每行像素所占的字节数（必须是4的倍数）
	RGBQUAD pColorTable[256];
	int lineByte = (bmpWidth * biBitCount / 8 + 3) / 4 * 4;//每行字节归四化
	if (biBitCount < 16 || (biBitCount == 16 && bmpinfo.biCompression == BI_BITFIELDS))//有调色板的非24位，有压缩的16位彩色图
	{
		//申请颜色表所需要的空间，读颜色表进内存
		fread(&pColorTable, sizeof(RGBQUAD), pow(2, biBitCount), fp);
	}
	//申请位图数据所需要的空间，读位图数据进内存
	unsigned char *pBmpBuf = new unsigned char[lineByte * bmpHeight];
	memset(pBmpBuf, 0, sizeof(pBmpBuf));
	auto p = pBmpBuf;
	fseek(fp, bmpheader.bfOffBits, SEEK_SET);//应为1078个字节（16位）/54个字节（24位）
	fread(pBmpBuf, lineByte * bmpHeight, 1, fp);
	cout << "Writing to file, please wait..." << endl;
	switch (biBitCount)
	{
	case 24:
		for (unsigned int j = 0; j < bmpHeight; j++)
		{
			for (unsigned int i = 0; i < lineByte; i += 3)
			{
				fprintf(file, "[%d ", (int)*p);
				p++;
				fprintf(file, "%d ", (int)*p);
				p++;
				fprintf(file, "%d]", (int)*p);
				p++;
			}
			fprintf(file, "\n");
		}
		break;
	case 16:
		if (bmpinfo.biCompression == BI_RGB)
		{
			fseek(fp, bmpheader.bfOffBits, SEEK_SET);//应为1078个字节（16位）/54个字节（24位）
			WORD color = 0;
			for (unsigned int j = 0; j < bmpHeight; j++)
			{
				for (unsigned int i = 0; i < lineByte; i += 2)
				{
					fread(&color, sizeof(WORD), 1, fp);
					int temp = color;
					int R = temp << 1;
					R = R & 0xf800;
					R = R >> 11;
					temp = color;
					int G = temp << 1;
					G = G & 0x07c0;
					G = G >> 6;
					temp = color;
					int B = temp << 1;
					B = B & 0x003e;
					B = B >> 1;
					fprintf(file, "[%d %d %d] ", R, G, B);
				}
				fprintf(file, "\n");
			}
		}
		else
		{
			printf("Not support RGB(565) format...\n");
			system("pause");
			exit(EXIT_FAILURE);
		}
		break;
	case 8:
	case 4:
	case 2:
		for (unsigned int j = 0; j < bmpHeight; j++)
		{
			for (unsigned int i = 0; i < bmpWidth; i++)
			{
				fprintf(file, "[%d %d %d] ", pColorTable[(int)*p].rgbRed, pColorTable[(int)*p].rgbGreen, pColorTable[(int)*p].rgbBlue);
			}
			fprintf(file, "\n");
		}
		break;
	default:
		fprintf(stderr, "Error!!\n");
		system("pause");
		exit(EXIT_FAILURE);
		break;
	}
	delete pBmpBuf;
	fclose(fp);//关闭文件
	fclose(file);//关闭文件
	cout << "Done!" << endl;
	cout << "The file name is: " << filename << endl;
	cout << "Press enter to exit." << endl;
	system("pause");
	return true;//写入文件成功
}