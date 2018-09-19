/*��ͷ�ļ�û���ṩ16λRGB555��RGB565�Ķ�ȡ��ʽ*/
#pragma once
#pragma warning(disable:4996)
#include <iostream>
#include <Windows.h>
#include <malloc.h>
#include <opencv2\highgui\highgui.hpp>
#define _CHANNEL_OF_COLOR 3
using namespace std;
using namespace cv;
class Image
{
public:
	Image(const char *filename)
	{
		ImgFilename = filename;
		memset(&bmpheader, 0, sizeof(BITMAPFILEHEADER));//����
		memset(&bmpinfo, 0, sizeof(BITMAPINFOHEADER));//����
	};
	Image(string filename)
	{
		ImgFilename = filename;
		memset(&bmpheader, 0, sizeof(BITMAPFILEHEADER));//����
		memset(&bmpinfo, 0, sizeof(BITMAPINFOHEADER));//����
	};
	bool BMPRead();
	void BMPCVShow();
	bool BMPWrite(const char *filename);
	bool WritePixelToText();
private:
	string ImgFilename;
	BITMAPFILEHEADER bmpheader;
	BITMAPINFOHEADER bmpinfo;
	bool BMPNone24Write(const char *filename);
};
bool Image::BMPRead()
{
	/*bmpfpΪ��ȡ�ļ�ָ��*/
	FILE *bmpfp;
	bmpfp = fopen(ImgFilename.c_str(), "rb");
	if (!bmpfp)
		return false;
	if (!(fread(&bmpheader, sizeof(BITMAPFILEHEADER), 1, bmpfp) && fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, bmpfp)))//��ȡ�ļ�ͷ����Ϣͷ
	{
		cerr << "Can't read image file, press any key to exit." << endl;
		cin.ignore();
		getchar();
		fclose(bmpfp);
		return false;
	}
	else if (bmpheader.bfType == 0x4d42)//ǰ������ͷ��4d=='B',42=='M'
	{
		cout << "The type of image is: BMP" << endl;
		cout << "���СΪ��" << bmpheader.bfSize / 1024 << "KB" << endl;
		cout << "��Ҫƫ�Ƶ�ƫ����Ϊ��" << bmpheader.bfOffBits <<"�ֽ�"<< endl;
		cout << "����Ϊ" << bmpinfo.biWidth << "����" << endl;
		cout << "��߶�Ϊ" << bmpinfo.biHeight << "����" << endl;
		cout << "��λ��Ϊ��" << bmpinfo.biBitCount << endl;
		switch (bmpinfo.biCompression)//ѹ����ʽ������16λBMPͼ��
		{
		case 0:
			cout << "û��ѹ��" << endl;
			break;
		case 1:
			cout << "����8�����γ̱���" << endl;
			break;
		case 2:
			cout << "����4�����γ̱���" << endl;
			break;
		default:
			cout << "��������ѹ�����뷽ʽ" << endl;
			break;
		}
		cout << "��ɫ���е���ɫ������Ϊ��" << bmpinfo.biClrUsed<<endl;
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
bool Image::BMPWrite(const char *filename)//�������ɫͼ���д��
{
	if (bmpinfo.biBitCount != 24)//��24λ���ɫͼ��
	{
		return BMPNone24Write(filename);
	}
	FILE *file = fopen(filename, "wb");
	if (!file)
		return false;
	//initalize the FILE pointer
	FILE *bmpfp;
	bmpfp = fopen(ImgFilename.c_str(), "rb");
	if (!bmpfp)
		return false;
	//Get width and height
	DWORD w, h;
	w = bmpinfo.biWidth;
	h = bmpinfo.biHeight;
	//allocate the buffer memory
	unsigned char *buff = new unsigned char[bmpinfo.biHeight*bmpinfo.biWidth*_CHANNEL_OF_COLOR];//�൱��malloc
	memset(buff, 0, sizeof(buff));
	auto p = buff;
	if (!(fread(&bmpheader, sizeof(BITMAPFILEHEADER), 1, bmpfp) && fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, bmpfp)))
	{
		cerr << "Can't read image file, press any key to exit." << endl;
		cin.ignore();
		getchar();
		free(buff);
		return false;
	}
	fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, file);
	fwrite(&bmpinfo, sizeof(BITMAPINFO), 1, file);
 	fseek(file,bmpheader.bfOffBits, SEEK_SET);//��0��ʼ��Ư��bfOffBits�ֽڣ������bfOffBits��Ա�����0��ʼƯ�Ƶ��ֽ����������ص�
	fread(buff,w*h*_CHANNEL_OF_COLOR,1, bmpfp);//��ȡ��*��*��ɫ����BGR����������buff��
	cout << "Writing to file, please wait..." << endl;
	for (unsigned int j = 0; j<h; j++)
	{
		for (unsigned int i = 0; i<w * 3; i++)
		{
			fwrite(p++, 1, 1, file);//ͨ��buff���ֽ�д��
		}
	}
	delete buff;//�൱��free
	fclose(bmpfp);
	fclose(file);
	cout << "Done!" << endl;
	cout << "The file name is: " << filename << endl;
	cout << "Press enter to exit." << endl;
	system("pause");
	return true;
}
void Image::BMPCVShow()
{
	//This function can only be executed in Release mode. I don't know why.
#ifndef _DEBUG
	Mat Img = imread(ImgFilename, CV_LOAD_IMAGE_UNCHANGED);//����ͼ��Mat��������
	if (Img.data == NULL)
		return;
	cout << "Showing picture..." << endl;
	namedWindow("MyWindow", CV_WINDOW_AUTOSIZE);
	//��MyWindow�Ĵ�������ʾ�洢��img�е�ͼƬ
	startWindowThread();
	imshow("MyWindow", Img);
	//�ȴ�ֱ���м�����
	waitKey(0);
	destroyAllWindows();
#endif
}
bool Image::BMPNone24Write(const char *filename)
{
	FILE *file = fopen(filename, "wb");
	if (!file)
		return false;
	FILE *fp = fopen(ImgFilename.c_str(), "rb");//�����ƶ���ʽ��ָ����ͼ���ļ�
	if (fp == 0)
		return 0;
	//fread(&bmpheader, sizeof(BITMAPFILEHEADER), 1, fp);
	////����λͼ��Ϣͷ�ṹ��������ȡλͼ��Ϣͷ���ڴ棬����ڱ���head��
	//fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, fp); //��ȡͼ����ߡ�ÿ������ռλ������Ϣ
	fseek(fp, 54, SEEK_SET);//Ư��54�ֽ�,54ΪBITMAPFILEHEADER(14�ֽ�)+BITMAPINFOHEADER(40�ֽ�)�Ĵ�С
	LONG bmpWidth = bmpinfo.biWidth;
	LONG bmpHeight = bmpinfo.biHeight;
	WORD biBitCount = bmpinfo.biBitCount;//�������������ͼ��ÿ��������ռ���ֽ�����������4�ı�����
	RGBQUAD pColorTable[256];
	int lineByte = (bmpWidth * biBitCount / 8 + 3) / 4 * 4;//ÿ���ֽڹ��Ļ����ڴ���룩
	//�����ɫ������Ҫ�Ŀռ䣬����ɫ����ڴ�
	fread(&pColorTable, sizeof(RGBQUAD), pow(2,bmpinfo.biBitCount), fp);
	//����λͼ��������Ҫ�Ŀռ䣬��λͼ���ݽ��ڴ�
	unsigned char *pBmpBuf = new unsigned char[lineByte * bmpHeight];
	memset(pBmpBuf, 0, sizeof(pBmpBuf));
	auto p = pBmpBuf;
	//д���ļ�ͷ����Ϣͷ
	fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, file);
	fwrite(&bmpinfo, sizeof(BITMAPINFO), 1, file);
	fseek(file,54, SEEK_SET);//��Ҫ��sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER),��Ī���������ĸ��ֽڣ�Ӧʹ��24λ�µ�bmpheader.biOffset����Ư�ƣ�����ã�Ϊ54
	//д���ɫ��
	for (unsigned int i = 0; i < pow(2, bmpinfo.biBitCount); i++)
	{
		fwrite(&pColorTable[i].rgbBlue, sizeof(BYTE), 1, file);
		fwrite(&pColorTable[i].rgbGreen, sizeof(BYTE), 1, file);
		fwrite(&pColorTable[i].rgbRed, sizeof(BYTE), 1, file);
		fwrite(&pColorTable[i].rgbReserved, sizeof(BYTE), 1, file);
	}
	fseek(file, bmpheader.bfOffBits,SEEK_SET);//bfOffBitsӦΪ1078���ֽڣ���Ϊ1078-54=1024=2^10,2^10/4=2^8=256����Ӧpow(2,bmpinfo.biBitCount)
	fread(pBmpBuf,lineByte * bmpHeight,1, fp);//���ֽ���*����=�������ֽ���
	cout << "Writing to file, please wait..." << endl;
	for (unsigned int j = 0; j<bmpHeight; j++)
	{
		for (unsigned int i = 0; i<lineByte; i++)
		{
			fwrite(p++, 1, 1, file);//����p++�����ú�������һ���ġ�
		}
	}
	delete pBmpBuf;
	fclose(fp);//�ر��ļ�
	fclose(file);//�ر��ļ�
	cout << "Done!" << endl;
	cout << "The file name is: " << filename << endl;
	cout << "Press enter to exit." << endl;
	system("pause");
	return true;//д���ļ��ɹ�
}
bool Image::WritePixelToText()
{
	cout << "Enter your filename that you want to write: ";
	char filename[100] = { 0 };
	cin >> filename;
	FILE *file = fopen(filename, "wb");
	if (!file)
		return false;
	FILE *fp = fopen(ImgFilename.c_str(), "rb");//�����ƶ���ʽ��ָ����ͼ���ļ�
	if (fp == 0)
		return 0;
	//fread(&bmpheader, sizeof(BITMAPFILEHEADER), 1, fp);
	////����λͼ��Ϣͷ�ṹ��������ȡλͼ��Ϣͷ���ڴ棬����ڱ���head��
	//fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, fp); //��ȡͼ����ߡ�ÿ������ռλ������Ϣ
	fseek(fp, 54, SEEK_SET);//Ư��54�ֽ�
	LONG bmpWidth = bmpinfo.biWidth;
	LONG bmpHeight = bmpinfo.biHeight;
	WORD biBitCount = bmpinfo.biBitCount;//�������������ͼ��ÿ��������ռ���ֽ�����������4�ı�����
	RGBQUAD pColorTable[256];
	int lineByte = (bmpWidth * biBitCount / 8 + 3) / 4 * 4;//ÿ���ֽڹ��Ļ�
	if (biBitCount == 8)
	{
		//������ɫ������Ҫ�Ŀռ䣬����ɫ����ڴ�
		fread(&pColorTable, sizeof(RGBQUAD), 256, fp);
	}
	//����λͼ��������Ҫ�Ŀռ䣬��λͼ���ݽ��ڴ�
	unsigned char *pBmpBuf = new unsigned char[lineByte * bmpHeight];
	memset(pBmpBuf, 0, sizeof(pBmpBuf));
	auto p = pBmpBuf;
	fseek(fp, bmpheader.bfOffBits, SEEK_SET);//ӦΪ1078���ֽڣ�16λ��/54���ֽڣ�24λ��
	fread(pBmpBuf, lineByte * bmpHeight, 1, fp);
	cout << "Writing to file, please wait..." << endl;
	for (unsigned int j = 0; j<bmpHeight; j++)
	{
		for (unsigned int i = 0; i<lineByte; i++)
		{
			fwrite(p++, 1, 1, file);//����p++�����ú�������һ���ġ�
		}
	}
	delete pBmpBuf;
	fclose(fp);//�ر��ļ�
	fclose(file);//�ر��ļ�
	cout << "Done!" << endl;
	cout << "The file name is: " << filename << endl;
	cout << "Press enter to exit." << endl;
	system("pause");
	return true;//д���ļ��ɹ�
}