#pragma warning(disable:4996)
#include "ImgHeader.h"
using namespace std;
using namespace cv;
int main()
{
	string str;
	cout << "Enter the filename you want to read:";
	cin >> str;
	Image bmp(str.c_str());
	bmp.BMPRead();
	bmp.BMPCVShow();
	//bmp.BMPWrite("result.bmp");
	bmp.WriteInfoToText();
	getchar();
	return 0;
}