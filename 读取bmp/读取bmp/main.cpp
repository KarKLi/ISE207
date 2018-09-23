#pragma warning(disable:4996)
#include "ImgHeader.h"
using namespace std;
using namespace cv;
int main()
{
	Image bmp("C:\\Users\\KarK\\Pictures\\test.bmp");
	bmp.BMPRead();
	bmp.BMPCVShow();
	bmp.BMPWrite("result.bmp");
	bmp.WritePixelToText();
	getchar();
	return 0;
}