#pragma warning(disable:4996)
#include "ImgHeader.h"
using namespace std;
using namespace cv;
int main()
{
	Image bmp("16.bmp");
	bmp.BMPRead();
	bmp.BMPCVShow();
	//bmp.BMPWrite("result.bmp");
	bmp.WriteInfoToText();
	getchar();
	return 0;
}