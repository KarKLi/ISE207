#include <opencv2\highgui\highgui.hpp>
#include <iostream>
using namespace cv;
using namespace  std;
int main(int argc, const char** argv)
{
	Mat img = imread("C:\\Users\\KarK\\Pictures\\test.bmp", CV_LOAD_IMAGE_UNCHANGED);
	if(img.empty())
	{
		cout << "图像加载失败！"<< endl;
		//system("pause");
		return -1;
	}
	//创建一个名字为MyWindow的窗口
	namedWindow("MyWindow", CV_WINDOW_AUTOSIZE);
	//在MyWindow的窗中中显示存储在img中的图片
	imshow("MyWindow", img);
	//等待直到有键按下
	waitKey(1);
	//销毁MyWindow的窗口
	destroyWindow("MyWindow");
	return 0;
}