#include <opencv2\highgui\highgui.hpp>
#include <iostream>
using namespace cv;
using namespace  std;
int main(int argc, const char** argv)
{
	Mat img = imread("C:\\Users\\KarK\\Pictures\\test.bmp", CV_LOAD_IMAGE_UNCHANGED);
	if(img.empty())
	{
		cout << "ͼ�����ʧ�ܣ�"<< endl;
		//system("pause");
		return -1;
	}
	//����һ������ΪMyWindow�Ĵ���
	namedWindow("MyWindow", CV_WINDOW_AUTOSIZE);
	//��MyWindow�Ĵ�������ʾ�洢��img�е�ͼƬ
	imshow("MyWindow", img);
	//�ȴ�ֱ���м�����
	waitKey(1);
	//����MyWindow�Ĵ���
	destroyWindow("MyWindow");
	return 0;
}