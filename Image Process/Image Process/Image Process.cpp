// Image Process.cpp : 定义应用程序的入口点。
//
#pragma comment(lib,"gdiplus.lib")
#include "stdafx.h"
#include "Image Process.h"

#define MAX_LOADSTRING 100

using namespace Gdiplus;
// 全局变量:
HINSTANCE hInst;                                // 当前实例
WCHAR szTitle[MAX_LOADSTRING];                  // 标题栏文本
WCHAR szWindowClass[MAX_LOADSTRING];            // 主窗口类名
static OPENFILENAME ofn = { 0 };
static wchar_t szFilePath[MAX_PATH] = L"";
static BITMAPFILEHEADER bmpheader;
static BITMAPINFOHEADER bmpinfo;
ULONG_PTR m_gdiplusToken;
int ImageOpen = 0;
int ImageChanged = 0;
int HasThreshold = 0;
int HasGrayScaled = 0;
int OnlyGrayScaled = 0;
int index = 0;//灰度图的像素点个数
int lineByte = 0;
wchar_t temp[100];//temp.bmp的路径
wchar_t temp2[100];//temp2.bmp的路径
wchar_t miao[100];//miao.bmp的路径
wchar_t about[100];//color1_small.jpg的路径
unsigned char *buff=nullptr;//图片缓冲区
unsigned char *graybuff=nullptr;//灰度图缓冲区
// 此代码模块中包含的函数的前向声明:
BOOL    CALLBACK    IsRunasAdmin();
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    About(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    LOG(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    THRESHOLD(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK    POWER(HWND, UINT, WPARAM, LPARAM);
errno_t CALLBACK    Reverse(BITMAPFILEHEADER, BITMAPINFOHEADER);
errno_t CALLBACK    Threshold(int);
errno_t CALLBACK    GrayScale(const wchar_t *);
errno_t CALLBACK    Power(double, double);
errno_t CALLBACK    Log(double);
errno_t CALLBACK    He(void);
void    CALLBACK    ProcessPath();
errno_t CALLBACK    CheckWhetherNumber(const wchar_t *);
int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: 在此处放置代码。

	// 初始化GDI+  
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	Gdiplus::GdiplusStartup(&m_gdiplusToken, &gdiplusStartupInput, NULL);

    // 初始化全局字符串
	BOOL result = IsRunasAdmin();
	LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
	if (result)
	{
		wchar_t *temp=new wchar_t[MAX_LOADSTRING];
		int x = 0;
		while (szTitle[x] != '\0')//wstrcpy
		{
			temp[x] = szTitle[x];
			x++;
		}
		temp[x] = '\0';
		swprintf(szTitle, MAX_LOADSTRING, L"%s(管理员)",temp);
		delete temp;
		temp = nullptr;
	}
    LoadStringW(hInstance, IDC_IMAGEPROCESS, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);
	ProcessPath();
    // 执行应用程序初始化:
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_IMAGEPROCESS));

    MSG msg;
	//路径字符串初始化
    // 主消息循环:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    return (int) msg.wParam;
}

BOOL IsRunasAdmin()
{
	BOOL bElevated = FALSE;
	HANDLE hToken = NULL;

	// Get current process token
	if (!OpenProcessToken(GetCurrentProcess(), TOKEN_QUERY, &hToken))
		return FALSE;

	TOKEN_ELEVATION tokenEle;
	DWORD dwRetLen = 0;

	// Retrieve token elevation information
	if (GetTokenInformation(hToken, TokenElevation, &tokenEle, sizeof(tokenEle), &dwRetLen))
	{
		if (dwRetLen == sizeof(tokenEle))
		{
			bElevated = tokenEle.TokenIsElevated;
		}
	}

	CloseHandle(hToken);
	return bElevated;
}


//
//  函数: MyRegisterClass()
//
//  目标: 注册窗口类。
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_IMAGEPROCESS));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = MAKEINTRESOURCEW(IDC_IMAGEPROCESS);
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

//
//   函数: InitInstance(HINSTANCE, int)
//
//   目标: 保存实例句柄并创建主窗口
//
//   注释:
//
//        在此函数中，我们在全局变量中保存实例句柄并
//        创建和显示主程序窗口。
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // 将实例句柄存储在全局变量中
   FILE *fp;
   if (_wfopen_s(&fp, miao, L"rb") != 0)
   {
	   MessageBox(NULL, L"没有找到miao.bmp，请检查你的设置！", L"警告", MB_OK | MB_ICONWARNING);
	   quick_exit(1);
   }
   fseek(fp, sizeof(BITMAPFILEHEADER), SEEK_SET);
   fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, fp);
   HWND hWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW^WS_THICKFRAME,
      0,0,bmpinfo.biWidth,bmpinfo.biHeight+59, nullptr, nullptr, hInstance, nullptr);
   fclose(fp);
   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, nCmdShow);
   UpdateWindow(hWnd);

   return TRUE;
}

//
//  函数: WndProc(HWND, UINT, WPARAM, LPARAM)
//
//  目标: 处理主窗口的消息。
//
//  WM_COMMAND  - 处理应用程序菜单
//  WM_PAINT    - 绘制主窗口
//  WM_DESTROY  - 发送退出消息并返回
//
//
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	//PAINTSTRUCT ps;
	static HBITMAP hbg=NULL;
	FILE *fp=NULL;
	static Gdiplus::Status status;
    switch (message)
    {
    case WM_COMMAND:
        {
            int wmId = LOWORD(wParam);
            // 分析菜单选择:
            switch (wmId)
            {
			case ID_OPEN: //打开
			{
				ImageChanged = 0;
				ofn.lStructSize = sizeof(ofn);
				ofn.hwndOwner = hWnd;
				ofn.lpstrFilter = L"24位真彩色bmp文件(*.bmp)\0*.bmp\0";//要选择的文件后缀   
				ofn.lpstrInitialDir = L"./";//默认的文件路径   
				ofn.lpstrFile = szFilePath;//存放文件的缓冲区   
				ofn.nMaxFile = sizeof(szFilePath) / sizeof(*szFilePath);
				ofn.nFilterIndex = 0;
				ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_EXPLORER; //标志如果是多选要加上OFN_ALLOWMULTISELECT
				if (!GetOpenFileName(&ofn))
				{
					break;
				}
				//文件名储存在szFilePath中
				if (lstrcmpW(szFilePath, L"") != 0)
				{
					ImageOpen = 1;
					_wfopen_s(&fp, szFilePath, L"rb");
					if (fp == NULL)
						break;
					fread(&bmpheader, sizeof(BITMAPFILEHEADER), 1, fp);
					fread(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, fp);
					if (bmpheader.bfType != 0x4d42 || bmpinfo.biBitCount<24)
					{
						MessageBox(hWnd, L"这不是一个有效的位图文件，该程序目前只支持24位bmp图像哦~", L"警告", MB_ICONWARNING | MB_OK);
						fclose(fp);
						break;
					}
					HasGrayScaled = 0;
					InvalidateRect(hWnd, NULL, TRUE);
				}
			}
				break;
			case ID_SAVE: //保存
			{
				if (ImageOpen == 0)
				{
					MessageBox(hWnd, L"请先打开文件！", L"警告！", MB_OK | MB_ICONWARNING);
					break;
				}
				OPENFILENAME sfn = { 0 };
				wchar_t SaveFileName[100] = { 0 };
				wchar_t SourceFileName[100] = { 0 };
				wsprintf(SaveFileName, L"%s", szFilePath);
				if (!HasGrayScaled)
					wsprintf(SourceFileName, L"%s", szFilePath);
				else
				{
					wchar_t temp[100];
					GetModuleFileName(NULL, temp, 100);
					int len = lstrlen(temp);
					while (temp[len - 1] != L'\\')
					{
						temp[len - 1] = L'\0';
						len--;	
					}
					temp[len - 1] = L'\0';
					if(OnlyGrayScaled==1)
					    wsprintf(SourceFileName, L"%s\\temp.bmp", temp);
					else
						wsprintf(SourceFileName, L"%s\\temp2.bmp", temp);
				}
				sfn.lStructSize = sizeof(ofn);
				sfn.hwndOwner = hWnd;
				sfn.lpstrFilter = L"bmp文件(*.bmp)\0*.bmp\0";//要选择的文件后缀   
				sfn.lpstrInitialDir = L"./";//默认的文件路径   
				sfn.lpstrFile = SaveFileName;//存放文件名的缓冲区   
				sfn.nMaxFile = sizeof(szFilePath) / sizeof(*szFilePath);
				sfn.nFilterIndex = 0;
				sfn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_EXPLORER | OFN_OVERWRITEPROMPT; //标志如果是多选要加上OFN_ALLOWMULTISELECT
				if (!GetSaveFileName(&sfn))
					break;
				
				if (CopyFile(SourceFileName, SaveFileName, FALSE) == FALSE)
					MessageBox(hWnd, L"保存失败！", L"警告", MB_OK | MB_ICONWARNING);
				else
					MessageBox(hWnd, L"保存成功！", L"提示", MB_OK | MB_ICONINFORMATION);
			}
			break;
			case ID_GRAYSCALE: //灰度变换
				if (ImageOpen != 0)
				{
						if (GrayScale(szFilePath) == 0)
						{
							ImageChanged = 1;
							HasGrayScaled = 1;
						}
					InvalidateRect(hWnd, NULL, TRUE);
				}
				else
				MessageBox(hWnd, L"请先打开图像文件！", L"警告", MB_OK | MB_ICONWARNING);
				break;
			case ID_LOG://对数变换
				if (ImageOpen != 0 && HasGrayScaled != 0)
				{
					INT_PTR result=DialogBox(hInst, MAKEINTRESOURCE(IDD_DLOG), hWnd, LOG);
					if (result == TRUE)
						InvalidateRect(hWnd, NULL, TRUE);
				}
				else if(ImageOpen==0)
					MessageBox(hWnd, L"请先打开图像文件！", L"警告", MB_OK | MB_ICONWARNING);
				else if (HasGrayScaled == 0)
				{
						MessageBox(hWnd, L"请先将图片进行灰度化！", L"警告", MB_OK | MB_ICONWARNING);
						return (INT_PTR)FALSE;
				}
				break;
			case ID_THRESHOLD: //二值化处理
				if (ImageOpen != 0 && HasGrayScaled!=0)
				{
					INT_PTR result=DialogBox(hInst, MAKEINTRESOURCE(IDD_THRESHOLD), hWnd, THRESHOLD);
					if(result==TRUE)
					    InvalidateRect(hWnd, NULL, TRUE);
				}
				else if (ImageOpen == 0)
					MessageBox(hWnd, L"请先打开图像文件！", L"警告", MB_OK | MB_ICONWARNING);
				else if (HasGrayScaled == 0)
				{
					MessageBox(hWnd, L"请先将图片进行灰度化！", L"警告", MB_OK | MB_ICONWARNING);
					return (INT_PTR)FALSE;
				}
				break;
			case ID_REVERSE: //图像反转
				if (ImageOpen != 0 && HasGrayScaled != 0)
				{
					Reverse(bmpheader, bmpinfo);
					InvalidateRect(hWnd, NULL, TRUE);
				}
				else if (ImageOpen == 0)
					MessageBox(hWnd, L"请先打开图像文件！", L"警告", MB_OK | MB_ICONWARNING);
				else if (HasGrayScaled == 0)
				{
					MessageBox(hWnd, L"请先将图片进行灰度化！", L"警告", MB_OK | MB_ICONWARNING);
					return (INT_PTR)FALSE;
				}
				break;
			case ID_POWER: //幂次变换
				if (ImageOpen != 0 && HasGrayScaled != 0)
				{
					INT_PTR result= DialogBox(hInst, MAKEINTRESOURCE(IDD_POWER), hWnd, POWER);
					if(result==TRUE)
					    InvalidateRect(hWnd, NULL, TRUE);
				}
				else if (ImageOpen == 0)
					MessageBox(hWnd, L"请先打开图像文件！", L"警告", MB_OK | MB_ICONWARNING);
				else if (HasGrayScaled == 0)
				{
					MessageBox(hWnd, L"请先将图片进行灰度化！", L"警告", MB_OK | MB_ICONWARNING);
					return (INT_PTR)FALSE;
				}
				break;
			case ID_HE: //直方图变换及均衡化
				if (ImageOpen != 0 && HasGrayScaled != 0)
				{
					He();
					InvalidateRect(hWnd, NULL, TRUE);
				}
				else if (ImageOpen == 0)
					MessageBox(hWnd, L"请先打开图像文件！", L"警告", MB_OK | MB_ICONWARNING);
				else if (HasGrayScaled == 0)
				{
					MessageBox(hWnd, L"请先将图片进行灰度化！", L"警告", MB_OK | MB_ICONWARNING);
					return (INT_PTR)FALSE;
				}
				break;
            case IDM_ABOUT:
                DialogBox(hInst, MAKEINTRESOURCE(IDD_ABOUTBOX), hWnd, About);
                break;
            case IDM_EXIT:
			{
				int result = MessageBox(hWnd, L"真的要离开吗？请确认您已经保存您所需要的图片", L"提示", MB_OKCANCEL | MB_ICONINFORMATION);
				if (result == IDOK)
				{
					if (status == Status::Ok)
						Gdiplus::GdiplusShutdown(m_gdiplusToken);
					DeleteFile(temp);
					DeleteFile(temp2);
					DestroyWindow(hWnd);
				}
			}
                break;
            default:
                return DefWindowProc(hWnd, message, wParam, lParam);
            }
        }
        break;
    case WM_PAINT:
        {
            PAINTSTRUCT ps;
            HDC hdc = BeginPaint(hWnd, &ps);
            // TODO: 在此处添加使用 hdc 的任何绘图代码...
			if (ImageOpen!=0)
			{
				RECT rect;
				GetWindowRect(hWnd, &rect);
				if(!MoveWindow(hWnd,rect.left,rect.top,bmpinfo.biWidth,bmpinfo.biHeight+59,TRUE))
					break;
				if (ImageChanged == 0)
				{
					Gdiplus::Image image(szFilePath);
					//把状态储存到全局变量
					status = image.GetLastStatus();
					if (image.GetLastStatus() != Status::Ok)
					{
						MessageBox(hWnd, L"加载图片失败!", L"提示", MB_OK);
						return -1;
					}
					//取得宽度和高度					 
					int width = image.GetWidth();
					int height = image.GetHeight();
					//绘图					 
					Graphics graphics(hdc);
					graphics.DrawImage(&image, 0, 0, width, height);
				}
				else
				{
					if (OnlyGrayScaled == 1)
					{
						Gdiplus::Image image(temp);
						status = image.GetLastStatus();
						if (image.GetLastStatus() != Status::Ok)
						{
							MessageBox(hWnd, L"加载图片失败!", L"提示", MB_OK);
							return -1;
						}
						//取得宽度和高度					 
						int width = image.GetWidth();
						int height = image.GetHeight();
						//绘图					 
						Graphics graphics(hdc);
						graphics.DrawImage(&image, 0, 0, width, height);
					}
					else
					{
						Gdiplus::Image image(temp2);
						status = image.GetLastStatus();
						if (image.GetLastStatus() != Status::Ok)
						{
							MessageBox(hWnd, L"加载图片失败!", L"提示", MB_OK);
							return -1;
						}
						//取得宽度和高度					 
						int width = image.GetWidth();
						int height = image.GetHeight();
						//绘图					 
						Graphics graphics(hdc);
						graphics.DrawImage(&image, 0, 0, width, height);
					}
				}
			}
			else
			{
				Gdiplus::Image image(miao);
				status = image.GetLastStatus();
				if (image.GetLastStatus() != Status::Ok)
				{
					MessageBox(hWnd, L"加载图片失败!", L"提示", MB_OK);
					return -1;
				}
				//取得宽度和高度					 
				int width = image.GetWidth();
				int height = image.GetHeight();
				//绘图					 
				Graphics graphics(hdc);
				graphics.DrawImage(&image, 0, 0, width, height);
			}
            EndPaint(hWnd, &ps);
        }
        break;
	case WM_CLOSE:
	{
		int result = MessageBox(hWnd, L"真的要离开吗？请确认您已经保存您所需要的图片", L"提示", MB_OKCANCEL | MB_ICONINFORMATION);
		if (result == IDOK)
		{
			if (status == Status::Ok)
				Gdiplus::GdiplusShutdown(m_gdiplusToken);
			DeleteFile(temp);
			DeleteFile(temp2);
			delete buff;
			delete graybuff;
			buff = nullptr;
			graybuff = nullptr;//野指针置零
			DestroyWindow(hWnd);
		}
	}
		break;
    case WM_DESTROY:
		MessageBox(NULL, L"感谢您的使用！\n再见！", L"再见~", MB_OK);
		PostQuitMessage(0);
		break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}

// “关于”框的消息处理程序。
INT_PTR CALLBACK About(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	HDC hdc;
	PAINTSTRUCT ps;
	Gdiplus::Image image(about);
	//取得宽度和高度					 
	int width = image.GetWidth();
	int height = image.GetHeight();
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
        return (INT_PTR)TRUE;
		break;
	case WM_PAINT:
	{
		hdc = BeginPaint(hDlg, &ps);
		Graphics graphics(hdc);
		if (image.GetLastStatus() != Status::Ok)
		{
			MessageBox(hDlg, L"加载图片失败!", L"提示", MB_OK);
			return -1;
		}
		//绘图					 
		graphics.DrawImage(&image, 0, 0, width, height);
		EndPaint(hDlg, &ps);
		break;
	}
    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK || LOWORD(wParam) == IDCANCEL)
        {
            EndDialog(hDlg, LOWORD(wParam));
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}
// “对数变换”框的消息处理程序
INT_PTR CALLBACK LOG(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	UNREFERENCED_PARAMETER(lParam);
	switch (message)
	{
	case WM_INITDIALOG:
		return (INT_PTR)TRUE;

	case WM_COMMAND:
	{
		if (LOWORD(wParam) == IDOK)
		{
			wchar_t buff1[10];
			double c=0;
			HWND editWnd1;
			editWnd1 = GetDlgItem(hDlg, IDC_EDIT1);
			GetWindowText(editWnd1, buff1, 10);
			if (!CheckWhetherNumber(buff1))
			{
				MessageBox(hDlg, L"请输入正确的数字！", L"警告", MB_OK | MB_ICONWARNING);
				return (INT_PTR)FALSE;
			}
			else
			    c = _wtof(buff1);
			if (c <= (double)0)
			{
				MessageBox(hDlg, L"输入的数据非法！", L"警告", MB_OK | MB_ICONWARNING);
				return (INT_PTR)FALSE;
			}
			EndDialog(hDlg, LOWORD(wParam));
			Log(c);
			return (INT_PTR)TRUE;
		}
		else if (LOWORD(wParam) == IDCANCEL)
		{
			EndDialog(hDlg, LOWORD(wParam));
			return (INT_PTR)FALSE;
		}
	}
	break;
	}
	return (INT_PTR)FALSE;
}
//二值化处理的消息处理程序
INT_PTR CALLBACK    THRESHOLD(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	HWND editWnd;
	wchar_t buff[100];
	int threshold = -1;
	UNREFERENCED_PARAMETER(lParam);
	switch (message)
	{
	case WM_INITDIALOG:
		return (INT_PTR)TRUE;

	case WM_COMMAND:
		if (LOWORD(wParam) == IDOK)
		{
			editWnd=GetDlgItem(hDlg, IDC_EDIT1);
			GetWindowText(editWnd, buff, 100);
			threshold = _wtof(buff);
			if (threshold < 0 || threshold>255)
			{
				MessageBox(hDlg, L"请输入正确的阈值！(0~255)", L"警告", MB_OK | MB_ICONWARNING);
				return (INT_PTR)FALSE;
			}
			EndDialog(hDlg, LOWORD(wParam));
			if (Threshold(threshold) == 0)
			{
				ImageChanged = 1;
				HasThreshold = 1;
			}
			else
				ImageChanged = 0;
			return (INT_PTR)TRUE;
		}
		else if (LOWORD(wParam) == IDCANCEL)
		{
			EndDialog(hDlg, LOWORD(wParam));
			return (INT_PTR)FALSE;
		}
		break;
	}
	return (INT_PTR)FALSE;
}
errno_t             Threshold(int Threshold)
{
	OnlyGrayScaled = 0;
	FILE *write;
	auto p1 = &buff[0];
	auto p2 = &buff[0];
	unsigned char * graybufftmp = new unsigned char[index];
	for (unsigned int x = 0; x < index; x++) {
		graybufftmp[x] = graybuff[x];
	}
	for (unsigned int x = 0; x < index; x++)//二值化
	{
		if (graybufftmp[x] >= Threshold)
			graybufftmp[x] = 255;
		else
			graybufftmp[x] = 0;
	}
	int indexs = 0;
	for (unsigned int y = 0; y < bmpinfo.biHeight; y++)//灰度生成
	{
		for (unsigned int x = 0; x < lineByte; x += 3)
		{
			int R, G, B;
			if (graybufftmp[indexs] == 255)
			{
				R = 255, G = 255, B = 255;
			}
			else
			{
				R = 0, G = 0, B = 0;
			}
			*p1 = (char)R;
			p1++;
			*p1 = (char)G;
			p1++;
			*p1 = (char)B;
			p1++;
			indexs++;
		}
	}
		errno_t status = _wfopen_s(&write,temp2, L"wb");
		if (status != 0)
		{
			return 1;
		}
		fseek(write, 0, SEEK_SET);
		fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, write);
		fwrite(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, write);
		fwrite(p2, lineByte*bmpinfo.biHeight, 1, write);
		fclose(write);
		return 0;
}
errno_t             GrayScale(const wchar_t *Filename)
{
	FILE *fp;
	FILE *write;
	if (_wfopen_s(&fp, Filename, L"rb") == 0)//成功返回0
	{
		lineByte = (bmpinfo.biWidth * bmpinfo.biBitCount / 8 + 3) / 4 * 4;//每行字节归四化
		buff = new unsigned char[lineByte*bmpinfo.biHeight];
		graybuff = new unsigned char[lineByte*bmpinfo.biHeight / 3];
		memset(buff, 0, sizeof(buff));
		auto p = buff;
		auto p1 = buff;
		auto p2 = buff;
		memset(graybuff, 0, sizeof(graybuff));
		fseek(fp, bmpheader.bfOffBits, SEEK_SET);
		fread(buff, lineByte*bmpinfo.biHeight, 1, fp);
		errno_t status = _wfopen_s(&write, temp, L"wb");
		if (status != 0)
		{
			fclose(fp);
			return 1;
		}
		index = 0;
		for (unsigned int y = 0; y < bmpinfo.biHeight; y++)//灰度生成
		{
			for (unsigned int x = 0; x < lineByte; x += 3)
			{
				int R, G, B;
				R = (int)*p;
				p++;
				G = (int)*p;
				p++;
				B = (int)*p;
				p++;
				graybuff[index] = (char)((R * 306 + G * 591 + B * 117) / 1024);
				index++;
			}
		}
		index = 0;
		for (unsigned int y = 0; y < bmpinfo.biHeight; y++)//灰度写入
		{
			for (unsigned int x = 0; x < lineByte; x += 3)
			{
				int R, G, B;
				R = G = B = graybuff[index];
				*p1 = (char)R;
				p1++;
				*p1 = (char)G;
				p1++;
				*p1 = (char)B;
				p1++;
				index++;
			}
		}
		fseek(write, 0, SEEK_SET);
		fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, write);
		fwrite(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, write);
		fwrite(p2, lineByte*bmpinfo.biHeight, 1, write);
		fclose(write);
		fclose(fp);
		OnlyGrayScaled = 1;
		HasGrayScaled = 1;
	}
	else return 1;
	return 0;
}
//幂次变换的消息处理程序
INT_PTR CALLBACK    POWER(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
	UNREFERENCED_PARAMETER(lParam);
	switch (message)
	{
	case WM_INITDIALOG:
		return (INT_PTR)TRUE;

	case WM_COMMAND:
	{
		if (LOWORD(wParam) == IDOK)
		{
			wchar_t buff1[10];
			wchar_t buff2[10];
			HWND editWnd1;
			HWND editWnd2;
			editWnd1 = GetDlgItem(hDlg, IDC_EDIT2);
			editWnd2 = GetDlgItem(hDlg, IDC_EDIT3);
			GetWindowText(editWnd1, buff1, 10);
			GetWindowText(editWnd2, buff2, 10);
			double c=0;
			double gamma=0;
			if ((!CheckWhetherNumber(buff1))||(!CheckWhetherNumber(buff2)))
			{
				MessageBox(hDlg, L"请输入正确的数字！", L"警告", MB_OK | MB_ICONWARNING);
				return (INT_PTR)FALSE;
			}
			else
			{
				c = _wtof(buff1);
				gamma = _wtof(buff2);
			}
			if (gamma >= 5)
			{
				MessageBox(hDlg, L"您输入的幂次太大，请将幂次控制在5以内", L"警告", MB_OK | MB_ICONWARNING);
				return (INT_PTR)FALSE;
			}
			else
			{
				EndDialog(hDlg, LOWORD(wParam));
				Power(c, gamma);
			}
			return (INT_PTR)TRUE;
		}
		else if (LOWORD(wParam) == IDCANCEL)
		{
			EndDialog(hDlg, LOWORD(wParam));
			return (INT_PTR)FALSE;
		}
	}
		break;
	}
	return (INT_PTR)FALSE;
}
errno_t             Power(double c, double gamma)
{
	OnlyGrayScaled = 0;
	FILE *write;
	auto p1 = &buff[0];
	auto p2 = &buff[0];
	int indexs = 0;
	for (unsigned int y = 0; y < bmpinfo.biHeight; y++)//灰度写入
	{
		for (unsigned int x = 0; x < lineByte; x += 3)
		{
			int R, G, B;
			R = G = B = (int)round((c*pow((double)graybuff[indexs]/(double)255, gamma))*255);
			if (R > 255)R = G = B = 255;
			*p1 = (char)R;
			p1++;
			*p1 = (char)G;
			p1++;
			*p1 = (char)B;
			p1++;
			indexs++;
		}
	}
	errno_t status = _wfopen_s(&write,temp2, L"wb");
	if (status != 0)
	{
		return 1;
	}
	fseek(write, 0, SEEK_SET);
	fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, write);
	fwrite(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, write);
	fwrite(p2, lineByte*bmpinfo.biHeight, 1, write);
	fclose(write);
	return 0;
}
errno_t             Reverse(BITMAPFILEHEADER bmpheader, BITMAPINFOHEADER bmpinfo)
{
	OnlyGrayScaled = 0;
	FILE *write;
	int indexs = 0;
	auto p1 = &buff[0];
	auto p2 = &buff[0];
	for (unsigned int y = 0; y < bmpinfo.biHeight; y++)//灰度写入
	{
		for (unsigned int x = 0; x < lineByte; x += 3)
		{
			int R, G, B;
			R = G = B = 255-(int)graybuff[indexs];
			*p1 = (char)R;
			p1++;
			*p1 = (char)G;
			p1++;
			*p1 = (char)B;
			p1++;
			indexs++;
		}
	}
	errno_t status = _wfopen_s(&write,temp2, L"wb");
	if (status != 0)
	{
		return 1;
	}
	fseek(write, 0, SEEK_SET);
	fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, write);
	fwrite(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, write);
	fwrite(p2, lineByte*bmpinfo.biHeight, 1, write);
	fclose(write);
	return 0;
}
errno_t             Log(double c)
{
	OnlyGrayScaled = 0;
	FILE *write;
	auto p1 = &buff[0];
	auto p2 = &buff[0];
	int index = 0;
	for (unsigned int y = 0; y < bmpinfo.biHeight; y++)//灰度写入
	{
		for (unsigned int x = 0; x < lineByte; x += 3)
		{
			int R, G, B;
			double tmp= c*(log2((double)1 + (double)graybuff[index] / (double)255));
			R = G = B = (int)round(tmp*(double)255);
			if (R > 255)R = G = B = 255;//溢出处理
			*p1 = (char)R;
			p1++;
			*p1 = (char)G;
			p1++;
			*p1 = (char)B;
			p1++;
			index++;
		}
	}
	errno_t status = _wfopen_s(&write,temp2, L"wb");
	if (status != 0)
	{
		return 1;
	}
	fseek(write, 0, SEEK_SET);
	fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, write);
	fwrite(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, write);
	fwrite(p2, lineByte*bmpinfo.biHeight, 1, write);
	fclose(write);
	return 0;
}
errno_t             He(void)
{

	OnlyGrayScaled = 0;
	FILE *write;
	auto p1 = &buff[0];
	auto p2 = &buff[0];
	int gray[256] = { 0 };
	for (int x = 0; x < index; x++)
	{
		gray[graybuff[x]]++;
	}
	double gray_prob[256] = { 0 };
	for (int i = 0; i < 256; i++)
	{
		gray_prob[i] = ((double)gray[i] / index);
	}
	double gray_distr[256] = { 0 };
	gray_distr[0] = gray_prob[0];
	for (int i = 1; i < 256; i++)
	{
		gray_distr[i] = gray_distr[i - 1] + gray_prob[i];
	}
	int gray_equal[256] = { 0 };
	for (int i = 0; i < 256; i++)
	{
		gray_equal[i] = (char)(255 * gray_distr[i] + 0.5);
	}
	for (int x = 0; x < index; x++)
	{
		graybuff[x] = gray_equal[graybuff[x]];
	}
    int indexs = 0;
	for (unsigned int y = 0; y < bmpinfo.biHeight; y++)//灰度写入
	{
		for (unsigned int x = 0; x < lineByte; x += 3)
		{
			int R, G, B;
			R = G = B = (int)graybuff[indexs];
			*p1 = (char)R;
			p1++;
			*p1 = (char)G;
			p1++;
			*p1 = (char)B;
			p1++;
			indexs++;
		}
	}
	errno_t status = _wfopen_s(&write,temp2, L"wb");
	if (status != 0)
	{
		return 1;
	}
	fseek(write, 0, SEEK_SET);
	fwrite(&bmpheader, sizeof(BITMAPFILEHEADER), 1, write);
	fwrite(&bmpinfo, sizeof(BITMAPINFOHEADER), 1, write);
	fwrite(p2, lineByte*bmpinfo.biHeight, 1, write);
	fclose(write);
	return 0;
}
void                ProcessPath()
{
	wchar_t temp1[100];
	GetModuleFileName(NULL, temp1, 100);
	int len = lstrlen(temp1);
	while (temp1[len - 1] != L'\\')
	{
		temp1[len - 1] = L'\0';
		len--;
	}
	temp1[len - 1] = L'\0';
	wsprintf(temp, L"%s\\temp.bmp", temp1);
	wsprintf(temp2, L"%s\\temp2.bmp", temp1);
	wsprintf(miao, L"%s\\miao.bmp", temp1);
	wsprintf(about, L"%s\\color1_small.jpg", temp1);
}
errno_t                CheckWhetherNumber(const wchar_t *str)
{
	const wchar_t *indexs = L"0123456789";
	int DotExist = 0;
	int len = 0;
	int x = 0;
	while (str[x] != '\0')
	{
		x++;
	}
	len = x;
	if (len > 1000)
		return FALSE;
	else
	{
		for (int i = 0; i < len; i++)
		{
			int flags = 0;
			for (int j = 0; j < 11; j++)
			{
				if (str[i] == indexs[j])
				{
					flags = 1;//matched
					break;
				}
				else if (str[i] == L'.')
				{
					DotExist++;
					flags = 1;
					break;
				}
				if (DotExist > 1)
					return FALSE;
			}
			if (flags == 0)
				return FALSE;
		}
		return TRUE;
	}
}