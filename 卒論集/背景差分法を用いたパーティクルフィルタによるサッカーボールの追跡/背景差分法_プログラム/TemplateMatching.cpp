#include <stdio.h>
#include <cv.h>
#include <cxcore.h>		//	フレーム画像の取得に必要	
#include <highgui.h>	//	画像の表示、及び読み込みに必要

#define THRESHOLD_MAX_VALUE		255
#define THRESHOLD				60

#define WIDTH	960		//	キャプチャ画像の横幅
#define HEIGHT	540		//	キャプチャ画像の縦幅

int main( int argc, char **argv ) { 
	int key;							//	キー入力用の変数
	CvCapture *capture = NULL;			//	カメラキャプチャ用の構造体
	IplImage *frameImage;				//	キャプチャ画像用IplqImage
	CvVideoWriter *writer = NULL;		//	動画書き出し用の構造体

	//int levels = 165;

	//	画像を生成する	IPL_DEPTH_8U, 1←8 bit 1 ch(256濃淡階調画像)
	IplImage *backgroundImage = cvCreateImage( cvSize( WIDTH, HEIGHT ), IPL_DEPTH_8U, 1 );	//	背景画像用IplImage
	IplImage *grayImage       = cvCreateImage( cvSize( WIDTH, HEIGHT ), IPL_DEPTH_8U, 1 );	//	グレースケール画像用IplImage
	IplImage *differenceImage = cvCreateImage( cvSize( WIDTH, HEIGHT ), IPL_DEPTH_8U, 1 );	//	差分画像用IplImage
	IplImage *binaryImage     = cvCreateImage( cvSize( WIDTH, HEIGHT ), IPL_DEPTH_8U, 1 );	//	2値化画像用IplImage

	char windowNameCapture[]    = "Capture"; 			//	キャプチャした画像を表示するウィンドウの名前
	char windowNameDifference[] = "Difference";		//	背景差分結果を表示するウィンドウの名前
	char windowNameBinary[]     = "Binary";				//	2値化を表示するウィンドウの名前
	char windowNameGray[]     = "Gray";				//	2値化を表示するウィンドウの名前

	//	aviファイル設定
	double fps = 29.0;
	//	Xvidのコーデックが必要、もしくはaviに再変換の必要有り
	writer = cvCreateVideoWriter( "data2_binary.avi", CV_FOURCC( 'X','V','I','D' ), fps, cvSize( WIDTH, HEIGHT ), 0 );

	//	カメラを初期化する
	if ( ( capture = cvCreateFileCapture( "data2.avi" ) ) == NULL ) {
		//	カメラが見つからなかった場合
		printf( "File Not Found\n" );
		return -1;
	}

	//	ウィンドウを生成する
	cvNamedWindow( windowNameCapture, CV_WINDOW_AUTOSIZE );
	cvNamedWindow( windowNameDifference, CV_WINDOW_AUTOSIZE );
	cvNamedWindow( windowNameBinary, CV_WINDOW_AUTOSIZE );
	cvNamedWindow( windowNameGray, CV_WINDOW_AUTOSIZE );

	//	初期背景を設定するための画像取得

	frameImage = cvLoadImage( "bg3.bmp" );
	if(frameImage == NULL){
		printf( "Can't Get\n" );
		return -1;
	}

	//	frameImageをグレースケール化し、背景画像とする
	//	8 bit 1 ch (256階調濃淡) 画像を作る
	cvCvtColor( frameImage, backgroundImage, CV_BGR2GRAY );

	//	メインループ
	while ( 1 ) {

		//	captureの入力画像1フレームをframeImageに格納する
		frameImage = cvQueryFrame( capture );

		if ( frameImage == NULL ){
			break;
		}

		//	frameImageをグレースケール化したものを、grayImageに格納する
		cvCvtColor( frameImage, grayImage, CV_BGR2GRAY );
		//	grayImageと背景画像との差分をとる
		cvSub( grayImage, backgroundImage, differenceImage );

		//cvCvtColor( frameImage, grayImage, CV_BGR2GRAY );

		//	cvTreshold( 元画像、処理画像、閾値、max値、変換方法 )
		//	CV_THRESH_BINARY : 閾値より上はmax値(白)，それ以外は0(黒)。
		cvThreshold( differenceImage, binaryImage, THRESHOLD, THRESHOLD_MAX_VALUE, CV_THRESH_BINARY );

		/*
		if ( differenceImage->origin == 0 ) {
		//　左上が原点の場合
		cvFlip( differenceImage, differenceImage, 0 );
		}
		*/

		//	1画面分の書込
		cvWriteFrame( writer, binaryImage ); 


		//	画像を表示する
		cvShowImage( windowNameCapture, frameImage );
		cvShowImage( windowNameDifference, differenceImage );
		cvShowImage( windowNameBinary, binaryImage );
		cvShowImage( windowNameGray, grayImage );

		//	キー入力判定
		key = cvWaitKey( 33 );
		if ( key == 'q' ) {
			//	'q'キーが押されたらループを抜ける
			break;
		} 
	}

	//	キャプチャを解放する
	cvReleaseCapture( &capture );
	cvReleaseVideoWriter( &writer );
	//	メモリを解放する
	cvReleaseImage( &backgroundImage );
	cvReleaseImage( &grayImage );
	cvReleaseImage( &differenceImage );
	cvReleaseImage( &binaryImage );
	//	ウィンドウを破棄する
	cvDestroyWindow( windowNameCapture );
	cvDestroyWindow( windowNameDifference );
	cvDestroyWindow( windowNameBinary );
	cvDestroyWindow( windowNameGray );

	return 0;
}