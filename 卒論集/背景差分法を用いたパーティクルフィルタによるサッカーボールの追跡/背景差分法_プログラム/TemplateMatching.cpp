#include <stdio.h>
#include <cv.h>
#include <cxcore.h>		//	�t���[���摜�̎擾�ɕK�v	
#include <highgui.h>	//	�摜�̕\���A�y�ѓǂݍ��݂ɕK�v

#define THRESHOLD_MAX_VALUE		255
#define THRESHOLD				60

#define WIDTH	960		//	�L���v�`���摜�̉���
#define HEIGHT	540		//	�L���v�`���摜�̏c��

int main( int argc, char **argv ) { 
	int key;							//	�L�[���͗p�̕ϐ�
	CvCapture *capture = NULL;			//	�J�����L���v�`���p�̍\����
	IplImage *frameImage;				//	�L���v�`���摜�pIplqImage
	CvVideoWriter *writer = NULL;		//	���揑���o���p�̍\����

	//int levels = 165;

	//	�摜�𐶐�����	IPL_DEPTH_8U, 1��8 bit 1 ch(256�Z�W�K���摜)
	IplImage *backgroundImage = cvCreateImage( cvSize( WIDTH, HEIGHT ), IPL_DEPTH_8U, 1 );	//	�w�i�摜�pIplImage
	IplImage *grayImage       = cvCreateImage( cvSize( WIDTH, HEIGHT ), IPL_DEPTH_8U, 1 );	//	�O���[�X�P�[���摜�pIplImage
	IplImage *differenceImage = cvCreateImage( cvSize( WIDTH, HEIGHT ), IPL_DEPTH_8U, 1 );	//	�����摜�pIplImage
	IplImage *binaryImage     = cvCreateImage( cvSize( WIDTH, HEIGHT ), IPL_DEPTH_8U, 1 );	//	2�l���摜�pIplImage

	char windowNameCapture[]    = "Capture"; 			//	�L���v�`�������摜��\������E�B���h�E�̖��O
	char windowNameDifference[] = "Difference";		//	�w�i�������ʂ�\������E�B���h�E�̖��O
	char windowNameBinary[]     = "Binary";				//	2�l����\������E�B���h�E�̖��O
	char windowNameGray[]     = "Gray";				//	2�l����\������E�B���h�E�̖��O

	//	avi�t�@�C���ݒ�
	double fps = 29.0;
	//	Xvid�̃R�[�f�b�N���K�v�A��������avi�ɍĕϊ��̕K�v�L��
	writer = cvCreateVideoWriter( "data2_binary.avi", CV_FOURCC( 'X','V','I','D' ), fps, cvSize( WIDTH, HEIGHT ), 0 );

	//	�J����������������
	if ( ( capture = cvCreateFileCapture( "data2.avi" ) ) == NULL ) {
		//	�J������������Ȃ������ꍇ
		printf( "File Not Found\n" );
		return -1;
	}

	//	�E�B���h�E�𐶐�����
	cvNamedWindow( windowNameCapture, CV_WINDOW_AUTOSIZE );
	cvNamedWindow( windowNameDifference, CV_WINDOW_AUTOSIZE );
	cvNamedWindow( windowNameBinary, CV_WINDOW_AUTOSIZE );
	cvNamedWindow( windowNameGray, CV_WINDOW_AUTOSIZE );

	//	�����w�i��ݒ肷�邽�߂̉摜�擾

	frameImage = cvLoadImage( "bg3.bmp" );
	if(frameImage == NULL){
		printf( "Can't Get\n" );
		return -1;
	}

	//	frameImage���O���[�X�P�[�������A�w�i�摜�Ƃ���
	//	8 bit 1 ch (256�K���Z�W) �摜�����
	cvCvtColor( frameImage, backgroundImage, CV_BGR2GRAY );

	//	���C�����[�v
	while ( 1 ) {

		//	capture�̓��͉摜1�t���[����frameImage�Ɋi�[����
		frameImage = cvQueryFrame( capture );

		if ( frameImage == NULL ){
			break;
		}

		//	frameImage���O���[�X�P�[�����������̂��AgrayImage�Ɋi�[����
		cvCvtColor( frameImage, grayImage, CV_BGR2GRAY );
		//	grayImage�Ɣw�i�摜�Ƃ̍������Ƃ�
		cvSub( grayImage, backgroundImage, differenceImage );

		//cvCvtColor( frameImage, grayImage, CV_BGR2GRAY );

		//	cvTreshold( ���摜�A�����摜�A臒l�Amax�l�A�ϊ����@ )
		//	CV_THRESH_BINARY : 臒l�����max�l(��)�C����ȊO��0(��)�B
		cvThreshold( differenceImage, binaryImage, THRESHOLD, THRESHOLD_MAX_VALUE, CV_THRESH_BINARY );

		/*
		if ( differenceImage->origin == 0 ) {
		//�@���オ���_�̏ꍇ
		cvFlip( differenceImage, differenceImage, 0 );
		}
		*/

		//	1��ʕ��̏���
		cvWriteFrame( writer, binaryImage ); 


		//	�摜��\������
		cvShowImage( windowNameCapture, frameImage );
		cvShowImage( windowNameDifference, differenceImage );
		cvShowImage( windowNameBinary, binaryImage );
		cvShowImage( windowNameGray, grayImage );

		//	�L�[���͔���
		key = cvWaitKey( 33 );
		if ( key == 'q' ) {
			//	'q'�L�[�������ꂽ�烋�[�v�𔲂���
			break;
		} 
	}

	//	�L���v�`�����������
	cvReleaseCapture( &capture );
	cvReleaseVideoWriter( &writer );
	//	���������������
	cvReleaseImage( &backgroundImage );
	cvReleaseImage( &grayImage );
	cvReleaseImage( &differenceImage );
	cvReleaseImage( &binaryImage );
	//	�E�B���h�E��j������
	cvDestroyWindow( windowNameCapture );
	cvDestroyWindow( windowNameDifference );
	cvDestroyWindow( windowNameBinary );
	cvDestroyWindow( windowNameGray );

	return 0;
}