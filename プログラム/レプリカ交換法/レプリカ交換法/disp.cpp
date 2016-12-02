#include "tsp.h"
#include <time.h>

void idle(void)
{
	glutPostRedisplay();			//�ĕ`��C�x���g�𔭐�������
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 's':
		glutIdleFunc(idle);		//idel���Ăё�����
		break;
	case 'w':
		glutPostRedisplay();
		glutIdleFunc(0);
		break;
	case 'q':
		exit(1);
		break;
	}
}

void resize(int w, int h)
{
	double margin = 0.05;

	glViewport(0, 0, w, h);
	glLoadIdentity();
	glOrtho(0.0 - margin, 1.0 + margin, 0.0 - margin, 1.0 + margin, -1.0, 1.0);
	//	glOrtho(0, 1.0, 0, 1.0, -1.0, 1.0); 
	glOrtho(-w*W / (WSIZE*1.0), w*W / (WSIZE*1.0), -h*W / (WSIZE*1.0), h*W / (WSIZE*1.0), -1.0, 1.0);
	//	glOrtho(-w*60/(WSIZE*1.0), w*60/(WSIZE*1.0), -h*80/(WSIZE*1.0), h*80/(WSIZE*1.0), -1.0, 1.0); //att48
	//	glOrtho(-w*4.6/(WSIZE*1.0), w*1.6/(WSIZE*1.0), -h*5.5/(WSIZE*1.0), h*2.3/(WSIZE*1.0), -1.0, 1.0); //eil101	
}


void drawBitmapString(void *font, char *string)
{
	glPushAttrib(GL_CURRENT_BIT);

	/* �r�b�g�}�b�v������̕`�� */
	while (*string)
		glutBitmapCharacter(font, *string++);

	glPopAttrib();
}

void draw_solution(int rt[N], double position[N][2], double x, double y)
{
	int i;
	char string[100];

	glClear(GL_COLOR_BUFFER_BIT);

	double dis = 0;
	for (i = 0; i<N; i++) {
		dis += hypot(position[rt[i]][0] - position[rt[(i + 1) % N]][0], position[rt[i]][1] - position[rt[(i + 1) % N]][1]);
	}

	glColor3d(1.0, 0.0, 0.0);
	//	glRasterPos2d(10000, 12500);
	glRasterPos2d(0.6*W, 1.0*W);
	//	glRasterPos2d(-1, 2.3); //eil101
	sprintf(string, "cost=%f", dis);
	drawBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, string);

	glColor3d(1.0, 0.0, 0.0);
	//	glRasterPos2d(4000, 12500);
	glRasterPos2d(0.3*W, 1.0*W); //att48
								 //	glRasterPos2d(0.3, 2.3); //eil101
	sprintf(string, "time=%f\n", (double)(y - x) / CLOCKS_PER_SEC);
	drawBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, string);


	/*	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINE_LOOP);
	glVertex2d(0, 0);
	glVertex2d(0, 1);
	glVertex2d(1, 1);
	glVertex2d(1, 0);
	glEnd();*/

	glColor3d(1.0, 1.0, 0.0); //yellow
	glBegin(GL_LINE_LOOP);

	for (i = 0; i<N; i++)
	{
		glVertex2dv(position[rt[i]]);
	}

	glEnd();

	glColor3d(1.0, 1.0, 0.0); //blue
	glPointSize(15);
	glBegin(GL_POINTS);
	for (i = 0; i<N; i++)
	{
		glVertex2dv(position[i]);
	}
	glEnd();

	glutSwapBuffers();
}

void random_route(int rt[N], int seed)//�����_���Ɍo�H�����߂Ă���B�����Ă������[�g�̔z��������_���ɓ���ւ��Ă���B
{
	int i, j, v[N], max, cid;

	srand((unsigned int)time(NULL));
	for (i = 0; i<N; i++)
	{
		v[i] = rand();		//��������
	}

	for (i = 0; i<N; i++)
	{
		max = -1;
		for (j = 0; j<N; j++)
		{
			if (v[j]>max)	//v[j]��max���傫���Ƃ����s
			{
				max = v[j];	//max���X�V
				cid = j;		//���̎��̔ԍ���cid�Ɋi�[
			}
		}
		rt[i] = cid;			//cid�̂Ƃ��̔ԍ���i(����)�Ŋi�[
		v[cid] = -1;			//�Q�x�I�΂�Ȃ��悤��v[cid]��-1�ŏ�����
	}
}