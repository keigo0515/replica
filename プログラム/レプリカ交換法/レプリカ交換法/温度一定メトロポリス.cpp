#include "tsp.h"
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdio>
using namespace std;

#define SIKOU 1 //���s��
#define M 7                 //�T���v�����Ƃ鉷�x��
double SAMPLING[M] = { 50000,50000,50000,50000,50000,50000,50000 };
//double SAMPLING[M] = { 1000000,1000000,1000000,50000,1000,1000,1000 };
//double SAMPLING[M] = { 3000000,3000000,1000000,10000,1000,1000,1000 };      //�T���v�����O��
#define SAMPLING2 50000     //���ϒl�v�Z��
 
double T[M] = { 0 };         //���x
double ave[M] = { 0 };       //���ϒl�p�z��
double var[M] = { 0 };       //�W���΍��p�z��
int route[N];                //��(�K�⏇��)
int pos[N][2];               //�X�̍��W
double stdpos[N][2];         //�ϊ���̊X�̍��W
double x[4] = { 0 };         //2-opt�p
double epsilon = 0.00001;     //�ŋ}�~���@�̊w�K�W��
double theta = 6;       //�ŋ}�~���@�̌덷����a��臒l

double cost(int rt[N], double position[N][2])
{
	double fx = 0;

	for (int i = 0; i < N; i++) {
		fx += hypot(position[rt[i]][0] - position[rt[(i + 1) % N]][0], position[rt[i]][1] - position[rt[(i + 1) % N]][1]);
	}

	return fx;
}

double fT(double T, double a, double b, double c, double d) {

	double fT = ((b*pow(T, a)) / (pow(T, a) + exp(d))) + c;

	return fT;
}

void display()
{
	for (int k = 0; k < SIKOU; k++) {
		printf("%d���\n", k+1);
		
		srand((unsigned int)time(NULL));

		FILE *fp;

		fp = fopen("att48_sampling.txt", "w");


		/////////////////////////////���W�̕ϊ�////////////////////////////////

		int maxx = pos[0][0];
		int minx = pos[0][0];
		int maxy = pos[0][1];
		int miny = pos[0][1];
		int max;

		for (int i = 1; i < N; i++) {
			if (pos[i][0] > maxx)maxx = pos[i][0];
			if (pos[i][0] > maxy)maxy = pos[i][1];
		}

		for (int i = 0; i < N; i++) {
			if (pos[i][0] < minx)minx = pos[i][0];
			if (pos[i][0] < miny)miny = pos[i][0];
		}

		if (maxx > maxy)max = maxx;
		if (maxx < maxy)max = maxy;

		for (int i = 0; i < N; i++) {
			stdpos[i][0] = (double)(pos[i][0] - minx) / (double)max;
			stdpos[i][1] = (double)(pos[i][1] - miny) / (double)max;
		}

		////////////////////////////////////////////////////////////////////////


		//////////////////////////////���x�ݒ�//////////////////////////////////

		T[M - 1] = 100; //�ō����x
		T[0] = 0.0001;  //�Œቷ�x

		for (int i = M - 1; i > 0; i--) {
			T[i - 1] = T[i] / pow(T[M - 1] / T[0], (double)1 / (M - 1));
		}

		////////////////////////////////////////////////////////////////////////


		clock_t start, end;
		start = clock();

		for (int l = 0; l < M; l++) {

			int count = 1;


			/////////////////////////////�����𐶐�/////////////////////////////////

			for (int i = 0; i < N; i++) {
				route[i] = i;
			}

			for (int i = N; i > 1; i--) {
				int a = i - 1;
				int b = rand() % i;
				int temp = route[a];
				route[a] = route[b];
				route[b] = temp;
			}

			////////////////////////////////////////////////////////////////////////


			////////////////////////////�����R�X�g�v�Z//////////////////////////////

			double *sum = (double *)malloc(sizeof(double)*SAMPLING2);

			for (int i = 0; i < SAMPLING2; i++) {
				sum[i] = 0;
			}

			///////////////////////////////////////////////////////////////////////


			draw_solution(route, stdpos, 0, 0);

			while (count <= (SAMPLING[l] + SAMPLING2)) {
				if ((count % 500000) == 0)printf("%d���\n", count);

				int t;
				int i = rand() % N;
				int j = rand() % N;

				if (i > j) { //�K��i<j��
					int temp3 = i;
					i = j;
					j = temp3;
				}

				if ((i + 1) != j) {
					x[0] = hypot(stdpos[route[i]][0] - stdpos[route[i + 1]][0], stdpos[route[i]][1] - stdpos[route[i + 1]][1]); //���݂̉�
					x[1] = hypot(stdpos[route[j - 1]][0] - stdpos[route[j]][0], stdpos[route[j - 1]][1] - stdpos[route[j]][1]);

					x[2] = hypot(stdpos[route[i]][0] - stdpos[route[j - 1]][0], stdpos[route[i]][1] - stdpos[route[j - 1]][1]); //�ߖT��
					x[3] = hypot(stdpos[route[i + 1]][0] - stdpos[route[j]][0], stdpos[route[i + 1]][1] - stdpos[route[j]][1]);
				}

				if ((i + 1) == j) {
					int i0 = i - 1;
					if (i0 < 0)i0 = N - 1;

					x[0] = hypot(stdpos[route[i0]][0] - stdpos[route[i]][0], stdpos[route[i0]][1] - stdpos[route[i]][1]);           //���݂̉�
					x[1] = hypot(stdpos[route[j]][0] - stdpos[route[(j + 1) % N]][0], stdpos[route[j]][1] - stdpos[route[(j + 1) % N]][1]);

					x[2] = hypot(stdpos[route[i0]][0] - stdpos[route[j]][0], stdpos[route[i0]][1] - stdpos[route[j]][1]);           //�ߖT��
					x[3] = hypot(stdpos[route[i]][0] - stdpos[route[(j + 1) % N]][0], stdpos[route[i]][1] - stdpos[route[(j + 1) % N]][1]);
				}


				double delta = x[3] + x[2] - (x[1] + x[0]);


				if (delta <= 0) { //���P�� 
					if ((i + 1) != j) {
						for (int c = 0; c < ((j - i - 1) / 2); c++) {
							t = route[i + 1 + c];
							route[i + 1 + c] = route[j - 1 - c];
							route[j - 1 - c] = t;
						}
					}


					if ((i + 1) == j) {
						t = route[i];
						route[i] = route[j];
						route[j] = t;
					}

					draw_solution(route, stdpos, 0, 0);
				}


				if (delta > 0) { //������ 
					double P = exp((-1 * delta) / T[l]); //�󗝊m��
					double p = (double)rand() / RAND_MAX; //�J�ڊm��
					if (p < P) {
						if ((i + 1) != j) {
							for (int c = 0; c < ((j - i - 1) / 2); c++) {
								t = route[i + 1 + c];
								route[i + 1 + c] = route[j - 1 - c];
								route[j - 1 - c] = t;
							}
						}


						if ((i + 1) == j) {
							t = route[i];
							route[i] = route[j];
							route[j] = t;
						}

						draw_solution(route, stdpos, 0, 0);
					}
				}


				if (count > SAMPLING[l]) {
					sum[count - (int)SAMPLING[l] - 1] = cost(route, stdpos);
					ave[l] += sum[count - (int)SAMPLING[l] - 1];
				}

				//			sum[count] = cost(route, stdpos);
				//			ave += sum[count];

				count++;
			}

			ave[l] /= SAMPLING2;  //���ϒl�v�Z


			for (int i = 0; i < SAMPLING2; i++) {  //�W���΍��v�Z
				var[l] += ((ave[l] - sum[i])*(ave[l] - sum[i]));
			}
			var[l] = sqrt(var[l] / SAMPLING2);

			printf("���xT=%f\n", T[l]);
			printf("���ϒl=%f\n", ave[l]);
			printf("�W���΍�=%f\n\n", var[l]);

			free(sum);

			draw_solution(route, stdpos, 0, 0);

			fprintf(fp, "%f %f\n", log(T[l]), ave[l]);

		}
		fclose(fp);

		////////////////////////////////�ŋ}�~���@/////////////////////////////////////

		double a = 1;
		double d = 0;

		double c = ave[0];
		for (int i = 1; i < M; i++) {
			if (c > ave[i])c = ave[i];
		}

		double b = ave[0];
		for (int i = 1; i < M; i++) {
			if (b < ave[i])b = ave[i];
		}
		b = b - c;

		printf("�ŋ}�~���@�p�����[�^�����ݒ�\n");
		printf("�w�K�W����=%f\n", epsilon);
		printf("�I�������F�덷���aE<=%f\n", theta);
		printf("a=%f\nb=%f\nc=%f\nd=%f\n\n",a,b,c,d);

		double E = 0;
		double da = 0;
		double db = 0;
		double dc = 0;
		double dd = 0;


		for (int i = 0; i < M; i++) {
			E += (fT(T[i], a, b, c, d) - ave[i])*(fT(T[i], a, b, c, d) - ave[i]);
		}

		printf("�덷���aE=%f\n\n", E);

		while (E > theta) {
			for (int i = 0; i < M; i++) {
				da += ((b*pow(T[i], a)*exp(d)*log(T[i])) / pow(pow(T[i], a) + exp(d), 2))*(fT(T[i], a, b, c, d) - ave[i]);

				db += (pow(T[i], a) / (pow(T[i], a) + exp(d)))*(fT(T[i], a, b, c, d) - ave[i]);

				dc += (fT(T[i], a, b, c, d) - ave[i]);

				dd += -(((b*pow(T[i], a)*exp(d)) / pow(pow(T[i], a) + exp(d),2))*(fT(T[i], a, b, c, d) - ave[i]));
			}

			a = a - epsilon*da;
			b = b - epsilon*db;
			c = c - epsilon*dc;
			d = d - epsilon*dd;

//			printf("a=%f\nd=%f\n\n", a, d);

			E = 0;
			for (int i = 0; i < M; i++) {
				E += (fT(T[i], a, b, c, d) - ave[i])*(fT(T[i], a, b, c, d) - ave[i]);
			}
//			printf("�덷���aE=%f\n", E);
		}

		printf("�ŏI�덷���aE=%f\n\n", E);
		printf("�ŏI�p�����[�^\n");
		printf("a=%f\nb=%f\nc=%f\nd=%f\n",a,b,c,d);

		FILE *fp2 = _popen("gnuplot.exe", "w");
		fprintf(fp2, "f(a,b,c,d,x)=(b/(1+exp(-a*x+d)))+c\n");
		fprintf(fp2, "plot f(%d,%d,%d,%d,x)\n", a, b, c, d);
		fflush(fp);
		cin.get();
		_pclose(fp2);


		///////////////////////////////////////////////////////////////////////////////
			
		end = clock();
		draw_solution(route, stdpos, start, end);
		printf("\n%f�b\n", (double)(end - start) / CLOCKS_PER_SEC);
	}

	system("pause");

}


void main(int argc, char *argv[])
{
	FILE *fp;
	int i;



	//�O���t�B�b�N�p�֐��D�폜����ȁI
	glutInitWindowPosition(500, 100);
	glutInitWindowSize(700, WSIZE);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutCreateWindow("TSP");
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keyboard);//s:start; w:wait; q:quit
	glClearColor(0.0, 0.0, 0.0, 1.0);


	//�t�@�C�����J��
	if ((fp = fopen("att48.txt", "r")) == NULL)
	{
		printf("no file\n");
		exit(0);
	}


	//�z�B����W���Q�����z��pos�ɓǂ݂��ށD
	for (i = 0; i < N; i++)
	{
		fscanf(fp, "%d,%d,%d\n", &route[i], &(pos[i][0]), &(pos[i][1]));//pos[i][0]��x���W�@pos[i][1]���W��y���W
	}
	fclose(fp);

	//�z�B�������P�����z��route1�ɓ����

	/*for(i=0; i<N; i++)	//�X�̔ԍ������킩���Ă���̂�route1��(�Ƃ肠����)����Ă���B�Ȃ̂ł܂��e���W�̏��͂Ȃ�(�͂�)�B
	{						//&route1�ɕύX��͊X�ƍ��W���z��̔ԍ��őΉ����Ă���͂�
	route1[i]=i;
	}
	*/
	//�����ɂ����������_�������\�b�h��diplay���\�b�h�̒���

	glutMainLoop();//�O���t�B�b�N�p�֐��D�폜����ȁI
}