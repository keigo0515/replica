#include "tsp.h"
#include <time.h>
#include <math.h>


#define L 1                           //�T���v����
#define M 4                           //���v���J��
#define LOOP 300000                   //�T����
#define TEMP 50000                    //���K���z����̃T���v�����O��
#define PI 3.141592653589793738462
#define EPS 1e-6
#define P0 0.4                        //�d�Ȃ藦
#define gosa 0.001                    //�d�Ȃ藦�̌덷
#define deltax 0.0001                 //���ϒl�ω��ʇ�x


int route[M][N];                      //���i�K�⏇���j
int exc[N];                           //�������p�z��
int exccount[L][M - 1] = { 0 };       //�������񐔃J�E���g�p�z��
double sum[M] = { 0 };                //�R�X�g�p�z��
double sum1 = 0, sum2 = 0;            //���������̃R�X�g�v�Z�p
double sum3 = 0;                      //��ۑ��p
int pos[N][2];                        //�X�̍��W
double stdpos[N][2];                  //�ϊ���̊X�̍��W
double x[4] = { 0 };                  //2-opt�p
double T[M] = { 0 };                  //���x�z��
int count = 1;                        //�T���񐔃J�E���g
int ans[N] = { 0 };                   //�o�͉�p�z��
double rep[L] = { 0 };                //�f�[�^���p�z��
double sum4 = 0;                      //�f�[�^���v�p
double avg = 0;                       //�f�[�^���ϗp
double var = 0;                       //�f�[�^���U�p
int routecopy[M][N] = { 0 };          //���x�������̉�ۑ�(�R�s�[)�p
double sumcopy[M] = { 0 };            //���x�������̐��K���z�v�Z�̃T���v�����O�R�X�g�ۑ��p





double cost(int rt[N], double position[N][2])
{
	double fx = 0;

	for (int i = 0; i < N; i++) {
		fx += hypot(position[rt[i]][0] - position[rt[(i + 1) % N]][0], position[rt[i]][1] - position[rt[(i + 1) % N]][1]);
	}

	return fx;
}

double kouten1(double myu1, double siguma1, double myu2, double siguma2)
{
	return (myu2*pow(siguma1, 2) - myu1*pow(siguma2, 2) + siguma1*siguma2*sqrt(pow(myu1 - myu2, 2) + 2 * (siguma1*siguma1 - siguma2*siguma2)*log(siguma1 / siguma2))) / (siguma1*siguma1 - siguma2*siguma2);
}

double kouten2(double myu1, double siguma1, double myu2, double siguma2)
{
	return (myu2*pow(siguma1, 2) - myu1*pow(siguma2, 2) - siguma1*siguma2*sqrt(pow(myu1 - myu2, 2) + 2 * (siguma1*siguma1 - siguma2*siguma2)*log(siguma1 / siguma2))) / (siguma1*siguma1 - siguma2*siguma2);
}

double fact(long n)
{
	double f = n;
	int i;

	if (n == 0)return 1.0;

	for (i = n - 1; i > 0; i--) {
		f *= i;
	}

	return f;
}

double prob_inte(double x)
{
	double p, pp, pn, q;
	double a, b;
	int sign;
	int i, j, k;

	p = a = x;

	pp = pn = 0.0;

	i = 1; j = 3; k = 1;

	sign = 1;

	do {

		a *= (-x*x);                             /* ���q */

		b = j*fact(k);                          /* ���� */

		sign *= (-1);

		if (sign > 0)

			pp += (a / b);                         /* �����̉��Z */

		if (sign < 0)

			pn += (a / b);                         /* �����̉��Z */

		j += 2;

		k++;

	} while (fabs(a / b) > EPS);                 /* EPS�F�ł��؂�덷 */

	p += (pp + pn);

	q = 2.0*p / sqrt(PI);

	return q;
}

double menseki(double u, double s, double a1, double b1)
{
	double a, b;

	double pa, pb, pab;

	double wk1, wk2;

	a = (a1 - u) / s;                                                 /* �W�����K���z�ϊ� */

	b = (b1 - u) / s;                                                 /* �W�����K���z�ϊ� */



	wk1 = fabs(a) / sqrt(2.0);                                      /* ���K���z���덷�֐��ɕϊ� */

	wk2 = fabs(b) / sqrt(2.0);                                      /* ���K���z���덷�֐��ɕϊ� */

	if (a == 0) {

		pa = 0;

	}

	else {

		pa = prob_inte(wk1);

	}

	if (b == 0) {

		pb = 0;

	}

	else {

		pb = prob_inte(wk2);

	}

	if (a*b >= 0)

		pab = fabs(pb - pa) / 2.0;

	else

		pab = fabs(pb + pa) / 2.0;


	a = a1;
	b = b1;

	//  printf("��� %g �` %g �̊m���� %7.5lf �D\n",a,b,pab);

	return pab;
}

double kasanariritu(double myu1, double siguma1, double myu2, double siguma2)
{
	double Ex1, Ex2, P;

	Ex1 = kouten1(myu1, siguma1, myu2, siguma2);
	Ex2 = kouten2(myu1, siguma1, myu2, siguma2);

	if (siguma1 < siguma2) {
		if (Ex1 > Ex2) {
			P = (0.5 - menseki(myu1, siguma1, myu1, Ex1)) + (0.5 - menseki(myu2, siguma2, Ex1, myu2));
		}

		if (Ex1 < Ex2) {
			P = (0.5 - menseki(myu1, siguma1, myu1, Ex2)) + (0.5 - menseki(myu2, siguma2, Ex2, myu2));
		}
	}

	if (siguma1 > siguma2) {
		if (Ex1 > Ex2) {
			P = (0.5 - menseki(myu1, siguma1, myu1, Ex2)) + (0.5 - menseki(myu2, siguma2, Ex2, myu2));
		}

		if (Ex1 < Ex2) {
			P = (0.5 - menseki(myu1, siguma1, myu1, Ex1)) + (0.5 - menseki(myu2, siguma2, Ex1, myu2));
		}
	}

	return P;
}

double ougonbunkatu(double a, double b, double myu1, double siguma1, double siguma2)
{
	//������
	double gamma = (3 - sqrt(5)) / 2;
	double x1 = a + gamma*(b - a);
	double x2 = a + (1 - gamma)*(b - a);

	//��ԍX�V�A�����_�X�V�A�I������
	while ((b - a) > 0.00000001) {
		double Px1 = kasanariritu(myu1, siguma1, x1, siguma2);
		double Px2 = kasanariritu(myu1, siguma1, x2, siguma2);

		if (pow(Px1 - P0, 2) > pow(Px2 - P0, 2)) {
			a = x1;
			x1 = x2;
			x2 = a + (1 - gamma)*(b - a);
		}

		if (pow(Px1 - P0, 2) < pow(Px2 - P0, 2)) {
			b = x2;
			x2 = x1;
			x1 = a + gamma*(b - a);
		}
	}

//	printf("�ŏI�I�ȋ�� %f �` %f \n\n", a, b);
//	printf("������̕��ϒl=%f\n\n", (a + b) / 2);

	return (a + b) / 2;
}

void display()
{
	/*	FILE *file = _popen("gnuplot.exe", "w");
	fputs("set yrange[0:1]\n", file);
	fputs("set xrange[-10:30]\n",file);
	fputs("plot 0\n", file);
	fflush(file);  */


	clock_t start1, end1;
	start1 = clock();
	for (int sample = 0; sample < L; sample++) {

		clock_t start, end;
		start = clock();

		count = 1;
		srand((unsigned int)time(NULL));


		/////////////////////////���W�̕ϊ�///////////////////////////

		int maxx = pos[0][0];
		int minx = pos[0][0];
		int maxy = pos[0][1];
		int miny = pos[0][1];
		int max;

		for (int i = 1; i < N; i++) {
			if (pos[i][0] > maxx)maxx = pos[i][0];
			if (pos[i][1] > maxy)maxy = pos[i][1];
		}

		for (int i = 1; i < N; i++) {
			if (pos[i][0] < minx)minx = pos[i][0];
			if (pos[i][1] < miny)miny = pos[i][1];
		}

		if (maxy < maxx)max = maxx;
		if (maxx < maxy)max = maxy;

		for (int i = 0; i < N; i++) {
			stdpos[i][0] = (double)(pos[i][0] - minx) / (double)max;
			stdpos[i][1] = (double)(pos[i][1] - miny) / (double)max;
		}

		////////////////////////////////////////////////////



		/////////////////////////�����𐶐�////////////////////////

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {//0����(N-1)�܂ł����Ɋi�[
				route[i][j] = j;
			}
		}

		for (int i = 0; i < M; i++) {
			for (int j = N; j > 1; j--) {//�����_���ɕ��ёւ��ď����𐶐�
				int a = j - 1;
				int b = rand() % j;
				int temp = route[i][a];
				route[i][a] = route[i][b];
				route[i][b] = temp;
			}
		}

		///////////////////////////////////////////////////////////



		//////////////////////////���x�ݒ�////////////////////////

		T[M - 1] = 4;
		T[0] = 0.0001;
		for (int i = M - 1; i > 0; i--) {
			T[i - 1] = T[i] / pow(T[M - 1] / T[0], (double)1 / (M - 1));
		}

		//////////////////////////////////////////////////////////


		for (int i = 0; i < M; i++) {
			printf("���v���J%d�̏������x=%f\n", i, T[i]);
		}

		printf("\n");


		/////////////////////�������R�X�g�v�Z///////////////////////

		for (int i = 0; i < M; i++) {
			sum[i] = cost(route[i], stdpos);
		}

		//////////////////////////////////////////////////////////


		for (int i = 0; i < M; i++) {
			printf("���v���J%d�̏����R�X�g=%f\n", i, sum[i]);
		}

		draw_solution(route[0], stdpos, 0, 0);


		double *sumcopy; //���x�������̃T���v�����O�R�X�g�ۑ��p
		sumcopy = (double *)malloc(sizeof(double) * TEMP);


//		printf("\n�T���v��%d���\n\n", (sample + 1));

		while (count <= LOOP) {
			//			printf("%d���\n", count);
			//////////////////////�ߖT�����A��//////////////////////////

			for (int a = 0; a < M; a++) {

				int t;
				int i = rand() % N;
				int j = rand() % N;

				while (i == j) {
					j = rand() % N;
				}

				if (i > j) { //�K��i<j�ɂ��� 
					int temp3 = i;
					i = j;
					j = temp3;
				}


				if ((i + 1) != j) {
					x[0] = hypot(stdpos[route[a][i]][0] - stdpos[route[a][i + 1]][0], stdpos[route[a][i]][1] - stdpos[route[a][i + 1]][1]); //���݂̉�
					x[1] = hypot(stdpos[route[a][j - 1]][0] - stdpos[route[a][j]][0], stdpos[route[a][j - 1]][1] - stdpos[route[a][j]][1]);

					x[2] = hypot(stdpos[route[a][i]][0] - stdpos[route[a][j - 1]][0], stdpos[route[a][i]][1] - stdpos[route[a][j - 1]][1]); //�ߖT��
					x[3] = hypot(stdpos[route[a][i + 1]][0] - stdpos[route[a][j]][0], stdpos[route[a][i + 1]][1] - stdpos[route[a][j]][1]);
				}

				if ((i + 1) == j) {
					int i0 = i - 1;
					if (i0 < 0)i0 = N - 1;

					x[0] = hypot(stdpos[route[a][i0]][0] - stdpos[route[a][i]][0], stdpos[route[a][i0]][1] - stdpos[route[a][i]][1]);           //���݂̉�
					x[1] = hypot(stdpos[route[a][j]][0] - stdpos[route[a][(j + 1) % N]][0], stdpos[route[a][j]][1] - stdpos[route[a][(j + 1) % N]][1]);

					x[2] = hypot(stdpos[route[a][i0]][0] - stdpos[route[a][j]][0], stdpos[route[a][i0]][1] - stdpos[route[a][j]][1]);           //�ߖT��
					x[3] = hypot(stdpos[route[a][i]][0] - stdpos[route[a][(j + 1) % N]][0], stdpos[route[a][i]][1] - stdpos[route[a][(j + 1) % N]][1]);
				}

				double delta = x[3] + x[2] - (x[1] + x[0]);

				if (delta <= 0) { //���P�� 
					if ((i + 1) != j) {
						for (int c = 0; c < ((j - i - 1) / 2); c++) {
							t = route[a][i + 1 + c];
							route[a][i + 1 + c] = route[a][j - 1 - c];
							route[a][j - 1 - c] = t;
						}
					}


					if ((i + 1) == j) {
						t = route[a][i];
						route[a][i] = route[a][j];
						route[a][j] = t;
					}

					draw_solution(route[0], stdpos, 0, 0);
				}

				if (delta > 0) { //������ 
					double P = exp((-1 * delta) / T[a]); //�󗝊m��
					double p = (double)rand() / RAND_MAX; //�J�ڊm��
					if (p < P) {
						if ((i + 1) != j) {
							for (int c = 0; c < ((j - i - 1) / 2); c++) {
								t = route[a][i + 1 + c];
								route[a][i + 1 + c] = route[a][j - 1 - c];
								route[a][j - 1 - c] = t;
							}
						}


						if ((i + 1) == j) {
							t = route[a][i];
							route[a][i] = route[a][j];
							route[a][j] = t;
						}

						draw_solution(route[0], stdpos, 0, 0);
					}
				}
			}

			////////////////////////////////////////////////////////////



			///////////////////////�T�����ł��ǂ����̋L��////////////////////////////

			double y[M] = { 0 }; //�R�X�g�L���p�z��
			for (int i = 0; i < M; i++) {
				y[i] = cost(route[i], stdpos);
			}

			if (count == 1) {
				sum3 = y[0];
				for (int f = 1; f < M; f++) {
					if (sum3 > y[f]) {
						sum3 = y[f];
						for (int i = 0; i < N; i++) {
							ans[i] = route[f][i];
						}
					}
				}
			}

			if (count != 1) {
				for (int f = 0; f < M; f++) {
					if (sum3 > y[f]) {
						sum3 = y[f];
						for (int i = 0; i < N; i++) {
							ans[i] = route[f][i];
						}
					}
				}
			}

			//////////////////////////////////////////////////////////////////////



			//////////////////////////////���v���J����////////////////////////////////

			if ((count % 5000) == 0) {
				int count2 = 0;
				while (count2 < (M - 1)) {
					if (count2 == (M - 1))break;

					sum1 = cost(route[count2], stdpos);
					sum2 = cost(route[count2 + 1], stdpos);

					if (((T[count2 + 1] - T[count2])*(sum2 - sum1)) <= 0) {
						for (int d = 0; d < N; d++) { //������ 
							exc[d] = route[count2][d];
							route[count2][d] = route[count2 + 1][d];
							route[count2 + 1][d] = exc[d];
						}
						exccount[sample][count2]++;
						count2++;
					}

					if (((T[count2 + 1] - T[count2])*(sum2 - sum1)) > 0) {
						double Q = exp(-(T[count2 + 1] - T[count2])*(sum2 - sum1) / (T[count2] * T[count2 + 1]));
						double q = (double)rand() / RAND_MAX;
						if (q < Q) {
							for (int d = 0; d < N; d++) {//������ 
								exc[d] = route[count2][d];
								route[count2][d] = route[count2 + 1][d];
								route[count2 + 1][d] = exc[d];
							}
							exccount[sample][count2]++;
							count2++;
						}
					}

					sum1 = 0;
					sum2 = 0;
					count2++;

				}
			}

			///////////////////////////////////////////////////////


			//////////////////////���x����/////////////////////////

			if ((count % 40000) == 0) {
//				printf("���x�����J�n\n");

				for (int i = 0; i < M; i++) { //���݂̉��̕ۑ�(�R�s�[)
					for (int j = 0; j < N; j++) {
						routecopy[i][j] = route[i][j];
					}
				}

				int i = 0;
				while (i < M - 1) {
//					printf("���v���J%d��%d\n", i, i + 1);

					double myu[M] = { 0 };           //���ϒl��
					double siguma[M] = { 0 };        //�W���΍���
					double oldmyu = 0;               //�ω��O�̕��ϒl


					for (int j = i; j <= i + 1; j++) {         //���v���J���Ƃ�

						for (int k = 0; k < TEMP; k++) { //�T���v�����O�R�X�g�ۑ��p�̏�����
							sumcopy[k] = 0;
						}

						for (int k = 0; k < TEMP; k++) { //���x��胁�g���|���X�ŃT���v�����O
							int t;
							int l = rand() % N;
							int m = rand() % N;

							while (l == m) {
								m = rand() % N;
							}

							if (l > m) { //�K��i<j��
								int q = l;
								l = m;
								m = q;
							}

							if ((l + 1) != m) {
								x[0] = hypot(stdpos[routecopy[j][l]][0] - stdpos[routecopy[j][l + 1]][0], stdpos[routecopy[j][l]][1] - stdpos[routecopy[j][l + 1]][1]); //���݂̉�
								x[1] = hypot(stdpos[routecopy[j][m - 1]][0] - stdpos[routecopy[j][m]][0], stdpos[routecopy[j][m - 1]][1] - stdpos[routecopy[j][m]][1]);

								x[2] = hypot(stdpos[routecopy[j][l]][0] - stdpos[routecopy[j][m - 1]][0], stdpos[routecopy[j][l]][1] - stdpos[routecopy[j][m - 1]][1]); //�ߖT��
								x[3] = hypot(stdpos[routecopy[j][l + 1]][0] - stdpos[routecopy[j][m]][0], stdpos[routecopy[j][l + 1]][1] - stdpos[routecopy[j][m]][1]);
							}

							double delta = x[3] + x[2] - x[1] - x[0];

							if (delta <= 0) { //���P�� 
								if ((l + 1) != m) {
									for (int c = 0; c < ((m - l - 1) / 2); c++) {
										t = routecopy[j][l + 1 + c];
										routecopy[j][l + 1 + c] = routecopy[j][m - 1 - c];
										routecopy[j][m - 1 - c] = t;
									}
								}


								if ((l + 1) == m) {
									t = routecopy[j][l];
									routecopy[j][l] = routecopy[i][m];
									routecopy[j][m] = t;
								}
							}


							if (delta > 0) { //������ 
								double P = exp((-1 * delta) / T[j]); //�󗝊m��
								double p = (double)rand() / RAND_MAX; //�J�ڊm��
								if (p < P) {
									if ((l + 1) != m) {
										for (int c = 0; c < ((m - l - 1) / 2); c++) {
											t = routecopy[j][l + 1 + c];
											routecopy[j][l + 1 + c] = routecopy[j][m - 1 - c];
											routecopy[j][m - 1 - c] = t;
										}
									}


									if ((l + 1) == m) {
										t = routecopy[j][l];
										routecopy[j][l] = routecopy[j][m];
										routecopy[j][m] = t;
									}

								}
							}

							sumcopy[k] = cost(routecopy[j], stdpos);

						}


						for (int u = 0; u < TEMP; u++) { //���ϒl�ʂ̌v�Z
							myu[j] += sumcopy[u];
						}
						myu[j] /= TEMP;


						for (int u = 0; u < TEMP; u++) { //�W���΍��Ђ̌v�Z
							siguma[j] += (myu[j] - sumcopy[u])*(myu[j] - sumcopy[u]);
						}
						siguma[j] = sqrt(siguma[j] / TEMP);


//						printf("���x%f�̃��v���J%d\n", T[j], j);
//						printf("���ϒl�ʁF%f\n", myu[j]);
//						printf("�W���΍��ЁF%f\n\n", siguma[j]);
						//   					fputs("replot exp(-(x-myu[i])**2/2*siguma[i]**2)/sqrt(2*pi*siguma**2)\n", file);
					}

//                  printf("\n");

					double Ex1[M - 1] = { 0 };
					double Ex2[M - 1] = { 0 };
					double kasanari[M - 1] = { 0 };

					oldmyu = myu[i + 1];
					if (myu[i] > myu[i + 1]) {
						myu[i + 1] += (2 * (myu[i] - myu[i + 1]));
						oldmyu = myu[i + 1];
//						printf("�ω���̍������ϒl=%f\n\n", myu[i + 1]);
					}

					//��x���g���p�^�[��
/*					do {

						kasanari[i] = kasanariritu(myu[i], siguma[i], myu[i + 1], siguma[i + 1]);

						//						printf("\n");

						if (kasanari[i] > (P0 + gosa)) { //�d�Ȃ藦���傫���ꍇ�͍������v���J�̕��ϒl��傫��
							myu[i + 1] += deltax;
						}

						if (kasanari[i] < (P0 - gosa)) { //�d�Ȃ藦���������ꍇ�͍������v���J�̕��ϒl��������
							myu[i + 1] -= deltax;
						}

					} while ((kasanari[i] < (P0 - gosa)) || (kasanari[i] > (P0 + gosa)));*/


					//���������@���g���p�^�[��
					kasanari[i] = kasanariritu(myu[i], siguma[i], myu[i + 1], siguma[i + 1]);

//					printf("�����O�̏d�Ȃ藦=%f\n\n", kasanari[i]);

					if (kasanari[i] > (P0 + gosa))myu[i + 1] = ougonbunkatu(myu[i + 1], myu[i + 1] + 1, myu[i], siguma[i], siguma[i + 1]);

					if (kasanari[i] < (P0 - gosa))myu[i + 1] = ougonbunkatu(myu[i] + (myu[i + 1] - myu[i]) / 100, myu[i + 1], myu[i], siguma[i], siguma[i + 1]);

					kasanari[i] = kasanariritu(myu[i], siguma[i], myu[i + 1], siguma[i + 1]);

/*					printf("������̏d�Ȃ藦=%f\n\n", kasanari[i]);

					printf("�ω��O�̕��ϒl=%f\n", oldmyu);
					printf("�ω���̕��ϒl=%f\n\n", myu[i + 1]);


					printf("�����O�̉��x=%f\n", T[i + 1]);*/

					//�������v���J�̉��x�𒲐�
					T[i + 1] = T[i] + ((T[i + 1] - T[i])*((myu[i + 1] - myu[i]) / (oldmyu - myu[i])));
						//					if (myu[i] < oldmyu)T[i + 1] = T[i] + ((T[i + 1] - T[i])*((myu[i + 1] - myu[i]) / (oldmyu - myu[i])));
						//					if (myu[i] > oldmyu)T[i + 1] = T[i] + ((T[i + 1] - T[i])*((myu[i + 1] - oldmyu) / (myu[i] - oldmyu)));

/*						printf("������̉��x=%f\n", T[i + 1]);
					printf("\n");

					system("pause");*/

					i++;
				}


//				printf("���x�����I��\n\n");
			}



			//////////////////////////////////////////////////////

			count++;
		}

		free(sumcopy);


		for (int i = 0; i < M; i++) {
			sum[i] = 0;
		}
		for (int j = 0; j < M; j++) {//�ŏI�R�X�g�v�Z
			sum[j] = cost(route[j], stdpos);
		}

		printf("\n");

		for (int i = 0; i < M; i++) {
			printf("���v���J%d�̍ŏI���x=%f\n", i, T[i]);
		}

		printf("\n");

		for (int i = 0; i < M; i++) {
			printf("���v���J%d�̍ŏI�R�X�g=%f\n", i, sum[i]);
		}

		rep[sample] = sum3; //�f�[�^���p�z��ɍŗǃR�X�g����

		end = clock();

		printf("\n�ŗǉ��̃R�X�g=%f\n\n", sum3);
		draw_solution(ans, stdpos, start, end);
		//	draw_solution(ans,pos,0,0);
	}

	printf("\n");

	for (int i = 0; i < L; i++) {
		sum4 += rep[i];
	}
	avg = sum4 / L;
	for (int i = 0; i < L; i++) {
		var += (rep[i] - avg)*(rep[i] - avg);
	}

	for (int i = 0; i < L; i++) {
		printf("%f\n", rep[i]);
	}

	printf("�f�[�^���F%d\n", L);
	printf("����=%f\n", avg);
	printf("�W���΍�=%f\n", sqrt(var / L));

	end1 = clock();
	printf("time=%.4f�b\n", (double)(end1 - start1) / CLOCKS_PER_SEC);

	for (int i = 0; i < L; i++) {
		for (int j = 0; j < (M - 1); j++) {
			printf("%d ", exccount[i][j]);
		}
		printf("\n");
	}

	system("pause");
	//	_pclose(file);

}


void main(int argc, char *argv[])
{
	FILE *fp;
	int i, j;



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
	for (j = 0; j < 10; j++) {
		for (i = 0; i < N; i++)
		{
			fscanf(fp, "%d,%d,%d\n", &route[j][i], &(pos[i][0]), &(pos[i][1]));//pos[i][0]��x���W�@pos[i][1]���W��y���W
		}
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
