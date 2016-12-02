#include "tsp.h"
#include <time.h>
#include <math.h>
#include <iostream>
#include <cstdio>
using namespace std;

#define SIKOU 1 //試行回数
#define M 7                 //サンプルをとる温度数
double SAMPLING[M] = { 50000,50000,50000,50000,50000,50000,50000 };
//double SAMPLING[M] = { 1000000,1000000,1000000,50000,1000,1000,1000 };
//double SAMPLING[M] = { 3000000,3000000,1000000,10000,1000,1000,1000 };      //サンプリング回数
#define SAMPLING2 50000     //平均値計算回数
 
double T[M] = { 0 };         //温度
double ave[M] = { 0 };       //平均値用配列
double var[M] = { 0 };       //標準偏差用配列
int route[N];                //解(訪問順序)
int pos[N][2];               //街の座標
double stdpos[N][2];         //変換後の街の座標
double x[4] = { 0 };         //2-opt用
double epsilon = 0.00001;     //最急降下法の学習係数
double theta = 6;       //最急降下法の誤差自乗和の閾値

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
		printf("%d回目\n", k+1);
		
		srand((unsigned int)time(NULL));

		FILE *fp;

		fp = fopen("att48_sampling.txt", "w");


		/////////////////////////////座標の変換////////////////////////////////

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


		//////////////////////////////温度設定//////////////////////////////////

		T[M - 1] = 100; //最高温度
		T[0] = 0.0001;  //最低温度

		for (int i = M - 1; i > 0; i--) {
			T[i - 1] = T[i] / pow(T[M - 1] / T[0], (double)1 / (M - 1));
		}

		////////////////////////////////////////////////////////////////////////


		clock_t start, end;
		start = clock();

		for (int l = 0; l < M; l++) {

			int count = 1;


			/////////////////////////////初期解生成/////////////////////////////////

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


			////////////////////////////初期コスト計算//////////////////////////////

			double *sum = (double *)malloc(sizeof(double)*SAMPLING2);

			for (int i = 0; i < SAMPLING2; i++) {
				sum[i] = 0;
			}

			///////////////////////////////////////////////////////////////////////


			draw_solution(route, stdpos, 0, 0);

			while (count <= (SAMPLING[l] + SAMPLING2)) {
				if ((count % 500000) == 0)printf("%d回目\n", count);

				int t;
				int i = rand() % N;
				int j = rand() % N;

				if (i > j) { //必ずi<jに
					int temp3 = i;
					i = j;
					j = temp3;
				}

				if ((i + 1) != j) {
					x[0] = hypot(stdpos[route[i]][0] - stdpos[route[i + 1]][0], stdpos[route[i]][1] - stdpos[route[i + 1]][1]); //現在の解
					x[1] = hypot(stdpos[route[j - 1]][0] - stdpos[route[j]][0], stdpos[route[j - 1]][1] - stdpos[route[j]][1]);

					x[2] = hypot(stdpos[route[i]][0] - stdpos[route[j - 1]][0], stdpos[route[i]][1] - stdpos[route[j - 1]][1]); //近傍解
					x[3] = hypot(stdpos[route[i + 1]][0] - stdpos[route[j]][0], stdpos[route[i + 1]][1] - stdpos[route[j]][1]);
				}

				if ((i + 1) == j) {
					int i0 = i - 1;
					if (i0 < 0)i0 = N - 1;

					x[0] = hypot(stdpos[route[i0]][0] - stdpos[route[i]][0], stdpos[route[i0]][1] - stdpos[route[i]][1]);           //現在の解
					x[1] = hypot(stdpos[route[j]][0] - stdpos[route[(j + 1) % N]][0], stdpos[route[j]][1] - stdpos[route[(j + 1) % N]][1]);

					x[2] = hypot(stdpos[route[i0]][0] - stdpos[route[j]][0], stdpos[route[i0]][1] - stdpos[route[j]][1]);           //近傍解
					x[3] = hypot(stdpos[route[i]][0] - stdpos[route[(j + 1) % N]][0], stdpos[route[i]][1] - stdpos[route[(j + 1) % N]][1]);
				}


				double delta = x[3] + x[2] - (x[1] + x[0]);


				if (delta <= 0) { //改善解 
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


				if (delta > 0) { //改悪解 
					double P = exp((-1 * delta) / T[l]); //受理確率
					double p = (double)rand() / RAND_MAX; //遷移確率
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

			ave[l] /= SAMPLING2;  //平均値計算


			for (int i = 0; i < SAMPLING2; i++) {  //標準偏差計算
				var[l] += ((ave[l] - sum[i])*(ave[l] - sum[i]));
			}
			var[l] = sqrt(var[l] / SAMPLING2);

			printf("温度T=%f\n", T[l]);
			printf("平均値=%f\n", ave[l]);
			printf("標準偏差=%f\n\n", var[l]);

			free(sum);

			draw_solution(route, stdpos, 0, 0);

			fprintf(fp, "%f %f\n", log(T[l]), ave[l]);

		}
		fclose(fp);

		////////////////////////////////最急降下法/////////////////////////////////////

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

		printf("最急降下法パラメータ初期設定\n");
		printf("学習係数ε=%f\n", epsilon);
		printf("終了条件：誤差二乗和E<=%f\n", theta);
		printf("a=%f\nb=%f\nc=%f\nd=%f\n\n",a,b,c,d);

		double E = 0;
		double da = 0;
		double db = 0;
		double dc = 0;
		double dd = 0;


		for (int i = 0; i < M; i++) {
			E += (fT(T[i], a, b, c, d) - ave[i])*(fT(T[i], a, b, c, d) - ave[i]);
		}

		printf("誤差二乗和E=%f\n\n", E);

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
//			printf("誤差二乗和E=%f\n", E);
		}

		printf("最終誤差二乗和E=%f\n\n", E);
		printf("最終パラメータ\n");
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
		printf("\n%f秒\n", (double)(end - start) / CLOCKS_PER_SEC);
	}

	system("pause");

}


void main(int argc, char *argv[])
{
	FILE *fp;
	int i;



	//グラフィック用関数．削除するな！
	glutInitWindowPosition(500, 100);
	glutInitWindowSize(700, WSIZE);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutCreateWindow("TSP");
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutKeyboardFunc(keyboard);//s:start; w:wait; q:quit
	glClearColor(0.0, 0.0, 0.0, 1.0);


	//ファイルを開く
	if ((fp = fopen("att48.txt", "r")) == NULL)
	{
		printf("no file\n");
		exit(0);
	}


	//配達先座標を２次元配列posに読みこむ．
	for (i = 0; i < N; i++)
	{
		fscanf(fp, "%d,%d,%d\n", &route[i], &(pos[i][0]), &(pos[i][1]));//pos[i][0]にx座標　pos[i][1]座標にy座標
	}
	fclose(fp);

	//配達順序を１次元配列route1に入れる

	/*for(i=0; i<N; i++)	//街の番号だけわかっているのでroute1に(とりあえず)入れてある。なのでまだ各座標の情報はない(はず)。
	{						//&route1に変更後は街と座標が配列の番号で対応しているはず
	route1[i]=i;
	}
	*/
	//ここにあったランダム化メソッドをdiplayメソッドの中へ

	glutMainLoop();//グラフィック用関数．削除するな！
}