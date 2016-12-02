#include "tsp.h"
#include <time.h>
#include <math.h>


#define L 1                           //サンプル回数
#define M 4                           //レプリカ数
#define LOOP 300000                   //探索回数
#define TEMP 50000                    //正規分布推定のサンプリング回数
#define PI 3.141592653589793738462
#define EPS 1e-6
#define P0 0.4                        //重なり率
#define gosa 0.001                    //重なり率の誤差
#define deltax 0.0001                 //平均値変化量⊿x


int route[M][N];                      //解（訪問順序）
int exc[N];                           //解交換用配列
int exccount[L][M - 1] = { 0 };       //解交換回数カウント用配列
double sum[M] = { 0 };                //コスト用配列
double sum1 = 0, sum2 = 0;            //解交換時のコスト計算用
double sum3 = 0;                      //解保存用
int pos[N][2];                        //街の座標
double stdpos[N][2];                  //変換後の街の座標
double x[4] = { 0 };                  //2-opt用
double T[M] = { 0 };                  //温度配列
int count = 1;                        //探索回数カウント
int ans[N] = { 0 };                   //出力解用配列
double rep[L] = { 0 };                //データ取る用配列
double sum4 = 0;                      //データ合計用
double avg = 0;                       //データ平均用
double var = 0;                       //データ分散用
int routecopy[M][N] = { 0 };          //温度調整時の解保存(コピー)用
double sumcopy[M] = { 0 };            //温度調整時の正規分布計算のサンプリングコスト保存用





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

		a *= (-x*x);                             /* 分子 */

		b = j*fact(k);                          /* 分母 */

		sign *= (-1);

		if (sign > 0)

			pp += (a / b);                         /* 正項の加算 */

		if (sign < 0)

			pn += (a / b);                         /* 負項の加算 */

		j += 2;

		k++;

	} while (fabs(a / b) > EPS);                 /* EPS：打ち切り誤差 */

	p += (pp + pn);

	q = 2.0*p / sqrt(PI);

	return q;
}

double menseki(double u, double s, double a1, double b1)
{
	double a, b;

	double pa, pb, pab;

	double wk1, wk2;

	a = (a1 - u) / s;                                                 /* 標準正規分布変換 */

	b = (b1 - u) / s;                                                 /* 標準正規分布変換 */



	wk1 = fabs(a) / sqrt(2.0);                                      /* 正規分布を誤差関数に変換 */

	wk2 = fabs(b) / sqrt(2.0);                                      /* 正規分布を誤差関数に変換 */

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

	//  printf("区間 %g ～ %g の確率は %7.5lf ．\n",a,b,pab);

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
	//初期化
	double gamma = (3 - sqrt(5)) / 2;
	double x1 = a + gamma*(b - a);
	double x2 = a + (1 - gamma)*(b - a);

	//区間更新、分割点更新、終了判定
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

//	printf("最終的な区間 %f ～ %f \n\n", a, b);
//	printf("調整後の平均値=%f\n\n", (a + b) / 2);

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


		/////////////////////////座標の変換///////////////////////////

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



		/////////////////////////初期解生成////////////////////////

		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {//0から(N-1)までを順に格納
				route[i][j] = j;
			}
		}

		for (int i = 0; i < M; i++) {
			for (int j = N; j > 1; j--) {//ランダムに並び替えて初期解生成
				int a = j - 1;
				int b = rand() % j;
				int temp = route[i][a];
				route[i][a] = route[i][b];
				route[i][b] = temp;
			}
		}

		///////////////////////////////////////////////////////////



		//////////////////////////温度設定////////////////////////

		T[M - 1] = 4;
		T[0] = 0.0001;
		for (int i = M - 1; i > 0; i--) {
			T[i - 1] = T[i] / pow(T[M - 1] / T[0], (double)1 / (M - 1));
		}

		//////////////////////////////////////////////////////////


		for (int i = 0; i < M; i++) {
			printf("レプリカ%dの初期温度=%f\n", i, T[i]);
		}

		printf("\n");


		/////////////////////初期解コスト計算///////////////////////

		for (int i = 0; i < M; i++) {
			sum[i] = cost(route[i], stdpos);
		}

		//////////////////////////////////////////////////////////


		for (int i = 0; i < M; i++) {
			printf("レプリカ%dの初期コスト=%f\n", i, sum[i]);
		}

		draw_solution(route[0], stdpos, 0, 0);


		double *sumcopy; //温度調整時のサンプリングコスト保存用
		sumcopy = (double *)malloc(sizeof(double) * TEMP);


//		printf("\nサンプル%d回目\n\n", (sample + 1));

		while (count <= LOOP) {
			//			printf("%d回目\n", count);
			//////////////////////近傍生成、受理//////////////////////////

			for (int a = 0; a < M; a++) {

				int t;
				int i = rand() % N;
				int j = rand() % N;

				while (i == j) {
					j = rand() % N;
				}

				if (i > j) { //必ずi<jにする 
					int temp3 = i;
					i = j;
					j = temp3;
				}


				if ((i + 1) != j) {
					x[0] = hypot(stdpos[route[a][i]][0] - stdpos[route[a][i + 1]][0], stdpos[route[a][i]][1] - stdpos[route[a][i + 1]][1]); //現在の解
					x[1] = hypot(stdpos[route[a][j - 1]][0] - stdpos[route[a][j]][0], stdpos[route[a][j - 1]][1] - stdpos[route[a][j]][1]);

					x[2] = hypot(stdpos[route[a][i]][0] - stdpos[route[a][j - 1]][0], stdpos[route[a][i]][1] - stdpos[route[a][j - 1]][1]); //近傍解
					x[3] = hypot(stdpos[route[a][i + 1]][0] - stdpos[route[a][j]][0], stdpos[route[a][i + 1]][1] - stdpos[route[a][j]][1]);
				}

				if ((i + 1) == j) {
					int i0 = i - 1;
					if (i0 < 0)i0 = N - 1;

					x[0] = hypot(stdpos[route[a][i0]][0] - stdpos[route[a][i]][0], stdpos[route[a][i0]][1] - stdpos[route[a][i]][1]);           //現在の解
					x[1] = hypot(stdpos[route[a][j]][0] - stdpos[route[a][(j + 1) % N]][0], stdpos[route[a][j]][1] - stdpos[route[a][(j + 1) % N]][1]);

					x[2] = hypot(stdpos[route[a][i0]][0] - stdpos[route[a][j]][0], stdpos[route[a][i0]][1] - stdpos[route[a][j]][1]);           //近傍解
					x[3] = hypot(stdpos[route[a][i]][0] - stdpos[route[a][(j + 1) % N]][0], stdpos[route[a][i]][1] - stdpos[route[a][(j + 1) % N]][1]);
				}

				double delta = x[3] + x[2] - (x[1] + x[0]);

				if (delta <= 0) { //改善解 
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

				if (delta > 0) { //改悪解 
					double P = exp((-1 * delta) / T[a]); //受理確率
					double p = (double)rand() / RAND_MAX; //遷移確率
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



			///////////////////////探索中最も良い解の記憶////////////////////////////

			double y[M] = { 0 }; //コスト記憶用配列
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



			//////////////////////////////レプリカ交換////////////////////////////////

			if ((count % 5000) == 0) {
				int count2 = 0;
				while (count2 < (M - 1)) {
					if (count2 == (M - 1))break;

					sum1 = cost(route[count2], stdpos);
					sum2 = cost(route[count2 + 1], stdpos);

					if (((T[count2 + 1] - T[count2])*(sum2 - sum1)) <= 0) {
						for (int d = 0; d < N; d++) { //解交換 
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
							for (int d = 0; d < N; d++) {//解交換 
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


			//////////////////////温度調整/////////////////////////

			if ((count % 40000) == 0) {
//				printf("温度調整開始\n");

				for (int i = 0; i < M; i++) { //現在の解の保存(コピー)
					for (int j = 0; j < N; j++) {
						routecopy[i][j] = route[i][j];
					}
				}

				int i = 0;
				while (i < M - 1) {
//					printf("レプリカ%dと%d\n", i, i + 1);

					double myu[M] = { 0 };           //平均値μ
					double siguma[M] = { 0 };        //標準偏差σ
					double oldmyu = 0;               //変化前の平均値


					for (int j = i; j <= i + 1; j++) {         //レプリカごとに

						for (int k = 0; k < TEMP; k++) { //サンプリングコスト保存用の初期化
							sumcopy[k] = 0;
						}

						for (int k = 0; k < TEMP; k++) { //温度一定メトロポリスでサンプリング
							int t;
							int l = rand() % N;
							int m = rand() % N;

							while (l == m) {
								m = rand() % N;
							}

							if (l > m) { //必ずi<jに
								int q = l;
								l = m;
								m = q;
							}

							if ((l + 1) != m) {
								x[0] = hypot(stdpos[routecopy[j][l]][0] - stdpos[routecopy[j][l + 1]][0], stdpos[routecopy[j][l]][1] - stdpos[routecopy[j][l + 1]][1]); //現在の解
								x[1] = hypot(stdpos[routecopy[j][m - 1]][0] - stdpos[routecopy[j][m]][0], stdpos[routecopy[j][m - 1]][1] - stdpos[routecopy[j][m]][1]);

								x[2] = hypot(stdpos[routecopy[j][l]][0] - stdpos[routecopy[j][m - 1]][0], stdpos[routecopy[j][l]][1] - stdpos[routecopy[j][m - 1]][1]); //近傍解
								x[3] = hypot(stdpos[routecopy[j][l + 1]][0] - stdpos[routecopy[j][m]][0], stdpos[routecopy[j][l + 1]][1] - stdpos[routecopy[j][m]][1]);
							}

							double delta = x[3] + x[2] - x[1] - x[0];

							if (delta <= 0) { //改善解 
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


							if (delta > 0) { //改悪解 
								double P = exp((-1 * delta) / T[j]); //受理確率
								double p = (double)rand() / RAND_MAX; //遷移確率
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


						for (int u = 0; u < TEMP; u++) { //平均値μの計算
							myu[j] += sumcopy[u];
						}
						myu[j] /= TEMP;


						for (int u = 0; u < TEMP; u++) { //標準偏差σの計算
							siguma[j] += (myu[j] - sumcopy[u])*(myu[j] - sumcopy[u]);
						}
						siguma[j] = sqrt(siguma[j] / TEMP);


//						printf("温度%fのレプリカ%d\n", T[j], j);
//						printf("平均値μ：%f\n", myu[j]);
//						printf("標準偏差σ：%f\n\n", siguma[j]);
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
//						printf("変化後の高温平均値=%f\n\n", myu[i + 1]);
					}

					//⊿xを使うパターン
/*					do {

						kasanari[i] = kasanariritu(myu[i], siguma[i], myu[i + 1], siguma[i + 1]);

						//						printf("\n");

						if (kasanari[i] > (P0 + gosa)) { //重なり率が大きい場合は高温レプリカの平均値を大きく
							myu[i + 1] += deltax;
						}

						if (kasanari[i] < (P0 - gosa)) { //重なり率が小さい場合は高温レプリカの平均値を小さく
							myu[i + 1] -= deltax;
						}

					} while ((kasanari[i] < (P0 - gosa)) || (kasanari[i] > (P0 + gosa)));*/


					//黄金分割法を使うパターン
					kasanari[i] = kasanariritu(myu[i], siguma[i], myu[i + 1], siguma[i + 1]);

//					printf("調整前の重なり率=%f\n\n", kasanari[i]);

					if (kasanari[i] > (P0 + gosa))myu[i + 1] = ougonbunkatu(myu[i + 1], myu[i + 1] + 1, myu[i], siguma[i], siguma[i + 1]);

					if (kasanari[i] < (P0 - gosa))myu[i + 1] = ougonbunkatu(myu[i] + (myu[i + 1] - myu[i]) / 100, myu[i + 1], myu[i], siguma[i], siguma[i + 1]);

					kasanari[i] = kasanariritu(myu[i], siguma[i], myu[i + 1], siguma[i + 1]);

/*					printf("調整後の重なり率=%f\n\n", kasanari[i]);

					printf("変化前の平均値=%f\n", oldmyu);
					printf("変化後の平均値=%f\n\n", myu[i + 1]);


					printf("調整前の温度=%f\n", T[i + 1]);*/

					//高温レプリカの温度を調整
					T[i + 1] = T[i] + ((T[i + 1] - T[i])*((myu[i + 1] - myu[i]) / (oldmyu - myu[i])));
						//					if (myu[i] < oldmyu)T[i + 1] = T[i] + ((T[i + 1] - T[i])*((myu[i + 1] - myu[i]) / (oldmyu - myu[i])));
						//					if (myu[i] > oldmyu)T[i + 1] = T[i] + ((T[i + 1] - T[i])*((myu[i + 1] - oldmyu) / (myu[i] - oldmyu)));

/*						printf("調整後の温度=%f\n", T[i + 1]);
					printf("\n");

					system("pause");*/

					i++;
				}


//				printf("温度調整終了\n\n");
			}



			//////////////////////////////////////////////////////

			count++;
		}

		free(sumcopy);


		for (int i = 0; i < M; i++) {
			sum[i] = 0;
		}
		for (int j = 0; j < M; j++) {//最終コスト計算
			sum[j] = cost(route[j], stdpos);
		}

		printf("\n");

		for (int i = 0; i < M; i++) {
			printf("レプリカ%dの最終温度=%f\n", i, T[i]);
		}

		printf("\n");

		for (int i = 0; i < M; i++) {
			printf("レプリカ%dの最終コスト=%f\n", i, sum[i]);
		}

		rep[sample] = sum3; //データ取る用配列に最良コストを代入

		end = clock();

		printf("\n最良解のコスト=%f\n\n", sum3);
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

	printf("データ数：%d\n", L);
	printf("平均=%f\n", avg);
	printf("標準偏差=%f\n", sqrt(var / L));

	end1 = clock();
	printf("time=%.4f秒\n", (double)(end1 - start1) / CLOCKS_PER_SEC);

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
	for (j = 0; j < 10; j++) {
		for (i = 0; i < N; i++)
		{
			fscanf(fp, "%d,%d,%d\n", &route[j][i], &(pos[i][0]), &(pos[i][1]));//pos[i][0]にx座標　pos[i][1]座標にy座標
		}
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
