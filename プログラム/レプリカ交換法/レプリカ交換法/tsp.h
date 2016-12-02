#define _CRT_SECURE_NO_DEPRECATE 1
#include <stdio.h>
#include <stdlib.h> //RAND_MAX
#include <math.h>
#include <GL/glut.h>

#define N 48          //“ss”
#define WSIZE 630
#define W 1
//•ÏŠ·‘O att48 7762   eil101 77   rat575 499   pr1002 16850
//•ÏŠ·Œã 1

void display();
void draw_solution(int rt[N], double position[N][2], double x, double y);
void keyboard(unsigned char key, int x, int y);
void resize(int w, int h);
void random_route(int rt[N], int seed);
