// ISODATA.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "math.h"
#include "malloc.h"

//存储样本点
struct sample{
	double x;
	double y;
};

//参数声明
/*
 * C		初始聚类中心
 * θN		一个类别至少应具有的样本数目
 * θS		一个类别样本标准差阈值
 * θC		聚类中心之间距离的阈值，即归并系数
 * L		在一次迭代中可以归并的类别的最多对数
 * I		允许迭代的最多次数
 */
int C, θN, θS, θC, L, I;

int couter;		//迭代次数

#define K 2		//所要求聚类中心数
#define N 8		//样本数
#define n 2		//样本维数
#define M 0.5	//不知道叫啥的参数

//求最小值的下标
int GetMin(double *temp) {
	double Min = 1000;
	int k=1000;
	for (int i = 0; i < C; i++) {
		if (temp[i] <= Min) {
			Min = temp[i];
			k = i;
		}		
	}
	return k;
}

//求样本点与聚类中心的距离
double GetDis(struct sample *x, sample *z, int i, int j) {
	double distance = 0;
	distance = sqrt(pow((x[i].x - z[j].x), 2)
		+ pow((x[i].y - z[j].y), 2));
	return distance;
}

//将N个样本按最小距离进行分类
void Classify(sample *x, sample *z, double *d, int *f, int *_N) {
	double *temp = (double*)malloc(C * sizeof(double));

	for (int j = 0; j < N; j++) {
		for (int i = 0; i < C; i++) {
			temp[i] = GetDis(x, z, j, i);
		}
		d[j] = temp[GetMin(temp)];
		f[GetMin(temp)*N + j] = j;	//分类结果
		_N[GetMin(temp)]++;
		//printf("GetMin(temp)=%d,j=%d\n", GetMin(temp), j);
	}

	for (int i = 0; i < C; i++) {
		for (int j = 0; j < N; j++) {
			//if (f[i*N + j] < 1000)
				//printf("%d\t", f[i*N + j]);
		}
		printf("\n");
	}
	//printf("\n_N 数目：\t");
	for (int i = 0; i < C; i++) {
		//printf("%d\t", _N[i]);
	}

	free(temp);
}

//修正聚类中心
void Correction(sample *x, sample *z, int *f, int *_N) {
	double sum_x = 0, sum_y = 0;

	for (int i = 0; i < C; i++) {
		for (int j = 0; j < N; j++) {
			sum_x += x[f[i*N + j]].x;
			sum_y += x[f[i*N + j]].y;
		}
		z[i].x = sum_x / _N[i];
		z[i].y = sum_y / _N[i];
		sum_x = 0;	sum_y = 0;
		
		//printf("%.2lf\t%.2lf\n", sum_x, sum_y);
		//printf("%lf\t%lf\n", z[i].x, z[i].y);
	}
}

//对每一聚类域f，计算器所有样本到其聚类中心距离的平均值
//计算所有样本到其相应聚类中心的距离的平均值
void GetMean(sample *x, sample *z, int *f, int *_N, double *mean_d, double m_d) {
	double sum = 0;
	
	for (int i = 0; i < C; i++) {
		for (int j = 0; j < N; j++) {
			if (f[i*N + j] < 1000)
				sum += GetDis(x, z, f[i*N + j], i);
		}
		mean_d[i] = sum / _N[i];
		m_d += mean_d[i] * _N[i];
		sum = 0;

		printf("%.2lf\n", mean_d[i]);
	}
	m_d = m_d / N;
	//printf("m_d = %.2lf\n", m_d);
}

//计算每一类别中样本与聚类中心距离的标准差向量
void GetVector(sample *x, sample *z, int *_N, double *sigma) {
	sample *temp_s = (sample*)malloc(C * sizeof(sample));
	double sum_x = 0, sum_y = 0;

	for (int j = 0; j < C; j++) {
		for (int i = 0; i < n; i++) {
			for (int _I = 0; _I < N; _I++) {
				if (i == 0) {
					sum_x += pow((x[_I].x - z[j].x), 2);
				}
				else {
					sum_y += pow((x[_I].y - z[j].y), 2);
				}
			}
		}
		sum_x = sum_x / _N[j];
		sum_y = sum_y / _N[j];
		temp_s[j].x = sqrt(sum_x);
		temp_s[j].y = sqrt(sum_y);
	}

	//赋值于sigma
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < C; j++) {
			if (i == 0)
				sigma[i*C + j] = temp_s[j].x;
			else
				sigma[i*C + j] = temp_s[j].y;
		}
	}

	free(temp_s);
}

//求每一标准差向量的最大分量，并用S返回其下标
void GetMax(double *temp, double *s_max, int *S) {
	double max = 0;

	for (int j = 0; j < C; j++) {
		for (int i = 0; i < n; i++) {
			if (temp[j*n + i] >= max) {
				max = temp[j*n + i];
				s_max[j] = max;			//储存最大值
				S[j] = i;				//储存最大值为第几分量
				//printf("%.3lf\t%d\n", max, i);
			}
		}
		max = 0;
	}
}

//聚类中心分裂
sample *Split(sample z[], int *S, double *s_max, int k) {

	sample zp, zm;			//z+和z-，即zplus和zminus
	sample *temp = (sample*)malloc(C * sizeof(sample));
	double r = 0;

	//赋值zp和zm
	zp.x = z[k].x;	zp.y = z[k].y;
	zm.x = z[k].x;	zm.y = z[k].y;
	r = M*s_max[k];
	//printf("%lf\n", r);
	if (S[k] == 0) {
		zp.x = z[k].x + r;
		zm.x = z[k].x - r;
	}
	else {
		zp.y = z[k].y + r;
		zm.y = z[k].y - r;
	}

	//将z传递于temp中
	for (int i = 0; i < C; i++) {
		temp[i].x = z[i].x;
		temp[i].y = z[i].y;
	}
	z = (sample*)malloc((C + 1) * sizeof(sample));	//扩大z内存

	//更新z
	z[k].x = zp.x;		z[k].y = zp.y;
	z[k + 1].x = zm.x;	z[k + 1].y = zm.y;

	for (int i = 0; i < C + 1; i++) {
		if (i < k) {
			z[i].x = temp[i].x;
			z[i].y = temp[i].y;
		}
		else if (i > k + 1) {
			z[i].x = temp[i - 1].x;
			z[i].y = temp[i - 1].y;
		}
	}
	C++;

	for (int i = 0; i < C; i++) {
		//printf("%.2lf\t%.2lf\n", z[i].x, z[i].y);
	}
	return z;
}

//计算每两个聚类中心间的距离
//	sample *z
double *GetDisCC(sample *z) {

	double *dij = (double *)malloc(C*C * sizeof(double));

	for (int i = 0; i < C; i++) {
		for (int j = 0; j < C; j++) {
			dij[i*C + j] = GetDis(z, z, i, j);
		}
	}
	return dij;
}

//dij与θc的比较，若小于则取出组成集合
//double *dij
double *Compare(double *dij,int &ok) {
	double *dt = (double*)malloc(C*C * sizeof(double));
	for (int i = 0; i < C; i++) {
		for (int j = 0; j < C; j++) {
			if ((i < j) && (dij[i*C + j] < θC)) {
				dt[i*C + j] = 1;
				ok = 1;
			}		
			else
				dt[i*C + j] = 0;
		}
	}
	return dt;
}

//合并聚类中心
void Combine() {}

//算法主要部分
void ISOTADA(struct sample x[], sample *z) {

	double *d = (double*)malloc(N * sizeof(double));
	double *mean_d = (double*)malloc(N * sizeof(double));
	double m_d = 0;
	int *f = (int*)malloc(C*N * sizeof(int));
	int *_N = (int*)malloc(C * sizeof(int));		//存储每一聚类中心样本数
	double *sigma = (double*)malloc(C*n * sizeof(double));
	double *s_max = (double*)malloc(C * sizeof(double));
	int *S = (int*)malloc(C * sizeof(int));
	double *dij = (double *)malloc(C*C * sizeof(double));
	double *dt;
	int ok = 0;

	//初始化f和_N，f也需要初始化
	for (int i = 0; i < C; i++) {
		for (int j = 0; j < N; j++) {
			f[i*N + j] = 1000;	
		}
		_N[i] = 0;
	}

	//以下为正式迭代内容

	couter++;	//迭代次数+1

	//第三步
	Classify(x, z, d, f, _N);

	//第五步
	Correction(x, z, f, _N);

	//第六、七步
	GetMean(x, z, f, _N, mean_d, m_d);

	
	//第八步 判断
	if (couter == I) {
		θC = 0;
		dij = GetDisCC(z);	//第十二步
		dt = Compare(dij, ok);	//第十三步
		Combine();	//第十四步

		printf("sss\n");

		for (int i = 0; i < C; i++) {
			for (int j = 0; j < N; j++) {
				//if (f[i*N + j] < 1000)
				printf("%d\t", f[i*N + j]);
			}
			printf("\n");
		}

		return;		//第十五步			
	}
	else if (C <= K / 2) {
		
		//第九步
		GetVector(x, z, _N, sigma);

		//第十步
		GetMax(sigma, s_max, S);

		//第十一步
		int temp_C = C;
		for (int i = 0; i < temp_C; i++) {
			if (s_max[i] > θS) {
				if (((mean_d[i] > m_d) && (_N[i] > 2 * (θN + 1))) || (C <= K / 2)) {
					//进行分裂
					z = Split(z, S, s_max, i);

					ISOTADA(x, z);
				}
				else {
					GetDisCC(z);	//第十二步
					Compare(dij,ok);	//第十三步
					Combine();	//第十四步

					//第十五步
					if (couter != I) {
						ISOTADA(x, z);
					}
					else
						return;
				}
			}
			
		}

	}
	else if (couter % 2 == 0 || C >= 2 * K) {
		dij = GetDisCC(z);	
		dt = Compare(dij, ok);
		Combine();	//第十四步

		if (!ok) {
			//第十五步
			if (couter != I) {
				ISOTADA(x, z);
			}
			else
				return;
		}
	}

	free(d);
	free(mean_d);
	free(f);
	free(_N);
	free(sigma);
	free(s_max);
	free(S);
}

int main()
{
	//原始样本点
	int X[16] = {0,0,
				1,1,
				2,2,
				4,3,
				5,3,
				4,4,
				5,4,
				6,5};
	//输入样本点
	sample x[8];
	int k = 0;
	for (int i = 0; i < 8; i++) {
		x[i].x = X[k++];
		x[i].y = X[k++];
		//printf("x[%d].x=%d\tx[%d].y=%d\n", i, x[i].x, i, x[i].y);
	}

	sample *z;

	couter = 0;//迭代次数为0

	//第一、二步
	//设置初始参数
	θN = 1, θS = 1, θC = 4, L = 0, I = 4, C = 1;
	
	//输入初始聚类中心
	z = (sample*)malloc(C * sizeof(sample));
	z[0].x = 0; z[0].y = 0;
	
	//开始迭代
	ISOTADA(x, z);


	free(z);
	printf("\n\nOver!\n");
    return 0;
}

