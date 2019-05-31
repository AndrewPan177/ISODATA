// ISODATA.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "math.h"
#include "malloc.h"

//�洢������
struct sample{
	double x;
	double y;
};

//��������
/*
 * C		��ʼ��������
 * ��N		һ���������Ӧ���е�������Ŀ
 * ��S		һ�����������׼����ֵ
 * ��C		��������֮��������ֵ�����鲢ϵ��
 * L		��һ�ε����п��Թ鲢������������
 * I		���������������
 */
int C, ��N, ��S, ��C, L, I;

int couter;		//��������

#define K 2		//��Ҫ�����������
#define N 8		//������
#define n 2		//����ά��
#define M 0.5	//��֪����ɶ�Ĳ���

//����Сֵ���±�
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

//����������������ĵľ���
double GetDis(struct sample *x, sample *z, int i, int j) {
	double distance = 0;
	distance = sqrt(pow((x[i].x - z[j].x), 2)
		+ pow((x[i].y - z[j].y), 2));
	return distance;
}

//��N����������С������з���
void Classify(sample *x, sample *z, double *d, int *f, int *_N) {
	double *temp = (double*)malloc(C * sizeof(double));

	for (int j = 0; j < N; j++) {
		for (int i = 0; i < C; i++) {
			temp[i] = GetDis(x, z, j, i);
		}
		d[j] = temp[GetMin(temp)];
		f[GetMin(temp)*N + j] = j;	//������
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
	//printf("\n_N ��Ŀ��\t");
	for (int i = 0; i < C; i++) {
		//printf("%d\t", _N[i]);
	}

	free(temp);
}

//������������
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

//��ÿһ������f��������������������������ľ����ƽ��ֵ
//������������������Ӧ�������ĵľ����ƽ��ֵ
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

//����ÿһ�����������������ľ���ı�׼������
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

	//��ֵ��sigma
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

//��ÿһ��׼��������������������S�������±�
void GetMax(double *temp, double *s_max, int *S) {
	double max = 0;

	for (int j = 0; j < C; j++) {
		for (int i = 0; i < n; i++) {
			if (temp[j*n + i] >= max) {
				max = temp[j*n + i];
				s_max[j] = max;			//�������ֵ
				S[j] = i;				//�������ֵΪ�ڼ�����
				//printf("%.3lf\t%d\n", max, i);
			}
		}
		max = 0;
	}
}

//�������ķ���
sample *Split(sample z[], int *S, double *s_max, int k) {

	sample zp, zm;			//z+��z-����zplus��zminus
	sample *temp = (sample*)malloc(C * sizeof(sample));
	double r = 0;

	//��ֵzp��zm
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

	//��z������temp��
	for (int i = 0; i < C; i++) {
		temp[i].x = z[i].x;
		temp[i].y = z[i].y;
	}
	z = (sample*)malloc((C + 1) * sizeof(sample));	//����z�ڴ�

	//����z
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

//����ÿ�����������ļ�ľ���
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

//dij���c�ıȽϣ���С����ȡ����ɼ���
//double *dij
double *Compare(double *dij,int &ok) {
	double *dt = (double*)malloc(C*C * sizeof(double));
	for (int i = 0; i < C; i++) {
		for (int j = 0; j < C; j++) {
			if ((i < j) && (dij[i*C + j] < ��C)) {
				dt[i*C + j] = 1;
				ok = 1;
			}		
			else
				dt[i*C + j] = 0;
		}
	}
	return dt;
}

//�ϲ���������
void Combine() {}

//�㷨��Ҫ����
void ISOTADA(struct sample x[], sample *z) {

	double *d = (double*)malloc(N * sizeof(double));
	double *mean_d = (double*)malloc(N * sizeof(double));
	double m_d = 0;
	int *f = (int*)malloc(C*N * sizeof(int));
	int *_N = (int*)malloc(C * sizeof(int));		//�洢ÿһ��������������
	double *sigma = (double*)malloc(C*n * sizeof(double));
	double *s_max = (double*)malloc(C * sizeof(double));
	int *S = (int*)malloc(C * sizeof(int));
	double *dij = (double *)malloc(C*C * sizeof(double));
	double *dt;
	int ok = 0;

	//��ʼ��f��_N��fҲ��Ҫ��ʼ��
	for (int i = 0; i < C; i++) {
		for (int j = 0; j < N; j++) {
			f[i*N + j] = 1000;	
		}
		_N[i] = 0;
	}

	//����Ϊ��ʽ��������

	couter++;	//��������+1

	//������
	Classify(x, z, d, f, _N);

	//���岽
	Correction(x, z, f, _N);

	//�������߲�
	GetMean(x, z, f, _N, mean_d, m_d);

	
	//�ڰ˲� �ж�
	if (couter == I) {
		��C = 0;
		dij = GetDisCC(z);	//��ʮ����
		dt = Compare(dij, ok);	//��ʮ����
		Combine();	//��ʮ�Ĳ�

		printf("sss\n");

		for (int i = 0; i < C; i++) {
			for (int j = 0; j < N; j++) {
				//if (f[i*N + j] < 1000)
				printf("%d\t", f[i*N + j]);
			}
			printf("\n");
		}

		return;		//��ʮ�岽			
	}
	else if (C <= K / 2) {
		
		//�ھŲ�
		GetVector(x, z, _N, sigma);

		//��ʮ��
		GetMax(sigma, s_max, S);

		//��ʮһ��
		int temp_C = C;
		for (int i = 0; i < temp_C; i++) {
			if (s_max[i] > ��S) {
				if (((mean_d[i] > m_d) && (_N[i] > 2 * (��N + 1))) || (C <= K / 2)) {
					//���з���
					z = Split(z, S, s_max, i);

					ISOTADA(x, z);
				}
				else {
					GetDisCC(z);	//��ʮ����
					Compare(dij,ok);	//��ʮ����
					Combine();	//��ʮ�Ĳ�

					//��ʮ�岽
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
		Combine();	//��ʮ�Ĳ�

		if (!ok) {
			//��ʮ�岽
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
	//ԭʼ������
	int X[16] = {0,0,
				1,1,
				2,2,
				4,3,
				5,3,
				4,4,
				5,4,
				6,5};
	//����������
	sample x[8];
	int k = 0;
	for (int i = 0; i < 8; i++) {
		x[i].x = X[k++];
		x[i].y = X[k++];
		//printf("x[%d].x=%d\tx[%d].y=%d\n", i, x[i].x, i, x[i].y);
	}

	sample *z;

	couter = 0;//��������Ϊ0

	//��һ������
	//���ó�ʼ����
	��N = 1, ��S = 1, ��C = 4, L = 0, I = 4, C = 1;
	
	//�����ʼ��������
	z = (sample*)malloc(C * sizeof(sample));
	z[0].x = 0; z[0].y = 0;
	
	//��ʼ����
	ISOTADA(x, z);


	free(z);
	printf("\n\nOver!\n");
    return 0;
}

