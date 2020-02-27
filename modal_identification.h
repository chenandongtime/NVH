#pragma once
#pragma comment(lib,"v4501v.lib")
#include"matlib.h"
#include<vector>

using std::vector;
Mm load_strain(int n, i_o_t, Mm& t, Mm& strain, double& f);//��txt�ļ�����Ӧ������
void SIMO_ITD(int N, Mm strain, double f,i_o_t,Mm& F,Mm& D,Mm& S);//ITD��ģ̬�����Ǳ�

void steady_diagram(vector<vector<Mm>>& data, i_o_t,Mm& F,Mm& D,Mm& S);//��̬ͼģ̬����ʶ��
inline void mytabulate(Mm m,i_o_t,Mm& table);//ͳ��ÿ������Ԫ�ظ�������
