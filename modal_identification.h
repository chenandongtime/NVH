#pragma once
#pragma comment(lib,"v4501v.lib")
#include"matlib.h"
#include<vector>

using std::vector;
Mm load_strain(int n, i_o_t, Mm& t, Mm& strain, double& f);//从txt文件加载应变数据
void SIMO_ITD(int N, Mm strain, double f,i_o_t,Mm& F,Mm& D,Mm& S);//ITD法模态参数是被

void steady_diagram(vector<vector<Mm>>& data, i_o_t,Mm& F,Mm& D,Mm& S);//稳态图模态参数识别
inline void mytabulate(Mm m,i_o_t,Mm& table);//统计每个独立元素个数函数
