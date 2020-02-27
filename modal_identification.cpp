#include"modal_identification.h"

Mm load_strain(int n, i_o_t, Mm& t, Mm& strain, double& f)
{
	initM(MATCOM_VERSION);//启动
	format(TM("long"));//设置格式
	
	int N = 2 * n + 1;//读取列数
	Mm filename, pathname;//路径与文件名
	uigetfile(TM("*.txt"), TM("Pick a file"), i_o, filename, pathname);//打开文件读取对话框
	Mm path;//完整路径
	Mm fid;//txt读取头
	path = strcat((CL(pathname), filename));//组成完整路径
	fid = fopen(path, TM("r"));//读取文件
							   //display(fid);
	Mm note;//临时读取变量
//读取采样频率，快速过度到数据部分
	while (1)
	{
		note = fscanf(fid, TM("%s"), 1);//读取一个临时字符串
		if (*strcmp(note, TM("Interleave:")).addr())
		{
			f = 2000 / eval(fscanf(fid, TM("%s"), 1)).r(1);//保存采样频率
		}
//过渡到数据部分
		if (*strcmp(note, TM("Time")).addr())
		{
			while (1)
			{
				note = fscanf(fid, TM("%s"), 1);
				if (*strcmp(note, TM("FBG_D1")).addr() || *strcmp(note, TM("FBG_C4")).addr())
				{
					break;
				}
			}
			break;
		}
	}
	Mm type;//数据读取类型
	type = TM("%f");
	for (int i = 2; i <= N; i++)
	{
		type = strcat((CL(type), TM(" %f")));//每次读取N个数据
	}
	Mm data,a;
	a = (BR(N), Inf);//N行无穷列
	for (int i = 1; i <= N * 10; i++)
	{
		note = fscanf(fid, TM("%s"), 1);//舍去前11列的数据
	}
	data = fscanf(fid, type, a);//将可用数据存入矩阵
	fclose(fid);//关闭文件打开
	double m= fix(data.cols() / 2);
	Mm dt;
	dt = 1 / f;
	t = linspace(0, (2 * m - 1)*dt, 2 * m);
 //提取应变数据
	for (int i = 1; i <= n; i++)
	{
		strain(i, c_p) = data(i + 1, c_p);
	}
	exitM();//关闭计算
	return data;
}

void SIMO_ITD(int N, Mm strain, double f, i_o_t, Mm& F, Mm& D, Mm& S)
{
	initM(MATCOM_VERSION);
	double n = fix(strain.cols() / 2);
	double Nn = strain.rows();//测点个数
	strain = strain(c_p, colon(1, 1, 2 * n));//取2n个时间段数据
	double dt = 1 / f;//计算时间间隔
	double M = 2 * N;//特征值个数
	double L = 2 * n - M - 10;
	Mm x0;
	x0 = strain;
	Mm med;
	//延时矩阵组装
	for (int i = Nn + 1; i <= M; i++)
	{
		med = x0(i - Nn, linspace(2, 2 * n,2*n-1));
		med(med.size() + 1) = 0;
		x0(i, c_p) = med;
	}
	Mm x1;
	Mm x2;
	x1 = (BR(x0(c_p, colon(1, 1, L))), semi, x0(c_p, colon(2, 1, L + 1)));
	x2 = (BR(x0(c_p, colon(2, 1, L+1))), semi, x0(c_p, colon(3, 1, L + 2)));
	
	//计算相关矩阵
	Mm B;
	B = x2*ctranspose(x1)*pinv(x1*ctranspose(x1));
	Mm A,V;
	eig(B, i_o, A, V);//求解特征值与特征向量
	Mm U;
	for (int k = 1; k <= M;k++)
	{
		U(k) = V(k, k);
	}
	U = ctranspose(U);
	Mm F1;
	F1 = abs(log(U)) / (2 * pi*dt);//计算固频率
	Mm D1,re,im,d;
	re = real(log(U));
	im = imag(log(U));
	d = rdivide(im, re);
	D1 = msqrt(rdivide(1,power(d,2)+1));//计算阻尼
	int l = 1;
	Mm V0;
	//滤除
	for (int k = 1; k <= M; k++)
	{
		re = real(U(k));
		im = imag(U(k));
		if (*abs(re).addr() <= 1 && *abs(im).addr() <= 1)
		{
			V0(l) = U(k);
			l = l + 1;
		}
	}
    //计算留数
	Mm Va,cV0;
	cV0 = conj(V0);
	for (int k = 0; k <= 2 * n - 1; k++)
	{
		Va(k + 1, c_p) = power(cV0, k);
	}
	//计算模态振型
	Mm S1;
	S1 = inv(conj(ctranspose(Va))*Va)*conj(ctranspose(Va))*ctranspose(strain);

	Mm F2, I;
	sort(F1, i_o, F2, I);//从小到大排序
	//剔除非共轭项
	int m = 0;
for (int k = 1; k <= M - 1; k++)
{
	if (F2.r(k) != F2.r(k + 1))
	{
		continue;
	}
	m = m + 1;
	l = I.r(k);
	F(m) = F1(l);//固有频率
	D(m) = D1(l);//阻尼比
	S(m, c_p) = S1(l, c_p);//模态振型 行为阶次，列为测点
}
for (int k = 1; k <= m; k++)
{
	S(k, c_p) = imag(S(k, c_p));//应变振型
}
exitM();
}

void steady_diagram(vector<vector<Mm>>& data, i_o_t, Mm& F, Mm& D, Mm& S)
{
	initM(MATCOM_VERSION);
	int n = data[1].size();//测点个数
	Mm SF;//稳定极点向量
	figure(CL(1));
	hold(TM("on"));
	Mm F0, D0, S0, N0, F1, D1, S1, N1, N;
	for (int i = 0; i <= n - 2; i++)
	{
		//取出第i阶和i+1阶计算结果
		F0 = data[0][i];
		D0 = data[1][i];
		S0 = data[2][i];
		N0 = length(F0);
		F1 = data[0][i + 1];
		D1 = data[1][i + 1];
		S1 = data[2][i + 1];
		N1 = length(F1);
		N = data[3][i];
		//绘制稳定图
		Mm mac;
		for (int j = 1; j <= N0.r(); j++)
		{
			for (int k = 1; k <= N1.r(); k++)
			{
				if ((abs(F0(j) - F1(k)) / abs(F1(k))).r() <= 0.01)
				{
					mac = mpower(S0(j, c_p)*ctranspose(S1(k, c_p)), 2) / ((S0(j, c_p)*transpose(S0(j, c_p)))*(S1(j, c_p)*transpose(S1(j, c_p))));//计算mac值
					if ((abs(D0(j) - D1(k)) / abs(D1(k))).r() <= 0.05)
					{
						if (1 - mac.r() <= 0.02)
						{
							plot((CL(F0(j)), N, TM("rv"), TM("linewidth"), 1));//绘出稳定极点，红色倒三角
							SF = (BR(SF), F0(j));
						}
						else
						{
							plot((CL(F0(j)), N, TM("b+"), TM("linewidth"), 1));//会出频率和阻尼稳定极点，蓝色+
						}
					}
					else
					{
						if (1 - mac.r() <= 0.02)
						{
							plot((CL(F0(j)), N, TM("gs"), TM("linewidth"), 1));//绘出频率稳定振型稳定极点，绿色方块
						}
						else
						{
							plot((CL(F0(j)), N, TM("c*"), TM("linewidth"), 1));//绘出频率稳定极点，蓝绿色*
						}
					}
				}
			}
		}
	}
	//title(TM("频率稳定图"));
	//xlabel(TM("频率/Hz"));
	//ylabel(TM("求解模型阶次"));
	grid(TM("on"));
	Mm F_all;
	F_all = rdivide(round(times(SF, 10)), 10);//将稳定极点保留一位有效数字

	Mm table,kong;
	mytabulate(F_all,i_o,table);//统计独立元素个数
	int i = 1;
	//合并相近元素
	while (i <= table.rows()-1)
	{
		if ((abs(table(i, 1) - table(i + 1, 1))/table(i,1) <= 0.01).r())
		{
			table(i, 2) = table(i, 2) + table(i + 1, 2);
			table(i + 1, c_p) = kong;
		}
		else
		{
			i = i + 1;
		}
	}
	Mm y;
	F = y;
	F0 = y;
	D = y;
	S = y;
	y = colon(0, 1, N);
	//将出现稳定次数超过4次的极点作为最终解的近似
	for (int i = 1; i <= table.rows(); i++)
	{
		if ((table(i, 2) >= 4).r())
		{
			F0 = (BR(F0), table(i, 1));
		}
	}
	display(F0);
    /*
	i = 1;
	//删除相同模态
	while (1)
	{
		if (i < length(F0))
		{
			if ((abs(F0(i) - F0(i + 1)) / abs(F0(i))).r() <= 0.01)
			{
				F0(i + 1) = kong;
			}
		}
		else
			break;
		i = i + 1;
	}
	*/

    //取出现稳定最高阶计算结果作为最终结果
	int a;
	for (int i = 1; i <= F0.size(); i++)
	{
		a = n - 1;
		F1 = data[0][a];
		D1 = data[1][a];
		S1 = data[2][a];
		for (int j = 1; j <= F1.size(); j++)
		{
			if ((abs(F1(j) - F0(i)) / abs(F0(i))).r() <= 0.01)
			{
				F = (BR(F), F1(j));
				D = (BR(D), D1(j));
				S = (BR(S),semi, S1(j, c_p));
				break;
			}
			if (j == F1.size())
			{
				a = a - 1;
				j = 0;
				F1 = data[0][a];
				D1 = data[1][a];
				S1 = data[2][a];
			}
		}

	}
	Mm x;
    hold(TM("on"));
	//将稳定解在图中以红色竖线标记
	for (int i = 1; i <= length(F); i++)
	{
		x = times(F(i), ones(y.size()));
		//line((CL(x), y, TM("Color"), TM("r"), TM("LineWidth"), 1));
		plot((CL(x), y, TM("r")));
	}
	exitM();
}

inline void mytabulate(Mm m,i_o_t,Mm& table)
{
	Mm um;
	um = unique(m);//取出独立元素
	Mm count;
	count = zeros(1, um.size());
	//统计每个独立元素个数
	for (int i = 1; i <= m.size(); i++)
	{
		for (int j = 1; j <= um.size(); j++)
		{
			if ((m(i) == um(j)).r())
			{
				count(j) = count(j) + 1;
				break;
			}
		}
	}
	//将统计结果保留在table内，第一列为独立元素，第二列为对应个数
	table(c_p, 1) = ctranspose(um);
	table(c_p, 2) = ctranspose(count);
}