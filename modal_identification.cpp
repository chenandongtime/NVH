#include"modal_identification.h"

Mm load_strain(int n, i_o_t, Mm& t, Mm& strain, double& f)
{
	initM(MATCOM_VERSION);//����
	format(TM("long"));//���ø�ʽ
	
	int N = 2 * n + 1;//��ȡ����
	Mm filename, pathname;//·�����ļ���
	uigetfile(TM("*.txt"), TM("Pick a file"), i_o, filename, pathname);//���ļ���ȡ�Ի���
	Mm path;//����·��
	Mm fid;//txt��ȡͷ
	path = strcat((CL(pathname), filename));//�������·��
	fid = fopen(path, TM("r"));//��ȡ�ļ�
							   //display(fid);
	Mm note;//��ʱ��ȡ����
//��ȡ����Ƶ�ʣ����ٹ��ȵ����ݲ���
	while (1)
	{
		note = fscanf(fid, TM("%s"), 1);//��ȡһ����ʱ�ַ���
		if (*strcmp(note, TM("Interleave:")).addr())
		{
			f = 2000 / eval(fscanf(fid, TM("%s"), 1)).r(1);//�������Ƶ��
		}
//���ɵ����ݲ���
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
	Mm type;//���ݶ�ȡ����
	type = TM("%f");
	for (int i = 2; i <= N; i++)
	{
		type = strcat((CL(type), TM(" %f")));//ÿ�ζ�ȡN������
	}
	Mm data,a;
	a = (BR(N), Inf);//N��������
	for (int i = 1; i <= N * 10; i++)
	{
		note = fscanf(fid, TM("%s"), 1);//��ȥǰ11�е�����
	}
	data = fscanf(fid, type, a);//���������ݴ������
	fclose(fid);//�ر��ļ���
	double m= fix(data.cols() / 2);
	Mm dt;
	dt = 1 / f;
	t = linspace(0, (2 * m - 1)*dt, 2 * m);
 //��ȡӦ������
	for (int i = 1; i <= n; i++)
	{
		strain(i, c_p) = data(i + 1, c_p);
	}
	exitM();//�رռ���
	return data;
}

void SIMO_ITD(int N, Mm strain, double f, i_o_t, Mm& F, Mm& D, Mm& S)
{
	initM(MATCOM_VERSION);
	double n = fix(strain.cols() / 2);
	double Nn = strain.rows();//������
	strain = strain(c_p, colon(1, 1, 2 * n));//ȡ2n��ʱ�������
	double dt = 1 / f;//����ʱ����
	double M = 2 * N;//����ֵ����
	double L = 2 * n - M - 10;
	Mm x0;
	x0 = strain;
	Mm med;
	//��ʱ������װ
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
	
	//������ؾ���
	Mm B;
	B = x2*ctranspose(x1)*pinv(x1*ctranspose(x1));
	Mm A,V;
	eig(B, i_o, A, V);//�������ֵ����������
	Mm U;
	for (int k = 1; k <= M;k++)
	{
		U(k) = V(k, k);
	}
	U = ctranspose(U);
	Mm F1;
	F1 = abs(log(U)) / (2 * pi*dt);//�����Ƶ��
	Mm D1,re,im,d;
	re = real(log(U));
	im = imag(log(U));
	d = rdivide(im, re);
	D1 = msqrt(rdivide(1,power(d,2)+1));//��������
	int l = 1;
	Mm V0;
	//�˳�
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
    //��������
	Mm Va,cV0;
	cV0 = conj(V0);
	for (int k = 0; k <= 2 * n - 1; k++)
	{
		Va(k + 1, c_p) = power(cV0, k);
	}
	//����ģ̬����
	Mm S1;
	S1 = inv(conj(ctranspose(Va))*Va)*conj(ctranspose(Va))*ctranspose(strain);

	Mm F2, I;
	sort(F1, i_o, F2, I);//��С��������
	//�޳��ǹ�����
	int m = 0;
for (int k = 1; k <= M - 1; k++)
{
	if (F2.r(k) != F2.r(k + 1))
	{
		continue;
	}
	m = m + 1;
	l = I.r(k);
	F(m) = F1(l);//����Ƶ��
	D(m) = D1(l);//�����
	S(m, c_p) = S1(l, c_p);//ģ̬���� ��Ϊ�״Σ���Ϊ���
}
for (int k = 1; k <= m; k++)
{
	S(k, c_p) = imag(S(k, c_p));//Ӧ������
}
exitM();
}

void steady_diagram(vector<vector<Mm>>& data, i_o_t, Mm& F, Mm& D, Mm& S)
{
	initM(MATCOM_VERSION);
	int n = data[1].size();//������
	Mm SF;//�ȶ���������
	figure(CL(1));
	hold(TM("on"));
	Mm F0, D0, S0, N0, F1, D1, S1, N1, N;
	for (int i = 0; i <= n - 2; i++)
	{
		//ȡ����i�׺�i+1�׼�����
		F0 = data[0][i];
		D0 = data[1][i];
		S0 = data[2][i];
		N0 = length(F0);
		F1 = data[0][i + 1];
		D1 = data[1][i + 1];
		S1 = data[2][i + 1];
		N1 = length(F1);
		N = data[3][i];
		//�����ȶ�ͼ
		Mm mac;
		for (int j = 1; j <= N0.r(); j++)
		{
			for (int k = 1; k <= N1.r(); k++)
			{
				if ((abs(F0(j) - F1(k)) / abs(F1(k))).r() <= 0.01)
				{
					mac = mpower(S0(j, c_p)*ctranspose(S1(k, c_p)), 2) / ((S0(j, c_p)*transpose(S0(j, c_p)))*(S1(j, c_p)*transpose(S1(j, c_p))));//����macֵ
					if ((abs(D0(j) - D1(k)) / abs(D1(k))).r() <= 0.05)
					{
						if (1 - mac.r() <= 0.02)
						{
							plot((CL(F0(j)), N, TM("rv"), TM("linewidth"), 1));//����ȶ����㣬��ɫ������
							SF = (BR(SF), F0(j));
						}
						else
						{
							plot((CL(F0(j)), N, TM("b+"), TM("linewidth"), 1));//���Ƶ�ʺ������ȶ����㣬��ɫ+
						}
					}
					else
					{
						if (1 - mac.r() <= 0.02)
						{
							plot((CL(F0(j)), N, TM("gs"), TM("linewidth"), 1));//���Ƶ���ȶ������ȶ����㣬��ɫ����
						}
						else
						{
							plot((CL(F0(j)), N, TM("c*"), TM("linewidth"), 1));//���Ƶ���ȶ����㣬����ɫ*
						}
					}
				}
			}
		}
	}
	//title(TM("Ƶ���ȶ�ͼ"));
	//xlabel(TM("Ƶ��/Hz"));
	//ylabel(TM("���ģ�ͽ״�"));
	grid(TM("on"));
	Mm F_all;
	F_all = rdivide(round(times(SF, 10)), 10);//���ȶ����㱣��һλ��Ч����

	Mm table,kong;
	mytabulate(F_all,i_o,table);//ͳ�ƶ���Ԫ�ظ���
	int i = 1;
	//�ϲ����Ԫ��
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
	//�������ȶ���������4�εļ�����Ϊ���ս�Ľ���
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
	//ɾ����ͬģ̬
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

    //ȡ�����ȶ���߽׼�������Ϊ���ս��
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
	//���ȶ�����ͼ���Ժ�ɫ���߱��
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
	um = unique(m);//ȡ������Ԫ��
	Mm count;
	count = zeros(1, um.size());
	//ͳ��ÿ������Ԫ�ظ���
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
	//��ͳ�ƽ��������table�ڣ���һ��Ϊ����Ԫ�أ��ڶ���Ϊ��Ӧ����
	table(c_p, 1) = ctranspose(um);
	table(c_p, 2) = ctranspose(count);
}