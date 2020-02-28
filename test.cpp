#include"modal_identification.h"
#include<iostream>
void main()
{
	initM(MATCOM_VERSION);//启动

	format(TM("long"));
	
	int n;//采样点个数
	n = 12;
	int N = 2 * n + 1;
	Mm strain, t;
	double f;
	load_strain(n, i_o, t, strain, f);//读取数据
	//plot((CL(t), strain(12, c_p)));
	Mm F, D, S;
	vector<vector<Mm>> data;
	vector<Mm> VF, VD, VS, VN;
	//ITD法计算各阶计算模态，将结果保留在data内
	for (int num = 1; num <=50; num++)
	{
		SIMO_ITD(num, strain, f, i_o, F, D, S);
		std::cout << num << std::endl;
		VF.push_back(F);
		VD.push_back(D);
		VS.push_back(S);
		VN.push_back(num);
	}
	data.push_back(VF);
	data.push_back(VD);
	data.push_back(VS);
	data.push_back(VN);

	steady_diagram(data,i_o,F,D,S);//稳态图模态识别

	//输出结果到文本
	Mm fno,fid;
	fno = TM("out_ITD.txt");
	fid = fopen(fno, TM("w"));
	fprintf(fid, TM("    阶数    频率（Hz)    阻尼比(%%)    振型系数\r\n"));
	Mm type;
	type = TM("%10.4f %10.4f %10.4f");
	for (int i = 1; i <= n;i++)
	{
		type = strcat((CL(type), TM(" %10.6f")));//数据输出格式

	}
	type = strcat((CL(type), TM("\r\n")));
	Mm note;
	//输出数据
	for (int i = 1; i <= length(F); i++)
	{
		note = (CL(i), F(i), D(i), S(i, c_p));
		fprintf(fid,type, note);
	}
	fclose(fid);
	exitM();
}