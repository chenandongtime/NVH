#include"modal_identification.h"
#include<iostream>
void main()
{
	initM(MATCOM_VERSION);//����

	format(TM("long"));
	
	int n;//���������
	n = 12;
	int N = 2 * n + 1;
	Mm strain, t;
	double f;
	load_strain(n, i_o, t, strain, f);//��ȡ����
	//plot((CL(t), strain(12, c_p)));
	Mm F, D, S;
	vector<vector<Mm>> data;
	vector<Mm> VF, VD, VS, VN;
	//ITD��������׼���ģ̬�������������data��
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

	steady_diagram(data,i_o,F,D,S);//��̬ͼģ̬ʶ��

	//���������ı�
	Mm fno,fid;
	fno = TM("out_ITD.txt");
	fid = fopen(fno, TM("w"));
	fprintf(fid, TM("    ����    Ƶ�ʣ�Hz)    �����(%%)    ����ϵ��\r\n"));
	Mm type;
	type = TM("%10.4f %10.4f %10.4f");
	for (int i = 1; i <= n;i++)
	{
		type = strcat((CL(type), TM(" %10.6f")));//���������ʽ

	}
	type = strcat((CL(type), TM("\r\n")));
	Mm note;
	//�������
	for (int i = 1; i <= length(F); i++)
	{
		note = (CL(i), F(i), D(i), S(i, c_p));
		fprintf(fid,type, note);
	}
	fclose(fid);
	exitM();
}