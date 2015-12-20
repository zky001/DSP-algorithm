/****************************************************************************/
/*                                                                          */
/*              快速傅里叶变换 / 快速傅里叶逆变换测试                       */
/*                                                                          */
/*              2014年04月20日                                              */
/*                                                                          */
/****************************************************************************/
#include <stdio.h>                  // C 语言标准输入输出函数库
#include <math.h>                   // C 数学函数库

#include "mathlib.h"                // DSP 数学函数库
#include "dsplib.h"                 // DSP 函数库

/****************************************************************************/
/*                                                                          */
/*              宏定义                                                      */
/*                                                                          */
/****************************************************************************/
// 软件断点
#define SW_BREAKPOINT     asm(" SWBP 0 ");

// 快速傅里叶变换
// π 及 浮点数极小值
#define PI                3.14159
#define F_TOL             (1e-06)

/****************************************************************************/
/*                                                                          */
/*              全局变量                                                    */
/*                                                                          */
/****************************************************************************/
// 快速傅里叶变换测试
// 测试快速傅里叶变换点数
// 注意:TI DSP库 最大支持一次性计算 128K 个点的 FFT
#define Tn  1024
// 采样频率
#define Fs  1000.0
// 信号
float Input[2*Tn+4];
// FFT 输入信号
#pragma DATA_ALIGN(CFFT_In, 8);
float CFFT_In[2*Tn+4];
// FFT 输入信号 副本
float CFFT_InOrig[2*Tn+4];
// FFT 输出
#pragma DATA_ALIGN(CFFT_Out, 8);
float CFFT_Out[2*Tn+4];
// IFFT 输出
#pragma DATA_ALIGN(CFFT_InvOut, 8);
float CFFT_InvOut[2*Tn+4];
// 中间运算临时变量
float CTemp[2*Tn+4];
// 存储旋转因子
float Cw[2*Tn];
// 模
float Cmo[Tn+2];

// 二进制位翻转
#pragma DATA_ALIGN (brev, 8);
unsigned char brev[64]=
{
	0x0, 0x20, 0x10, 0x30, 0x8, 0x28, 0x18, 0x38,
	0x4, 0x24, 0x14, 0x34, 0xc, 0x2c, 0x1c, 0x3c,
	0x2, 0x22, 0x12, 0x32, 0xa, 0x2a, 0x1a, 0x3a,
	0x6, 0x26, 0x16, 0x36, 0xe, 0x2e, 0x1e, 0x3e,
	0x1, 0x21, 0x11, 0x31, 0x9, 0x29, 0x19, 0x39,
	0x5, 0x25, 0x15, 0x35, 0xd, 0x2d, 0x1d, 0x3d,
	0x3, 0x23, 0x13, 0x33, 0xb, 0x2b, 0x1b, 0x3b,
	0x7, 0x27, 0x17, 0x37, 0xf, 0x2f, 0x1f, 0x3f
};

/****************************************************************************/
/*                                                                          */
/*              函数声明                                                    */
/*                                                                          */
/****************************************************************************/
// 产生旋转因子
void tw_gen(float *w, int n);
// FFT 测试
void FFTTest();

/****************************************************************************/
/*                                                                          */
/*              主函数                                                      */
/*                                                                          */
/****************************************************************************/
int main(void)
{
	// FFT 测试
	FFTTest();

	// 断点
	SW_BREAKPOINT;
}

/****************************************************************************/
/*                                                                          */
/*              快速傅里叶变换测试                                          */
/*                                                                          */
/****************************************************************************/
// 产生旋转因子
void tw_gen(float *w, int n)
{
	int i,j,k;
	double x_t,y_t,theta1,theta2,theta3;

	for(j=1,k=0;j<=n>>2;j=j<<2)
	{
		for(i=0;i<n>>2;i += j)
		{
			theta1=2*PI*i/n;
			x_t=cos(theta1);
			y_t=sin(theta1);
			w[k]=(float)x_t;
			w[k+1]=(float)y_t;

			theta2=4*PI*i/n;
			x_t=cos(theta2);
			y_t=sin(theta2);
			w[k+2]=(float)x_t;
			w[k+3]=(float)y_t;

			theta3=6*PI*i/n;
			x_t=cos(theta3);
			y_t=sin(theta3);
			w[k+4]=(float)x_t;
			w[k+5]=(float)y_t;
			k+=6;
		}
	}
}

// 快速傅里叶变换
void FFTTest(void)
{
	// 产生待测试信号
	unsigned int i;
	for (i=0;i<Tn;i++)
		Input[i]=5*sin(2*PI*150*(i/Fs))+15*sin(2*PI*350*(i/Fs));

	// 确定快速傅里叶变换基
	unsigned char rad;
	if(Tn==16 || Tn==64 || Tn==256 || Tn==1024 || Tn==4096 || Tn==16384 || Tn==65536)
		rad=4;
	else if(Tn==8 || Tn==32 || Tn==128 || Tn==512 || Tn==2048 || Tn==8192 || Tn==32768)
		rad=2;
	else
	{
		printf ("不支持 计算 %d 点快速傅里叶变换！\n",Tn);
		return;
	}

	// 复数 FFT
	for (i=0;i<2*Tn;i++)
		CFFT_In[i]=0.0;
	for (i=0;i<Tn;i++)
	{
		CFFT_In[2*i]=Input[i];		// 实部
		CFFT_In[2*i+1]=0;     		// 虚部为 0
	}

	// 保留一份输入信号副本
	memcpy(CFFT_InOrig,CFFT_In,2*Tn*sizeof(float));

	// 产生旋转因子
	tw_gen(Cw,Tn);

	// FFT 计算
	DSPF_sp_fftSPxSP(Tn,CFFT_In,Cw,CFFT_Out,brev,rad,0,Tn);

	// 计算振幅
	for(i=0;i<Tn;i++)
		Cmo[i]=0.0;
	for(i=0;i<Tn+2;i++)
	{
		Cmo[i]=sqrtsp(CFFT_Out[2*i]*CFFT_Out[2*i]+CFFT_Out[2*i+1]*CFFT_Out[2*i+1]);
		Cmo[i]=Cmo[i]*2/Tn;
	}

	// 保留一份 FFT 结果副本
	memcpy(CTemp,CFFT_Out,2*Tn*sizeof(float));

	// IFFT 计算
	DSPF_sp_ifftSPxSP(Tn,CFFT_Out,Cw,CFFT_InvOut,brev,rad,0,Tn);

	// 恢复 FFT 结果
	memcpy(CFFT_Out,CTemp,2*Tn*sizeof(float));

	printf("\n复数 FFT 测试结果:");

	unsigned char Flag;
	for(i=0;i<Tn;i++)
		if(abs(CFFT_InOrig[i]-CFFT_InvOut[i])>F_TOL)
			Flag=1;

	if(Flag==1)
			printf ("失败！\n");
	else
			printf ("成功！\n");
}
