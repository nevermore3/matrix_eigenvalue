/********************************
本程序输入矩阵号从零开始
考虑了盖尔圆为一点的情况
*********************************/
#include "type.h"
void DV_GetCeig(double *ceig, SPMultiply SPMul, Info info);
void DV_GetInterval(double *ceig, Info info);
int DV_GetNeig(SPMultiply SPMul, LinearSystem LSystem, Info info);
void DV_GetTag(int *tag, int neig_loc, int times, Info info);
void DV_Modify(int tag, double *ceig, Info info);

void Division(SPMultiply SPMul, LinearSystem LSystem, Info info)
{
    int i;

    //根据矩阵盖尔圆以及最小特征值得初始的估计特征值分布
    double *ceig;
    if (info.intervalid == 0)
        ceig = calloc(info.N, sizeof(double));
    DV_GetCeig(ceig, SPMul, info);
    //迭代循环，划分区间
    int tag_loc = 0;
    int neig_loc;
    int times = 1;

    while (tag_loc == 0)
    {
        times++;
        info.times=times;

        //分割区间
        DV_GetInterval(ceig, info);

        //计算每个子区间内特征值的个数
        neig_loc = DV_GetNeig(SPMul, LSystem, info);
        //判断终止条件
        DV_GetTag(&tag_loc, neig_loc, times, info);


        //根据子区间特征值个数进行修正
        DV_Modify(tag_loc, ceig, info);
    }

    DV_End(info);
}
