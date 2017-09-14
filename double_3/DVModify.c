#include "type.h"
void Modify_Info(Info info);
void Modify_Ceig(double *ceig, Info info);
void Modify_Scal(Info info);

void DV_Modify(int tag, double *ceig, Info info)
{
    //根据子区间特征值个数进行修正
    if (tag == 0)
    {
        Modify_Info(info);

        Modify_Ceig(ceig, info);

        Modify_Scal(info);

    }
}

/******************************************************
history
*******************************************************/
void Modify_Info(Info info)
{
    if (info.myid < info.ninterval)
    {
        int *neig = info.neig;
        int *sneig = info.sneig;
        int *history = info.history;
        double *interval = info.interval;
        double *modifyregion = info.modifyregion;

        int N = info.N;
        int ninterval = info.ninterval;
        int myid = info.myid;
        int M0 = info.M0;
        int i;
        int maxneig = neig[0];
        printf("maxneig  %d  M0  %d   \n", maxneig, info.M0);
        double dist;
        if (myid == 0)
        {
            dist = N * info.propotion / ninterval;//每个子区间目标特征值的个数
            printf("sneig  %d   dist %lf history  %d  %d\n", sneig[0], dist, history[0], history[1]  );

            if (maxneig < M0)
            {
                if (sneig[0] < dist && sneig[0] >= history[0] - 1)
                {
                    history[0] = sneig[0];
                    modifyregion[0] = interval[1];
                }
                if (sneig[0] > dist && sneig[0] <= history[1] - 1)
                {
                    history[1] = sneig[0];
                    modifyregion[1] = interval[1];
                }
            }
            // else
            // {
            //     history[1] = sneig[myid];
            //     modifyregion[1] = interval[1];
            // }

        }
        else
        {
            for (i = 1; i < myid; i++)
                if (maxneig < neig[i])
                    maxneig = neig[i];
            dist = (myid + 1) * N * info.propotion / ninterval;
            if (maxneig < M0)
            {
                if (sneig[myid - 1] < dist && sneig[myid - 1] >= history[0] - 1)
                {
                    history[0] = sneig[myid - 1];
                    modifyregion[0] = interval[0];
                }
                if (sneig[myid - 1] > dist && sneig[myid - 1] <= history[1] - 1)
                {
                    history[1] = sneig[myid - 1];
                    modifyregion[1] = interval[0];
                }
            }

            if (maxneig < neig[myid])
                maxneig = neig[myid];
            if (maxneig < M0)
            {
                if (sneig[myid] < dist && sneig[myid] >= history[0] - 1)
                {
                    history[0] = sneig[myid];
                    modifyregion[0] = interval[1];
                }
                if (sneig[myid] > dist && sneig[myid] <= history[1] - 1)
                {
                    history[1] = sneig[myid];
                    modifyregion[1] = interval[1];
                }
            }
            // else
            // {
            //     history[1] = sneig[myid];
            //     modifyregion[1] = interval[1];

            // }

        }
        if (history[0] == 0)
            history[0] = 1;

        printf("modyfyinfo:  myid  %d \n\t\t history: %d  %d  %lf \n\t\t modifyregion %lf %lf\n", myid, history[0], history[1], dist, modifyregion[0], modifyregion[1]);

    }

}


/********************************************************
根据每一个子区间内特征值的个数对特征值的分布进行修正
************************************************************/
void Modify_Ceig(double *ceig, Info info)
{
    if (info.myid < info.ninterval)
    {
        int N = info.N;

        int *history = info.history;
        double *interval = info.interval;
        double *modifyregion = info.modifyregion;

        int i, j;
        int beginid = history[0] - 1;
        int endid = history[1];
        double a = (modifyregion[1] - modifyregion[0]) / (ceig[endid - 1] - ceig[beginid]);
        double b = modifyregion[0] - a * ceig[beginid];

        for (i = beginid; i < endid; i++)
        {
            ceig[i] = a * ceig[i] + b;
            if (ceig[i] < ceig[i - 1] && i > beginid)
            {
                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~i %d a  %lf  b %lf  ceig[i]  %lf\n", i, a, b, ceig[i]);
            }

        }
    }
}


/******************************************************
Scal
*******************************************************/
void Modify_Scal(Info info)
{
    if (info.myid < info.ninterval)
    {
        int times = info.times;
        double *scal = info.scal;
        int *sneig = info.sneig;
        int *neig = info.neig;
        int *history = info.history;

        int M0 = info.M0;
        int myid = info.myid;

        int i;

        int tag = 0;
        int maxneig = neig[0];
        for (i = 1; i <= myid; i++)
            if (maxneig < neig[i])
                maxneig = neig[i];
        if (maxneig == M0)
            tag = 1;

        double dist = (myid + 1.0) * info.N * info.propotion / info.ninterval;
        int sneig_loc = sneig[myid];

        if (sneig_loc >= dist)
        {
            if (sneig_loc / dist > 1.5 || history[1] == info.N)
            {
                if (history[0] / dist < 0.6)
                    scal[0] = 0.6;
                else
                    scal[0] = 0.5 + 0.5 * history[0] / dist;
            }
            else if (sneig_loc / dist > 1.3)
            {
                if (history[0] / dist < 0.8)
                    scal[0] = 0.8;
                else
                    scal[0] = 0.5 + 0.5 * history[0] / dist;
            }
            else if (sneig_loc / dist > 1.2)
            {
                if (history[0] / dist < 0.9)
                    scal[0] = 0.9;
                else
                    scal[0] = 0.5 + 0.5 * history[0] / dist;
            }

            else if (sneig_loc / dist > 1.1)
            {
                if (history[0] / dist < 0.95)
                    scal[0] = 0.95;
                else
                    scal[0] = 0.5 + 0.5 * history[0] / dist;
            }
            else
            {
                if (history[1] / dist < 1.15)
                    scal[0] = 1;
                else
                    scal[0] = 1.1;
            }

        }
        else
        {
            if (history[1] / dist > 4)
            {
                scal[0] = 1.2 * dist / sneig_loc ;
                if (scal[0] > 2)
                {
                    scal[0] = times + 1;
                    if (scal[0] > history[1] / dist)
                        scal[0] = (int)(history[1] / dist) - 2;
                }
            }
            else if (sneig_loc / dist < 0.7)
            {
                if (history[1] / dist > 1.3)
                    scal[0] = 1.3;
                else
                    scal[0] = 0.5 + 0.5 * history[1] / dist;
            }
            else if (sneig_loc / dist < 0.8)
            {
                if (history[1] / dist > 1.2)
                    scal[0] = 1.2;
                else
                    scal[0] = 0.5 + 0.5 * history[1] / dist;
            }
            else if (sneig_loc / dist < 0.9)
            {
                if (history[1] / dist > 1.1)
                    scal[0] = 1.1;
                else
                    scal[0] = 0.5 + 0.5 * history[1] / dist;
            }
            else
            {
                if (history[1] / dist < 1.11)
                    scal[0] = 1.01;
                else
                    scal[0] = 1.05;
            }


        }


    }

}
