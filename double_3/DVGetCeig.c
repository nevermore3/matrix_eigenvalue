#include "type.h"
/****************************************************
通过矩阵的盖尔圆信息以及最小特征值的信息，得初始的估计的特征值分布
*************************************************/
void Ger_setGerinfo(Gerinfo *Gerinfo, SPMultiply SPMul, Info info);
void Ger_getGerinfo(Vector *Gerdot, Density *Gerdensity, Gerinfo Gerinfo);
void Ger_getCdensity(Density *cdensity, Density Gerdensity, Vector  Gerdot);
void Ger_getCeig(double *ceig, int N, Density cdensity);
void Ger_ModifyCeig(double *ceig, int N, double mineig);

void DV_GetCeig(double *ceig, SPMultiply SPMul, Info info)
{
    double *modifyregion = info.modifyregion;
    Gerinfo Gerinfo;
    Ger_setGerinfo(&Gerinfo, SPMul, info);
    if (info.intervalid == 0)
    {
        Vector Gerdot;
        Density Gerdensity;

        //统计盖尔圆的信息
        Ger_getGerinfo(&Gerdot, &Gerdensity, Gerinfo);

        Density cdensity;
        //计算特征值分布概率
        Ger_getCdensity(&cdensity, Gerdensity, Gerdot);

        //特征值分布的重构
        Ger_getCeig(ceig, info.N, cdensity);

        //最小特征值修正
        Ger_ModifyCeig(ceig, info.N, info.mineig);
        modifyregion[0] = ceig[0];
        modifyregion[1] = ceig[info.N - 1];


    }
}

/***************************************
统计盖尔圆信息，根据矩阵A，得到盖尔圆的断点信息以及其对应的概率密度变化量
******************************************/
void MPI_Gatherv_double(double *sendbuf, int sendcount, double *recvbuf, MPI_Comm comm);
void Ger_setGerinfo(Gerinfo *Gerinfo, SPMultiply SPMul, Info info)
{
    //获取输入信息
    int fstrow = SPMul.fstrow;
    int N_loc = SPMul.N_loc;
    int nnz_in = SPMul.nnz_in;
    int *ia_in = SPMul.ia_in;
    int *ja_in = SPMul.ja_in;
    double *a_in = SPMul.a_in;
    int nnz_out = SPMul.nnz_out;
    int *ia_out = SPMul.ia_out;
    int *ja_out = SPMul.ja_out;
    double *a_out = SPMul.a_out;
    int *rowid_out = SPMul.rowid_out;

    int myid = info.intervalid;
    int N = info.N;

    int i, j, from, to, id, rowid, colid;
    double aa;
    //获取盖尔圆信息——圆心以及半径
    /**************  *mycenter, *myradius  *********************/
    double *mycenter = calloc(N_loc, sizeof(double));
    double *myradius = calloc (N_loc, sizeof(double));
    for (i = 0; i < N_loc; i++)
    {
        rowid = fstrow + 1;
        from = ia_in[i] - 1;
        to = ia_in[i + 1] - 1;
        for (j = from; j < to; j++)
        {
            colid = ja_in[j];
            aa = a_in[j];
            if (rowid == colid)
                mycenter[i] = aa;
            else
                myradius[i] = myradius[i] + fabs(aa);
        }
    }

    for (i = 0; i < nnz_out; i++)
    {
        id = ia_out[i];
        rowid = fstrow + ia_out[i];
        j = ja_out[i];
        colid = rowid_out[j];
        aa = a_out[i];

        if (rowid == colid)
            mycenter[id] = aa;
        else
            myradius[id] = myradius[id] + fabs(aa);
    }

    //搜集各个进程的数据
    /*********************  *center, *radius  *****************************/
    double *center, *radius;
    if (myid == 0)
    {
        center = calloc(N, sizeof(double));
        radius = calloc(N, sizeof(double));
    }

    MPI_Gatherv_double(mycenter, N_loc, center, info.intervalcomm);
    MPI_Gatherv_double(myradius, N_loc, radius, info.intervalcomm);
    free(myradius);
    free(mycenter);

    //out
    (*Gerinfo).N = N;
    (*Gerinfo).center = center;
    (*Gerinfo).radius = radius;
}

/***************************************
统计盖尔圆信息，根据矩阵A，得到盖尔圆的断点信息以及其对应的概率密度变化量
******************************************/
void Ger_Sort(int N, double *A, double *density);
void Ger_getGerinfo(Vector *Gerdot, Density *Gerdensity, Gerinfo Gerinfo)
{
    int N = Gerinfo.N;
    double *center = Gerinfo.center;
    double *radius = Gerinfo.radius;

    int i, j;
    double aa;
    //统计盖尔圆内分布概率————节点，以及概率变化值
    /********************* Gerdot_n  *Gerdensity_node, *Gerdensity_density  *****************************/
    int Gerdot_n;
    double *Gerdensity_node, *Gerdensity_density;
    Gerdensity_node = calloc(2 * N, sizeof(double));
    Gerdensity_density = calloc(2 * N, sizeof(double));
    Gerdot_n = 0;
    for (i = 0; i < N; i++)
    {
        if (radius[i] == 0)
        {
            Gerdot_n = Gerdot_n + 1;
            Gerdensity_node[2 * i] = center[i] - radius[i];
            Gerdensity_node[2 * i + 1] = center[i] + radius[i];
        }
        else
        {
            Gerdensity_node[2 * i] = center[i] - radius[i];
            Gerdensity_node[2 * i + 1] = center[i] + radius[i];
            Gerdensity_density[2 * i] = 0.5 / radius[i];
            Gerdensity_density[2 * i + 1] = -0.5 / radius[i];
        }
    }

    for (i = 0; i < 2 * N; i++)
    {
        Gerdensity_density[i] = Gerdensity_density[i] / (N - Gerdot_n);
    }


    //统计盖尔点的信息，并且对其进行排序
    /********************* Gerdot_data  *****************************/
    double *Gerdot_data;
    if (Gerdot_n != 0)
    {
        Gerdot_data = calloc(Gerdot_n, sizeof(double));
        for (i = 0; i < N; i++)
        {
            if (radius[i] == 0)
            {
                Gerdot_data[i] = Gerdensity_node[2 * i];
            }
        }

        //排序
        for (i = 0; i < Gerdot_n; i++)
        {
            for (j = i + 1; j < Gerdot_n; j++)
            {
                if (*(Gerdot_data + i) > *(Gerdot_data + j))
                {
                    aa = *(Gerdot_data + j);
                    *(Gerdot_data + j) = *(Gerdot_data + i);
                    *(Gerdot_data + i) = aa;
                }
            }
        }

    }

    //根据分布节点进行排序
    Ger_Sort(2 * N, Gerdensity_node, Gerdensity_density);

    //out
    free(center);
    free(radius);

    (*Gerdot).n = Gerdot_n;
    (*Gerdot).data = Gerdot_data;

    (*Gerdensity).nnode = 2 * N;
    (*Gerdensity).node = Gerdensity_node;
    (*Gerdensity).density = Gerdensity_density;

}


/***************************************************
根据盖尔圆信息估算特征值的分布
****************************************************/
void Ger_getCdensity(Density *cdensity, Density Gerdensity, Vector  Gerdot)
{
    //输入
    int Gerdensity_nnode = Gerdensity.nnode;
    double *Gerdensity_node = Gerdensity.node;
    double *Gerdensity_density = Gerdensity.density;
    int Gerdot_n = Gerdot.n;
    double *Gerdot_data = Gerdot.data;

    int N = Gerdensity_nnode / 2;

    int i, id;
    //统计节点的个数，并分配空间
    /************************  cdensity_nnode *********************/
    int cdensity_nnode;
    cdensity_nnode = 1;
    for (i = 1; i < Gerdensity_nnode; i++)
    {
        if (Gerdensity_node[i] > Gerdensity_node[i - 1])
            cdensity_nnode++;
    }

    //统计节点的个数，并分配空间
    /************************ *cdensity_node *********************/
    double *cdensity_node, *density;
    cdensity_node = calloc(cdensity_nnode, sizeof(double));
    density = calloc(cdensity_nnode, sizeof(double));
    id = 0;
    density[id] = Gerdensity_density[0];
    cdensity_node[id] = Gerdensity_node[0];
    for (i = 1; i < Gerdensity_nnode; i++)
    {
        if (Gerdensity_node[i] > Gerdensity_node[i - 1])
        {
            id++;
            density[id] = Gerdensity_density[i];
            cdensity_node[id] = Gerdensity_node[i];
        }
        else
        {
            density[id] = density[id] + Gerdensity_density[i];
        }
    }
    if (id + 1 != cdensity_nnode)
        printf("error1 in getcdensity  id %d cdensity_nnode %d\n", id, cdensity_nnode);

    for (i = 1; i < cdensity_nnode; i++)
        density[i] = density[i] + density[i - 1];

    /************************** *cdensity_density ***********************/
    double *cdensity_density;
    cdensity_density = calloc(cdensity_nnode, sizeof(double));
    for (i = 1; i < cdensity_nnode; i++)
        cdensity_density[i] = density[i - 1] * (cdensity_node[i] - cdensity_node[i - 1]);
    if (Gerdot_n != 0)
    {
        int idnode = 0;
        int iddot = 0;
        while (iddot < Gerdot_n)
        {
            if (cdensity_node[idnode] == Gerdot_data[iddot])
            {
                iddot++;
                cdensity_density[idnode] = cdensity_density[idnode] + 1.0 / N;
            }
            else
            {
                idnode++;
            }
        }
        free(Gerdot_data);
    }

    for (i = 1; i < cdensity_nnode; i++)
    {
        cdensity_density[i] = cdensity_density[i] + cdensity_density[i - 1];
        if (cdensity_density[i] > 1)
            cdensity_density[i] = 1;
    }


    (*cdensity).nnode = cdensity_nnode;
    (*cdensity).node = cdensity_node;
    (*cdensity).density = cdensity_density;


    free(Gerdensity_node);
    free(Gerdensity_density);
    free(density);
}




/******************************************
特征值分布的重构，根据盖尔圆得到的特征值分布概率进行插值，得到估计的每一个特征值的位置
********************************************/
void Ger_getCeig(double *ceig, int N, Density cdensity)
{
    int nnode = cdensity.nnode;
    double *density_node = cdensity.node;
    double *density_density = cdensity.density;

    int i;

    int idceig, iddensity;
    double c, c1, c2, x1, x2;
    idceig = 0;
    iddensity = 1;
    c = 1.0 / N;
    while (idceig < N)
    {
        if (iddensity < nnode)
        {
            if (density_density[iddensity] < c )
            {
                iddensity++;
            }
            else
            {
                c1 = density_density[iddensity - 1];
                c2 = density_density[iddensity];

                x1 = density_node[iddensity - 1];
                x2 = density_node[iddensity];

                ceig[idceig] = x1 + (c - c1) * (x2 - x1) / (c2 - c1);

                idceig++;
                c = c + 1.0 / N;
            }
        }
        else
        {
            ceig[idceig] = density_node[nnode - 1];
            idceig++;
        }

    }

    free(density_node);
    free(density_density);
}


/*****************************************
根据最小的特征值修正估计的特征值
****************************************/
void Ger_ModifyCeig(double *ceig, int N, double mineig)
{
    //构造density
    double a = (mineig - ceig[N - 1]) / (ceig[0] - ceig[N - 1]);
    double b = mineig - a * ceig[0];
    int i;
    for (i = 0; i < N; i++)
        ceig[i] = a * ceig[i] + b;
}


/**************************************************************************************
**
**                                    little function
**
**
***************************************************************************************/


/********************************************************
排序程序
*******************************************************/

void Ger_MiniSort(int N, int block, double *A, double *WA, double *density, double *Wdensity)
{
    int begin, id, ida, idb, endida, endidb;

    id = 0;
    for (begin = 0; begin < N; begin = begin + 2 * block)
    {
        ida = begin;
        idb = ida + block;
        endida = ida + block;
        endidb = idb + block;
        if (idb > N) idb = N;
        if (endida > N) endida = N;
        if (endidb > N) endidb = N;

        while (ida < endida || idb < endidb)
        {
            if (ida < endida && idb < endidb)
            {
                if (A[ida] < A[idb])
                {
                    WA[id] = A[ida];
                    Wdensity[id] = density[ida];
                    ida++;
                }
                else
                {
                    WA[id] = A[idb];
                    Wdensity[id] = density[idb];
                    idb++;
                }
            }
            else
            {
                if (idb == endidb)
                {
                    WA[id] = A[ida];
                    Wdensity[id] = density[ida];
                    ida++;
                }
                else
                {
                    WA[id] = A[idb];
                    Wdensity[id] = density[idb];
                    idb++;
                }
            }

            id++;
        }

    }

    if (id != N)
        printf("error in minisort  id: %d\n ", id);

    if (endidb != N)
        printf("error minisort endidb: %d\n", endidb);

}
void Ger_Sort(int N, double *A, double *density)
{
    int i, block, tag;
    double *WA, *Wdensity;

    WA = calloc(N, sizeof(double));
    Wdensity = calloc(N, sizeof(double));

    tag = 0;
    for (block = 1; block < N; block = block * 2)
    {
        if ( tag == 1)
        {
            Ger_MiniSort(N, block, WA, A, Wdensity, density);
            tag = 0;
        }
        else
        {
            Ger_MiniSort(N, block, A, WA, density, Wdensity);
            tag = 1;
        }
    }

    if (tag == 1)
    {
        for (i = 0; i < N; i++)
        {
            A[i] = WA[i];
            density[i] = Wdensity[i];
        }
    }

    free(Wdensity);
    free(WA);
}

void MPI_Gatherv_double(double *sendbuf, int sendcount, double *recvbuf, MPI_Comm comm)
{
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int *recvcnts = calloc(np, sizeof(int));
    int *displs = calloc(np, sizeof(int));

    MPI_Allgather(&sendcount, 1, MPI_INT, recvcnts, 1, MPI_INT, comm);
    int i;
    for (i = 1; i < np; i++)
    {
        displs[i] = displs[i - 1] + recvcnts[i - 1];
    }

    MPI_Gatherv(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcnts, displs, MPI_DOUBLE, 0, comm);

    free(recvcnts);
    free(displs);
}