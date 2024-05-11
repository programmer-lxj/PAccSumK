#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//此代码是实现2008年 A parallel algorithm for accurate dot product中并行求和算法PSumK
//也可以进一步完善将剩余的静态数组全改为动态数组
//int L=3;

//double p[15]={99,100,2.3,10,21,1,2,3,4,5,1.2,244,312,40.4,58};

double *TwoSum(double a,double b)
{
    double x,z,y,*d;
    x=a+b;
    z=x-a;
    y=(a-(x-z))+(b-z);
	d=(double *)malloc(2*sizeof(double));
    d[0]=x;
    d[1]=y;
    return d;
}


double *VecSum(double *p,int lv)
{
    double *t,*q;
    int i;
	q=(double *)malloc(lv*sizeof(double));
    q[0]=p[0];
    for(i=1;i<lv;i++)//此处不应该是n
    {
        t=TwoSum(p[i],q[i-1]);
        q[i]=t[0];
        q[i-1]=t[1];
    }
    return q;
}

double *SumL(double *p,int L,int lp)
{
    double *temp,*q2;
    int k,i;
	q2=(double *)malloc(L*sizeof(double));
    for(k=1;k<L;k++)
    {
        temp=VecSum(p,lp-k+1);
        for(i=0;i<lp-k+1;i++)
        {
            p[i]=*(temp+i);
        }
        q2[k-1]=p[lp-k];
    }
    for(i=0;i<lp-L+1;i++)
    {
        q2[L-1]+=p[i];
    }
    return q2;
}

double SumK(double *p,int K,int lq)
{
    int i,k;
    double res,*t;
	t=(double *)malloc(2*sizeof(double));
    for(k=1;k<K;k++)
    {
        for(i=1;i<lq;i++)
		{
			t=TwoSum(p[i],p[i-1]);
			p[i]=t[0];
			p[i-1]=t[1];
		}
    }
    for(i=0;i<lq-1;i++)
    {
        res+=p[i];
    }
	res+=p[lq-1];
    return res;
}


double maxabs(double *p,int n)
{
    int i;
    double maxp;
    maxp=fabs(*(p+0));
    for(i=0;i<n;i++)
    {
        if(fabs(*(p+i))>maxp)
        {
            maxp=fabs(*(p+i));
        }
    }
    return maxp;
}

double NextPowerTwo(double p)
{
    double u,q,L;
    u=DBL_EPSILON/2;
    if(fabs(p)<((1-u)*u*DBL_MAX))
    {
        q=(1/u)*p;
        L=fabs((q-p)-q);
    }
    else
        L=2*(1/u)*NextPowerTwo(0.5*u*p);
        //L=2*pow(u,-1)*NextPowerTwo(0.5*u*p);
    return L;
}

double *Transform(double *p,double rho,int n)
{
    double *tf,u,mu,tau1,tau2,Phi,*q,tau=0,sigma,Ms,phi,t;
    int i;
    u=0.5*DBL_EPSILON;
    mu=maxabs(p,n);
    tf=(double *)malloc((n+3)*sizeof(double));
    if(n==0 || mu==0)
    {
        tau1=rho;
        tau2=0;
        sigma=0;
        tf[0]=tau1;
        tf[1]=tau2;
        tf[2]=sigma;
        for(i=3;i<n+3;i++)
        {
            tf[i]=0;
        }
        return tf;
    }
    Ms=NextPowerTwo(n+2);
	//printf("%lf  %lf\n",Ms,mu);
    sigma=Ms*NextPowerTwo(mu);
    phi=Ms*u;
    Phi=Ms*Ms*u;
    t=rho;
    q=(double *)malloc(n*sizeof(double));
    while(1)
    {
		tau=0;
        for(i=0;i<n;i++)
        {
            q[i]=(sigma+p[i])-sigma;
            tau+=q[i];
            p[i]=p[i]-q[i];
        }
        tau1=t+tau;
        if((fabs(tau1)>=Phi*sigma) || (sigma<=DBL_MIN))
        {
            tau2=tau-(tau1-t);
            tf[0]=tau1;
            tf[1]=tau2;
            tf[2]=sigma;
            for(i=3;i<n+3;i++)
            {
                tf[i]=p[i-3];
            }
            return tf;
        }
        t=tau1;
        if(t==0)
        {
            //应该用一个指针来接收
            tf=Transform(p,0,n);
            return tf;
        }
        sigma=phi*sigma;
    }
}

double *TransformK(double *p,double rho,int n)//好像sigma和Ms也是需要的
{
    //还要进行动态内存分配与释放
    double *tk,*tfk,s=0,res,R;
    int i;
    //matlab中if还没写，要弄清楚sigma和Ms的作用(不重要)
    tk=Transform(p,rho,n);
    for(i=3;i<n+3;i++)
    {
        s=*(tk+i)+s;
    }
    s=*(tk+1)+s;
    res=*(tk+0)+s;
    R=*(tk+1)-(res-(*(tk+0)));
    //给tfk赋值
    tfk=(double *)malloc((n+2)*sizeof(double));
    tfk[0]=res;
    tfk[1]=R;
    for(i=2;i<n+2;i++)
    {
        tfk[i]=*(tk+i+1);
    }
    return tfk;
}

double *AccSumK(double *p,int K,int n)
{
    double *temp,*tsk,R=0,*Res;
    int i,j,k;
    Res=(double *)malloc(K*sizeof(double));
    //temp=TransformK(p,R,n);//matlab里面这里有if-else，应该要实现，c语言好像没有可变参数，难道只能写2个函数？
    for(i=0;i<K;i++)
    {
        temp=TransformK(p,R,n);//res也就是每一次的第一个值都不同
        Res[i]=*(temp+0);
        R=*(temp+1);
        for(j=0;j<n;j++)
        {
            *(p+j)=*(temp+j+2);
        }
        if(fabs(Res[i])<=DBL_MIN)//不是Resk!
        {
            for(k=i+1;k<K;k++)
            {
                Res[k]=0;
            }
            return Res;
        }
        //temp=NULL;
    }
    return Res;
}



int main(int argc,char **argv)
{
	int n,c,c1,i,j,id,nprocs,M=6,K=7,*scounts;
	double local_int,total_int;
	double *RecPosX; //receive buffer
	int reminder; //reminder if it is not divisible
	double startwtime,endwtime;
    c=(n-1)/M+1;  //ÏòÉÏÈ¡Õû
    c1=n-c*(M-1);
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    int index1[c],index2[K];	
	double *a=NULL,*local_a=NULL,*local_q=NULL,*t,*total_q=NULL,*local_ti=NULL;
	int* sendcounts = NULL ;
    int* displs = NULL ;
/*	
    if(id==0)
    {
        scanf("%d",&n);
    }
*/
n=2002;
//    MPI_Bcast(&n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//startwtime=MPI_Wtime();
	c=(n-1)/nprocs+1;  //other cpu
    c1=n-c*(nprocs-1); //0 cpu
    
	a = (double *)malloc(n * sizeof(double));
	FILE *fpRead=fopen("20021e80e3.txt","r");
	if(fpRead==NULL)
	{
		return 0;
	}
    for(i = 0; i < n; i++)
    {
        fscanf(fpRead,"%lf ",&a[i]);
    }
	
	startwtime=MPI_Wtime();
	sendcounts = (int *)malloc(nprocs * sizeof(int)); //allocate memory for array storing number of elements
	//reminder = n%nprocs; //remineder
	
	//calculate number of elements need to be scattered in each process
	//目前最后一个进程分配的数据大于等于其他进程分配的数据
    for(i=1;i<nprocs;i++)
    {
        sendcounts [i] = c; 
    }
    //number of elements in the last process
    sendcounts [0] = c1;
	
	//calculate corresponding displacement  

    //displs = new int [ nProcs ];
	displs = (int *)malloc(nprocs * sizeof(int));
    for(i=0;i<nprocs;i++)
        {
            displs[i] = 0;
            for(j=0;j<i;j++)
            {
                displs[i] = displs[i] + sendcounts[j];
            }
        }
    //allocate the receive buffer   
    local_a = (double*)malloc(sendcounts[id]*sizeof(double));
    //now everything is ready and we can start MPI_Scatterv operation.  
    MPI_Scatterv(a, sendcounts , displs, MPI_DOUBLE, local_a, sendcounts[id], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	

	local_q=(double *)malloc(K*sizeof(double));
	local_q=SumL(local_a,K,sendcounts[id]);//L=K;

    free(a);
    free(local_a);
//	free(t);
	total_q=(double *)malloc((K*nprocs)*sizeof(double));
	MPI_Gather(local_q,K,MPI_DOUBLE,total_q,K,MPI_DOUBLE,0,MPI_COMM_WORLD); 
//	MPI_Gather(q2,K,MPI_DOUBLE,total_q,K,MPI_DOUBLE,0,MPI_COMM_WORLD);  //显示q2没有定义，明天再看看！
	free(local_q);
//if(id==0)
//{
    //total_int=SumK(total_q,K,K*nprocs);
	local_ti=(double *)malloc(K*sizeof(double));
	local_ti=AccSumK(total_q,K,K*nprocs);
//}
	for(i=0;i<K;i++)
	{
		total_int+=local_ti[i];
	}
	free(local_ti);
	free(total_q);
	free(sendcounts);
	free(displs);
	if(id==0)
    {
        printf("%d number thread, every part length is %d, result is %.100e\n",nprocs,n/nprocs,total_int);
		endwtime=MPI_Wtime();
		printf("wall clock time = %lf\n", endwtime-startwtime);
		fflush(stdout);
    }
	MPI_Finalize();
	return 0;
}
