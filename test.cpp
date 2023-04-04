#include<iostream>
#include<cmath>
#include <ctime>
#include <cstdlib>
#include <stdlib.h>
#include "test.hpp"

using namespace std;

const float pi = atan(1) * 4;

int sum(int **x, int p, int q)
{
int sum=0;
for(int i=0; i<q ; i++)
{
//sum+=x[p*i+i];
sum+=x[p][i];
}
return sum;
}


double sinfun1(double q) //test function
{
double test=0;
if (q>0&&q<1)
{
test=sin(pi*q);
}
else
{
test=0;
}

return test;
}



double sinfunn(double * Q) //test function
{
double test=0;
if (Q[0]>0&&Q[0]<1&&Q[1]>0&&Q[1]<1&&Q[2]>0&&Q[2]<1)
{
test=sin(pi*Q[0])*sin(pi*Q[1])*sin(pi*Q[2]);
}
else
{
test=0;
}

return test;
}



double normal(double mu , double sigma, double x) //normal distribution
{
double n=0;
n=(1/(sqrt(2*pi)*sigma))*exp(-pow(x-mu,2)/(2*pow(sigma,2)));
return n;
}


double uniform(double a, double b) //uniform distribution
{
double u=0;
u=1/(b-a);
return u;
}



double lagrangeBasis(vector<double> X, double x, int i)
{
          int len=X.size();
          double l_i = 1;//拉格朗日基函数li(x)
	  if(i<0 || i >=len)
	  {
	    cout<<"This point does not exist"<<endl;
	  }
	  else if(len==1)
	  {
	    l_i=1;
	  }
          else
          {
	    for (int j = 0; j < len; j++) 
	    {
	         if (j != i) 
	         {
	         l_i *= (x - X[j]) / (X[i] - X[j]);
	         }
	    }
	  }
return l_i;
}



int C(int n, int m)  //combinatorial number
{
    int ans = 1;
    for(int i = 1;i<= n;i++)
    {
        ans *= i;
    }
    for(int i = 1;i <= m;i++)
    {
        ans /=i;
    }
    for(int i = 1;i<= n-m;i++) 
    {
        ans /= i;
    }
    return ans;
}





int m(int i) //total number of interpolation points
{
    if(i < 0)
        return 0;
    if(i == 1)
        return 1;
    return pow(2,i-1)+1;
}



vector<double> ccpoints(int n, double a, double b)//Clenshaw-Curtis points in Arbitrary interval [a,b], n is total numbers of points
{
    vector<double> l(n);
    if(n > 1)
        {
          for (int j=1; j<n+1; j++)
          {
            l[j-1]=(1-cos(pi*(j-1)/(n-1)))/2.0*(b-a)+a;
          }
          return l;
        }
    else if(n==1)
        {
          l[0]=(a+b)/2;
          return l;
        }
    else
        {
          //cout<<"cc point error"<<end;
          exit(0);
        }
}



vector<double> weight(double a, double b, vector<double> x) //return weight
{
//double q=0;
double step=0.000001;
int n=x.size();
vector<double> w(n,1);

for (int k=0 ; k<n ; k++)
{
    double S = 0;
    double value1=lagrangeBasis(x, a , k)*uniform(a,b);
    double value2;
    double trapezoidalArea;
    for (double i = a; i < b; i = i + step)
    {
        value2 = lagrangeBasis(x, i + step , k)*uniform(a,b);
        trapezoidalArea = (value1 + value2)*step / 2;
        value1 = value2;
        S += trapezoidalArea;
    }
    
    w[k]=S;
    //cout<<S<<" "<<endl;
}

return w;

}




int index(const int d, int N, int X, int Y) //return order of interpolation
{

int mu;

if (N>=d)
{
mu=N-d;
}
else 
{
cout<<"N must larger than d"<<endl;
return 0;
}


if(d==1)
{
int in=d+mu;
return in;
}


int n=pow((mu+2),d); //preallocate grid
int **ind=new int*[n];
for (int i=0; i<n; i++)
{
ind[i]=new int [d];
}

for(int i=0; i<n; i++)
{
  for (int j=0; j<d; j++)
  {
  ind[i][j]=0;
  }
}
//cout<<mu<<" "<<n<<endl;

for (int k = 0; k < d ; k++)
{
ind[0][k]=1;
}

int cont=1;
int j=d-2;
int i=0;


while (cont==1)
{
i+=1;

for (int k = 0; k < d ; k++)
{
ind[i][k]=ind[i-1][k];
}


   while (j<d-2 && ind[i][j]<=mu+j+1-sum(ind , i , j)) 
   {
    j=j+1;
    //cout<<"1"<<" i ="<<i<<" j ="<<j<<endl;
   }

   if (j==d-2)
   {
   ind[i][j]+=1;
   //cout<<"2"<<" i ="<<i<<" j ="<<j<<endl;
   }
   
   while (j>0 && ind[i][j]>mu+j+1-sum(ind , i , j))
   {
   ind[i][j]=1;
   ind[i][j-1]+=1;
   j-=1;
   //cout<<"3"<<" i ="<<i<<" j ="<<j<<endl;
   }
  
  if (j==0 && ind[i][j]>mu+1)
  {
  i-=1;
  cont=0;
  //cout<<"4"<<" i ="<<i<<" j ="<<j<<endl;
  }
   
   
}

int y[i+1];
for (int p=0 ; p<i+1 ; p++)
{
int s=0;
for(int q=0 ; q<d-1 ; q++)
{
s+=ind[p][q];
}
//cout<<s<<endl;
y[p]=s;
}


for (int k=0; k<i+1; k++)
{
ind[k][d-1]=mu+d-y[k];
}


int le[i+1][d];

for (int a=0; a<i+1; a++)
{
for (int b=0; b<d; b++)
{
le[a][b]=ind[a][b];
}
}

return le[X][Y];

/*
cout<<"N= "<<N<<endl;
cout<<"d= "<<d<<endl;

for (int a=0; a<i+1; a++)
{
for (int b=0; b<d; b++)
{
cout<<ind[a][b]<<" ";
}
cout<<endl;
}
*/

}





int index_num(const int d, int N)//return number of index
{

int mu;

if (N>=d)
{
mu=N-d;
}
else 
{
cout<<"N must larger than d"<<endl;
return 0;
}


if(d==1)
{
int in=d+mu;
return in;
}


int n=pow((mu+2),d); //preallocate grid
int **ind=new int*[n];
for (int i=0; i<n; i++)
{
ind[i]=new int [d];
}

for(int i=0; i<n; i++)
{
  for (int j=0; j<d; j++)
  {
  ind[i][j]=0;
  }
}
//cout<<mu<<" "<<n<<endl;

for (int k = 0; k < d ; k++)
{
ind[0][k]=1;
}

int cont=1;
int j=d-2;
int i=0;


while (cont==1)
{
i+=1;

for (int k = 0; k < d ; k++)
{
ind[i][k]=ind[i-1][k];
}


   while (j<d-2 && ind[i][j]<=mu+j+1-sum(ind , i , j)) 
   {
    j=j+1;
    //cout<<"1"<<" i ="<<i<<" j ="<<j<<endl;
   }

   if (j==d-2)
   {
   ind[i][j]+=1;
   //cout<<"2"<<" i ="<<i<<" j ="<<j<<endl;
   }
   
   while (j>0 && ind[i][j]>mu+j+1-sum(ind , i , j))
   {
   ind[i][j]=1;
   ind[i][j-1]+=1;
   j-=1;
   //cout<<"3"<<" i ="<<i<<" j ="<<j<<endl;
   }
  
  if (j==0 && ind[i][j]>mu+1)
  {
  i-=1;
  cont=0;
  //cout<<"4"<<" i ="<<i<<" j ="<<j<<endl;
  }
   
   
}



return i+1;


}



















