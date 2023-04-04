#include <iostream>
#include "test.hpp"
#include<vector>
#include<cmath>
#include <fstream>
#include <iomanip>

using namespace std;

int main()
{
double a=0;
double b=1;
const int d=2; //dimension
ofstream SCmean;
SCmean.open ("meandata-2D.txt");
ofstream SCvar;
SCvar.open ("vardata-2D.txt");
for (int q=2; q<=8; q++)
{
//int q; //depth of sparse grid interpolation
int wn=0;

//point set
vector<vector<double> > B(1, vector<double>(d, (a+b)/2)); //grid 
vector<vector<double> > WE(1, vector<double>(d, 1)); 
vector<double> F(1, 1);

if (q<d || d<=0)
{
cout<<"q or d error"<<endl;
}
else if(d==1)
{

int n=m(q); //numbers of nodes
vector<double> nod=ccpoints(n, a, b);
vector<double> w=weight(a, b, nod);
vector<double> sa(n,1);
double m=0;
cout<<"The function value are"<<" ";
for (int i=0; i<n; i++)
{
  sa[i]=sinfun1(nod[i]);
  cout<<sa[i]<<" ";
  m+=sa[i]*w[i];
}
cout<<endl;
cout<<"The mean is "<<m<<endl;

double v=0;
for (int i=0; i<n; i++)
{
  v+=pow(sa[i]-m,2)*w[i];
}



cout<<"The variance is "<<v<<endl;
cout<<"The weight are "<<" ";

for (int ii=0; ii<w.size(); ii++)
{
  cout<<w[ii]<<" ";
}

cout<<endl;
cout<<"The points are "<<" ";
for (int ii=0; ii<nod.size(); ii++)
{
  cout<<nod[ii]<<" ";
}
cout<<endl;



}
else
{
double SGm=0; //mean
double SGv=0; //variance
int pp=0;

for (int K=q-d+1; K<=q; K++) //solve mean
{ 
  double U=0;
  double coef=pow(-1, q-K)*C(d-1, q-K);
  if (K<d)
  {
    U=0;
  }
  else
  {
  int p=0;
  //cout<<"K="<<K<<" coef="<<coef<<endl;
  int n=index_num(d,K);
  int x[n][d]; //dimensions of univariate grids
  vector<vector<vector<double> > > A(n, vector<vector<double>>(d, vector<double>(1, 1))); //univariate grids
  vector<vector<vector<double> > > W(n, vector<vector<double>>(d, vector<double>(1, 1))); //univariate weight
  
  for (int i=0 ; i<n ; i++)
  {
    for (int j=0 ; j<d; j++)
    {
      x[i][j]=m(index( d, K, i, j));
      A[i][j]=ccpoints( x[i][j],  a,  b);
      W[i][j]=weight(a, b, A[i][j]); //The integration interval in each direction is the same
    }

      vector<int> multi_ind(d, 0);  //multi-index for each Cartesian product
      int k=d;
      int cont=1;       
      while (cont==1)
      {
        p+=1;
        vector<double> aa(d, 1);
        vector<double> ww(d, 1);
        for (int qq=0; qq<d; qq++)
        {
          aa[qq]=A[i][qq][multi_ind[qq]];
          ww[qq]=W[i][qq][multi_ind[qq]];
        }
        
        //cout<<"p= "<<p<<" ";p<<" "<<ww[0]<<" "<<ww[1]<<" "<<aa[0]<<" "<<aa[1]<<" "<<sinfunn(aa[0],aa[1])<<endl;
        // construct points set
          
          double *aaa = new double[d];
          for(int kk=0; kk<d; kk++)
          {
          aaa[kk]=aa[kk];
          }
          
        
          if (p==1&&wn==0)
          { 
            pp+=1;
            //cout<<pp<<endl;
            double we=1;         
            B[0]=aa;
            F[0]=sinfunn(aaa);
            WE[0]=ww;
            for(int ii=0; ii<d; ii++)
            {
              we*=ww[ii];
            }
            U+=F[0]*we;
            wn+=1;
            //cout<<U<<endl;
          }
          else
          { 
            int re=0;
            for (int ii=0; ii<B.size(); ii++) //find duplicate points
            {
              if (aa!=B[ii])
              {
               re=re+1;
              } 
              else if(aa==B[ii])
              {
               break;
              }
            }
         
            if (re==B.size()) //new point
            { 
              pp+=1;
              //cout<<pp<<endl;
              double we=1;
              B.push_back(aa);
              WE.push_back(ww);
              F.push_back(sinfunn(aaa));
              for(int ii=0; ii<d; ii++)
              {
                we*=ww[ii];
              }
              U+=F[re]*we;
              //cout<<U<<endl;
            } 
            else if(re!=B.size()) //old point
            {
              //cout<<endl;
              double we=1;
              for(int ii=0; ii<d; ii++)
              {
                we*=ww[ii];
              }
              U+=F[re]*we;
             // cout<<U<<endl;
            }
           
           }
     
        
        int upd_ind=1; //update multi-index
        while (upd_ind==1 && k>0)
        {
          if (multi_ind[k-1]< x[i][k-1]-1)
          {
            multi_ind[k-1]+=1;
            upd_ind=0;
            k=d;
          }
          else
          {
            multi_ind[k-1]=0;
            k-=1;
          }  
        }
      
        if(k==0)
        {
          cont=0;
        }
       
       delete[] aaa;
      }

  }
  
  }
 //cout<<U<<" "<<coef<<endl;
 SGm+=U*coef;

}

SCmean<<B.size()<<'\t'<<setprecision(15)<<SGm<<endl;

cout<<"The numbers of points= "<<B.size()<<endl;
cout<<"The mean is "<<SGm<<endl;


for (int K=q-d+1; K<=q; K++) //solve variance
{ 
  double V=0;
  double coef=pow(-1, q-K)*C(d-1, q-K);
  if (K<d)
  {
    V=0;
  }
  else
  {
  int p=0;
  //cout<<"K="<<K<<" coef="<<coef<<endl;
  int n=index_num(d,K);
  int x[n][d]; //dimensions of univariate grids
  vector<vector<vector<double> > > A(n, vector<vector<double>>(d, vector<double>(1, 1))); //univariate grids
  vector<vector<vector<double> > > W(n, vector<vector<double>>(d, vector<double>(1, 1))); //univariate weight
  
  for (int i=0 ; i<n ; i++)
  {

    for (int j=0 ; j<d; j++)
    {
      x[i][j]=m(index( d, K, i, j));
      A[i][j]=ccpoints( x[i][j],  a,  b);
      W[i][j]=weight(a, b, A[i][j]); //The integration interval in each direction is the same
    }
    
      vector<int> multi_ind(d, 0);  //multi-index for each Cartesian product
      int k=d;
      int cont=1;       
      while (cont==1)
      {
        p+=1;
        vector<double> aa(d, 1);
        vector<double> ww(d, 1);
        for (int qq=0; qq<d; qq++)
        {
          aa[qq]=A[i][qq][multi_ind[qq]];
          ww[qq]=W[i][qq][multi_ind[qq]];
        }
        
        //cout<<"p= "<<p<<" "<<ww[0]<<" "<<ww[1]<<" "<<aa[0]<<" "<<aa[1]<<" "<<sinfunn(aa[0],aa[1])<<endl;
        //construct points set
 
        int re=0;
        for (int ii=0; ii<B.size(); ii++) //find duplicate points
        {
              if (aa!=B[ii])
              {
               re=re+1;
              } 
              else if(aa==B[ii])
              {
               break;
              }
         }
         
        double we=1;
        for(int ii=0; ii<d; ii++)
        {
          we*=ww[ii];
        }
        
        V+=pow(F[re]-SGm,2)*we;
           
        int upd_ind=1; //update multi-index
        while (upd_ind==1 && k>0)
        {
          if (multi_ind[k-1]< x[i][k-1]-1)
          {
            multi_ind[k-1]+=1;
            upd_ind=0;
            k=d;
          }
          else
          {
            multi_ind[k-1]=0;
            k-=1;
          }  
        }
      
        if(k==0)
        {
          cont=0;
        }
    
      }

  }
 }
 //cout<<U<<" "<<coef<<endl;
 SGv+=V*coef;
}

SCvar<<B.size()<<'\t'<<setprecision(15)<<SGv<<endl;
cout<<"The variance is "<<SGv<<endl;

/////////////////


}



}

SCmean.close();
SCvar.close();


/*

for (int ii=0; ii<B.size(); ii++)
{
for (int jj=0; jj<B[ii].size(); jj++)
{
cout<<B[ii][jj]<<" ";
}
cout<<endl;
}

*/


}
