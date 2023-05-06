package math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import math.SpVect;
import static java.lang.Math.*;
import math.Vect;

public class BandMat  {

	public double[][] el;
	public int nRow,width,kL;
	public boolean lower=true;
	
	public BandMat(){}

	public BandMat(SpMat A){
		int dx=0;
		 lower=true;
		for(int i=0;i<A.nRow;i++){
			int m=A.row[i].nzLength-1;
			int n=abs(A.row[i].index[m]-i);
			if(n>dx)
				dx=n;	
			if(lower)
			{
				if(A.row[i].index[m]>i)
					lower=false;
			}
		}
		

	 
/*		 if(!lower){
			
			 width=dx+1;
			 this.kL=(dx)/2;

		 }
		 else{
			 
			 width=dx+1;
			 this.kL=dx;

		}*/
		 width=2*dx+1;
		 this.kL=dx;

		 this.nRow=A.nRow;
		el=new double[nRow][width];

		if(this.lower){
			for(int i=0;i<A.nRow;i++){
				int p=A.row[i].nzLength-1;
				for(int k=0;k<=p;k++){
					int m=(i-A.row[i].index[k]);
					el[i][kL-m]=A.row[i].el[k];
				}
					
			}
			
		}
		else{

			for(int i=0;i<A.nRow;i++){
			int p=A.row[i].nzLength-1;

			for(int k=0;k<=p;k++){
				int m=(A.row[i].index[k]-i)+kL;
				el[i][m]=A.row[i].el[k];
			}
		}
				
		}
		

	}
	
	

	public BandMat(int nRow,int kL){
		this.nRow=nRow;
		el=new double[nRow][kL];
	}


	public BandMat (Mat M,boolean sym){
		
		this.nRow=M.nRow;
		int kL=0;
		for(int i=0;i<nRow;i++){
			int ix=0;
			for(int j=0;j<=i;j++)
				if(M.el[i][j]!=0)
					ix++;
			if(ix>kL)
				kL=ix;
		}

		this.width=kL;
		
		if(lower)
		 this.kL=	this.width-1;
		else
			 this.kL=(this.width-1)/2;
		
		for(int i=0;i<nRow;i++){
			int ix=0;
			for(int j=0;j<=i;j++)
				if(M.el[i][j]!=0){
					el[i][kL-1-j]=M.el[i][j];
					ix++;
				}
		}
	}
	


	public BandMat deepCopy(){
		BandMat M=new BandMat(nRow,this.width);
		for(int i=0;i<nRow;i++)
			for(int k=0;k<width;k++){
			M.el[i][k]=this.el[i][k];
		}
		M.lower=this.lower;
		M.width=this.width;
		M.kL=this.kL;

		return M;
	}
	
	public BandMat eye(int L){
		BandMat M=new BandMat(L,1);
		for(int i=0;i<L;i++){
			M.el[i][0]=1;
		}

		M.width=1;
		M.kL=0;
		return M;
	}
	public Vect mul(Vect u){
		if(nRow!=u.length) throw new IllegalArgumentException("Dimensions do not agree.");

		if(this.lower)
			return smul(u);
		
		Vect v=new Vect(this.nRow);
		
		int dx=(width-1)/2;

		for(int i=0;i<nRow;i++){

			for(int k=0;k<this.width;k++){
				if(el[i][k]!=0) {
					int j=i-dx+k;
				v.el[i]+=el[i][k]*u.el[j];		
							}
			}
		}
		

			
	return v;
	}
	
	public Vect smul(Vect u){
		Vect v=new Vect(this.nRow);
		Vect w=new Vect(this.nRow);
		int dx=(width-1);
		
		for(int i=0;i<nRow;i++){
			for(int k=0;k<this.width;k++){
				int kx=dx-k;
				if(el[i][kx]!=0) {
					int j=i-k;
					v.el[i]+=el[i][kx]*u.el[j];	
					
					if(k!=0){
						w.el[j]+=el[i][kx]*u.el[i];
					}
					
				}
			}
		}

	return v.add(w);
	}
	
	public double get(int i, int j){

		if(abs(i-j)>kL) return 0;
		
		if(lower){
			if(j>i) return get(j,i);
			
			return el[i][j-i+kL];
		}
		else{
			
			return el[i][j-i+kL];
		}

	}
	
	public void setE(double a,int i, int j){
		
		
		if(lower){
			if(j>i) return;
			
			 el[i][j-i+kL]=a;
		}
		else{
			
		el[i][j-i+kL]=a;
		}

	}
	

	
public void lu(){
		


		
		for(int i=0;i<this.nRow-1;i++){
			
			
		
			double pw=this.get(i,i);
				double c;
				for(int j=i+1;j<min(nRow,i+kL+1);j++){
		
					double x=this.get(j,i);
					if(x==0) continue;
					c=x/pw;
					
					this.setE(c,j,i);
					for(int k=i+1;k<min(nRow,i+kL+1);k++){
						
						double y=this.get(i,k);
						if(y!=0){
							double p=this.get(j,k)-c*y;
							this.setE(p,j,k);
						}
					}
				}
			
		}

	
	}
	
	

	
	public void show(){

		this.matForm().show();
	}
	
	
	public Mat matForm(){
		int I,J;
		I=this.nRow;
		Mat A=new Mat(I,I);

		for(int i=0;i<I;i++){
			if(this.lower) J=i+1;
			else J=min(I,i+kL+1);
			for(int j=max(0,i-kL);j<J;j++)
				A.el[i][j]=get(i,j);
		}
		
		return A;
	}
	

	

}
