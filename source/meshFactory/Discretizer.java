package meshFactory;
import math.*;

import static java.lang.Math.*;
import java.util.Arrays;
public class Discretizer {

	private double[][] sortedBoundary;
	private int[][] sortedIndex;
	public double[][] blockBoundary;
	public double[][] minMeshLeft,minMeshRight,baseLeft,baseRight;
	public boolean[][] discLeft,discRight;
	public double[] X,Y,Z;
	public double minMesh,factor;
	int dim=3;
	public int nBlocks,decimal,numberOfRegions,numberOfElements,numberOfNodes;
	private int maxMeshInOneDirection;

	public Discretizer(Geometry mg){
		
		int nB=mg.blockBoundary[0].length;
		factor=1.0/mg.scaleFactor;
		this.nBlocks=mg.nBlocks;
		dim=nB/2;
		this.numberOfRegions=this.nBlocks;
		this.blockBoundary=new double[this.numberOfRegions][nB];
		this.minMeshLeft=new double[this.numberOfRegions][nB];
		this.minMeshRight=new double[this.numberOfRegions][nB];
		this.baseLeft=new double[this.numberOfRegions][nB];
		this.baseRight=new double[this.numberOfRegions][nB];
		this.discLeft=new boolean[this.numberOfRegions][nB];
		this.discRight=new boolean[this.numberOfRegions][nB];

		
		for(int i=0;i<this.nBlocks;i++)
			for(int j=0;j<nB;j++){
				this.blockBoundary[i][j]=mg.blockBoundary[i][j];
				this.minMeshLeft[i][j]=mg.minMeshLeft[i][j];
				this.minMeshRight[i][j]=mg.minMeshRight[i][j];
				this.baseLeft[i][j]=mg.baseLeft[i][j];
				this.baseRight[i][j]=mg.baseRight[i][j];
				this.discLeft[i][j]=mg.discLeft[i][j];
				this.discRight[i][j]=mg.discRight[i][j];
			}	

		sortBoundary();
		setXYZ();

	}


	private void setXYZ(){

		double maxDimension;
		maxDimension=max(this.blockBoundary[0][1]-this.blockBoundary[0][0],
				this.blockBoundary[0][3]-this.blockBoundary[0][2]);
		
		if(dim==3) maxDimension=max(max(this.blockBoundary[0][1]-this.blockBoundary[0][0],
				this.blockBoundary[0][3]-this.blockBoundary[0][2]),
				this.blockBoundary[0][5]-this.blockBoundary[0][4]);
		this.minMesh=maxDimension;
		for(int i=0;i<this.nBlocks;i++)
			for(int j=0;j<this.blockBoundary[0].length;j++){
				if(this.minMeshLeft[i][j]!=0 && this.minMeshLeft[i][j]<this.minMesh) this.minMesh=this.minMeshLeft[i][j];
				if(this.minMeshRight[i][j]!=0 && this.minMeshRight[i][j]<this.minMesh) this.minMesh=this.minMeshRight[i][j];
			}

		this.maxMeshInOneDirection=(int)(maxDimension/this.minMesh)+1000;

		double[][] XYZ=new double[this.maxMeshInOneDirection][this.dim];
		double[] v;

		double d1=1e10,d2=d1,b1,b2;
		int N= this.sortedBoundary.length;
	
		int blockLeft=0,blockRight=0;
		int[] nXYZ=new int[this.dim];
		int j1,j2;
		for(int j=0;j<this.dim;j++){

			for(int i=0;i<N-1;i++) {
				if(this.sortedIndex[i][j]%2==0)
					j1=2*j;
				else
					j1=2*j+1;
				if(this.sortedIndex[i+1][j]%2==0)
					j2=2*j;
				else
					j2=2*j+1;

				if(this.sortedBoundary[i][j]==this.sortedBoundary[i+1][j])
					XYZ[nXYZ[j]++][j]=this.sortedBoundary[i+1][j];
				else{
					if(i==-1){
						b1=0;
						blockLeft=0;
					
					}
					else{
						
						blockLeft=(this.sortedIndex[i][j])/2;
						
						
						d1=this.minMeshRight[blockLeft][j1];

						if(this.discRight[blockLeft][j1])
							b1=this.baseRight[blockLeft][j1];
						else
							b1=0;
					}		
					if(i==N-1){
						b2=0;
						blockRight=0;
					}
					else{
						blockRight=(this.sortedIndex[i+1][j])/2;
						

						d2=this.minMeshLeft[blockRight][j2];

						if(this.discLeft[blockRight][j2])
							b2=this.baseLeft[blockRight][j2];	
						else b2=0;
					}							 

					v=meshExp(this.sortedBoundary[i][j],this.sortedBoundary[i+1][j],d1,d2,b1,b2);
					for(int p=0;p<v.length;p++){
						XYZ[nXYZ[j]++][j]=v[p];
					}
				}
			}
			
	
		}
		
		double[] tempX=new double[nXYZ[0]],	tempY=new double[nXYZ[1]],tempZ=new double[1];
		if(dim==3)tempZ=new double[nXYZ[2]];
		tempX[0]=XYZ[0][0];
		tempY[0]=XYZ[0][1];
		if(dim==3)
		tempZ[0]=XYZ[0][2];
		int ix=1,iy=1,iz=1;

		//double epsilon=1e-6/pow(10,decimal);
		double epsilon=1e-8;
		for(int i=1;i<tempX.length;i++)
			if(abs(XYZ[i][0]-XYZ[i-1][0])>epsilon)
				tempX[ix++]=XYZ[i][0];


		this.X=new double[ix];
		this.X=Arrays.copyOf(tempX,ix);

		for(int i=1;i<tempY.length;i++)
			if(abs(XYZ[i][1]-XYZ[i-1][1])>epsilon)
				tempY[iy++]=XYZ[i][1];


		this.Y=new double[iy];
		this.Y=Arrays.copyOf(tempY,iy);

		for(int i=1;i<tempZ.length;i++)
			if(abs(XYZ[i][2]-XYZ[i-1][2])>epsilon)
				tempZ[iz++]=XYZ[i][2];
		
		

if(dim==3){
		this.Z=new double[iz];
		this.Z=Arrays.copyOf(tempZ,iz);
}
		
		for(int i=0;i<X.length;i++)
			X[i]*=factor;
		
		for(int i=0;i<Y.length;i++)
			Y[i]*=factor;
		
if(dim==3){
		for(int i=0;i<Z.length;i++)
			Z[i]*=factor;
		
}

		for(int i=0;i<this.numberOfRegions;i++)
			for(int j=0;j<2*dim;j++){
				this.blockBoundary[i][j]*=factor;
				
			}	
		
	}

	public void sortBoundary(){
		int I,J;
		I=2*this.blockBoundary.length;

		J=this.dim;
		this.sortedBoundary=new double[I][J];	
		this.sortedIndex=new int[I][J];				
		Vect v=new Vect(I);
		for(int j=0;j<J;j++){
			for(int i=0;i<I;i++){
				if(i%2==0)  v.el[i]=this.blockBoundary[i/2][2*j];
				else
					v.el[i]=this.blockBoundary[(i-1)/2][2*j+1];

			}
			
		
			int[] indice=v.bubble();

			for(int i=0;i<I;i++){				
				this.sortedBoundary[i][j]=v.el[i];
				this.sortedIndex[i][j]=indice[i];

			}
			

		
		}		

	}



	private static double[] meshExp(double a,double b, double d1,double d2,double base1,double base2){

		if(d1*d2==0)  throw new IllegalArgumentException("Initial mesh size is zero!");
		if(base1==1 || base2==1){
			if(base1==1 && base2==1){
				Vect v=new Vect().linspace(a, b, min(d1,d2));
				return v.el;
			}
			else if(base1==1 &&( base2>1 && d1<d2) || base2==0){
				Vect v=new Vect().linspace(a, b, d1);
				return v.el;
			}
			else if(base2==1 && (base1>1 && d2<d1) || base1==0){
				Vect v=new Vect().linspace(a, b, d2);
				return v.el;
			}
		}
		if(base1==0 && base2==0) {double[] meshed=new double[2];
		meshed[0]=a;
		meshed[1]=b;
		return meshed;
		}
		if(a==b){ double[] meshed=new double[1];
		meshed[0]=a;
		return meshed;
		}
		if((b-a)<=min(d1,d2)){ double[] meshed=new double[2];
		meshed[0]=a;
		meshed[1]=b;
		return meshed;
		}
		if((b-a)<=2*min(d1,d2)){ double[] meshed=new double[3];
		meshed[0]=a;
		meshed[1]=(a+b)/2;
		meshed[2]=b;
		return meshed;
		}
		if((b-a)<=3*min(d1,d2)){ double[] meshed=new double[4];
		meshed[0]=a;
		meshed[1]=a+(b-a)/3;
		meshed[2]=a+2*(b-a)/3;
		meshed[3]=b;
		return meshed;
		}

		int n1max=(int)(Math.abs(b-a)/d1);
		int n2max=(int)(Math.abs(b-a)/d2);
		int Nmax=2*(4+Math.max(n1max,n2max)); 

		if(Nmax>50000) Nmax=50000;
		double[] temp1=new double[Nmax]; 
		double[] temp2=new double[Nmax]; 

		int n1,n2;
		temp1[0]=a;
		n1=1;
		if(base1!=0){
			while(temp1[n1-1]+2*pow(base1,n1-1)*d1<b){
				temp1[n1]=temp1[n1-1]+pow(base1,n1-1)*d1;
				n1++;
			}

			temp1[n1++]=b;}

		temp2[0]=b;
		n2=1;
		if(base2!=0){
			while(temp2[n2-1]-2*pow(base2,n2-1)*d2>a){
				temp2[n2]=temp2[n2-1]-pow(base2,n2-1)*d2;
				n2++;
			}
			temp2[n2++]=a;
		}
		if(n2==1){
			double[] discret=new double[n1];

			discret=Arrays.copyOf(temp1,n1);

			return discret; 
		}
		else if(n1==1){
			double[] discret=new double[n2];
			for(int i=0;i<n2;i++)
				discret[i]=temp2[n2-1-i];
			return discret; 
		}


		if(d1<=d2){
			int k=0;
			boolean juncture=true;
			double[] discret1=new double[n1]; 	
			double[] discret2=new double[n2];
			discret1=Arrays.copyOf(temp1,n1);
			for(int i=0;i<n2;i++)
				discret2[i]=temp2[n2-1-i];

			double[] temp=new double[n1+n2];
			temp[k++]=discret1[0];
			temp[k++]=discret1[1];
			for(int i=2;i<n1;i++)
				if(discret1[i]-discret1[i-1]<=intervalAsc(discret2,discret1[i]))
					temp[k++]=discret1[i];

			int k1=k-1;

			for(int i=1;i<n2-1;i++){
				if(discret2[i]<=temp[k1]) continue;

				if(juncture){
					juncture=false;
					if(discret2[i]-temp[k1]<.9*(temp[k1]-temp[k1-1])){
						temp[k1]=temp[k1-1]+(discret2[i+1]-temp[k1-1])/3;	
						temp[k++]=temp[k1]+(discret2[i+1]-temp[k1-1])/3;	
						temp[k++]=discret2[i+1];
						i++;
					}
					else{

						temp[k++]=discret2[i];
					}
				}
				else
					temp[k++]=discret2[i];
			}


			temp[k++]=b;



			double[] discret=new double[k];
			discret=Arrays.copyOf(temp,k);

			return discret; 
		}
		else{
			int k=0;
			boolean juncture=true;
			double[] discret1=new double[n1]; 	
			double[] discret2=new double[n2];
			discret2=Arrays.copyOf(temp2,n2);
			for(int i=0;i<n1;i++)
				discret1[i]=temp1[n1-1-i];

			double[] temp=new double[n1+n2];
			temp[k++]=discret2[0];
			for(int i=1;i<n2;i++)			
				if(discret2[i-1]-discret2[i]<=intervalDesc(discret1,discret2[i]))
					temp[k++]=discret2[i];

			int k1=k-1;
			for(int i=1;i<n1-1;i++){
				if(discret1[i]>=temp[k1]) continue;
				if(juncture){
					juncture=false;
					if(i<n1-2)
						temp[k1]=(discret1[i]+temp[k1])/2;
					else{
						temp[k1]=(discret1[i]+temp[k1-1])/2;
						temp[k++]=discret1[i];

					}
				}

				else
					temp[k++]=discret1[i];
			}
			temp[k++]=a;

			double[] discret=new double[k];
			for(int i=0;i<k;i++)
				discret[i]=temp[k-1-i];

			return discret; 

		}

	}

	private static double intervalAsc(double[] A,double a){

		int N=A.length;
		double d=0;
		for(int i=0;i<N-1;i++)
			if(A[i]<=a && a<A[i+1]){
				d=A[i+1]-A[i];
				break;
			}
		return d;		
	}
	private static double intervalDesc(double[] A,double a){

		int N=A.length;
		double d=0;
		for(int i=0;i<N-1;i++)
			if(A[i]>=a && a>A[i+1]){
				d=A[i]-A[i+1];
				break;
			}
		return d;		
	}



}
