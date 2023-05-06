package math;

import static java.lang.Math.*;

import java.awt.Color;

import io.Loader;

import java.util.Arrays;
import java.util.Random;

import javax.swing.JFrame;

import org.math.plot.Plot2DPanel;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import Jama.QRDecomposition;

public class Mat {
	public double[][] el=null;
	public int nRow;
	public int nCol;
	public static void main(String[] Args){
		
		String file=System.getProperty("user.dir") + "\\mat.txt";
		
		double Jz=1e5;
		double Lx=2,Ly=2;
		double Area=Lx*Ly;
		double R=.1;
		
		int n=400;

		double[][] FJ=new double[n][n];
		Random r=new Random(1000);		
		
		for(int i=0;i<n;i++)
			for(int j=0;j<=i;j++){
		
			if(i==j)
				FJ[i][j]=0e0*abs(2*i)+i+10;
			else
			{
				//if(abs(i-j)==1){
				if(r.nextDouble()<.1){
			FJ[i][j]=(i+j+1)/(2*n);
		//FJ[j][i]=FJ[i][j];
			//FJ[j][i]=i+j+1;
				}
			}

		}


		      
		Mat K=new Mat(FJ);
		SpMat Ks=new SpMat(K);
		
		int L=100;
		Vect[] F=new Vect[L];
		for(int i=0;i<L;i++){
			F[i]=new Vect(n);
			
			for(int j=0;j<n;j++)
				F[i].el[j]=100*(Math.cos(i*10.0/L)+.3*Math.cos(i*20.0/L))*(sin(3.0/n*j)+.5*sin(8.0/n*j))*(1+.0*r.nextDouble());
		}
		
		SpMatSolver ms=new SpMatSolver();
		int p=4;
		if(p>1){
		Mat W=new Mat(n,p);
		for(int i=0;i<p;i++)
		{
			Vect x=ms.CG(Ks,F[i*4],1e-6,1000);
			
			W.setCol(x, i);
			
		}
		
		Mat C=W.transp().mul(W).times(1.0/p);
		
		 Eigen eg=new Eigen(C);
		
		 eg.lam.hshow();
		 
		 Mat Q=eg.V;
		 
	
	//	 Vect lam2=C.mul(Q).normCol();
		// lam2.hshow();
		
		 
		 Mat Phi=W.mul(Q);

Mat PhiT=Phi.transp();
		 
		Mat  Kr=PhiT.mul(K.mul(Phi));
	
	//	Kr.show();
	
		Vect[] Fr=new Vect[L];
		for(int i=0;i<L;i++){
			Fr[i]=PhiT.mul(F[i]);
		
		}
		
		MatSolver ss=new MatSolver();
		
		Vect outr=new Vect(L);
		for(int i=0;i<L;i++)
		{
		Vect xr=ss.CG(Kr,Fr[i],1e-6,1000);
		/*	Vect xr=new Vect(1);
			xr.el[0]=Fr[i].el[0]/Kr.el[0][0];*/
			Vect x=Phi.mul(xr);
			outr.el[i]=x.el[n/2];

		}
	util.plot(outr);
		

		outr.hshow();

		}
		 
		Vect out=new Vect(L);
		for(int i=0;i<L;i++)
		{
			Vect x=ms.CG(Ks,F[i],1e-6,1000);
			
			out.el[i]=x.el[n/2];
			
		}
		util.plot(out);
		

		out.hshow();
	
		
		 
		/*		Vect out=new Vect(L);
				for(int i=0;i<L;i++)
				{
					Vect x=ms.gaussel(Mr, F[i]);
					
					out.el[i]=x.el[n/2];
					
				}*/
		
	//	util.plot(out);
		
		//Mr.show();
	
		
/*		SpMat Ms=new SpMat(Mr);
		//SpMat Mt=Ms.reflect(300);
	//Ms.show();
		//Ms.lower=false;
		
		//Ms.plot();
		
	//	Ms=Ms.RCM();
		
		
		//Ms.showLower();

		//Ms.plot();
		BandMat B=new BandMat(Ms);
	//	B.show();
	//	B.lower=false;
	//	Mat B1=B.matForm();
		//B.lu();
		//B.lower=false;
		
		Mat B1=B.matForm();

		B1.lu();
		B1.show();
		
		Vect v=new Vect().linspace(1,n,n);
		Vect b=B.mul(v);
		
		B.lu();
				B.show();
				
				B.lower=false;
		Vect x=new MatSolver().solvelu(B1, b);

		Vect x2=new SpMatSolver().solvelu(B, b);
		x.hshow();
x2.hshow();
*/
	}



	public Mat(){};


	public  Mat(double[][] array){
		nRow=array.length;
		nCol=array[0].length;
		el=array;
	}

	public  Mat(int nRow, int nCol){
		this.nRow=nRow;
		this.nCol=nCol;
		this.el=new double[nRow][nCol];
	}

	public  Mat(int[] dim){
		this.nRow=dim[0];
		this.nCol=dim[1];
		this.el=new double[dim[0]][dim[1]];
	}

	public  Mat(int nCol){
		this.nRow=1;
		this.nCol=nCol;
		this.el=new double[nRow][nCol];
	}

	public  Mat (int I, int J,double a, double b){
		nRow=I;
		nCol=J;
		el=new double[I][J];
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				el[i][j]=a+(b-a)*Math.random();
	}

	public void set(double[][] A){
		this.el=A;
		this.nRow=A.length;;
		this.nCol=A[0].length;
	}

	public void set(int[][] A){

		nRow=A.length;;
		nCol=A[0].length;
		el=new double[nRow][nCol];
		for(int i=0;i<nRow;i++)
			for(int j=0;j<nCol;j++)
				el[i][j]=A[i][j];
	}	

	public void set(double[] v){
		int J;
		J=v.length;
		this.el[0]=v;
		this.nRow=1;
		this.nCol=J;
	}

	public int[] size(){
		int[] size=new int[2];
		size[0]=this.nRow;
		size[1]=this.nCol;
		return size;
	}

	public void sz(){
		util.show(this.size());
	}

	public Mat mat1D(double[] v){
		int J;
		J=v.length;
		Mat M=new Mat(1,J);

		M.el[0]=v;
		return M;
	}	

	public Mat deepCopy(){
		int I=this.nRow;
		int J=this.nCol;
		Mat A=new Mat(I,J);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				A.el[i][j]=this.el[i][j];
		return A;
	}	




	public void ones(int I,int J){
		el=new double[I][J];
		this.nRow=I;
		this.nCol=J;
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				el[i][j]=1;

	}

	public void eye(int I){
		this.nRow=I;
		this.nCol=I;
		el=new double[I][I];
		for(int i=0;i<I;i++)
			el[i][i]=1;	
	}
	public void eye(){
		if(this.nRow!=this.nCol) throw new IllegalArgumentException("Matrix is not square.");

		el=new double[nRow][nRow];
		for(int i=0;i<nRow;i++)
			el[i][i]=1;	
	}

	public Mat add(Mat A){
		int I1,J1,I2,J2;
		I1=this.nRow;
		J1=this.nCol;
		I2=A.nRow;
		J2=A.nCol;

		if(I1!=I2 || J1!=J2) throw new IllegalArgumentException("Array dimensions do not agree");

		Mat C=new Mat(I1,J1);

		for(int i=0;i<I1;i++)
			for(int j=0;j<J1;j++)
				C.el[i][j]=this.el[i][j]+A.el[i][j];
		return C;
	}

	public void addToDiag(Vect D){
		int I=this.nRow;
		int J=this.nCol;
		int L=D.length;
			if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
			if(I!=L) throw new IllegalArgumentException("Arrays dimensions do not agree.");

	for(int i=0;i<I;i++)
				this.el[i][i]+=D.el[i];
	}
	
	public Mat sub(Mat A){
		int I1,J1,I2,J2;
		I1=this.nRow;
		J1=this.nCol;
		I2=A.nRow;
		J2=A.nCol;

		if(I1!=I2 || J1!=J2) throw new IllegalArgumentException("Array dimensions do not agree");

		Mat C=new Mat(I1,J1);

		for(int i=0;i<I1;i++)
			for(int j=0;j<J1;j++)
				C.el[i][j]=this.el[i][j]-A.el[i][j];
		return C;
	}
	


	public Mat mul(Mat A){

		return mul(A,A.nCol);
	}

	public Mat mul(Mat A,int K){
		int I1,J1,I2,J2;
		I1=this.nRow;
		J1=this.nCol;
		I2=A.nRow;
		J2=A.nCol;
		if(J1!=I2 || K>J2) throw new IllegalArgumentException("Matrix dimensions do not agree");
		double s;
		Mat C=new Mat(I1,K);
		for(int i=0;i<I1;i++)
			for(int j=0;j<K;j++){
				s=0;
				for(int n=0;n<J1;n++)
				
					s=s+this.el[i][n]*A.el[n][j];
				C.el[i][j]=s;
			}
		return C;						
	}

	public Vect mul(Vect u){
		int I,J,L;
		I=this.nRow;
		J=this.nCol;
		L=u.length;

		if(J!=L) throw new IllegalArgumentException("Matrix dimensions do not agree");
		double s;
		Vect v=new Vect(I);
		for(int i=0;i<I;i++){		
			s=0;
			for(int n=0;n<L;n++)
				if(this.el[i][n]!=0 && u.el[n]!=0)
				s=s+this.el[i][n]*u.el[n];
			v.el[i]=s;
		}
		return v;						
	}
	
	public Complex[] mul(Complex[] u){
		int I,J,L;
		I=this.nRow;
		J=this.nCol;
		L=u.length;

		if(J!=L) throw new IllegalArgumentException("Matrix dimensions do not agree; "+J+" and "+L);
		double s,t;
		Complex[] v=new Complex[I];
		for(int i=0;i<I;i++){		
			s=0;
			t=0;
			for(int n=0;n<L;n++){
				if(this.el[i][n]!=0 && u[n].re!=0)
				s=s+this.el[i][n]*u[n].re;
				
				if(this.el[i][n]!=0 && u[n].im!=0)
					t=t+this.el[i][n]*u[n].im;
			}
			v[i]=new Complex(s,t);
		}
		return v;						
	}

	public Vect getColVect(int j){

		Vect v=new Vect(this.nRow);
		for(int i=0;i<this.nRow;i++){			
			v.el[i]=this.el[i][j];
		}
		return v;						
	}


	public Vect getColVect(int j,int i1,int i2){

		Vect v=new Vect(i2-i1+1);
		for(int i=i1;i<i2;i++){	
			if(i>=0 && i<this.nRow)
				v.el[i-i1]=this.el[i][j];
		}
		return v;						
	}

	public Vect rowVect(int i){


		return new Vect(this.el[i]);						
	}

	public void setRow(Vect v,int i){
		if(v.length!=this.nCol) throw new IllegalArgumentException("Array dimensions do not agree");
		this.el[i]=Arrays.copyOf(v.el,v.length);
	}

	public void setCol(Vect v,int j){

		if(v.length!=this.nRow) throw new IllegalArgumentException("Array dimensions do not agree");

		for(int i=0;i<this.nRow;i++){			
			this.el[i][j]=v.el[i];
		}
	}


	public void setCol(Vect v,int j,int i1,int i2){

		if(v.length!=i2-i1+1) throw new IllegalArgumentException("Array dimensions do not agree");

		for(int i=i1;i<i2;i++){	
			if(i>=0 && i<this.nRow)
				this.el[i][j]=v.el[i-i1];
		}
	}

	public void setDiag(Vect v){

		if(this.nCol!=this.nRow) throw new IllegalArgumentException("Matrix is not square.");

		for(int i=0;i<this.nRow;i++){			
			this.el[i][i]=v.el[i];
		}
	}

	public Mat times(double p){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		Mat C=new Mat(I,J);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				C.el[i][j]=p*this.el[i][j];
		return C;						
	}


	public Mat times(int p){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		Mat C=new Mat(I,J);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				C.el[i][j]=p*this.el[i][j];
		return C;						
	}

	public Mat ddom(double a,double b){
		int I;
		I=this.nRow;
		Mat C=this;
		for(int i=0;i<I;i++)
			C.el[i][i]=b+a*this.el[i][i];
		return C;						
	}

	public Mat ddom(int a,int b){
		int I;
		I=this.nRow;
		Mat C=this;
		for(int i=0;i<I;i++)
			C.el[i][i]=a+b*this.el[i][i];
		return C;						
	}

	public Mat flr(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		Mat Af=new Mat(I,J);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				Af.el[i][j]=floor(this.el[i][j]);
		return Af;						
	}

	public Mat matForm(Mat[][] MB){
		int I,J;
		I=MB.length;
		J=MB[0].length;
		int m,n;
		m=MB[0][0].nRow;
		n=MB[0][0].nCol;
		
		Mat M=new Mat(I*m,J*n);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				for(int k=0;k<m;k++)
					for(int p=0;p<n;p++)
				M.el[i*m+k][j*n+p]=MB[i][j].el[k][p];
		
		return M;
	}
	
	public Mat transp(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		Mat AT=new Mat(J,I);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				AT.el[j][i]=this.el[i][j];
		return AT;
	}

	public Mat inv(){
		if(this.nRow!=this.nCol) throw new IllegalArgumentException("Matrix is not square.");
		if(this.nRow==2) return this.inv2();
		else if(this.nRow==3) return this.inv3();
		double[][] array=util.copy(el);

		return new Mat(Inverse.invert(array));
	}


	public Mat invU(){

		int I=this.nRow;
		int J=this.nCol;
		if(this.nRow!=this.nCol) throw new IllegalArgumentException("Matrix is not square.");
		double elmax=abs(this.maxElement());
		for(int i=0;i<I;i++)
			for(int j=0;j<i;j++){
				if(Math.abs(el[i][j])/elmax>1e-16) throw new IllegalArgumentException("Matrix is not upper triangular.");
			}
		MatSolver solver=new MatSolver();
		Mat U=new Mat(size());
		U.eye();
		U=solver.backSubMat(this.aug(U));

		return U;

	}

	public Mat invL(){

		int I=this.nRow;
		int J=this.nCol;
		if(this.nRow!=this.nCol) throw new IllegalArgumentException("Matrix is not square.");
		double elmax=abs(this.maxElement());
		for(int i=0;i<I;i++)
			for(int j=i+1;j<I;j++){
				if(Math.abs(el[i][j])/elmax>1e-16) throw new IllegalArgumentException("Matrix is not lower triangular.");
			}

		MatSolver solver=new MatSolver();
		Mat U=new Mat(size());
		U.eye();
		U=solver.forwardSub(this.aug(U));

		return U;

	}

	public void transpVoid(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		if(I!=J) throw new IllegalArgumentException("Only square natrices cab be transposed.");

		double temp;
		for(int i=0;i<I;i++)
			for(int j=0;j<i;j++){
				temp=el[i][j];
				el[i][j]=el[j][i];
				el[j][i]=temp;

			}

	}


	public  Mat rand(double a, double b){
		Mat r=new Mat(this.nRow,this.nCol);
		for(int i=0;i<this.nRow;i++)
			for(int j=0;j<this.nCol;j++)
				r.el[i][j]=a+(b-a)*Math.random();
		return r;
	}	

	public  Mat rand(){
		Mat r=new Mat(this.nRow,this.nCol);
		for(int i=0;i<this.nRow;i++)
			for(int j=0;j<this.nCol;j++)
				r.el[i][j]=Math.random();
		return r;
	}

	public  void rand(int I,int J){
		this.nRow=I;
		this.nCol=J;
		this.el=new double[I][J];
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				this.el[i][j]=Math.random();
	}

	public  void symrand(int I){
		this.nRow=I;
		this.nCol=I;
		this.el=new double[I][I];
		for(int i=0;i<I;i++)
			for(int j=0;j<i+1;j++)
				this.el[i][j]=Math.random();
		for(int i=0;i<I;i++)
			for(int j=i+1;j<I;j++)
				this.el[i][j]=this.el[j][i];
	}

	public  void rband(int I,int J,int b){
		this.nRow=I;
		this.nCol=J;
		this.el=new double[I][J];
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				if(abs(i-j)<0.5*b)
					this.el[i][j]=1+Math.random();
	}

	public  void symSprand(int I,int J,int n){
		this.nRow=I;
		this.nCol=J;
		this.el=new double[I][J];
		for(int i=0;i<I;i++){
			this.el[i][i]=Math.random();
			for(int k=0;k<(n-1)/2;k++)
				this.el[i][(int)(J*random())]=Math.random();
		}
		for(int i=0;i<I;i++)
			for(int j=i+1;j<J;j++)
				this.el[i][j]=this.el[j][i];
	}


	public Mat diag(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		if(I!=J) throw new IllegalArgumentException("Only square natrices cab be diagonalized");
		Mat D=new Mat(I,J);
		for(int i=0;i<I;i++)
			D.el[i][i]=this.el[i][i];
		return D;
	}

	public Vect diagonal(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		if(I!=J) throw new IllegalArgumentException("Only square natrices cab be diagonalized");
		Vect D=new Vect(I);
		for(int i=0;i<I;i++)
			D.el[i]=this.el[i][i];
		return D;
	}

	public Mat lowerst(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		if(I!=J) throw new IllegalArgumentException("Only square natrices cab be diagonalized");
		Mat D=new Mat(I,I);
		for(int i=0;i<I;i++)
			for(int j=0;j<I;j++)
				if(j<=i)
					D.el[i][j]=this.el[i][j];
		return D;
	}

	public Mat chol()
	{
		double er=1e-6*this.absMat().maxElement();
		return chol(er);
	}
	
	public Mat chol(double symcr){
		/** Performs the Cholesky decomposiretuns  of the matrix and the lower triangular matrix L;
		 * 
		 */
		util.pr(symcr);
		double eps=symcr;
		int I,J;
		I=nRow;
		J=nCol;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				if(abs(el[i][j]-el[j][i])>eps)  throw new IllegalArgumentException("Matrix is not symmetric");
		Mat L=new Mat(I,I);
		Vect D=new Vect(I);

		double s;
		for(int i=0;i<I;i++){
			s=0;
			for(int k=0;k<i;k++)
				s=s+L.el[i][k]*L.el[i][k]*D.el[k];
			D.el[i]=el[i][i]-s;
		
			for(int j=i+1;j<I;j++){
				s=0;
				for(int k=0;k<i;k++)
					s=s+L.el[j][k]*L.el[i][k]*D.el[k];
				L.el[j][i]=(el[j][i]-s)/D.el[i];		
			}
		}
		
		for(int i=0;i<I;i++)
			for(int j=0;j<i+1;j++)
				if(i==j)
					L.el[i][i]=sqrt(D.el[i]);
				else
					L.el[i][j]=L.el[i][j]*sqrt(D.el[j]);
		return L;
	}
	


	public int negOnDiag(){
		/** Performs the LU decomposition the number of negative enteris on the diagonal;
		 * 
		 */
		Mat A=this.deepCopy();
		A.gaussUpper();

		int nneg=0;
		for(int i=0;i<A.nRow;i++)
			if(A.el[i][i]<0) nneg++;

		return nneg;
	}

	public Mat ichol(){
		/** Performs an incomplete Cholesky decomposiretuns  of the matrix and the lower triangular matrix L;
		 * 
		 */
		int I,J;
		I=this.nRow;
		J=this.nCol;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				if(el[i][j]!=el[j][i])  throw new IllegalArgumentException("Matrix is not symmetric");
		Mat L=new Mat(I,I);
		Vect D=new Vect(I);
		double s;
		for(int i=0;i<I;i++){
			for(int j=0;j<=i;j++){
				if(el[i][j]==0) continue;
				s=0;
				for(int k=0;k<j;k++)
					s=s+L.el[j][k]*L.el[j][k]*D.el[k];

				D.el[j]=el[j][j]-s;

				if(j<i && el[i][j]!=0){
					s=0;
					for(int k=0;k<j;k++)
						s+=L.el[j][k]*L.el[i][k]*D.el[k];

					L.el[i][j]=(el[i][j]-s)/D.el[j];		
				}
			}
		}


		for(int i=0;i<I;i++)
			for(int j=0;j<=i;j++)
				if(i==j)
					L.el[i][i]=sqrt(D.el[j]);
				else
					L.el[i][j]=L.el[i][j]*sqrt(D.el[j]);
		return L;

	}

	public void ichol2(){
		int I,J;
		I=nRow;
		J=nCol;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				if(el[i][j]!=el[j][i])  throw new IllegalArgumentException("Matrix is not symmetric");
		//Mat L=new Mat(I,I);
		Vect D=new Vect(I);

		double s;
		for(int i=0;i<I;i++){
			s=0;
			for(int k=0;k<i;k++)
				s=s+el[i][k]*el[i][k]*D.el[k];
			D.el[i]=el[i][i]-s;
			for(int j=i+1;j<I;j++)
				if(el[i][j]!=0){
					s=0;
					for(int k=0;k<i;k++)
						s=s+el[j][k]*el[i][k]*D.el[k];
					el[j][i]=(el[j][i]-s)/D.el[i];		
				}
		}
		for(int i=0;i<I;i++)
			for(int j=0;j<=i;j++)
				if(i==j)
					el[i][i]=sqrt(D.el[i]);
				else
					el[i][j]=el[i][j]*sqrt(D.el[j]);

	}
/*	public Mat upperx(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		if(I!=J) throw new IllegalArgumentException("Only square natrices cab be diagonalized");
		Mat D=new Mat(I,J);
		for(int i=0;i<I;i++)
			for(int j=0;j<J;j++)
				if(j>i)
					D.el[i][j]=this.el[i][j];
		return D;
	}*/



	public Mat offDiag(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		if(I!=J) throw new IllegalArgumentException("Only square natrices cab be diagonalized");
		Mat D=this.deepCopy();
		for(int i=0;i<I;i++)
			D.el[i][i]=0;
		return D;
	}

	public Vect scale(Vect b){

		Vect Ci=this.diagonal().inv();


		for(int i=0;i<this.nRow;i++)
			for(int j=0;j<this.nCol;j++){
				el[i][j]*=Ci.el[i];
			}

		for(int i=0;i<this.nRow;i++)
			el[i][i]*=Ci.el[i];

		return Ci;
	}

	public int low0(){
		int flag;
		for(int i=0;i<this.nRow-1;i++){
			flag=bestpiv(i);
			if(flag==0) return 0;
			eLow(i);
		}
		return 1;
	}

	public int low0(int m){
		int flag;
		for(int i=0;i<m-1;i++){
			flag=bestpiv(i);
			if(flag==0) return 0;
			eLow(i);
		}
		return 1;
	}
	
	
	public void lu(){
		
		
		for(int i=0;i<this.nRow-1;i++){
		
			
				double c;
				for(int j=i+1;j<this.nRow;j++){
					
					if(this.el[j][i]==0) continue;
					c=this.el[j][i]/this.el[i][i];
					
					this.el[j][i]=c;
					for(int k=i+1;k<this.nCol;k++)
						if(this.el[i][k]!=0)
						this.el[j][k]=this.el[j][k]-c*this.el[i][k];
				}
			
		}
	
	}

	public void gaussUpper(){
		for(int i=0;i<this.nRow-1;i++)
			eLow(i);

	}

	public int bestpiv(int i){
		int im;
		im=util.indpiv(this,i);
		if(this.el[im][im]==0) return 0;
		if(im!=i)
			rswap(i,im);	

		return 1;

	}

	public void normalizeColumns(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		Vect v;
		double vn;
		for(int j=0;j<nCol;j++){
			v=this.getColVect(j);
			vn=v.norm();
			if(vn>0)
				v=v.times(1.0/vn);
			this.setCol(v,j);
		}



	}

	public void normalizeColumns(boolean[] locked){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		Vect v;
		double vn;
		for(int j=0;j<nCol;j++){
			if(locked[j]) continue;
			v=this.getColVect(j);
			vn=v.norm();
			if(vn>0)
				v=v.times(1.0/vn);
			this.setCol(v,j);
		}



	}


	public Vect normRow(){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		double nr=0;
		Vect v=new Vect(nRow);
		Vect vr=new Vect(nCol);
		for(int i=0;i<I;i++){
			vr.set(el[i]);
			v.el[i]=vr.norm();
		}

		return v;
	}

	public Vect normCol(){
		int I,J;
		I=this.nRow;
		J=this.nCol;

		Vect v=new Vect(J);
		for(int j=0;j<J;j++){
			v.el[j]=this.getColVect(j).norm();
		}

		return v;
	}


	public double norm(){

		return this.normRow().norm();
	}

	public void show(){
		
		show("%6.6e");
	}

	public void show(String format){
		int I,J;
		I=this.nRow;
		J=this.nCol;
		for(int i=0;i<I;i++){
			for(int j=0;j<J;j++)
				System.out.format(format+" \t",el[i][j]);
			System.out.println();
		}
		System.out.println();
	}
	
	public void plotx(){
		
		for(int i=0;i<this.nRow;i++){
			for(int j=0;j<this.nCol;j++){
				if(this.el[i][j]==0)
					System.out.print(".");
				else
					System.out.print("x");
			}
			System.out.println();
	}
	}
	
	public  void plot(){
		
		 Plot2DPanel plot = new Plot2DPanel();
	
		 for(int i=0;i<nRow;i++){
	
			 int L=nCol;
			 Vect x=new Vect(L);
			 Vect y=new Vect(L);
			 for(int j=0;j<L;j++){
				
				 if(el[i][j]==1){
				 x.el[j]=j;
				 y.el[j]=i;
				 }
				
			 }

			 
			 plot.addScatterPlot("",Color.red, x.el, y.el);

			 
		 }

		  JFrame frame = new JFrame("a plot panel");
		   frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		  frame.setSize(800,800);
		  frame.setContentPane(plot);
		  frame.setVisible(true);
		}
	
	public void rswap(int m, int n){
		double[] temp=new double[this.nCol];
		temp=this.el[m];
		this.el[m]=this.el[n];
		this.el[n]=temp;
	}
	public void cswap(int m, int n){
		double temp;
		for(int i=0;i<this.nRow;i++){
			temp=this.el[i][m];
			this.el[i][m]=this.el[i][n];
			this.el[i][n]=temp;}
	}

	public void eLow(int p){
		double c;
		for(int i=p+1;i<this.nRow;i++){
			c=this.el[i][p]/this.el[p][p];
			this.el[i][p]=0;
			for(int j=p+1;j<this.nCol;j++)
				this.el[i][j]=this.el[i][j]-c*this.el[p][j];
		}
	}


	public void eLowband(int p,int b){

		double c;
		int I=this.nRow;
		int J=this.nCol;
		for(int i=p+1;i<I;i++){
			if(this.el[i][p]!=0){
				c=this.el[i][p]/this.el[p][p];
				this.el[i][p]=0;
				for(int j=p+1;j<min(J-1,p+b+1);j++)
					this.el[i][j]=this.el[i][j]-c*this.el[p][j];

				this.el[i][J-1]=this.el[i][J-1]-c*this.el[p][J-1];
			}
		}

	}

	public Mat aug(Mat A){
		int I1,J1,I2,J2;
		I1=this.nRow;
		J1=this.nCol;
		I2=A.nRow;
		J2=A.nCol;
		if(I1!=I2) throw new IllegalArgumentException("Matrix domensions do not agree");
		Mat C=new Mat(I1,J1+J2);
		for(int i=0;i<I1;i++)
			for(int j=0;j<J1+J2;j++)
				if(j<J1)
					C.el[i][j]=this.el[i][j];
				else
					C.el[i][j]=A.el[i][j-J1];

		return C;
	}
	public Mat aug(Vect A){
		int I1,J1,I2,J2;
		I1=this.nRow;
		J1=this.nCol;
		I2=A.length;
		J2=1;
		if(I1!=I2) throw new IllegalArgumentException("Matrix domensions do not agree");
		Mat C=new Mat(I1,J1+J2);
		for(int i=0;i<I1;i++)
			for(int j=0;j<J1+J2;j++)
				if(j<J1)
					C.el[i][j]=this.el[i][j];
				else
					C.el[i][j]=A.el[i];

		return C;
	}

	public Mat augv(Mat A){
		int I1,J1,I2,J2;
		I1=this.nRow;
		J1=this.nCol;
		I2=A.nRow;
		J2=A.nCol;
		if(J1!=J2) throw new IllegalArgumentException("Matrix domensions do not agree");
		Mat C=new Mat(I1+I2,J1);
		for(int i=0;i<I1;i++)
			for(int j=0;j<J1;j++)
				C.el[i][j]=this.el[i][j];
		
		for(int i=0;i<I2;i++)
			for(int j=0;j<J1;j++)
				C.el[i+I1][j]=A.el[i][j];
				

		return C;
	}
	

	public Mat absMat() {
		Mat M=deepCopy();
		for(int i=0;i<nRow;i++)
			for(int j=0;j<nCol;j++)
				M.el[i][j]=abs(el[i][j]);

		return M;
	}

	public double maxElement() {

		double max=0;
		double tmp;
		for(int i=0;i<nRow;i++)
			for(int j=0;j<nCol;j++){
				tmp=abs(el[i][j]);
				if(tmp>max) max=tmp;
			}

		return max;
	}

	public double det3Diag(){
		double det=0;



		return det;

	}

	public double determinant() {
		int I,J;
		I=this.nRow;
		J=this.nCol;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");

		double result;

		if(nRow == 3) {
			result = el[0][0] * el[1][1]*el[2][2]+ el[0][1] * el[1][2]*el[2][0]+el[0][2] * el[1][0]*el[2][1]
			                                                                                              -el[0][2] * el[1][1]*el[2][0]- el[0][0] * el[1][2]*el[2][1]-el[0][1] * el[1][0]*el[2][2];
			return result;
		}
		else if(nRow == 1) {
			result = el[0][0];
			return result;
		}

		else if(nRow == 2) {
			result = el[0][0] * el[1][1] - el[0][1] * el[1][0];
			return result;
		}

		else
		{
			int flag;
			Mat M=new Mat(nRow,nCol);
			for(int i=0;i<nRow;i++)
				for(int j=0;j<nCol;j++)
					M.el[i][j]=el[i][j];
			flag=M.low0(); 
			if(flag==0) return 0;
			result=1;
			for(int i=0;i<nRow;i++)
				result*=M.el[i][i];
			return result;
		}




	}

	
	public  Mat inv2(){

		if(nRow!=nCol || nRow!=2)throw new IllegalArgumentException("Matrix is not 2 by 2");

		Mat M=new Mat(2,2); 

		double det=el[0][0]*el[1][1]-el[1][0]*el[0][1];
		
		if(det==0)throw new IllegalArgumentException("Matrix is singular");
		double detRev=1./det;
		M.el[0][0]=el[1][1]*detRev;
		M.el[0][1]=-el[0][1]*detRev;
		M.el[1][0]=-el[1][0]*detRev;
		M.el[1][1]=el[0][0]*detRev;

		return M;
	}

	public  Mat inv3(){

		if(nRow!=nCol || nRow!=3)throw new IllegalArgumentException("Matrix is not 3 by 3");

		Mat M=new Mat(nRow,nRow); 

		double det=determinant();
		
		if(det==0)throw new IllegalArgumentException("Matrix is singular");

		M.el[0][0]=(el[1][1]*el[2][2]-el[1][2]*el[2][1])/det;
		M.el[0][1]=(el[0][2]*el[2][1]-el[0][1]*el[2][2])/det;
		M.el[0][2]=(el[0][1]*el[1][2]-el[0][2]*el[1][1])/det;

		M.el[1][0]=(el[1][2]*el[2][0]-el[1][0]*el[2][2])/det;
		M.el[1][1]=(el[0][0]*el[2][2]-el[0][2]*el[2][0])/det;
		M.el[1][2]=(el[0][2]*el[1][0]-el[0][0]*el[1][2])/det;

		M.el[2][0]=(el[1][0]*el[2][1]-el[1][1]*el[2][0])/det;
		M.el[2][1]=(el[0][1]*el[2][0]-el[0][0]*el[2][1])/det;
		M.el[2][2]=(el[0][0]*el[1][1]-el[0][1]*el[1][0])/det;

		return M;
	}


	public Mat QRD(Mat R) {
		QRDecomposition qd=new QRDecomposition(new Matrix(this.el));
		double[][] R1=qd.getR().getArray();
		for(int i=0;i<R.nRow;i++)
			R.el[i]=Arrays.copyOf(R1[i], R.nCol);
		return new Mat(qd.getQ().getArray());

	}

	public Mat QRD() {
		QRDecomposition qd=new QRDecomposition(new Matrix(this.el));

		return new Mat(qd.getQ().getArray());

	}



	public Mat QR(Mat R) {
		Mat Q=this.transp();
		double c;
		Vect v=new Vect(),e=new Vect(nRow),vr=new Vect(nRow);
		v=Q.getColVect(0);
		c=v.norm();
		v=v.times(1/c);

		R.el[0][0]=c;
		Q.setCol(v, 0);


		for(int i=1;i<nRow;i++){
			v=Q.getColVect(i);
			vr=v;
			for(int j=0;j<i;j++){
				e=Q.getColVect(j);
				c=v.dot(e);
				R.el[j][i]=c;
				v=v.sub(e.times(c));
			}

			v.normalize();
			R.el[i][i]=vr.dot(v);;	
			Q.setCol(v, i);

		}


		return Q.transp();

	}

	public Mat QR() {
		Mat Q=new Mat(this.size());
		double c;
		Vect v=new Vect(),e=new Vect(nRow);
		v=this.getColVect(0);
		c=v.norm();
		if(c!=0)
			v=v.times(1/c);

		Q.setCol(v,0);


		for(int i=1;i<nCol;i++){
			v=this.getColVect(i);
			for(int j=0;j<i;j++){
				e=Q.getColVect(j);
				c=v.dot(e);
				v=v.sub(e.times(c));
			}

			v.normalize();
			Q.setCol(v,i);

		}

		return Q;

	}


	public Mat QRtri() {
		Mat Q=new Mat(this.size());
		double c;
		Vect v=new Vect(nRow),e=new Vect(nRow);
		//v=this.getColVect(0,0,3);
		v=this.getColVect(0);

		//c=v.norm();
		c=4;
		if(c!=0)
			v=v.times(1/c);

		//  Q.setCol(v,0,0,3);
		Q.setCol(v,0);

		for(int i=1;i<nCol;i++){
			//v=this.getColVect(i,i-2,i+1);
			v=this.getColVect(i);
			//v.hshow();
			for(int j=0;j<i;j++){
				//e=Q.getColVect(j,j-2,j+1);
				e=Q.getColVect(j);
				//e.hshow();
				//v.hshow();
				c=dot(v,e,i-j);
				if(i-j<=2) 
					v=sub(v,e,c,i);
				//v=v.sub(e.times(c));
			}

			//	v.normalize();
			v=normalize(v,i);
			//	Q.setCol(v,i,i-2,i+1);
			Q.setCol(v,i);


		}

		return Q;

	}

	public double dot(Vect v1, Vect v2,int i){
		if (i>2) return 0;
		else
			return v1.dot(v2);

	}

	public Vect sub(Vect v1, Vect v2,double c,int i){
		Vect v=v1.deepCopy();
		for(int j=max(0,i-3);j<=i;j++)

			v.el[j]-=c*v2.el[j];
		return v;

	}
	public Vect normalize(Vect v1, int i){
		double s=0;
		Vect v=v1.deepCopy();
		for(int j=max(0,i-3);j<=i;j++)
			s+=pow(v1.el[j],2);
		s=sqrt(s);
		for(int j=max(0,i-3);j<=i;j++)
			v.el[j]*=1/s;
		return v;

	}

	public Mat QRHous() {
		Mat A=times(-1);
		Vect v;
		double sum=0;
		double vn;
		for(int i=0;i<nRow-1;i++){
			v=new Vect(nRow-i);
			for(int j=0;j<v.length;j++)
				v.el[j]=A.el[i+j][i];

			vn=v.norm();
			v.el[0]-=signum(v.el[0])*vn;

			vn=v.norm();
			if(vn==0) continue;
			v=v.times(1/vn);

			for(int k=i;k<nCol;k++){
				sum=0;
				for(int m=0;m<v.length;m++)
					sum+=v.el[m]*A.el[i+m][k];

				for(int j=i;j<nRow;j++)
					A.el[j][k]-=2*v.el[j-i]*sum;
			}
		}	

		return A;

	}

	public Mat QRHous(Mat R) {

		Mat A=deepCopy();
		Mat P=new Mat();
		Mat Q=new Mat();
		P.eye(nRow);
		Q.eye(nRow);
		for(int i=0;i<nRow;i++){
			P=A.hsHold(i);
			A=P.mul(A);

			Q=Q.mul(P);
		}	
		for(int i=0;i<nRow;i++)
			for(int j=0;j<nRow;j++)
				R.el[i][j]=A.el[i][j];

		return Q;

	}

	public Mat QRHous2(Mat R) {
		Mat A=deepCopy();
		Mat P=new Mat();
		Mat Q=new Mat();
		P.eye(nRow);
		Q.eye(nRow);
		for(int i=0;i<nRow-1;i++){
			P=A.hsHold(i);
			A.mulPart(P,i);

			//A=P.mul(A);

			//Q=Q.mul(P);
			//Q=Q.transp();
			//Q.mulPart(P, i);
		}	
		for(int i=0;i<nRow;i++)
			for(int j=0;j<nRow;j++)
				R.el[i][j]=A.el[i][j];

		MatSolver solver =new MatSolver();

		Mat S=new Mat(size());
		Q=solver.forwardSub(R.transp().aug(this.transp())).transp();

		return Q;


	}


	/*  private void QR3d () {

			Mat A=deepCopy();
			Mat P=new Mat();
			Mat Q=new Mat();
		   //  This is derived from the Algol procedures tql2, by
		   //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		   //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		   //  Fortran subroutine in EISPACK.

		      for (int i = 1; i < nRow; i++) {
		         e[i-1] = e[i];
		      }
		      e[n-1] = 0.0;

		      double f = 0.0;
		      double tst1 = 0.0;
		      double eps = Math.pow(2.0,-52.0);
		      for (int l = 0; l < n; l++) {

		         // Find small subdiagonal element

		         tst1 = Math.max(tst1,Math.abs(d[l]) + Math.abs(e[l]));
		         int m = l;
		         while (m < n) {
		            if (Math.abs(e[m]) <= eps*tst1) {
		               break;
		            }
		            m++;
		         }

		         // If m == l, d[l] is an eigenvalue,
		         // otherwise, iterate.

		         if (m > l) {
		            int iter = 0;
		            do {
		               iter = iter + 1;  // (Could check iteration count here.)

		               // Compute implicit shift

		               double g = d[l];
		               double p = (d[l+1] - g) / (2.0 * e[l]);
		               double r = Maths.hypot(p,1.0);
		               if (p < 0) {
		                  r = -r;
		               }
		               d[l] = e[l] / (p + r);
		               d[l+1] = e[l] * (p + r);
		               double dl1 = d[l+1];
		               double h = g - d[l];
		               for (int i = l+2; i < n; i++) {
		                  d[i] -= h;
		               }
		               f = f + h;

		               // Implicit QL transformation.

		               p = d[m];
		               double c = 1.0;
		               double c2 = c;
		               double c3 = c;
		               double el1 = e[l+1];
		               double s = 0.0;
		               double s2 = 0.0;
		               for (int i = m-1; i >= l; i--) {
		                  c3 = c2;
		                  c2 = c;
		                  s2 = s;
		                  g = c * e[i];
		                  h = c * p;
		                  r = Maths.hypot(p,e[i]);
		                  e[i+1] = s * r;
		                  s = e[i] / r;
		                  c = p / r;
		                  p = c * d[i] - s * g;
		                  d[i+1] = h + s * (c * g + s * d[i]);

		                  // Accumulate transformation.

		                  for (int k = 0; k < n; k++) {
		                     h = V[k][i+1];
		                     V[k][i+1] = s * V[k][i] + c * h;
		                     V[k][i] = c * V[k][i] - s * h;
		                  }
		               }
		               p = -s * s2 * c3 * el1 * e[l] / dl1;
		               e[l] = s * p;
		               d[l] = c * p;

		               // Check for convergence.

		            } while (Math.abs(e[l]) > eps*tst1);
		         }
		         d[l] = d[l] + f;
		         e[l] = 0.0;
		      }

		      // Sort eigenvalues and corresponding vectors.

		      for (int i = 0; i < n-1; i++) {
		         int k = i;
		         double p = d[i];
		         for (int j = i+1; j < n; j++) {
		            if (d[j] < p) {
		               k = j;
		               p = d[j];
		            }
		         }
		         if (k != i) {
		            d[k] = d[i];
		            d[i] = p;
		            for (int j = 0; j < n; j++) {
		               p = V[j][i];
		               V[j][i] = V[j][k];
		               V[j][k] = p;
		            }
		         }
		      }
		   }*/

	public Mat QRHous3Diag(Mat R) {

		Mat A=deepCopy();
		Mat P=new Mat();
		Mat Q=new Mat();
		P.eye(nRow);
		Q.eye(nRow);
		for(int i=0;i<nRow;i++){
			P=A.hsHoldOf3Diag(i);
			A.mul(P,i,3);
			Q=P.mul(Q);
		}	

		for(int i=0;i<nRow;i++)
			for(int j=0;j<nRow;j++)
				R.el[i][j]=A.el[i][j];

		return Q.transp();

	}

	public void mul(Mat P,int i,int n){
		if(i>=nRow-2) n=nRow-i;
		Mat C=new Mat(n,n);
		double sum;
		for(int j=0;j<n;j++)
			for(int k=0;k<n;k++){
				sum=0;
				for(int p=0;p<n;p++)
					sum+=P.el[i+j][i+p]*el[i+p][i+k];
				C.el[j][k]=sum;
			}

		for(int j=0;j<n;j++)
			for(int k=0;k<n;k++){
				el[i+j][i+k]=C.el[j][k];
			}

	}

	public void mulPart(Mat P,int i){

		Mat C=new Mat(nRow-i,nCol-i);

		for(int j=0;j<C.nRow;j++)
			for(int k=0;k<C.nCol;k++){
				for(int p=0;p<C.nRow;p++)
					C.el[j][k]+=P.el[i+j][i+p]*this.el[i+p][i+k];
			}
		for(int j=0;j<C.nRow;j++)
			for(int k=0;k<C.nCol;k++)
				this.el[i+j][i+k]=C.el[j][k];

	}


	public Mat hsHoldOf3Diag(int k) {

		Mat P=new Mat();
		P.eye(nRow);
		Vect v=new Vect(2);
		v.el[0]=el[k][k];
		if(k<nRow-1)
			v.el[1]=el[k+1][k];

		double vn=v.norm();
		v.el[0]-=vn;
		vn=v.norm();
		if(vn==0) return P;
		v=v.times(1/vn);
		P.el[k][k]=1-2*v.el[0]*v.el[0];

		if(k<nRow-1)
		{
			P.el[k+1][k]=-2*v.el[0]*v.el[1];
			P.el[k][k+1]=P.el[k+1][k];
			P.el[k+1][k+1]=1-2*v.el[1]*v.el[1];
		}

		return P;
	}

	public Mat hsHold(int k) {

		Mat P=new Mat();
		P.eye(nRow);
		Vect v=new Vect(nRow-k);
		for(int i=k;i<nRow;i++)
			v.el[i-k]=el[i][k];
		double vn=v.norm();
		if(vn==v.el[0]) return P;
		v.el[0]-=vn;
		v.normalize();

		for(int i=0;i<v.length;i++)
			for(int j=i;j<v.length;j++){
				P.el[i+k][j+k]=-2*v.el[i]*v.el[j];
				if(i!=j)
					P.el[j+k][i+k]=P.el[i+k][j+k];
			}

		for(int i=k;i<P.nRow;i++){
			P.el[i][i]=1+P.el[i][i];
		}
		return P;
	}

	public Mat hsHold4Trid(int k) {

		Mat P=new Mat();
		P.eye(nRow);
		Vect v=new Vect(nRow-k);
		double alpha=0;
		for(int i=k+1;i<nRow;i++)
			alpha+=pow(el[i][k],2);

		if(el[k+1][k]>0)
			alpha=-sqrt(alpha);
		else
			alpha=sqrt(alpha);
		double s=.5*(pow(alpha,2)-el[k+1][k]*alpha);
		if(s<=0) return P;
		double rInv=1.0/sqrt(s);

		v.el[1]=.5*(el[k+1][k]-alpha)*rInv;
		for(int i=k+2;i<nRow;i++)
			v.el[i-k]=.5*el[i][k]*rInv;
		v.normalize();
		for(int i=0;i<v.length;i++)
			for(int j=i;j<v.length;j++){
				P.el[i+k][j+k]=-2*v.el[i]*v.el[j];
				if(i!=j)
					P.el[j+k][i+k]=P.el[i+k][j+k];
			}

		for(int i=k;i<P.nRow;i++){
			P.el[i][i]=1+P.el[i][i];
		}
		return P;
	}

	public Mat tridiagHous() {
		Mat A=deepCopy();
		Mat P;
		for(int i=0;i<A.nRow-2;i++){
			P=A.hsHold4Trid(i);
			A=P.mul(A).mul(P);

		}
		return A;

	}


	public Mat Jacpq(int p, int q) {
		double teta=(el[q][q]-el[p][p])/(2*el[p][q]);
		double t=signum(teta)/(abs(teta)+sqrt(1+pow(teta,2)));
		double c=1.0/sqrt(pow(t,2)+1);
		double s=c*t;
		Mat P=new Mat();
		P.eye(nRow);
		P.el[p][p]=c;
		P.el[q][q]=c;
		P.el[p][q]=s;
		P.el[q][p]=-s;

		return P;

	}

	public Mat eigVectLanc(int m){
		Mat Q=new Mat(nRow,m);
		Mat V=new Mat(nCol,m);
		Mat T=tridiagLanc(m,V);
		Q=V.mul(T.eigVect());

		return Q;
	}

	public Vect eigLanc(int m){
		Mat T=tridiagLanc(m);
		Vect lam= T.eigVal();
		lam.quickSort();
		return lam;


	}



	public Mat tridiagLanc(int m,Mat V) {


		Vect v=new Vect(nRow);
		Vect r=new Vect(nRow);
		Vect vp=new Vect(nRow);
		double[] a=new double[m];
		double[] b=new double[m+1];
		r.rand();
		b[0]=r.norm();
		for(int i=0;i<m;i++){
			vp=v;
			v=r.times(1.0/b[i]);
			r=mul(v).sub(vp.times(b[i]));
			a[i]=r.dot(v);
			r=r.sub(v.times(a[i]));
			b[i+1]=r.norm();

			for(int j=0;j<nRow;j++)
				V.el[j][i]=v.el[j];

		}

		Mat P=new Mat(m,m);
		P.el[0][0]=a[0];
		P.el[0][1]=b[1];
		P.el[m-1][m-2]=b[m-1];
		P.el[m-1][m-1]=a[m-1];
		for(int i=1;i<m-1;i++){
			P.el[i][i-1]=b[i];
			P.el[i][i]=a[i];
			P.el[i][i+1]=b[i+1];
		}

		return P;
	}

	public Mat tridiagLanc(int m) {


		Vect v=new Vect(nRow);
		Vect r=new Vect(nRow);
		Vect vp=new Vect(nRow);
		double[] a=new double[m];
		double[] b=new double[m+1];
		r.rand();
		b[0]=r.norm();
		for(int i=0;i<m;i++){
			vp=v;
			v=r.times(1.0/b[i]);
			r=mul(v).sub(vp.times(b[i]));
			a[i]=r.dot(v);
			r=r.sub(v.times(a[i]));
			b[i+1]=r.norm();


		}

		Mat P=new Mat(m,m);
		P.el[0][0]=a[0];
		P.el[0][1]=b[1];
		P.el[m-1][m-2]=b[m-1];
		P.el[m-1][m-1]=a[m-1];
		for(int i=1;i<m-1;i++){
			P.el[i][i-1]=b[i];
			P.el[i][i]=a[i];
			P.el[i][i+1]=b[i+1];
		}

		return P;
	}


	public double eigMax(double errc) {

		Vect v=new Vect(nRow);
		v.rand();
		v.normalize();
		Vect vp;
		double err=1;
		int k=0;
		while(err>errc){
			k++;
			if(k>1000) {  break;}				
			vp=v;
			v=mul(v);
			v.normalize();
			err=v.sub(vp).max();

		}
		System.out.println(" iterantions: "+k+" error: "+  err);

		return mul(v).norm();

	}

	public Vect eigMax(Vect v,double errc) {


		v.normalize();
		Vect vp;
		double err=1;
		int k=0;
		while(err>errc){
			k++;
			if(k>1000) {  break;}				
			vp=v;
			v=mul(v);
			v.normalize();
			err=v.sub(vp).max();

		}
		System.out.println(" iterantions: "+k+" error: "+  err);

		return v;

	}
	public Vect  eigGen(Mat B,Mat Q){

		Mat R=new Mat(size());
		Vect D=B.eigValVect(R);
		Mat A=R.transp().mul(this).mul(R);
		Mat Qa=new Mat(size());
		Vect lam=A.eigValVect(Qa);
		Qa=R.transp().mul(Qa);
		for(int i=0;i<Q.nRow;i++)
			Q.el[i]=Arrays.copyOf(Qa.el[i],Q.nCol);
		return lam.times(D.inv());

	}


	public Vect  eigGenSym(Mat B,Mat Q){
		if(B.norm()==0) return new Vect(this.nRow);
		Mat L=B.chol();
		Mat invL=L.invL();
		Mat A=invL.mul(this).mul(invL.transp());		
		Mat Q1=new Mat(A.size());
		Vect D=A.eigValVect(Q1);
		Mat Q2=invL.transp().mul(Q1);
		//Q2.normalizeColumns();
		for(int i=0;i<Q.nRow;i++)
			Q.el[i]=Arrays.copyOf(Q2.el[i],Q.nCol);

		return D;

	}

	public Vect  eigGenSym(Mat B){

		Mat L=B.chol();
		Mat invL=L.invL();
		Mat A=invL.mul(this).mul(invL.transp());		
		Vect D=A.eigVal();	
		return D;

	}

	public Vect  eigValGen(Mat B){

		Mat R=new Mat(size());
		Vect D=B.eigValVect(R);
		Mat A=R.transp().mul(this).mul(R);
		Vect lam=A.eigVal();

		return lam.times(D.inv());

	}

	public Vect eigInv(double errc) {
		MatSolver solver=new MatSolver();
		Vect v=new Vect(nRow);
		v.rand();
		v.normalize();
		Vect vp;
		double err=1;
		int k=0;
		while(err>errc){
			k++;
			if(k>1000) {  break;}			
			vp=v;
			v=solver.gaussel(this,v);
			v.normalize();
			err=v.sub(vp).max();

		}
		System.out.println(" iterantions: "+k+" error: "+  err);

		return v;

	}

	public double eigMin(double errc) {
		MatSolver solver=new MatSolver();
		Vect v=new Vect(nRow);
		v.rand();
		v.normalize();
		Vect vp;
		double err=1;
		int k=0;
		while(err>errc){
			k++;
			if(k>1000) {  break;}			
			vp=v;
			v=solver.gaussel(this,v);
			v.normalize();
			err=v.sub(vp).max();

		}
		System.out.println(" iterantions: "+k+" error: "+  err);

		return mul(v).norm();

	}

	public Vect eigShift(double a, double errc) {
		Mat B=this;
		for(int i=0;i<nRow;i++)
			B.el[i][i]-=a;

		return B.eigInv(errc);

	}

	public Vect eigRayShift() {
		MatSolver solver=new MatSolver();
		Mat B=deepCopy();
		Vect v=new Vect(nRow);
		v.rand();
		v.normalize();
		double a=v.dot(this.mul(v));
		System.out.println(" a: "+a);
		double ap=0;
		double err=1,errc=1e-6;
		int k=0;
		while(err>errc){
			k++;
			if(k>500) {  break;}

			for(int i=0;i<nRow;i++)
				B.el[i][i]=el[i][i]-a;

			v=solver.gaussel(B,v);
			v.normalize();
			ap=a;
			a=v.dot(this.mul(v));

			err=abs(a-ap);

		}
		System.out.println(" iterantions: "+k+" error: "+  err);

		return v;

	}

	public Vect eigQR(double errc) {

		Mat A=this.tridiagHous();
		//Mat A=this.deepCopy();
		Mat Ap;
		Mat Q=new Mat();
		Mat R=new Mat(A.size());
		int k=0;
		double err=1;
		while(err>errc){
			k++;
			if(k>1000) {  break;}

			//	Q=A.QR(R);
			Q=A.QRHous3Diag(R);
			Ap=A.deepCopy();
			A=R.mul(Q);	

			err=A.sub(Ap).maxElement()/A.maxElement();
		}
		System.out.println(" iterantions: "+k+" error: "+  err);
		Vect v=R.diagonal();
		v.bubble();
		return v;

	}

	public Vect eigSimul(double errc) {
		Vect lam=new Vect(this.nRow);
		Vect lamr=new Vect(this.nRow);
		Mat Q=new Mat(this.size());
		Q=Q.rand();
		int k=0;
		double err=1;
		while(err>errc){
			k++;
			if(k>1000) {  break;}
			Q=Q.QR();
			Q=this.mul(Q);
			Q.normalizeColumns();
			if(k%10==0){
				lamr=lam;
				lam=this.mul(Q).normCol();
				err=lam.sub(lamr).abs().max()/lam.abs().max();
			}
		}
		System.out.println(" iterantions: "+k+" error: "+  err);
		lam.bubble();
		return lam;

	}

	public Mat eigVectQR(double errc) {
		Mat A=deepCopy();
		Mat Ap=new Mat(A.size());
		Mat Q=new Mat();
		Mat Q2=new Mat();
		Q2.eye(A.nRow);
		Mat R=new Mat(A.size());
		int k=0;
		int I=A.nRow;
		double errMax=1;
		while(errMax>errc){
			k++;
			if(k%10==9) Ap=A;
			if(k>1000) {  break;}
			Q=A.QRHous(R);
			Q2=Q2.mul(Q);
			if(k%10==0){
				errMax=A.sub(Ap).maxElement()/A.diagonal().abs().min();

			}
			A=R.mul(Q);	
		}
		System.out.println(" iterantions: "+k+" error: "+  errMax);
		return Q2;

	}


	public Vect eigQR3Diag(double errc) {
		Mat A=deepCopy();
		Mat Q=new Mat();
		Mat R=new Mat(A.size());
		int k=0;
		int I=A.nRow;
		double err=0,errMax=1,diag, diagMin=1e20;
		while(errMax>errc){
			k++;
			if(k>1) {  break;}
			Q=A.QRHous3Diag(R);
			if(k%10==0){
				errMax=0;
				for(int i=0;i<I;i++)
					for(int j=i;j<I;j++)
						if(i!=j){
							err=abs(R.el[i][j]);		
							if(err>errMax)
								errMax=err;
						}
						else{
							diag=abs(R.el[i][i]);
							if(diag<diagMin)
								diagMin=diag;

						}
				errMax=errMax/diagMin;

			}
			A=R.mul(Q);	

		}
		System.out.println(" iterantions: "+k+" error: "+  errMax);
		Vect v=R.diagonal();
		return v;

	}

	public Vect eigVal(){

		EigenvalueDecomposition evd=new EigenvalueDecomposition(new Matrix(util.copy(el)));
		Vect eval= new Vect(evd.getRealEigenvalues());
		eval.bubble();
		return eval;
	}

	public Mat eigVect(){
		EigenvalueDecomposition evd=new EigenvalueDecomposition(new Matrix(util.copy(el)));
		Vect eval= new Vect(evd.getRealEigenvalues());
		Mat Q=new Mat(evd.getV().getArrayCopy());

		int[] index=eval.bubble();

		Mat X= new Mat(Q.size());
		for(int i=0;i<eval.length;i++)
			X.setCol(Q.getColVect(index[i]), i);

		return X;
	}

	public Vect eigValVect(Mat Q){

		EigenvalueDecomposition evd=new EigenvalueDecomposition(new Matrix(util.copy(el)));
		Vect eval= new Vect(evd.getRealEigenvalues());
		Mat Q1=new Mat(evd.getV().getArrayCopy());


		int[] index=eval.bubble();
		for(int i=0;i<Q.nRow;i++)
			for(int j=0;j<Q.nCol;j++)
				Q.el[i][j]=	Q1.el[i][index[j]];
		for(int i=0;i<Q.nRow;i++)
			Q.el[i]=Arrays.copyOf(Q1.el[i], Q1.nCol);

		return eval;
	}



	public Vect eigSubspace(int p,double errMax){

		MatSolver solver=new MatSolver();
		int n=this.nRow;
		int q=min(min(p+8,2*p),n);

		Mat X=new Mat(n,q);

		Mat X_hat=new Mat(n,q);
		Vect x=new Vect(n);
		Vect errv=new Vect(q);
		Vect lamqr=new Vect(q);
		Vect lamq=new Vect(q);

		int[] index=new int[q];
		for(int i=0;i<q;i++)
			index[i]=i;


		X.setCol(new Vect().ones(n), 0);
		for(int i=1;i<q-1;i++)
			X.el[i][i]=1;
		x.rand();
		X.setCol(x, q-1);

		boolean[] locked=new boolean[q];

		Mat Q=new Mat(q,q);
		Mat I=new Mat(q,q);
		I.eye();
		Mat K_hat=new Mat();
		Mat M_hat=new Mat();
		int iter;
		int nconv=0;

		for(iter=0;iter<200;iter++){

			X_hat=solver.gaussel(this,X);	

			K_hat=X_hat.transp().mul(X);
			M_hat=X_hat.transp().mul(X_hat);		
			K_hat.eigGenSym(M_hat,Q);
			X=X_hat.mul(Q);

			lamqr=lamq.deepCopy();
			lamq=this.mul(X).normCol();
			index=lamq.bubble();

			errv=lamq.sub(lamqr).abs().div0(lamq);


			for( int i=0;i<p;i++)
				if(!locked[i] && errv.el[i]<errMax){

					locked[i]=true;
					nconv++;

				}

			//errv.hshow();

			if( nconv==p || solver.terminate) {  break;}
		}


		Vect lam=new Vect(p);		
		for(int i=0;i<p;i++)
			lam.el[i]=lamq.el[i];

		System.out.println("Subspace Iteration: " +iter+"   error: "+errv.el[p-1]);

		lam.quickSort();

		return lam;	

	}


	public Vect eigSubspace(Mat M,int p,double errMax){
		return  eigSubspace(M,new Mat(this.size()),p,errMax);
	}

	public Vect eigSubspace(Mat M,Mat P,int p,double errMax){

		MatSolver solver=new MatSolver();
		int n=this.nRow;
		int q=min(min(p+8,2*p),n);

		Mat X=new Mat(n,q);

		Mat X_hat=new Mat(n,q);
		Vect x=new Vect(n);
		Vect errv=new Vect(q);
		Vect lamqr=new Vect(q);
		Vect lamq=new Vect(q);

		int[] index=new int[q];
		for(int i=0;i<q;i++)
			index[i]=i;


		X.setCol(new Vect().ones(n), 0);
		for(int i=1;i<q-1;i++)
			X.el[i][i]=1;
		x.rand();
		X.setCol(x, q-1);

		boolean[] locked=new boolean[q];

		Mat Q=new Mat(q,q);
		Mat I=new Mat(q,q);
		I.eye();
		Mat K_hat=new Mat();
		Mat M_hat=new Mat();
		int iter;
		int nconv=0;

		for(iter=0;iter<200;iter++){
			X=M.mul(X);
			X_hat=solver.gaussel(this,X);	

			K_hat=X_hat.transp().mul(X);
			M_hat=X_hat.transp().mul(M.mul(X_hat));		
			K_hat.eigGenSym(M_hat,Q);
			X=X_hat.mul(Q);

			lamqr=lamq.deepCopy();
			lamq=this.mul(X).normCol().div0(M.mul(X).normCol());
			index=lamq.bubble();

			errv=lamq.sub(lamqr).abs().div0(lamq);


			for( int i=0;i<p;i++)
				if(!locked[i] && errv.el[i]<errMax){

					locked[i]=true;
					nconv++;

				}

			//errv.hshow();

			if( nconv==p || solver.terminate) {  break;}
		}


		Vect lam=new Vect(p);		
		for(int i=0;i<p;i++)
			lam.el[i]=lamq.el[i];

		System.out.println("Subspace Iteration: " +iter+"   error: "+errv.el[p-1]);

		lam.quickSort();

		return lam;	

	}
	
	public Mat lower(){
		int I,J;
		I=this.nRow;
		J=this.nCol;

		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");

		Mat C=new Mat(I,J);

		for(int i=0;i<I;i++)
			for(int j=0;j<=i;j++)
				C.el[i][j]=this.el[i][j];
		return C;
	}
	public Mat upper(){
		int I,J;
		I=this.nRow;
		J=this.nCol;

		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");

		Mat C=new Mat(I,J);

		for(int i=0;i<I;i++)
			for(int j=i;j<J;j++)
				C.el[i][j]=this.el[i][j];
		return C;
	}



}
