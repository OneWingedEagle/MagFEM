package math;
import static java.lang.Math.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import Jama.QRDecomposition;
public class SpMatSolver {
	public boolean terminate;
	public int totalIter;
	public List<Double> errs = new ArrayList<Double>();
	public double resRef=0;


	public SpMatSolver(){}

	public static void main2(String[] args){
		
		int N=100;
		
		int Nm=1000;
		Mat A=new Mat();
/*		A.symSprand(N, N, 1);
		A.ddom(1, 1);

		
		SpMat Ms=new SpMat(A);*/
		SpMat Ms=new SpMat(N);
		Ms.lower=true;

	//Ms.show();
		for(int i=0;i<N;i++){
			Vect v=new Vect().rand(N, 0, 1);
			for(int j=0;j<N;j++)
				if(j>i)v.el[j]=0;
				else
				if(v.el[j]<.96) v.el[j]=0;
			v.el[i]=(1+Math.random())*2;
			Ms.row[i]=new SpVect(v);
		}
		//Ms.show();
		int M=2;
		Mat B=new Mat(N,M);
		for(int i=0;i<M;i++){
			Vect v=new Vect().rand(N, 0, 1);
			//Vect v=new Vect(N);
			//v.el[i]=1000;;
			Vect v2=Ms.smul(v);
			B.setCol(v2, i);
		}
		
		QRDecomposition qr=new QRDecomposition(new Matrix(B.el));
		Mat Q=new Mat(qr.getQ().getArray());
		B=Q.deepCopy();
		/*Mat R=new Mat(qr.getR().getArray());
		Q.show();
		R.show();*/
		
		
		Mat X=new Mat(N,M);
	
		
		SpMat Ls=Ms.ichol(1.05);
	//	Ls.show();
		
		
		SpMatSolver solver=new SpMatSolver();
		double t1=System.currentTimeMillis();
		for(int i=0;i<M;i++){
			Vect x=solver.ICCG(Ms,Ls,  B.getColVect(i), 1e-6, Nm);
			//Vect x=solver.CG(Ms,  B.getColVect(i), 1e-6, Nm);
			
		}
		double t2=System.currentTimeMillis();
		
		//Ms.show();
		//System.out.println("Main method is empty");
		
	  X=solver.blockICCG(Ms, Ls, B, 1e-6, Nm);
	//X=solver.blockCG(Ms, B, 1e-6, Nm);
	
	double t3=System.currentTimeMillis();

		System.out.println(t2-t1);
		System.out.println(t3-t2);
	//	util.pr(Ms.smul(X).sub(B).norm()/B.norm());	
		
		
	}

	
	public Vect forwardSub(SpMat L, Vect b,SpMat A){

		int I=L.getnRow();
		int J=L.getnCol();
		int K=b.length;
		if(I!=K) throw new IllegalArgumentException("Array dimensions do not agree.");
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		Vect x=new Vect(I);
		x.el[0]=b.el[0]/L.row[0].el[0];
		int j;
		for(int i=1;i<I;i++){
			double s=0;
			for(j=0;j<L.row[i].nzLength;j++){
			if(A.row[i].index[j]==i) break;
				s+=L.row[i].el[j]*x.el[A.row[i].index[j]];
			}
			
			if(j<L.row[i].nzLength)
			x.el[i]=(b.el[i]-s)/L.row[i].el[j];
			
			}
	
		return x;
	}
	

	public Vect forwardSub(SpMat L, Vect b){

		int I=L.getnRow();
		int J=L.getnCol();
		int K=b.length;
		if(I!=K) throw new IllegalArgumentException("Array dimensions do not agree.");
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		Vect x=new Vect(I);
		x.el[0]=b.el[0]/L.row[0].el[0];
		int j;
		for(int i=1;i<I;i++){
			double s=0;
			for(j=0;j<L.row[i].nzLength;j++){
			if(L.row[i].index[j]==i) break;
				s+=L.row[i].el[j]*x.el[L.row[i].index[j]];
			}
			
			if(j<L.row[i].nzLength)
			x.el[i]=(b.el[i]-s)/L.row[i].el[j];
			
			}
	
	
		return x;
	}
	
	public Vect backSub(SpMat L, Vect b,SpMat A){
		int I=L.getnRow();
		int J=L.getnCol();
		int K=b.length;
		if(I!=K) throw new IllegalArgumentException("Array dimensions do not agree.");
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		double[] su=new double[I];
		Vect x=new Vect(I);
		x.el[I-1]=b.el[I-1]/L.row[I-1].el[L.row[I-1].nzLength-1];
		
		for(int i=I-2;i>=0;i--){
			for(int k=0;k<A.row[i+1].nzLength;k++)
				su[A.row[i+1].index[k]]+=L.row[i+1].el[k]*x.el[i+1];		
				x.el[i]=(b.el[i]-su[i])/L.row[i].el[L.row[i].nzLength-1];
			}
	
		return x;
	}
	
	//============
	
	
	public Vect forwardSub(SpMatAsym L, Vect b,SpMatAsym A){

		int I=L.getnRow();
		int J=L.getnCol();
		int K=b.length;
		if(I!=K) throw new IllegalArgumentException("Array dimensions do not agree.");
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		Vect x=new Vect(I);
		x.el[0]=b.el[0]/L.row[0].el[0];
		int j;
		for(int i=1;i<I;i++){
			double s=0;
			for(j=0;j<L.row[i].nzLength;j++){
			if(A.row[i].index[j]==i) break;
				s+=L.row[i].el[j]*x.el[A.row[i].index[j]];
			}
			
			if(j<L.row[i].nzLength)
			x.el[i]=(b.el[i]-s)/L.row[i].el[j];
			
			}
	
		return x;
	}
	

	public Vect forwardSub(SpMatAsym L, Vect b){

		int I=L.getnRow();
		int J=L.getnCol();
		int K=b.length;
		if(I!=K) throw new IllegalArgumentException("Array dimensions do not agree.");
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		Vect x=new Vect(I);
		x.el[0]=b.el[0]/L.row[0].el[0];
		int j;
		for(int i=1;i<I;i++){
			double s=0;
			for(j=0;j<L.row[i].nzLength;j++){
			if(L.row[i].index[j]==i) break;
				s+=L.row[i].el[j]*x.el[L.row[i].index[j]];
			}
			
			if(j<L.row[i].nzLength)
			x.el[i]=(b.el[i]-s)/L.row[i].el[j];
			
			}
	
		return x;
	}
	
	public Vect backSub(SpMatAsym L, Vect b,SpMatAsym A){
		int I=L.getnRow();
		int J=L.getnCol();
		int K=b.length;
		if(I!=K) throw new IllegalArgumentException("Array dimensions do not agree.");
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		double[] su=new double[I];
		Vect x=new Vect(I);
		x.el[I-1]=b.el[I-1]/L.row[I-1].el[L.row[I-1].nzLength-1];
		for(int i=I-2;i>=0;i--){
			for(int k=0;k<A.row[i+1].nzLength;k++)
				su[A.row[i+1].index[k]]+=L.row[i+1].el[k]*x.el[i+1];		
				x.el[i]=(b.el[i]-su[i])/L.row[i].el[L.row[i].nzLength-1];
			}
	
		return x;
	}
	
	public Vect solveTriangular(SpMat L, Vect b, SpMat A){
		
		Vect y=forwardSub(L,b,A);
		
		return backSub(L,y,A);

	}
	
	public Vect solveTriangular(SpMatAsym L, Vect b, SpMatAsym A){
		
		Vect y=forwardSub(L,b,A);
		
		return backSub(L,y,A);

	}
	public VectComp solveTriangular(SpMatComp L, VectComp b, SpMatComp A,boolean hermit){
		
		VectComp y=forwardSub(L,b,A);

		VectComp x=backSub(L,y,A,hermit);

		return x;
		

	}
	
	public VectComp forwardSub(SpMatComp L, VectComp b,SpMatComp A){

		int I=L.getnRow();
		int J=L.getnCol();
		int K=b.length;
		if(I!=K) throw new IllegalArgumentException("Array dimensions do not agree.");
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		VectComp x=new VectComp(I);
		x.el[0]=b.el[0].times(L.row[0].el[0].inv());
		int j;
		for(int i=1;i<I;i++){
			Complex s=new Complex();
			for(j=0;j<L.row[i].nzLength;j++){
			if(A.row[i].index[j]==i) break;
				s=s.add(L.row[i].el[j].times(x.el[A.row[i].index[j]]));
			}
			
			if(j<L.row[i].nzLength)
			x.el[i]=b.el[i].sub(s).times(L.row[i].el[j].inv());
			
			}
		
		return x;
	}
	
	public VectComp backSub(SpMatComp L, VectComp b,SpMatComp A,boolean hermit){

		int I=L.getnRow();
		int J=L.getnCol();
		int K=b.length;
		if(I!=K) throw new IllegalArgumentException("Array dimensions do not agree.");
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		Complex[] su=new Complex[I];
		
		for(int i=0;i<I;i++)
			su[i]=new Complex();
		
		VectComp x=new VectComp(I);
		x.el[I-1]=b.el[I-1].times(L.row[I-1].el[L.row[I-1].nzLength-1].inv());
	
		if(hermit){

			for(int i=I-2;i>=0;i--){
			for(int k=0;k<A.row[i+1].nzLength;k++){
				int j=A.row[i+1].index[k];
				su[j]=su[j].add(L.row[i+1].el[k].conj().times(x.el[i+1]));	
				}


			x.el[i]=b.el[i].sub(su[i]).times(L.row[i].el[L.row[i].nzLength-1].inv());
			}
		}
		else{

			for(int i=I-2;i>=0;i--){
			for(int k=0;k<A.row[i+1].nzLength;k++)
				su[A.row[i+1].index[k]]=su[A.row[i+1].index[k]].add(L.row[i+1].el[k].times(x.el[i+1]));		
				x.el[i]=b.el[i].sub(su[i]).times(L.row[i].el[L.row[i].nzLength-1].inv());
			}
		}

		return x;
	}
	
	public Mat backElim(Mat A, Mat b){
		int I=A.nRow;
		int J=A.nCol;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		Mat x=new Mat(I,1);
		x.el[0][0]=b.el[0][0]/A.el[0][0];
	
		for(int i=1;i<I;i++){
			double s=0;
			for(int j=0;j<i;j++)
			s=s+A.el[i][j]*x.el[j][0];
			x.el[i][0]=(b.el[i][0]-s)/A.el[i][i];
		}
	
		return x;
	}
	
	public Vect forwardSub(Mat A, Vect b){
		int I=A.nRow;
		int J=A.nCol;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		Vect x=new Vect(I);
		x.el[0]=b.el[0]/A.el[0][0];
	
		for(int i=1;i<I;i++){
			double s=0;
			for(int j=0;j<i;j++)
			s=s+A.el[i][j]*x.el[j];
			x.el[i]=(b.el[i]-s)/A.el[i][i];
		}
	
		return x;
	}
	

	
	
	public Vect jacobi(SpMat A, Vect b,double errMax,int N){

		Vect x=new Vect(b.length);
		x.rand();
		Vect xp= new Vect(x.length);
		Vect D=A.diagSym();
		SpMat R=A.offDiagSym();
		
		Vect invD=D.inv();
	
		double err=1;
		double beta=.3;
		int k=0;
		for( k=0;(k<N+1 && err>errMax);k++){
			xp=x.deepCopy();
			x=b.sub(R.smul(x)).times(invD);
			x=x.times(beta).add(xp.times(1-beta));

			err=x.sub(xp).abs().max()/x.abs().max();

			if(k%20==0)
			System.out.println("k= "+k+"  residual= "+err);	
		}
		
		System.out.println("k= "+k+"  residual= "+err);	
		return x;
		
		
	}
	
	public Vect jacobi(SpMatAsym A, Vect b,double errMax,int N){

		Vect x=new Vect(b.length);
		x.rand();
		Vect xp= new Vect(x.length);
		Vect D=A.diag();
		SpMatAsym R=A.offDiag();
		
		Vect invD=D.inv();
	
		double err=1;
		double beta=.3;
		int k=0;
		for( k=0;(k<N+1 && err>errMax);k++){
			xp=x.deepCopy();
			x=b.sub(R.mul(x)).times(invD);
			x=x.times(beta).add(xp.times(1-beta));

			err=x.sub(xp).abs().max()/x.abs().max();

			if(k%500==0)
			System.out.println("k= "+k+"  residual= "+err);	
		}
		
		System.out.println("k= "+k+"  residual= "+err);	
		return x;
		
		
	}
	
	public Vect silentGaussIter(SpMat A, Vect b,double errMax,int N){
		 return gaussIter(A,b,errMax,N,false);
	}
	
	public Vect gaussIter(SpMat A, Vect b,double errMax,int N){
		 return gaussIter(A,b,errMax,N,true);
	}
	
	public Vect gaussIter(SpMat A, Vect b,double errMax,int N,boolean echo){

		Vect x=new Vect(b.length);
		x.rand();
		Vect xp= new Vect(x.length);
	
		double err=1;
		double beta=1.2;
		int k=0;
		for( k=0;(k<N+1 && err>errMax);k++){
			xp=x.deepCopy();
			x=b.sub(A.mulUpper(x));
			x=forwardSub(A,x);
			x=x.times(beta).add(xp.times(1-beta));
			err=x.sub(xp).abs().max()/x.abs().max();

			if(k%50==0 || echo)
			System.out.println("k= "+k+"  residual= "+err);	
		}
		//if(echo)
		System.out.println("k= "+k+"  residual= "+err);	
		return x;
		
		
	}
	
	
	public Vect gaussIter(SpMatAsym A, Vect b,double errMax,int N,boolean echo){

		Vect x=new Vect(b.length);
		x.rand();
		Vect xp= new Vect(x.length);
	
		double err=1;
		double beta=1.2;
		int k=0;
		for( k=0;(k<N+1 && err>errMax);k++){
			xp=x.deepCopy();
			x=b.sub(A.mul(x));
			x=forwardSub(A,x);
			x=x.times(beta).add(xp.times(1-beta));
			err=x.sub(xp).abs().max()/x.abs().max();

			if(k%50==0 || echo)
			System.out.println("k= "+k+"  residual= "+err);	
		}
		//if(echo)
		System.out.println("k= "+k+"  residual= "+err);	
		return x;
		
		
	}
	
	public Vect steepGrad(Mat A,Vect b, double erc,int N){
		int I=A.nRow;
		int J=A.nCol;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		Vect x=new Vect(b.length);
		x.rand();

		Vect r=new Vect(I);
		r=b.sub(A.mul(x));
		double resMax0=b.abs().max();
		double resMax=x.abs().max();
		int k=0;
		double alpha;
		double resRatio=resMax/resMax0;

		for(k=1;(k<=N && resRatio>erc) ;k++){

			if(k==1 || k%10==0) report("steepestDesent",k,resRatio, resMax);
			
			alpha=r.dot(r)/(r.dot(A.mul(r)));
			x=x.add(r.times(alpha));
			r=r.sub(A.mul(r).times(alpha));
			resMax=r.abs().max();
			resRatio=resMax/resMax0;	
		}
		report("steepestDesent",k,resRatio, resMax);
		return x;
	}
	
	public Vect steepDec(SpMat A,Vect b, double erc, int N){
		int I=A.getnRow();
		int J=A.getnCol();
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		if(I!=b.length) throw new IllegalArgumentException("Array dimensions do not agree.");
		Vect x=new Vect(b.length);
		x.rand();

		Vect r=new Vect(I);
		r=b.sub(A.smul(x));
		double resMax0=b.abs().max();
		double resMax=x.abs().max();
		int k=0;
		double alpha;
		double resRatio=resMax/resMax0;

		for(k=1;(k<=N && resMax>erc) ;k++){
		
			if(k==1 || k%10==0)report("steepestDesent",k,resRatio, resMax);
			
			alpha=r.dot(r)/(r.dot(A.smul(r)));
			x=x.add(r.times(alpha));
			r=r.sub(A.smul(r).times(alpha));
			resMax=r.abs().max();
			resRatio=resMax/resMax0;	
		}
		report("steepestDesent",k,resRatio, resMax);
		return x;
	}
	
	public Vect precSteepDec(SpMat A,Vect b, double erc, int N){
		int I=A.getnRow();
		int J=A.getnCol();
		if(I!=J) throw new IllegalArgumentException("Matrix is not square.");
		if(I!=b.length) throw new IllegalArgumentException("Array dimensions do not agree.");

		Vect x=new Vect(b.length);

		Vect Ci=A.scale(b);
		x=steepDec(A,b,erc,N);
		x.timesVoid(Ci);
	
		return x;
	}
	


	
	
	public Vect ICCG(SpMat A,SpMat L,Vect b,double errMax,int N){
		return ICCG(A,L, b, errMax,N,new Vect(b.length),false);
		}
	

	
	public Mat ICCG(SpMat A,SpMat L,Mat B,double errMax,int N,boolean echo){
		Mat X=new Mat(B.size());
			for(int i=0;i<B.nCol;i++)
		X.setCol(ICCG(A,L, B.getColVect(i), errMax,N,new Vect(B.nRow),echo),i);
		
			return X;
	}
	
	public Mat silentICCG(SpMat A,SpMat L,Mat B,double errMax,int N){
		
			return ICCG(A,L,B,errMax,N,false);
	}
	
	public Vect silentICCG(SpMat A,SpMat L,Vect b,double errMax,int N,Vect x){
		return ICCG( A, L, b, errMax, N, x, false);

		}
	
	public Vect ICCG(SpMat A,SpMat L,Vect b,double errMax,int N,Vect x){

		return ICCG( A, L, b, errMax, N, x, false);
		}
	public Vect ICCG(SpMat A,SpMat L,Vect b,double errMax,int N,Vect x, boolean echo){
		return err1ICCG( A, L, b, errMax, N, x, echo);
	}
	public Vect err1ICCG(SpMat A,SpMat L,Vect b,double errMax,int N,Vect x, boolean echo){
		return ICCG( A, L, b, errMax, N, x,1, echo);
	}
	public Vect err1ICCG(SpMat A,SpMat L,Vect b,double errMax,int N,Vect x){
		return ICCG( A, L, b, errMax, N, x,1, false);
	}
	public Vect err1ICCG(SpMat A,SpMat L,Vect b,double errMax,int N){
		return ICCG( A, L, b, errMax, N, new Vect(b.length),1, true);
	}
	public Vect err0ICCG(SpMat A,SpMat L,Vect b,double errMax,int N,Vect x){
		return ICCG( A, L, b, errMax, N, x,0, false);
	}
	public Vect err0ICCG(SpMat A,SpMat L,Vect b,double errMax,int N,Vect x, boolean echo){
		return ICCG( A, L, b, errMax, N, x,0, echo);
	}
	public Vect err0ICCG(SpMat A,SpMat L,Vect b,double errMax,int N){
		return err0ICCG( A, L, b, errMax, N, new Vect(b.length), false);
	}
	
	public Vect ICCG(SpMat A,SpMat L,Vect b,double errMax,int N,Vect x,int errType, boolean echo){

		
		int I=A.getnRow();
		int J=A.getnCol();
		int Lb=b.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=Lb) throw new IllegalArgumentException("Arrays do not agree");
	
		Vect r=b.sub(A.smul(x));
		Vect v;
		Vect z=solveTriangular(L,r,A);
		
		Vect p=z;
		double res=b.norm();
		double resIni=res;
		if(resRef>0)resIni=resRef;
		int k=0;
		double temp,alpha,beta;
		double  error;
		if(errType==0 || errType==2) error=resIni; 
		else error=res/resIni;
		
	
		report("ICCG",k,error, resIni);

		for(k=1;(k<=N &&  error>errMax && !this.terminate) ;k++){
				
			if(resRef>0)
				errs.add(log10(error*resIni/resRef));
			else
				errs.add(log10(error));

				totalIter++;
				

				v=A.smul(p);

			temp=z.dot(r);

			alpha=temp/(p.dot(v));
			x=x.add(p.times(alpha));
			r=r.sub(v.times(alpha));
		
			z=solveTriangular(L,r,A);
			beta=z.dot(r)/temp;
			p=z.add(p.times(beta));
			
			if(k>0) {
				
					if(errType==0){
						res=r.norm();
						error=res;
					}
					else if(errType==1){
						
						res=r.norm();
						error=res/resIni;
						
						}
					else if(errType==2) {
						res=r.abs().max();
						error=res;
					}
					else if(errType==3) {
						res=r.norm()/I;
						error=res;
					}
					else{
						res=r.abs().max();
						error=res/resIni;
					}
										
					if(k%50==0 || echo){
						report("ICCG",k,error, res);
					}
				}
		}
		
		if(resRef>0)
			errs.add(log10(error*resIni/resRef));
		else
			errs.add(log10(error));

			totalIter++;
		
		//if(echo)
		report("ICCG",k,error, res);


		return x;
	}
	
	public Mat blockCG(SpMat A,Mat B,double errMax,int N){
		return blockCG( A, B, errMax, N,new Mat(B.size()), true);
	}
	
	public Mat blockCG(SpMat A,Mat B,double errMax,int N,Mat X, boolean echo){
		// block CG
			int I=A.getnRow();
			int J=A.getnCol();
			MatSolver ms=new MatSolver();
			
			int L=B.nRow;
			if(I!=J) throw new IllegalArgumentException("Matrix is not square");
			if(I!=L) throw new IllegalArgumentException("Arrays dimensions do not agree");
			Mat R=B.sub(A.smul(X));

			Mat V;
			int m=B.nCol;
			Mat P=R.deepCopy();
			Mat T,G=new Mat(I,m),RtR=new Mat(m,m),C=new Mat(I,m);

			Vect  allRes=null;
			
			Vect resMax0=B.normCol();
			double resMax=resMax0.max();
			int k=0;
			double resRatio=1;

				for(k=1;(k<=N &&  resRatio>errMax) ;k++){
				
					
					allRes=R.normCol();
					
					resRatio=allRes.div(resMax0).max();
					
					if(k%50==1 && echo){
						resMax=allRes.max();
						report("BLCG",k,resRatio, resMax);

					}
			
					
				V=A.smul(P);
				
				 T=P.transp().mul(V);
				 
		
				 RtR=R.transp().mul(R);
		
				 
				G=ms.gaussel(T, RtR);
			
				X=X.add(P.mul(G));

				R=R.sub(V.mul(G));

				C=ms.gaussel(RtR,R.transp().mul(R));


				P=R.add(P.mul(C));
				

				
			
		
		}
	
		resMax=allRes.max();
		report("BLCG",k,resRatio, resMax);
		
		return X;
		}
	
	public Mat blockICCG(SpMat A,SpMat L,Mat B,double errMax,int N){
		
		return blockICCG( A, L, B, errMax, N,new Mat(B.size()), true);
	}
	
	public Mat blockICCG(SpMat A,SpMat L,Mat B,double errMax,int N,Mat X, boolean echo){
		// block CG
			int I=A.getnRow();
			int J=A.getnCol();
			MatSolver ms=new MatSolver();
			
			int K=B.nRow;
			if(I!=J) throw new IllegalArgumentException("Matrix is not square");
			if(I!=K) throw new IllegalArgumentException("Arrays dimensions do not agree");
			Mat R=B.sub(A.smul(X));
			
		//	QRDecomposition qr=new QRDecomposition(new Matrix(R.el));
		// R=new Mat(qr.getQ().getArray());
		/*	int rnk=qr.getQ().rank();
			util.show(;
			util.pr("Rank= "+rnk);*/
			
			Mat Z=new Mat(R.size());
			for(int j=0;j<Z.nCol;j++){
				Vect z=solveTriangular(L,R.getColVect(j),A);
				Z.setCol(z, j);
			}
			

			Mat V;
			int m=B.nCol;
			Mat P=Z.deepCopy();
			Mat T,G=new Mat(I,m),RtZ=new Mat(m,m),C=new Mat(I,m);

			//double resMax0=B.absMat().maxElement();
			Vect resMax0=B.normCol();
			double resMax=resMax0.max();
			int k=0;
			double resRatio=1;

				for(k=1;(k<=N &&  resRatio>errMax) ;k++){
					
				V=A.smul(P);
				
				 T=P.transp().mul(V);

		
				 RtZ=R.transp().mul(Z);

				G=ms.gaussel(T, RtZ);
				
				X=X.add(P.mul(G));

				R=R.sub(V.mul(G));
				
				for(int j=0;j<Z.nCol;j++){
					Vect z=solveTriangular(L,R.getColVect(j),A);
					Z.setCol(z, j);
				}

		
				C=ms.gaussel(RtZ.transp(),Z.transp().mul(R));
			
				
				P=Z.add(P.mul(C));
					
					//resMax=R.absMat().maxElement();
				//	resMax=R.absMat().maxElement();
				resRatio=R.normCol().div(resMax0).max();
					//resMax=R.getColVect(0).norm();
					//resRatio=resMax/resMax0;
					if(k%50==0 || echo){
						report("BLICCG",k,resRatio, resMax);

					}
			
		
		}
		//if(echo)
		report("BLICCG",k,resRatio, resMax);
		return X;
		}
	
	
	public Mat uneqConvBLCG(SpMat A,Mat B,double errMax,int N,Mat X,int errType, boolean echo){
		// block CG unequal converg. ; not fine yet
			int I=A.getnRow();
			int J=A.getnCol();
		//	MatSolver ms=new MatSolver();
			
			int L=B.nRow;
			if(I!=J) throw new IllegalArgumentException("Matrix is not square");
			if(I!=L) throw new IllegalArgumentException("Arrays dimensions do not agree");
			Mat R=B.sub(A.smul(X));
			Mat V;
			Mat P0=R.deepCopy();
			double resMax0=B.absMat().maxElement();
			double resMax=resMax0;
			int k=0;
			double resRatio=resMax/resMax0;
			
			int m=B.nCol;
			boolean[] conv=new boolean[m];

			int q=m;
			
				for(k=1;(k<=N &&  resRatio>errMax) ;k++){
				
				Mat P=new Mat(I,q);
				for(int j=0;j<q;j++)
					P.setCol(P.getColVect(j), j);
				
				V=A.smul(P);
				Mat T=P.transp().mul(V);
				Mat RtR=R.transp().mul(R);
				Mat G=T.inv().mul(RtR);
				X=X.add(P.mul(G));
				R=R.sub(V.mul(G));
				Mat C=RtR.inv().mul(R.transp().mul(R));
				P=R.add(P.mul(C));
				
					if(echo) {
					resMax=R.absMat().maxElement();
					resRatio=resMax/resMax0;
					if(k%50==0 || echo)
						report("CG",k,resRatio, resMax);
				}
		
		}
	//	if(echo)
		report("BLCG",k,resRatio, resMax);
		return X;
		}
	

	
public Vect AICCGx(SpMat A,SpMat L,Vect b,double errMax,int N,Vect x,int errType, boolean echo){
		
		int I=A.getnRow();
		int J=A.getnCol();
		int Lb=b.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=Lb) throw new IllegalArgumentException("Arrays do not agree");
	
		Vect r=b.sub(A.amul(x));
		Vect v;
		Vect z=solveTriangular(L,r,A);

		Vect p=z;
		double resIni=r.norm(),res=resIni;
		int k=0;
		double temp,alpha,beta;
		double  error;
		if(errType==0 || errType==2) error=resIni; 
		else error=1;

			for(k=1;(k<=N &&  error>errMax && !this.terminate) ;k++){

			temp=z.dot(r);
			v=A.amul(p);
			alpha=temp/(p.dot(v));
			x=x.add(p.times(alpha));
			r=r.sub(v.times(alpha));
			z=solveTriangular(L,r,A);
			beta=z.dot(r)/temp;
			p=z.add(p.times(beta));
			
			if(k>0) {
				
					if(errType==0){
						res=r.norm();
						error=res;
					}
					else if(errType==1){
						res=r.norm();
						error=res/resIni;
					}
					else if(errType==2) {
						res=r.abs().max();
						error=res;
					}
					else{
						res=r.abs().max();
						error=res/resIni;
					}
										
					if(k%50==0 || echo){
						report("ICCG",k,error, res);
					}
				}
		
		}
		
	//	if(echo)
		report("ICCG",k,error, res);


		return x;
	}
	
	public Vect ICCGr0max(SpMat A,SpMat L,Vect b,double errMax,int N){
		return ICCG(A,L,b,errMax,N,new Vect(b.length),2,true);
		
	}
	
	public Vect ICCGr1max(SpMat A,SpMat L,Vect b,double errMax,int N){
		return ICCG(A,L,b,errMax,N,new Vect(b.length),3,true);
		
	}

	
	public Vect silentCG(SpMat A,Vect b,double errMax,int N,Vect xini){
		return CG(A,b,errMax,N,xini,false);
	}
	
	public Vect CG(SpMat A,Vect b,double errMax,int N,Vect xini){
		return CG(A,b,errMax,N,xini,true);
	}
	
	public Vect CG(SpMat A,Vect b,double errMax,int N){
		return CG(A,b,errMax,N,new Vect(b.length),true);
	}
	
public Vect errMixICCG(SpMat A,SpMat L,Vect b,double errMax,double er2Max,int N,Vect x, boolean echo){
		
		int I=A.getnRow();
		int J=A.getnCol();
		int Lb=b.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=Lb) throw new IllegalArgumentException("Arrays do not agree");
	
		Vect r=b.sub(A.smul(x));
		Vect v;
		Vect z=solveTriangular(L,r,A);

		Vect p=z;
		double resIni=r.norm(),res=resIni;
		int k=0;
		double temp,alpha,beta;
		double  error;
		 error=1;

			for(k=1;(k<=N &&  (error>errMax && res>er2Max) && !this.terminate) ;k++){

			temp=z.dot(r);
			v=A.smul(p);
			alpha=temp/(p.dot(v));
			x=x.add(p.times(alpha));
			r=r.sub(v.times(alpha));
			z=solveTriangular(L,r,A);
			beta=z.dot(r)/temp;
			p=z.add(p.times(beta));
			
			if(k>0) {
				
					
						res=r.norm();
						error=res/resIni;
					
						
					if(k%50==0 || echo){
						report("ICCG",k,error, res);
					}
				}
		
		}
		
	//	if(echo)
		report("ICCG",k,error, res);


		return x;
	}
	
	
	public Vect CG(SpMat A,Vect b,double errMax,int N,Vect x,boolean echo){
		int I=A.getnRow();
		int J=A.getnCol();
		int L=b.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=L) throw new IllegalArgumentException("Arrays dimensions do not agree");
		Vect r=b.sub(A.smul(x));
		Vect p=new Vect(I);
		Vect v;
		p=r.deepCopy();
		double res0=b.norm();
		double res=r.norm();
		int k=0;
		double c,temp,alpha;
		double resRatio=res/res0;

			for(k=1;(k<=N &&  resRatio>errMax) ;k++){
				
				res=r.norm();
				resRatio=res/res0;
				
				if(k%50==1 && echo)
					report("CG",k,resRatio, res);
				
				
			v=A.smul(p);
	
			temp=r.dot(r);
			alpha=temp/(p.dot(v));
			
			 x=x.add(p.times(alpha));
		
			r=r.sub(v.times(alpha));
			c=r.dot(r)/temp;
			p=r.add(p.times(c));
	
				
	
	
	
	}
//	if(echo)
	report("CG",k,resRatio, res);
	return x;
	}
	
	
	public Vect BiCGSTAB(SpMatAsym A,Vect b,double errMax,int N,Vect x,boolean echo){
		int I=A.getnRow();
		int J=A.getnCol();
		int L=b.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=L) throw new IllegalArgumentException("Arrays dimensions do not agree");
		Vect rh=new Vect(I);
		Vect r=new Vect(I);
		Vect s, t;

		r=b.sub(A.mul(x));
		rh=r.deepCopy();
		double rho=1,alpha=1,w=1,wp=1,rhop,beta;
		Vect v=new Vect(I);
		Vect p=new Vect(I);
		double resMax0=b.abs().max();
		double resMax=r.abs().max();
		int k=0;

		double resRatio=resMax/resMax0;

			for(k=1;(k<=N &&  resRatio>errMax) ;k++){
				rhop=rho;
				rho=rh.dot(r);
				beta=(rho/rhop)*alpha/w;
				
				p=r.add(p.sub(v.times(w)).times(beta));
				v=A.mul(p);
				double rv=rh.dot(v);

				alpha=rho/rv;
				s=r.sub(v.times(alpha));
				t=A.mul(s);
				double t2=t.dot(t);

				w=t.dot(s)/t2;
			
				x=x.add(p.times(alpha)).add(s.times(w));	
				r=s.sub(t.times(w));
			
				if(echo) {
				resMax=r.abs().max();
				resRatio=resMax/resMax0;
				if(k%50==0 || echo)
					report("BiCGSTAB",k,resRatio, resMax);
			}
	
	}
//	if(echo)
	report("BiCGSTAB",k,resRatio, resMax);
	return x;
	}
	
	public Vect BiICCGSTAB(SpMatAsym A,SpMatAsym L,Vect b,double errMax,int N,Vect x,boolean echo){
		int I=A.getnRow();
		int J=A.getnCol();
		int Lb=b.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=Lb) throw new IllegalArgumentException("Arrays dimensions do not agree");
		Vect rh=new Vect(I);
		Vect r=new Vect(I);
		Vect s, t;

		r=b.sub(A.mul(x));
		rh=r.deepCopy();
		double rho=1,alpha=1,w=1,wp=1,rhop,beta;
		Vect v=new Vect(I);
		Vect p=new Vect(I);
		double resMax0=b.abs().max();
		double resMax=r.abs().max();
		int k=0;

		double resRatio=resMax/resMax0;

			for(k=1;(k<=N &&  resRatio>errMax) ;k++){
				rhop=rho;
				rho=rh.dot(r);
				beta=(rho/rhop)*alpha/w;
				
				p=r.add(p.sub(v.times(w)).times(beta));
				Vect y=solveTriangular(L,p,A);
				v=A.mul(y);
				
				double rv=rh.dot(v);

				alpha=rho/rv;
				s=r.sub(v.times(alpha));

				Vect z=solveTriangular(L,s,A);
				
				t=A.mul(z);

				Vect f=solveTriangular(L,t,A);
				double f2=f.dot(f);
				w=f.dot(z)/f2;
			
				x=x.add(y.times(alpha)).add(z.times(w));	
				r=s.sub(t.times(w));
			
				if(echo) {
				resMax=r.abs().max();
				resRatio=resMax/resMax0;
				if(k%50==0 || echo)
					report("BiICCGSTAB",k,resRatio, resMax);
			}
	
	}
	//if(echo)
	report("BiICCGSTAB",k,resRatio, resMax);
	return x;
	}
	
	public void report(String str,int k,double resRatio, double resMax){
		DecimalFormat formatter=new DecimalFormat("0.00E00");

		
		System.out.format("%s iteration: %5d\t error : %s\t res max: %s\n",str,k, formatter.format(resRatio),formatter.format(resMax));	
	}
	
	public Vect DCG(SpMat Ar,Vect br,double errMax,int N,Vect x){
		int I=Ar.getnRow();
		int J=Ar.getnCol();
		int L=br.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=L) throw new IllegalArgumentException("Arrays do not agree");
		Vect b=br.deepCopy();
		
		SpMat A=Ar.deepCopy();
		Vect Ci=A.scale(b);
		x.timesVoid(Ci);
		Vect r=b.sub(A.smul(x));
		Vect p=new Vect(I);
		Vect v;
		p=r.deepCopy();
		double resMax0=b.abs().max();
		double resMax=r.abs().max();
		int k=0;
		double c,temp,alpha;
		double resRatio=resMax/resMax0;

			for(k=1;(k<=N &&  resRatio>errMax) ;k++){
				if(this.terminate) break;
			v=A.smul(p);
			temp=r.dot(r);
			alpha=temp/(p.dot(v));
			x=x.add(p.times(alpha));
			r=r.sub(v.times(alpha));
			c=r.dot(r)/temp;
			p=r.add(p.times(c));
			resMax=r.abs().max();
			resRatio=resMax/resMax0;
			if(k>0) {
				resMax=r.abs().max();
				resRatio=resMax/resMax0;
				if(k%50==0)
					report("DCG",k,resRatio, resMax);
			}
	
	}
	
			report("DCG",k,resRatio, resMax);
	x.timesVoid(Ci);
	return x;
	
	}
	
	public Vect solveQuadratic(SpMat A,SpMat B,Vect c,double erc,int N){
		int I=A.getnRow();
		int J=B.getnRow();
		int K=c.length;
		if(I!=J) throw new IllegalArgumentException("Array dimensions do not agree.");
		
		
		Vect x=new Vect(c.length);
		Vect dx=new Vect(c.length);


		/*	Vect r=new Vect(I);
		r=c.sub(A.smul(x.times(x))).sub(B.smul(x));
	
		double err0=r.norm();
		SpMat A1=A.deepCopy();
		A1.times(0);
		A1.addToDiag(A.diagSym().times(2));
		A1.add(B);
double a=0;
		int k=0;
		double err=1;
		while (err>erc){
			k++;
			if(k%10==0) System.out.println("kkNNR= "+k+" res="+r.norm());
			if(k>N)  {
				System.out.println("kNNR= "+k+" res="+r.norm());
				break;
			}
			Mat A2=x.mul(A1).add(B);
			dx=CG(A2,r,1e-8,500);
			//dx=gaussel(A2,r);
			x=x.add(dx);
			r=c.sub(A.mul(x.times(x))).sub(B.mul(x));
			err=r.norm()/err0;
		}
		System.out.println("kNNR= "+k+" res="+r.norm());*/

		return x;
	}
	
	public VectComp CG(SpMatComp A,VectComp b,double errMax,int N,VectComp x,boolean echo){
		int I=A.getnRow();
		int J=A.getnCol();
		int L=b.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=L) throw new IllegalArgumentException("Arrays dimensions do not agree");
		VectComp r=b.sub(A.mul(x));
		VectComp p=new VectComp(I);
		VectComp v;
		p=r.deepCopy();
		double resMax0=b.abs().max();
		double resMax=r.abs().max();
		int k=0;
		double c,temp,alpha;
		double resRatio=resMax/resMax0;

			for(k=1;(k<=N &&  resRatio>errMax) ;k++){
			v=A.mul(p);
			temp=r.norm2();
			
			alpha=temp/(p.dot(v.conj())).re;
			
	

			x=x.add(p.times(alpha));
			r=r.sub(v.times(alpha));
			c=r.norm2()/temp;
			p=r.add(p.times(c));
			
				if(echo) {
				resMax=r.abs().max();
				resRatio=resMax/resMax0;
				if(k%50==0 || echo)
					report("CG",k,resRatio, resMax);
			}
	
	}
//	if(echo)
	report("CG",k,resRatio, resMax);
	
	return x;
	}
	
	
	public VectComp COICCG(SpMatComp A,SpMatComp L,VectComp b,double errMax,int N,VectComp x,int errType, boolean echo){
	     

		int I=A.getnRow();
		int J=A.getnCol();
		int Lb=b.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=Lb) throw new IllegalArgumentException("Arrays do not agree");

		VectComp r=b.sub(A.smul(x));
		VectComp v;
		VectComp z=solveTriangular(L,r,A,A.hermit);
		VectComp p=z.deepCopy();

		double resIni=r.norm(),res=resIni;
		int k=0;
		Complex temp,alpha,beta;
		double  error;
		if(errType==0 || errType==2) error=resIni; 
		else error=1;

		report("COICCG",k,error, resIni);
		
			for(k=1;(k<=N &&  error>errMax && !this.terminate) ;k++){
				
				if(resRef>0)
					errs.add(log10(error*resIni/resRef));
				else
					errs.add(log10(error));

					totalIter++;

					
			temp=z.dot(r);

			v=A.smul(p);
			Complex tt=p.dot(v);
			alpha=temp.times(tt.inv());

			x=x.add(p.times(alpha));
			r=r.sub(v.times(alpha));
		
			z=solveTriangular(L,r,A,A.hermit);
	
			beta=z.dot(r).times(temp.inv());
			p=z.add(p.times(beta));
	
			if(k>0) {
				
					if(errType==0){
						res=r.norm();
						error=res;
					}
					else if(errType==1){
						res=r.norm();
						error=res/resIni;
					}
					else if(errType==2) {
						res=r.abs().max();
						error=res;
					}
					else{
						res=r.abs().max();
						error=res/resIni;
					}
										
					if(k%50==0 || echo){
						report("COICCG",k,error, res);
					}
				}
		
		}
		
	//	if(echo)
		report("COICCG",k,error, res);


		return x;
	}
	
	public VectComp COCG(SpMatComp A,VectComp b,double errMax,int N,VectComp x,int errType, boolean echo){
	     

		int I=A.getnRow();
		int J=A.getnCol();
		int Lb=b.length;
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=Lb) throw new IllegalArgumentException("Arrays do not agree");

		VectComp r=b.sub(A.smul(x));
		VectComp v;

		VectComp p=r.deepCopy();

		double resIni=r.norm(),res=resIni;
		int k=0;
		Complex temp,alpha,beta;
		double  error;
		if(errType==0 || errType==2) error=resIni; 
		else error=1;

			for(k=1;(k<=N &&  error>errMax && !this.terminate) ;k++){
			temp=p.dot(r);

			v=A.smul(p);
			Complex tt=p.dot(v);
			alpha=temp.times(tt.inv());

			x=x.add(p.times(alpha));
			r=r.sub(v.times(alpha));
		
			beta=r.dot(r).times(temp.inv());
			p=r.add(p.times(beta));
	
			if(k>0) {
				
					if(errType==0){
						res=r.norm();
						error=res;
					}
					else if(errType==1){
						res=r.norm();
						error=res/resIni;
					}
					else if(errType==2) {
						res=r.abs().max();
						error=res;
					}
					else{
						res=r.abs().max();
						error=res/resIni;
					}
										
					if(k%50==0 || echo){
						report("COCG",k,error, res);
					}
				}
		
		}
		
	//	if(echo)
		report("COCG",k,error, res);


		return x;
	}
	
	
	
	public Vect CG(BandMat A,Vect b,double errMax,int N,Vect x,boolean echo){
		int I=A.nRow;
	//	int J=A.getnCol();
		int L=b.length;
		//if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		if(I!=L) throw new IllegalArgumentException("Arrays dimensions do not agree");
		Vect r=b.sub(A.mul(x));
		Vect p=new Vect(I);
		Vect v;
		p=r.deepCopy();
		double resMax0=b.abs().max();
		double resMax=r.abs().max();
		int k=0;
		double c,temp,alpha;
		double resRatio=resMax/resMax0;

			for(k=1;(k<=N &&  resRatio>errMax) ;k++){
				
			v=A.mul(p);
	
			temp=r.dot(r);
			alpha=temp/(p.dot(v));
			
			 x=x.add(p.times(alpha));
		
			r=r.sub(v.times(alpha));
			c=r.dot(r)/temp;
			p=r.add(p.times(c));
	
				resMax=r.abs().max();
				resRatio=resMax/resMax0;
				if(k%50==0 || echo)
					report("CG",k,resRatio, resMax);
	
	
	}
	//if(echo)
	report("CG",k,resRatio, resMax);
	return x;
	}
	
	
	public Vect solvelu(BandMat A, Vect b){
		Vect x=solveLowlu(A,b);
		x=solveUp(A,x);
		return x;
		
	}
	
	public Vect solveLowlu(BandMat A, Vect b){
		int I=A.nRow;
		Vect x=new Vect(I);
		x.el[0]=b.el[0];
	
		
		for(int i=1;i<I;i++){
			double s=0;
			for(int j=max(0,i-A.kL-1);j<i;j++)
			s=s+A.get(i,j)*x.el[j];
			x.el[i]=(b.el[i]-s);
		}
	
		return x;
	}
	
	public Vect solveUp(BandMat A, Vect b){
		int I=A.nRow;
		Vect x=new Vect(I);
		x.el[I-1]=b.el[I-1]/A.get(I-1,I-1);
	
		for(int i=I-2;i>=0;i--){
			double s=0;
			double dd=A.get(i,i);
			for(int j=i+1;j<min(I,i+A.kL+1);j++){
			s=s+A.get(i,j)*x.el[j];
			}
			x.el[i]=(b.el[i]-s)/dd;
		}
	
		return x;
			}
	
	
/*	public Vect gaussel(BandMat A, Vect b){
		
		
	}*/
	

	public void terminate(boolean b){
		this.terminate=b;

	}


	
}
