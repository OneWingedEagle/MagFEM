package math;

import static java.lang.Math.*;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import Jama.QRDecomposition;

public class Eigen {
	public Vect lam;
	public Mat V;
	public Vect vmax;
	public Vect vmin;

	EigenvalueDecomposition evd;

	public Eigen(){}

	public Eigen(Mat A){
		V=new Mat(A.size());
		Mat B=A.deepCopy();
		evd=new EigenvalueDecomposition(new Matrix(B.el));
		lam=new Vect(evd.getRealEigenvalues());
		Mat X=new Mat(evd.getV().getArray());

		int[] index=lam.bubble();
		for(int i=0;i<lam.length;i++)
			V.setCol(X.getColVect(index[i]),i);

		vmin=V.getColVect(0);
		vmax=V.getColVect(lam.length-1);

	}

	public double eigMax(Mat A){
		return eigMax(A,1e-8);
	}

	public double eigMax(Mat A,double errMax) {

		Vect v=new Vect(A.nRow);
		v.rand();
		v.normalize();
		Vect vp;
		double err=1;
		int k=0;
		while(err>errMax){
			k++;
			if(k>1000) {  break;}				
			vp=v;
			v=A.mul(v);
			v.normalize();
			err=v.sub(vp).max();

		}
		vmax=v;
		return A.mul(v).norm();

	}


	public double eigMin(Mat A){
		return eigMin(A,1e-8);
	}
	public double eigMin(Mat A,double errMax) {
		MatSolver solver=new MatSolver();
		Vect v=new Vect(A.nRow);
		v.rand();
		v.normalize();
		Vect vp;
		double err=1;
		int k=0;
		while(err>errMax){
			k++;
			if(k>1000) {  break;}			
			vp=v;
			v=solver.gaussel(A,v);
			v.normalize();
			err=v.sub(vp).max();

		}
		System.out.println(" iterantions: "+k+" error: "+  err);

		vmin=v;
		return A.mul(v).norm();

	}

	public void general(Mat A,Mat B){

		Eigen eg1=new Eigen(B.inv().mul(A));

		this.V=eg1.V;
		this.lam=eg1.lam;

	}

	public void generalSym(Mat A,Mat B){

		if(B.norm()!=0) {
			Mat L=B.chol();
			Mat invL=L.invL();
			Mat A_hat=invL.mul(A).mul(invL.transp());	
			Eigen eg1=new Eigen(A_hat);
			V=invL.transp().mul(eg1.V);
			this.lam=A.mul(V).normCol().div(B.mul(V).normCol());
			vmin=V.getColVect(0);
			vmax=V.getColVect(lam.length-1);
		}
		else
			throw new IllegalArgumentException("Right hand side Matrix is zero.");


	}

	public void subspace(Mat A,int p,double errMax){

		Eigen eg1=new Eigen();
		MatSolver solver=new MatSolver();
		int n=A.nRow;
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

			X_hat=solver.gaussel(A,X);	

			K_hat=X_hat.transp().mul(X);
			M_hat=X_hat.transp().mul(X_hat);		
			eg1.generalSym(K_hat, M_hat);
			X=X_hat.mul(eg1.V);

			lamqr=lamq.deepCopy();
			lamq=A.mul(X).normCol();
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

		V=new Mat(n,p);
		for(int j=0;j<p;j++)
			V.setCol(X.getColVect(j), j);
		this.lam=lam;	

	}

	public void eigSubspace(Mat K,Mat M,int p,double errMax){

		MatSolver solver=new MatSolver();
		int n=K.nRow;
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
			X_hat=solver.gaussel(K,X);	

			K_hat=X_hat.transp().mul(X);
			M_hat=X_hat.transp().mul(M.mul(X_hat));		
			K_hat.eigGenSym(M_hat,Q);
			X=X_hat.mul(Q);

			lamqr=lamq.deepCopy();
			lamq=K.mul(X).normCol().div0(M.mul(X).normCol());
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

		V=new Mat(n,p);
		for(int j=0;j<p;j++)
			V.setCol(X.getColVect(j), j);
		this.lam=lam;	


	}

	public void eigQR(Mat M, double errMax) {
		QRDecomposition qd;

		Mat A=M.deepCopy();
		Mat Ap=new Mat(A.size());
		Mat Q=new Mat(M.size());
		Mat P=new Mat(M.size());
		P.eye();
		Mat R=new Mat(A.size());
		int k=0;
		double err=1;
		while(err>errMax){
			k++;
			if(k%10==9) Ap=A;
			if(k>1000) {  break;}

			//Q=A.QRHous(R);
			Q=A.QR(R);
			//P=P.mul(Q);
			if(k%10==0){
				err=A.sub(Ap).absMat().maxElement()/A.diagonal().abs().min();

			}
			A=R.mul(Q);	
		}
		System.out.println(" iterantions: "+k+" error: "+  errMax);
		Vect v=R.diagonal();
		v.bubble();
		lam=v;

		V=P.mul(A);
		V.normalizeColumns();
	}

	public double eigMax(SpMat As,double errc) {



		Vect v=new Vect(As.nRow);
		v.rand();
		v.normalize();
		double lam=1,lamr;
		double err=1;
		int k=0;
		while(err>errc){
			k++;
			if(k>1000) {  break;}				
			v=As.smul(v);
			v.normalize();

			lamr=lam;
			lam=As.smul(v).norm();
			if(lam!=0)
				err=abs(lam-lamr)/abs(lam);
			else if(lam!=0)
				err=abs(lam-lamr)/abs(lamr);
			else
				err=abs(lam-lamr)/1e-20;

		}
		System.out.println(" iterantions: "+k+" error: "+  err);

		return lam;
	}

	public double eigMin(SpMat As,double errc, SpMatSolver solver) {


		Vect v=new Vect(As.nRow);
		Vect xini=new Vect(As.nRow);
		v.rand();
		v.normalize();
		Vect b;
		double lam=1,lamr;
		double err=1;
		int k=0;
		SpMat A=As.deepCopy();
		//Vect Ci=A.scale(xini);
		SpMat L=A.ichol();
		xini=new Vect(A.nRow);
		while(	!solver.terminate &&  err>errc){
			k++;
			if(k>1000) {  break;}			
			//	b=v.times(Ci);
			b=v.deepCopy();
			v=solver.silentICCG(A,L, b,1e-6,1000,v.times(lam));
			//	v=solver.silentICCG(Ks,L,x,1e-6,1000,x.times(1/lamq.el[j]).div(Ci));	
			xini=v.deepCopy();
			//v.timesVoid(Ci);
			v.normalize();
			lamr=lam;
			lam=As.smul(v).norm();
			if(lam!=0)
				err=abs(lam-lamr)/abs(lam);
			else if(lam!=0)
				err=abs(lam-lamr)/abs(lamr);
			else
				err=abs(lam-lamr)/1e-20;

			System.out.println(" iterantions: "+k+" error: "+  err);
		}

		System.out.println(" iterantions: "+k+" error: "+  err);
		return lam;

	}

	public double eigMin(SpMat As,SpMat Bs, double errc, SpMatSolver solver) {

		Vect v=new Vect(As.nRow);
		Vect xini=new Vect(As.nRow);
		v.rand();
		v.normalize();
		Vect vp,b;

		double err=1;
		int k=0;
		SpMat A=As.deepCopy();
		Vect Ci=A.scale(xini);
		SpMat L=A.ichol();
		xini=new Vect(As.nRow);
		while(	!solver.terminate &&  err>errc){
			k++;
			if(k>1000) {  break;}			
			vp=v;
			b=Bs.smul(v);		

			b.timesVoid(Ci);
			v=solver.silentICCG(A,L, b,1e-6,1000,xini);

			v.timesVoid(Ci);
			v.normalize();
			err=v.sub(vp).abs().max();
			System.out.println(" iterantions: "+k+" error: "+  err);
		}

		System.out.println(" iterantions: "+k+" error: "+  err);
		return As.smul(v).norm()/Bs.smul(v).norm();

	}

	public Vect subspace(SpMat As,int p,double errMax){

		return subspace(As,p,errMax,new SpMatSolver());

	}
	public Vect subspace(SpMat As,int p,double errMax, SpMatSolver solver){
		return subspace(As,p, new Mat(As.nRow,p),errMax,solver);
	}

	public Vect subspace(SpMat As,int p, Mat Q,double errMax, SpMatSolver solver){

		Eigen eg=new Eigen();
		int n=As.nRow;
		p=min(p,n);
		int q=min(min(p+8,2*p),n);
		Mat X=new Mat(n,q);
		Mat X_hat=new Mat(n,q);
		Vect x=new Vect(n);
		Vect errv=new Vect(q);
		Vect lamqr=new Vect(q);
		Vect lamq=new Vect(q);
		lamq=lamq.ones(q);

		int[] index=new int[q];
		for(int i=0;i<q;i++)
			index[i]=i;

		X.setCol(new Vect().ones(n), 0);
		for(int i=1;i<q-1;i++)
			X.el[i][i]=1;
		x.rand();
		X.setCol(x, q-1);

		boolean[] locked=new boolean[q];

		Mat K_hat=new Mat();
		Mat M_hat=new Mat();
		int iter;
		int nconv=0;
		SpMat Ks=As.deepCopy();
		Vect Ci=Ks.scale();
		SpMat L=Ks.ichol(1);

		for(iter=0;iter<200;iter++){
			for(int j=0;j<q;j++){				
				x=X.getColVect(index[j]);	
				if(!locked[j])
				{

					Vect xini=x.times(1/lamq.el[j]);
					x.timesVoid(Ci);
					x=solver.silentICCG(Ks,L,x,1e-6,1000,xini.div(Ci));	
					x.timesVoid(Ci);
				}
				else{
					if(lamq.el[j]!=0)
						x=x.times(1/lamq.el[j]);
				}

				if( solver.terminate) break;	

				X_hat.setCol(x,index[j]);		
			}
			K_hat=X_hat.transp().mul(X);
			M_hat=X_hat.transp().mul(X_hat);	
			eg.generalSym(K_hat, M_hat);
			X=X_hat.mul(eg.V);


			lamqr=lamq.deepCopy();
			lamq=As.smul(X).normCol();
			index=lamq.bubble();

			errv=lamq.sub(lamqr).abs().div0(lamq);


			for( int i=0;i<p;i++)
				if(!locked[i] && errv.el[i]<errMax){
					locked[i]=true;
					nconv++;
				}


			//errv.hshow();
			System.out.println("Subspace Iteration: " +iter+"   | error: "+errv.el[p-1]+   "  |  Numver of converged eigenvalues: "+nconv);

			if( nconv==p || solver.terminate) {  break;}

		}

		Vect lam=new Vect(p);		
		for(int i=0;i<p;i++)
			lam.el[i]=lamq.el[i];

		for(int i=0;i<p;i++)
			Q.setCol(X.getColVect(index[i]),i);


		return lam;	

	}


	public Vect subspace(SpMat As,int p,SpMat Ms, double errMax,int itMax, SpMatSolver solver)
	{
		return subspace(As,Ms,p,new Mat(As.nRow,p),errMax,itMax,solver);
	}
	public Vect subspace(SpMat As,SpMat Ms,int p,Mat Q,double errMax,int itMax, SpMatSolver solver){

		solver.terminate=false;

		Eigen eg=new Eigen();
		int n=As.nRow;

		int q=min(min(p+8,2*p),n);


		Mat X=new Mat(n,q);

		Mat X_hat=new Mat(n,q);
		Vect x=new Vect(n);
		Vect errv=new Vect(q);
		Vect lamqr=new Vect(q);
		Vect lamq=new Vect(q);
		lamq=lamq.ones(q);

		int[] index=new int[q];
		for(int i=0;i<q;i++)
			index[i]=i;

		X.setCol(new Vect().ones(n), 0);
		for(int i=1;i<q-1;i++)
			X.el[i][i]=1;
		x.rand();
		X.setCol(x, q-1);



		boolean[] locked=new boolean[q];

		Mat K_hat=new Mat();
		Mat M_hat=new Mat();
		int iter;
		int nconv=0;
		SpMat Ks=As.deepCopy();
		Vect Ci=Ks.scale();
		SpMat L=Ks.ichol();
		Mat X2;
		for(iter=0;iter<200;iter++){
			X2=Ms.smul(X);
			for(int j=0;j<q;j++){				

				if(!locked[j])
				{
					x=X2.getColVect(index[j]);		
					Vect xini=X.getColVect(index[j]).times(1/lamq.el[j]);
					x.timesVoid(Ci);

					//x=solver.errMixICCG(Ks,L,x,errMax,1e-10,itMax,xini.div(Ci),false);	
					x=solver.CG(Ks,x,1e-6,itMax,xini.div(Ci),false);	

					x.timesVoid(Ci);

				}

				else{
					x=X.getColVect(index[j]);	
					if(lamq.el[j]!=0)
						x=x.times(1/lamq.el[j]);

				}


				X_hat.setCol(x,index[j]);		



				if( solver.terminate) break;	
			}



			K_hat=X_hat.transp().mul(X2);
			M_hat=X_hat.transp().mul(Ms.smul(X_hat));	

			eg.generalSym(K_hat, M_hat);

			X=X_hat.mul(eg.V);


			lamqr=lamq.deepCopy();

			lamq.hshow();
			lamq=As.smul(X).normCol().div0(Ms.smul(X).normCol());
			index=lamq.bubble();

			errv=lamq.sub(lamqr).abs().div0(lamq);

			for( int i=0;i<p;i++)
				if(!locked[i] && errv.el[i]<errMax){	
					locked[i]=true;
					nconv++;
				}		

			errv.hshow();
			System.out.println("Subs.Iter.: " +iter+   " | numb. of conv. eigenvalues: "+nconv+"  | error: "+errv.el[p-1]);

			if( nconv==p || solver.terminate) {  break;}
		}

		Vect lam=new Vect(p);	

		for(int i=0;i<p;i++)
			lam.el[i]=lamq.el[i];

		System.out.println("Subspace Iteration: " +iter+"   error: "+errv.el[p-1]);

		for(int i=0;i<p;i++)
			Q.setCol(X.getColVect(index[i]),i);

		return lam;	
	}



}
