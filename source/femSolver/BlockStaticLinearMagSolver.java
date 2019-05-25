package femSolver;

import Jama.Matrix;
import Jama.QRDecomposition;
import fem.Model;
import math.Mat;
import math.SpMat;
import math.SpMatSolver;
import math.SpVect;
import math.Vect;
import math.util;


public class BlockStaticLinearMagSolver{
	int stepNumb;
	boolean usePrev=false;

	public BlockStaticLinearMagSolver(){	}

	public Mat solve(Model model, int step){


		int N=8;

		SpMat L=new SpMat();

		

		model.solver.terminate(false);

		Mat bkRHS=new Mat(model.numberOfUnknowns,N);
		
		Mat X=new Mat(bkRHS.size());
		


		double dt=.003;

		
		for(int j=0;j<N;j++){
			double a=Math.cos(j*314*dt);//*(1+.5*Math.random());
		
			for(int ir=1;ir<=model.numberOfRegions;ir++){
				if(model.region[ir].hasJ){
					Vect J=model.region[ir].getJ();
					J=J.times(a);
					for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


						model.element[i].setJ(J);
					}
				}

			}

			model.magMat.setRHS(model);

			bkRHS.setCol(model.RHS, j);
		}
		
		QRDecomposition qr=new QRDecomposition(new Matrix(bkRHS.el));
		Mat Q=new Mat(qr.getQ().getArray());
		Mat R=new Mat(qr.getR().getArray());
		
		int ix=0;
		for(int i=0;i<R.nRow;i++){
			if(R.rowVect(i).norm()>1e-6) ix++;
		}
		Mat B=new Mat(model.numberOfUnknowns,ix);
		
		 ix=0;
		for(int i=0;i<Q.nCol;i++){
			if(R.rowVect(i).norm()>1e-6){
				B.setCol(Q.getColVect(i),ix);
				ix++;
			}
		}

		if(step==0)
			model.setMagMat();


		SpMat Ks=model.Hs.deepCopy();

	//	Vect Ci=Ks.scale(model.RHS);



		L=Ks.ichol();

				//x=model.solver.ICCG(Ks,L, model.RHS,model.errCGmax,model.iterMax,model.xp);

				X=model.solver.blockICCG(Ks,L, B,model.errCGmax,model.iterMax);	
		
		for(int j=0;j<B.nCol;j++){

			Vect x=model.solver.ICCG(Ks,L, B.getColVect(j),model.errCGmax,model.iterMax);	
		//	Vect x=X.getColVect(j);
		model.setSolution(x);
		model.setB();
		X.setCol(x, j);
		}

			
	
		System.out.println("Bmax ( linear analysis): "+model.Bmax);

		return X;



	}
	
	public Mat solvexx(Model model, int step){

			model.setMagMat();
			
			int M=1;
			
			Mat bkRHS=new Mat(model.numberOfUnknowns,M);
			
			Mat X=new Mat(bkRHS.size());
/*			for(int j=0;j<M;j++){
				Vect x=new Vect(model.numberOfUnknowns);
				x.rand();
				//X.setCol(x, j);
			}*/

			double dt=.0033;

			for(int j=0;j<M;j++){
				double a=Math.cos(j*314*dt);//*(1+.5*Math.random());
				for(int ir=1;ir<=model.numberOfRegions;ir++){
					if(model.region[ir].hasJ){
						Vect J=model.region[ir].getJ();
						J=J.times(a);
						for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


							model.element[i].setJ(J);
						}
					}



				}
				model.magMat.setRHS(model);

				bkRHS.setCol(model.RHS, j);
			}
			

			

/*			QRDecomposition qr=new QRDecomposition(new Matrix(bkRHS.el));
			Mat Q=new Mat(qr.getQ().getArray());
			bkRHS=Q.deepCopy();*/
		

			int Nmax=model.iterMax;

		SpMat Ms=model.Hs.deepCopy();
		
		Vect Ci=Ms.scale();
		//Ms.lower=true;
int N=Ms.nRow;

		Mat B=new Mat(N,M);
		for(int i=0;i<M;i++){
			Vect v=new Vect().rand(N, 0, 1);
			Vect v2=Ms.smul(v);
			B.setCol(bkRHS.getColVect(i), i);
			//B.setCol(v2, i);
		}
		

		
		SpMat Ls=Ms.ichol(1.05);
	//	Ls.show();
		
		
		SpMatSolver solver=new SpMatSolver();
		double t1=System.currentTimeMillis();
		for(int i=0;i<M;i++){
			Vect x=solver.ICCG(Ms,Ls,  B.getColVect(i), 1e-6, Nmax);
		//	Vect x=solver.CG(Ms,  B.getColVect(i), 1e-6, Nmax);
			
		}
		double t2=System.currentTimeMillis();
		
		//Ms.show();
		//System.out.println("Main method is empty");
		
	X=solver.blockICCG(Ms, Ls, B, 1e-6, Nmax);
		//X=solver.blockCG(Ms, B, 1e-6, Nmax,X,1,true);
	
	double t3=System.currentTimeMillis();

		System.out.println("Futsuo : "+(t2-t1));
		System.out.println("Block : "+(t3-t2));
		System.out.println("speed up : "+(t2-t1)/(t3-t2));
	//	util.pr(Ms.smul(X).sub(B).norm()/B.norm());	
		
		
		return X;
		
	}
	




}
