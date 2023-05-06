package fem;


import static java.lang.Math.PI;

import java.io.File;
import java.text.DecimalFormat;

import Jama.Matrix;
import Jama.QRDecomposition;
import femSolver.ACMagSolver;
import femSolver.StaticMORNonlinearMagSolver;
import main.Main;
import math.Complex;
import math.DFT;
import math.Eigen;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatSolver;
import math.SpVect;
import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class POD {
	private DecimalFormat formatter=new DecimalFormat("0.00");

	public POD(){}

	
	public void setMagPOD(Model model, Main main){
		
		//model.POD=-1;
		model.snapShot=1;
		
		double tStart=System.currentTimeMillis();

		String fluxFolder="";
	
		if(model.saveFlux){
			fluxFolder = model.resultFolder+"\\fluxes";
		File dfolder = new File(fluxFolder);
		if(dfolder.exists())
			util.deleteDir(dfolder);
		dfolder.mkdir();

	}
		
		//model.writeMesh(fluxFolder+File.separator+"bun.txt");


		model.setMagBC();

		model.solveCoils();
		
		String solutionfile="\\D:\\Works and Studies\\POD2023\\inductor3D\\solutions10.txt";
		util.pr(solutionfile);
		//String solutionfile=System.getProperty("user.dir") +"\\solutionsA18LockedRotor.txt";
		
		//String deltasolutionfile=System.getProperty("user.dir") +"\\solutionsdeltaA180deg0.5Rough.txt";
	//	String solutionfile=System.getProperty("user.dir") +"\\solutionsA90deg1Rough.txt";
		//solutionfile=System.getProperty("user.dir") +"\\solutions90Cogging.txt";

		int L=model.numberOfUnknowns;

		int D=model.snapShot;


		util.pr(" POD Using "+D+" raw basis.");
		
		int fineStep=1;//model.nInc;

		Vect TB=new Vect(model.nEnd-model.nBegin+1);
		

		Vect dA=new Vect(model.numberOfUnknowns);
		Vect dAp=new Vect(model.numberOfUnknowns);
		MatSolver ms=new MatSolver();
		
		Mat	As=new Mat(model.loader.loadArrays(L,D, solutionfile));
	//	Mat	dAs=new Mat(model.loader.loadArrays(L,D1, dsolutionfile));

		

		int numCompBlocks=(model.nEnd)+1;
	
		int jx=0;
		
		model.magMat.setReactMat(model);
		
		if(model.analysisMode>0){
			model.magMat.setConductMat(model);
			model.Ss.times(model.nInc);

			}



		Mat W=new Mat(L,D);
		
		
		for(int j=0;j<D;j++)
			W.setCol(As.getColVect(j), j);
	
		W.normalizeColumns();
	
		
	//	util.pr("POD use Using "+D+" indepenent basis.");
		
		Mat C=W.transp().mul(W);
			
		 Eigen eg2=new Eigen(C);
		 
		 Mat Q=eg2.V;
		 
		 Mat Phi=W.mul(Q);

		 Mat PhiT=Phi.transp();

			int nx0=2352;
			int cmp=0;
			nx0=1558; cmp=1; // reactor
			//nx0=15;
			nx0=24697; cmp=0;// motor half
	
			double dts=model.dt/model.nInc;

		Mat Sr=null;

		if(model.analysisMode>0)
			Sr=PhiT.mul(model.Ss.smul(Phi));
		
		Mat Mr1=PhiT.mul(model.Hs.smul(Phi));
		
		
		Mat Kr=null;
		if(model.analysisMode>0){
			Kr=Mr1.add(Sr);
		}

		Vect dAr=new Vect(D);
		Vect dArp=new Vect(D);
		int nTsteps=model.nTsteps;
		Vect T=new Vect(nTsteps);
		int ix=0;
	//
			for(int i=	model.nBegin;i<=	model.nEnd;	i+=model.nInc){
				
			//	main.gui.tfX[0].setText(Integer.toString(j)+"/"+(model.nEnd));
			//	main.gui.tfX[1].setText(formatter.format(model.TrqZ));
				model.currentTimeStep=i+1;
				
			if(model.phiCoils==null){
				model.setJ0();	
			}
			
			model.magMat.setRHS(model);
			//model.magMat.setRHS(model,false);
			//model.RHS.show();
					
			Vect br=null;
			if(model.analysisMode==0){
				br=PhiT.mul(model.RHS);
			}
			else{
				//br=PhiT.mul(model.RHS.add(model.Ss.smul(dAp)));
				br=PhiT.mul(model.RHS).add(Sr.mul(dArp));
				//br=new Vect().ones(br.length);

			}

	
			if(br.norm()>1e-6){

				dAr=ms.gaussel(Kr, br);
			}
			else
				dAr=new Vect(br.length);	

			dArp=dAr.deepCopy();

			//dAp=dA.deepCopy();
		    dA=Phi.mul(dAr);
		     
		    // Vect Ap=A1.times(1-alpha).add(A2.times(alpha));

		    Vect A=dA;//Ap.add(dA); 
		    
		//model.setPODSolution(A,alpha);


			 model.setSolution(A);
				
			 model.setB();
				
				TB.el[ix++]=model.element[100].getB().el[1];

				int nx=Math.min(nx0,model.numberOfNodes);
				

				if(model.saveFlux)
					if(model.saveFlux){
						String fluxFile = fluxFolder+"\\flux"+i+".txt";
					
						model.writeB(fluxFile);
					}


				if(model.solver.terminate) break;

			}
			


		if(model.saveFlux)
			model.writeMesh(fluxFolder+"\\bun"+0+".txt");
		
		util.plot(TB);
		
		TB.show();
		
		double tEnd=System.currentTimeMillis();
		
		util.pr("-------------------");
		util.pr("Elaspsed Time (sec):");
		util.pr((tEnd-tStart)/1000.0);
		
/*		 Complex[] Y=DFT.dft(TB.el);
		 int N=Y.length;
		  for(int k=0;k<Y.length;k++){
				 util.pr(Y[k].norm()/N);
		}*/
	}
		
	}


















