package femSolver;

import java.text.DecimalFormat;

import fem.Model;
import math.Complex;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatComp;
import math.SpVect;
import math.Vect;
import math.util;


public class TransientNolinearMagSolver {

	boolean relaxedNR=true;
	int totalNonlinIter;
	double coilInduct=0e-2;
	double theta=.5;
	
	boolean usePrev=false;
	int stpNumb=0;
	
	public TransientNolinearMagSolver(){	}

	public TransientNolinearMagSolver(Model model)	{}


	public Vect solve(Model model,Vect x,boolean echo,int step){


		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;

		Vect Ci;
		SpMat L;
		Vect dA=new Vect(model.numberOfUnknowns);

		Vect[] B1,B2;

		model.solver.terminate(false);

		double errNR=1,resNR=0,resNR0=1;
		double errFlux=1;

		int nonLinIter=0;
		
		
		if(model.Ss==null)
			model.magMat.setConductMat(model);

		
		model.magMat.setRHS(model);

		for( nonLinIter=0; errNR>model.errNRmax && nonLinIter<nonLinIterMax;nonLinIter++)

			
		{
			totalNonlinIter++;

			if(echo)
			{
				System.out.println();
				System.out.println();
				System.out.println("    nonlinear Iteration: "+nonLinIter+
						"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
				System.out.println(" Flux    error: "+dfe.format(errFlux));
				System.out.println();
			}

			SpMat Ks=new SpMat();

			Vect b=new Vect();

						//model.setMagMat();
			model.magMat.setReactMat(model);
			
	
			Vect dv=new Vect();


				if(model.eddyTimeIntegMode==0)
				{


					Ks=model.Hs.addSmallerNew(model.Ss);
					Vect vp=model.getUnknownAp();

					Vect vi=model.getUnknownA();

					dv=vp.sub(vi);
					
				

					b=model.RHS.sub(model.HkAk).add(model.Ss.smul(dv));
				}

				else if(model.eddyTimeIntegMode==1){

//crank
					Vect vp=model.getUnknownAp();
					Vect vi=model.getUnknownA();
			

					b=model.RHS.add(model.bT).sub(model.HkAk.add(model.HpAp).times(0.5)).add(model.Ss.smul(vp.sub(vi)));

					Ks=model.Hs.timesNew(.5).addSmallerNew(model.Ss);

				}

	


			if(nonLinIter==0){
				resNR0=b.norm();
				model.solver.resRef=resNR0;
			}
			
			resNR=b.norm();
			
		
			errNR=resNR/resNR0;



			Ci=Ks.scale(b);
			L=Ks.ichol();

			//if(b.abs().max()>1e-6){
				//dA=model.solver.err0ICCG(Ks,L,b,1e-3*model.errCGmax,iterMax);	
				dA=model.solver.ICCG(Ks,L, b,model.errCGmax,iterMax);

			//}

			if(model.solver.terminate) break;


			dA.timesVoid(Ci);
			
			
			x=x.add(dA);
			

			model.up=x.deepCopy();

			
			
			B1=model.getAllB();

			model.setSolution(x);	

			model.setB();	
			
			B2=model.getAllB();
			
			errFlux=model.getDiffMax(B1,B2);

			
	

		}
		


		model.HpAp=model.HkAk.deepCopy();
		
		System.out.println();
		System.out.println("-----------------------------------------------------");
		System.out.println(" Nonlinear Iteration: "+nonLinIter+
				"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
		System.out.println("Flux    error: "+dfe.format(errFlux));
		System.out.println("-----------------------------------------------------");
		System.out.println();
		
		System.out.println();
		System.out.println("=======================================================");
		System.out.println("Total number of ICCG iterations: "+model.solver.totalIter);
		System.out.println("Total number of Nonlinear iterations: "+totalNonlinIter);
		
		
		System.out.println();
	
		System.out.println(" Final NR residual: "+resNR);
		
		System.out.println("=======================================================");
		System.out.println();
		


		return x;

	}
	
	public SpMat getCircuitHs(Model model,int step){


		SpMat Hs=model.Hs.deepCopy();


		
	double cf=1;
	
		if(model.eddyTimeIntegMode==-3) cf=this.theta;

		
		if(model.analysisMode>0){
			Hs.addSmaller(model.Ss);

		}
		
		return Hs;

	}
	



}
