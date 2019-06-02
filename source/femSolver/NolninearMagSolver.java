package femSolver;

import java.text.DecimalFormat;

import fem.Model;
import math.Complex;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.SpMatComp;
import math.SpVect;
import math.SpVectComp;
import math.Vect;
import math.VectComp;
import math.util;

public class NolninearMagSolver {

	boolean relaxedNR=true;
	int totalNonlinIterx;
	double coilInduct=0e-2;
	double theta=.5;
	
	boolean usePrev=false;
	int stpNumb=0;
	
	public NolninearMagSolver(){	}

	public NolninearMagSolver(Model model)	{}


	public Vect solve(Model model,Vect x,boolean echo,int step){

	//	if(relaxedNR)
		//return solveNonLinearRelaxed( model, x, echo, step);

		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");

		int mode=1;
		
		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;
		double errNRMax=model.errNRmax;
		double errFluxMx=model.errFluxMax;
		double errFlux=0;
		double errNR=1;
		double resNR=0;
		double resNR0=0;
		boolean fluxErr=false;
		boolean scaling=true;
		
		if(fluxErr) errNRMax=errFluxMx;
		Vect Ci=new Vect();
		SpMat L;
		Vect dA=new Vect(model.numberOfUnknowns);

		Vect[] B1=null,B2=null;

		model.solver.terminate(false);

		int nonLinIter=0;


		for( nonLinIter=0; errNR>errNRMax && nonLinIter<nonLinIterMax;nonLinIter++)
			
		{
			totalNonlinIterx++;

			if(echo)
			{
				System.out.println();
				System.out.println("-----------------------------------------------------");
				System.out.println(" Nonlinear Iteration: "+nonLinIter+"     Bmax: "+dfB.format(model. Bmax));
				System.out.println("     error: "+dfe.format(errNR)+ "     dB max: "+dfe.format(errFlux));
				//		System.out.println("Flux    error: "+dfe.format(err));
				System.out.println("-----------------------------------------------------");
				System.out.println();

			}

			SpMat Ks=new SpMat();

			Vect b=new Vect();


			model.setMagMat();


			Vect dv=new Vect();

			if(model.analysisMode>0)
			{


				if(model.eddyTimeIntegMode==0/* || step==0*/)
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

			}


			else
			{
				Ks=model.Hs.deepCopy();

				b=model.RHS.sub(model.HkAk);
			}

			resNR=b.norm();
			
			if(nonLinIter==0){
			resNR0=resNR;
			model.solver.resRef=resNR0;
			}

			errNR=resNR/resNR0;
			
			if(scaling)
			Ci=Ks.scale(b);

			
			L=Ks.ichol();
			
			double errCG=model.errCGmax*model.errCG_NR*errNR/model.errNRmax;
				
			if(mode==0)
				dA=model.solver.err0ICCG(Ks,L,b,model.errCGmax*model.errNRmax,iterMax);	
			else if(mode==1)
				dA=model.solver.ICCG(Ks,L, b,errCG,iterMax);

			if(model.solver.terminate) break;

	

			if(scaling)
			dA.timesVoid(Ci);

			x=x.add(dA);


			model.up=x.deepCopy();



				B1=model.getAllB();

				model.setSolution(x);	

				model.setB();	

				B2=model.getAllB();
				
			
			
				errFlux=model.getDiffMax(B1,B2);//+Math.abs(dA.el[dA.length-1]);



				if(fluxErr){
					errNR=errFlux;
					}
		
			


		}



		model.HpAp=model.HkAk.deepCopy();

		System.out.println();
		System.out.println("-----------------------------------------------------");
		System.out.println(" Nonlinear Iteration: "+nonLinIter+"     Bmax: "+dfB.format(model.Bmax));
		System.out.println("     error: "+dfe.format(errNR)+ "     dB max: "+dfe.format(errFlux));
		System.out.println("-----------------------------------------------------");
		System.out.println();
		
		System.out.println();
		System.out.println("=======================================================");
		System.out.println("Total number of ICCG iterations: "+model.solver.totalIter);
		System.out.println("Total number of Nonlinear iterations: "+totalNonlinIterx);
		
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
