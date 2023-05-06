package femSolver;

import java.text.DecimalFormat;

import fem.Model;
import math.SpMat;
import math.Vect;
import math.util;


public class StaticNonlinearMagSolver{
	int stepNumb;
	boolean usePrev=false;
	int totalNonlinIter;
	boolean relaxedNR=true;

	public StaticNonlinearMagSolver(){	}


	
	
	public Vect solve(Model model,Vect x,boolean echo,int step){

	
		boolean old=false;
		if(old) return solveOld(model,x,echo,step);

		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;
		

		Vect Ci;
		SpMat L;
		

		int kx=0;

		
		
		Vect dA=new Vect(kx+model.numberOfUnknowns);
		
		model.setSolution(x);

		 x.zero();
			
		Vect[] B1,B2;

		model.solver.terminate(false);

		model.magMat.setRHS(model);
		
		double errNR=1;
		double fluxErr=1;
		
		double residual;
		boolean useFluxErr=false;

		int nonLinIter=0;
		model.errNRmax=1e-4; // not working reliably
		double convRatio=1e-2;

		double updatedConvCrot=1;
		
		double residual0=model.RHS.norm();
			
	
		for( nonLinIter=0; errNR>model.errNRmax && nonLinIter<nonLinIterMax;nonLinIter++)

			
		{
			totalNonlinIter++;

			if(echo)
			{
				System.out.println();
				System.out.println();
				System.out.println("    nonlinear Iteration: "+nonLinIter+
						"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
				System.out.println(" Flux    error: "+dfe.format(fluxErr));
				System.out.println();
				System.out.println("    updatedConvCrot: "+updatedConvCrot);
			}

			SpMat Ks=new SpMat();

			Vect b=new Vect();
			

			model.setMagMat();
		

				Ks=model.Hs.deepCopy();
				//Ks.shownz();

				b=model.RHS.sub(model.HkAk);
				
				
				residual=b.norm();
				if(useFluxErr)
					errNR=fluxErr;
				else
					errNR=residual/residual0;
	
				updatedConvCrot=errNR*convRatio;
			
			
			Ci=Ks.scale(b);
			
			L=Ks.ichol();
			
			if(nonLinIter==0)
			model.solver.resRef=b.norm();

			
			if(b.abs().sum()>1e-11){
			//	dA=model.solver.ICCG(Ks,L, b,model.errCGmax,iterMax);

				dA=model.solver.ICCG(Ks,L, b,updatedConvCrot,iterMax);	
			//dA=model.solver.ICCG(Ks,L, b,model.errCGmax*1e-3,iterMax);	
				//if(b.abs().max()>1e-6)
			//dA=model.solver.err0ICCG(Ks,L,b,errNR*iccgConvRatio,iterMax);	

			}

			if(model.solver.terminate) break;


			dA.timesVoid(Ci);
		
			x=x.add(dA);	
				
			B1=model.getAllB();

			model.setSolution(x);	
			
		
			B2=model.getAllB();
			
			fluxErr=model.getFluxErrSquared(B1,B2);

			
	

		}
		

	
		
		System.out.println();
		System.out.println("-----------------------------------------------------");
		System.out.println(" Nonlinear Iteration: "+nonLinIter+
				"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
		System.out.println("Flux    error: "+dfe.format(fluxErr));
		System.out.println("-----------------------------------------------------");
		System.out.println();
		
		System.out.println();
		System.out.println("=======================================================");
		System.out.println("Total number of ICCG iterations: "+model.solver.totalIter);
		System.out.println("Total number of Nonlinear iterations: "+totalNonlinIter);
		
		
		System.out.println();
	
		System.out.println(" Final NR residual: "+residual0);
		
		System.out.println("=======================================================");
		System.out.println();
		
		return x;

	}



	public Vect solveOld(Model model,Vect x,boolean echo,int step){

		

		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int iterMax=model.iterMax;
		int nonLinIterMax=model.nonLinIterMax;

		Vect Ci;
		SpMat L;
		Vect dA=new Vect(model.numberOfUnknowns);

		Vect[] B1,B2;

		model.solver.terminate(false);
		
		double errNR=1;
		double fluxErr=1;
		
		double residual;
		boolean useFluxErr=false;

		int nonLinIter=0;
		model.errNRmax=1e-2; // not working reliably
		double convRatio=1e-6;

		double updatedConvCrot=1;
		
		double residual0=model.RHS.norm();
			model.solver.resRef=residual0;
		
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
				System.out.println(" Flux    error: "+dfe.format(fluxErr));
				System.out.println();
			}

			SpMat Ks=new SpMat();

			Vect b=new Vect();

			model.setMagMat();
		

				Ks=model.Hs.deepCopy();

				b=model.RHS.sub(model.HkAk);
		

				residual=b.norm();


				
				if(useFluxErr)
					errNR=fluxErr;
				else
					errNR=residual/residual0;
					
				

					updatedConvCrot=errNR*convRatio;
					
					Ci=Ks.scale(b);
					
					L=Ks.ichol();
					
					
					if(b.abs().sum()>1e-11){
					//	dA=model.solver.ICCG(Ks,L, b,model.errCGmax,iterMax);

						dA=model.solver.ICCG(Ks,L, b,updatedConvCrot,iterMax);	
					//dA=model.solver.ICCG(Ks,L, b,model.errCGmax*1e-3,iterMax);	
						//if(b.abs().max()>1e-6)
					//dA=model.solver.err0ICCG(Ks,L,b,errNR*iccgConvRatio,iterMax);	

					}



			if(model.solver.terminate) break;


			dA.timesVoid(Ci);
			util.pr(x.length+"  "+dA.length);

			x=x.add(dA);

			B1=model.getAllB();

			model.setSolution(x);	

			model.setB();	
			
			B2=model.getAllB();
			
					
			fluxErr=model.getDiffMax(B1,B2);


		}
		

	
		
		System.out.println();
		System.out.println("-----------------------------------------------------");
		System.out.println(" Nonlinear Iteration: "+nonLinIter+
				"     Bmax: "+dfB.format(model.Bmax)+ "     error: "+dfe.format(errNR));
		System.out.println("Flux    error: "+dfe.format(fluxErr));
		System.out.println("-----------------------------------------------------");
		System.out.println();
		
		System.out.println();
		System.out.println("=======================================================");
		System.out.println("Total number of ICCG iterations: "+model.solver.totalIter);
		System.out.println("Total number of Nonlinear iterations: "+totalNonlinIter);
		
		
		System.out.println();
	
		System.out.println(" Final NR residual: "+residual0);
		
		System.out.println("=======================================================");
		System.out.println();


		return x;

	}

}
