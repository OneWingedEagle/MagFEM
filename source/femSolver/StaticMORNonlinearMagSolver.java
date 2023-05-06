package femSolver;

import java.text.DecimalFormat;

import fem.Model;
import math.Mat;
import math.MatSolver;
import math.SpMat;
import math.Vect;
import math.util;


public class StaticMORNonlinearMagSolver{
	int stepNumb;
	boolean usePrev=false;
	int totalNonlinIter;
	boolean relaxedNR=true;

	public StaticMORNonlinearMagSolver(){	}


	
	
	public Vect solve(Model model,Mat Phi,Vect x,boolean echo,int step){

		Mat PhiT=Phi.transp();
		Vect[] B1,B2;
		MatSolver ms=new MatSolver();
		
		DecimalFormat dfB=new DecimalFormat("0.0000");
		DecimalFormat dfe=new DecimalFormat("0.0E00");


		int nonLinIterMax=1;//model.nonLinIterMax;

		model.solver.terminate(false);

		model.magMat.setRHS(model);
		
		double errNR=1,resNR=0,resNR0=1;
		double errFlux=1.0;

		int nonLinIter;
		model.errNRmax=1e-6;

		for( nonLinIter=0; (errFlux>1e-2 || errNR>model.errNRmax) && nonLinIter<nonLinIterMax;nonLinIter++)

			
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


			Vect b=new Vect();

			
			model.setMagMat();
		
				b=model.RHS.sub(model.HkAk);
				
				Mat Kr=PhiT.mul(model.Hs.smul(Phi));
	
				Vect br=PhiT.mul(b);
						
			if(nonLinIter==0){
				resNR0=b.norm();
				model.solver.resRef=resNR0;
			}

			resNR=br.norm();
			
		
			errNR=resNR/resNR0;

		
			Vect dAr;

			if(br.norm()>1e-10){

				dAr=ms.gaussel(Kr, br);
			}
			else
				dAr=new Vect(br.length);	
			


			if(model.solver.terminate) break;

		     Vect dA=Phi.mul(dAr);

			x=x.add(dA);
			

			
			
			
			B1=model.getAllB();

			model.setSolution(x);	

			model.setB();	
			
			B2=model.getAllB();
			
			errFlux=model.getDiffMax(B1,B2);

			
	

		}
		

	
		
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

}