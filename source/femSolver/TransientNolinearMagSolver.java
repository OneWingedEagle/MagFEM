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

				else  if(model.eddyTimeIntegMode<=-2){


					 Ks=getCircuitHs(model,step);


					int nNeut=model.nNeutral;
					int nUnCur=model.numberOfUnknownCurrents;


					Vect vk=model.getUnknownA();

					Vect vk2=new Vect(model.numberOfUnknowns);

					for(int j=0;j<vk.length;j++)
						vk2.el[j]=vk.el[j];

					
					Vect vp=model.getUnknownAp();

					Vect vp2=new Vect(model.numberOfUnknowns);

					for(int j=0;j<vp.length;j++)
						vp2.el[j]=vp.el[j];

					dv=vp2.sub(vk2);

					model.RHS=model.RHS.add(model.Ss.smul(dv));


					
					for(int i=0;i<nUnCur;i++){
						int nr=model.unCurRegNumb[i];

						double vprev=model.region[nr].terminalVoltagep;
						if(	this.stpNumb==0)
							vprev=model.region[nr].terminalVoltage;
						
				
						double ip=model.region[nr].currentp;

		
						if(model.eddyTimeIntegMode==-3){
														
							double cf=this.theta;
							
							if(model.HpAp!=null){
								model.RHS=model.RHS.sub(model.HpAp.times(1-cf));
									}
			
					
						model.RHS.el[model.Hs.nRow-nUnCur+i-nNeut]=cf*(model.lastRows[i].dot(vp2)
						+(((1-cf)*model.region[nr].getWireRes()*ip
								+(1-cf)*model.vNeutral
								-(cf*model.region[nr].terminalVoltage+(1-cf)*vprev)))*model.dt			
						-this.coilInduct*ip)/model.height;
											

						model.RHS=model.RHS.sub(model.lastRows[i].times((1-cf)*model.region[nr].currentp).vectForm());

						}
						else{
				
							model.RHS.el[model.Hs.nRow-nUnCur+i-nNeut]=model.lastRows[i].dot(vp2)
							+((model.vNeutral-model.region[nr].terminalVoltage)*model.dt
							-this.coilInduct*ip)/model.height;

						}

					}

					SpMat Q=new SpMat(model.RHS.length,model.RHS.length);
					for(int i=0;i<model.RHS.length;i++){
						if(i<vk.length)
							Q.row[i]=null;
						else
							Q.row[i]=Ks.row[i].deepCopy();
					}

			
					Vect v3=new Vect(model.numberOfUnknowns);
					for(int j=0;j<vk.length;j++){
						v3.el[j]=vk.el[j];
					}

					for(int i=0;i<nUnCur;i++){
						int nr=model.unCurRegNumb[i];
						v3.el[vk.length+i]=model.region[nr].current;

					}
					if(model.nNeutral>0)
					v3.el[vk.length+nUnCur]=model.vNeutral;
					
			
			
					model.RHS=model.RHS.sub(Q.smul(v3));


					if(model.eddyTimeIntegMode==-3){
				
						b=model.RHS.sub(model.HkAk.times((this.theta)));

					}
					else{

						b=model.RHS.sub(model.HkAk);
						
					}
			

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
			
			

			int nUnCur=model.numberOfUnknownCurrents;
			int nNeut=model.nNeutral;


			x=x.add(dA);
			

			model.up=x.deepCopy();


			if(model.eddyTimeIntegMode<=-2){

				
				Vect vk=model.getUnknownA();

				Vect vk2=new Vect(model.numberOfUnknowns);

				for(int j=0;j<vk.length;j++)
					vk2.el[j]=vk.el[j];

				Vect vp=model.getUnknownAp();

				Vect vp2=new Vect(model.numberOfUnknowns);

				for(int j=0;j<vp.length;j++)
					vp2.el[j]=vp.el[j];

				dv=vp2.sub(vk2);

				if(nNeut>0)
					model.vNeutral=x.el[x.length-nNeut];

				for(int k=0;k<nUnCur;k++){
					int nr=model.unCurRegNumb[k];
					model.region[nr].inducedVoltage=model.lastRows[k].dot(dv)/model.dt*model.height;

					model.region[nr].current=x.el[x.length-nUnCur-nNeut+k];

				}
			}
			
			
			
			
			B1=model.getAllB();

			model.setSolution(x);	

			model.setB();	
			
			B2=model.getAllB();
			
			errFlux=model.getDiffMax(B1,B2);

			
	

		}
		

		int nUnCur=model.numberOfUnknownCurrents;


		for(int k=0;k<nUnCur;k++){
			int nr=model.unCurRegNumb[k];

	 model.region[nr].currentp=model.region[nr].current;
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

		int nUnCur=model.numberOfUnknownCurrents;

		int nNeut=model.nNeutral;
		
	double cf=1;
	
		if(model.eddyTimeIntegMode==-3) cf=this.theta;
		
		
		int nr=Hs.nRow-nUnCur-nNeut;
		

		for(int i=0;i<nr;i++)
		 Hs.row[i]=Hs.row[i].times(cf);
		
		

		for(int i=0;i<nUnCur;i++){


			int jr=Hs.nRow-nUnCur-nNeut+i;

			Hs.row[jr]=model.lastRows[i].times(cf);
			int nz=Hs.row[jr].nzLength;
			Hs.row[jr].extend(1);

	
			double R=model.region[model.unCurRegNumb[i]].getWireRes();

			Hs.row[jr].el[nz]=-cf*(cf*model.dt*R+this.coilInduct)/model.height;

			Hs.row[jr].index[nz]=jr;
		}

		int jr=Hs.nRow-1;

		
		if(nNeut>0){
			Hs.row[jr]=new SpVect(Hs.nRow,1+nUnCur);		

			for(int i=0;i<nUnCur;i++){

				Hs.row[jr].el[i]=-cf*model.dt/model.height;

				Hs.row[jr].index[i]=jr-nUnCur+i;
			}

			Hs.row[jr].el[nUnCur]=cf*model.dt/model.Rg/model.height;

			Hs.row[jr].index[nUnCur]=jr;

		}

		if(model.analysisMode>0){
			Hs.addSmaller(model.Ss);

		}
		
		return Hs;

	}
	



}
