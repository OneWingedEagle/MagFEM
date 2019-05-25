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
				//B1=model.m2d.getAllB();

				model.setSolution(x);	

				model.setB();	

				B2=model.getAllB();
				
			
			
				errFlux=model.getDiffMax(B1,B2);//+Math.abs(dA.el[dA.length-1]);



				if(fluxErr){
					errNR=errFlux;
					}
		
			


		}
		
		int nUnCur=model.numberOfUnknownCurrents;


		for(int k=0;k<nUnCur;k++){
			int nr=model.unCurRegNumb[k];

	 model.region[nr].currentp=model.region[nr].current;
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
