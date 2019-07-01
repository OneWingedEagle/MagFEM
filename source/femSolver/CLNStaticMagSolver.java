package femSolver;

import static java.lang.Math.abs;

import fem.Calculator;
import fem.Model;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;


public class CLNStaticMagSolver{
	
	private Calculator calc;
	
	public SpMat magMat;
	public Vect RHS;

	public CLNStaticMagSolver(Model model){

		this.calc=new Calculator(model);
	}

	public void setMatrix(Model model){


		setMagMat(model);
		//setTmat(model);
		//annexTmat();


	}


	public Vect solve(Model model){
		
	SpMat L=new SpMat();

	Vect x=new Vect(model.numberOfUnknowns);

	model.solver.terminate(false);

	SpMat Ks=magMat.deepCopy();

	Vect Ci=Ks.scale(RHS);

		L=Ks.ichol(1.2);

		if(RHS.norm()/RHS.length>.01){
			x=model.solver.ICCG(Ks,L, RHS,model.errCGmax,model.iterMax);
			
		}else if(RHS.norm()>1e-6)
			//x=model.solver.ICCGr0max(Ks, L, RHS,model.errCGmax*1e-2,model.iterMax);
		x=model.solver.ICCG(Ks,L, RHS,model.errCGmax,model.iterMax);

		else
			x=new Vect(x.length);


		x.timesVoid(Ci);		

		return x;

}
	
	public void setMagMat(Model model){		


		double eps=1e-10,cPB=model.cpb; 
		double[][] H1=new double[model.nElEdge][model.nElEdge];

		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,ext=6;
	
		int[] nz=new int[model.numberOfUnknowns];

		magMat=new SpMat(model.numberOfUnknowns);

	

		for(int i=0;i<model.numberOfUnknownEdges;i++){

			magMat.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
		}


		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


				H1=this.calc.He(model,0,0,i,false,false,false,false);
			

				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;

					for(int k=0;k<model.nElEdge;k++){

						columnEdgeNumb=edgeNumb[k];



						//===== Periodic BC ================

						if(model.edge[columnEdgeNumb].aPBC ||model.edge[rowEdgeNumb].aPBC){

							cPB=model.edge[columnEdgeNumb].getnPBC()*model.edge[rowEdgeNumb].getnPBC();

							H1[j][k]*=cPB;					

						}	

						if(model.edge[columnEdgeNumb].edgeKnown  ){
						
							continue;
						}


						columnIndex=model.edgeUnknownIndex[columnEdgeNumb]-1;
						if(columnIndex>matrixRow) continue;
						m=util.search(magMat.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	
							if(abs(H1[j][k])>eps  ){	

								magMat.row[matrixRow].index[nz[matrixRow]]=columnIndex;

								magMat.row[matrixRow].el[nz[matrixRow]]=H1[j][k];


								nz[matrixRow]++;

								//===========================
								if(nz[matrixRow]==magMat.row[matrixRow].nzLength-1){
									magMat.row[matrixRow].extend(ext);
								}
								//===========================
							}
						}
						else{

							magMat.row[matrixRow].addToNz(H1[j][k],m);
						
						}
					}			
				}
			}
		}


		for(int i=0;i<model.numberOfUnknownEdges;i++){
			magMat.row[i].sortAndTrim(nz[i]);
		}


	}

	public void setRHS(Model model,double[] totalLoss){		

		RHS=new Vect(model.numberOfUnknowns);
		
		double[] loss=new double[1];

		int matrixRow=0, rowEdgeNumb;
		
		Vect elemV=new Vect(model.nElEdge);
		
		double heat=0;
		
		for(int ir=1;ir<=model.numberOfRegions;ir++){

			if(!model.region[ir].isConductor) continue;
			
			double conductivity=model.region[ir].getSigma().el[0];
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


				elemV=calc.elemVectPhiCoil(model,i,loss,1.);
				
				loss[0]*=conductivity;
				
				heat+=loss[0];
		
				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;

			
					//===========  right-hand side
				
					RHS.el[matrixRow]+=elemV.el[j];	

					
				}
			}
		}
		
		totalLoss[0]=heat;

	//	RHS.show();

	}

	public void setRHSCLN(Model model,double[] totalLoss){		

		RHS=new Vect(model.numberOfUnknowns);
		
		double[] loss=new double[1];

		int matrixRow=0, rowEdgeNumb;
		
		Vect elemV=new Vect(model.nElEdge);
		
		double heat=0;
		
		for(int ir=1;ir<=model.numberOfRegions;ir++){

			if(!model.region[ir].isConductor) continue;
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


				elemV=calc.elemVectCLN(model,i,loss);
		
				heat+=loss[0];
		
				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;

			
					//===========  right-hand side
				
					RHS.el[matrixRow]+=elemV.el[j];	

					
				}
			}
		}
		
		totalLoss[0]=heat;

	//	RHS.show();

	}




}
