package fem;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import fem.Network.Elem;
import fem.Network.ElemType;
import math.Complex;
import math.Mat;
import math.SpMat;
import math.SpMatComp;
import math.SpVect;
import math.SpVectComp;
import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class MagMatAssembler {

	private Calculator calc;

	public  MagMatAssembler(){}

	public  MagMatAssembler(Model model){

		this.calc=new Calculator(model);
	}




	public void setMagMat(Model model){

		setReactMat(model);


		if(model.analysisMode>0)
			setConductMat(model);

		if(model.analysisMode>1 && model.dim==3){

			setA_Phi_couplingMat( model);

			if(!model.AC){

				for(int i=0;i<model.numberOfVarNodes;i++){

					SpVect spv=model.Ps.row[i].deepCopy();

					spv=spv.augh(model.Qs.row[i].times(model.dt));
					model.Hs.row[i+model.numberOfUnknownEdges]=spv.deepCopy();
				}

			}
		}
		
		annexNetworkCoupling(model);

	}

	public void setReactMat(Model model){		

		double cPB=model.cpb; 
		boolean nonLinear;
		double[][] H1=new double[model.nElEdge][model.nElEdge];

		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,nBH=0,nLam=0,ext=6;

		int[] nz=new int[model.numberOfUnknowns];

		model.Hs=new SpMat(model.numberOfUnknowns);

		for(int i=0;i<model.numberOfUnknownEdges;i++){

			model.Hs.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
		}




		if(model.eddyTimeIntegMode==1){
			if(model.RHS!=null)
				model.bT=model.RHS.deepCopy();
			else
				model.bT=null;
		}

		//model.RHS=new Vect(model.numberOfUnknowns);
		model.HkAk= new Vect(model.numberOfUnknowns);

		model.RHS_boundary= new Vect(model.numberOfUnknowns);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			nBH=model.region[ir].BHnumber;
			nLam=model.region[ir].lamBNumber;
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				nonLinear=model.element[i].isNonlin();

				H1=this.calc.He(model,nBH,nLam,i,nonLinear,false,false,false);


				int[] edgeNumbs=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumbs[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;

					for(int k=0;k<model.nElEdge;k++){

						columnEdgeNumb=edgeNumbs[k];



						//===== Periodic BC ================

						if(model.edge[columnEdgeNumb].aPBC ||model.edge[rowEdgeNumb].aPBC){

							cPB=model.edge[columnEdgeNumb].getnPBC()*model.edge[rowEdgeNumb].getnPBC();

							H1[j][k]*=cPB;
							if(nonLinear){
								model.H2[j][k]*=cPB;
							}


						}	

						//===========================

						//===========  right-hand side 1

						double Ak=0;
						if(!model.edge[columnEdgeNumb].hasPBC()) Ak= model.edge[columnEdgeNumb].A;
						else {

							Ak= model.edge[model.edge[columnEdgeNumb].map].A;
						}

					

						if(model.edge[columnEdgeNumb].edgeKnown  ){

							model.RHS_boundary.el[matrixRow]+=H1[j][k]*Ak;
		
							continue;
						}else
							model.HkAk.el[matrixRow]+=H1[j][k]*Ak;


						//=======================

						if(nonLinear){

							H1[j][k]+= model.H2[j][k];
						}


						columnIndex=model.edgeUnknownIndex[columnEdgeNumb]-1;
						if(columnIndex>matrixRow) continue;

						m=util.search(model.Hs.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	



							model.Hs.row[matrixRow].index[nz[matrixRow]]=columnIndex;

							model.Hs.row[matrixRow].el[nz[matrixRow]]=H1[j][k];

							nz[matrixRow]++;

							//===========================
							if(nz[matrixRow]==model.Hs.row[matrixRow].nzLength-1){
								model.Hs.row[matrixRow].extend(ext);
							}
							//===========================

						}
						else{

							model.Hs.row[matrixRow].addToNz(H1[j][k],m);

						}


					}			
				}
			}
		}


		for(int i=0;i<model.numberOfUnknownEdges;i++){
			model.Hs.row[i].sortAndTrim(nz[i]);
		}


	}




	public void setConductMat(Model model){		

		double eps=1e-10,cPB=model.cpb; 
		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,ext=6;

		int[] nz=new int[model.numberOfUnknowns];

		model.Ss=new SpMat(model.numberOfUnknowns);
		for(int i=0;i<model.numberOfUnknownEdges;i++){

			model.Ss.row[i]=new SpVect(model.numberOfUnknowns,model.nEdEd);
		}





		if(model.eddyTimeIntegMode==1){
			if(model.RHS!=null)
				model.bT=model.RHS.deepCopy();
			else
				model.bT=null;
		}




		for(int ir=1;ir<=model.numberOfRegions;ir++){


			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


				if(!model.element[i].isConductor()) continue;


				this.calc.He(model,0,0,i,false,true,false,false);




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

							model.H3[j][k]*=cPB;

						}	

						//===========================

						if(model.edge[columnEdgeNumb].edgeKnown){

							continue;
						}


						columnIndex=model.edgeUnknownIndex[columnEdgeNumb]-1;


						if(columnIndex>matrixRow) continue;
						m=util.search(model.Ss.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	


							if(abs(model.H3[j][k])>eps  ){	

								model.Ss.row[matrixRow].index[nz[matrixRow]]=columnIndex;

								model.Ss.row[matrixRow].el[nz[matrixRow]]=model.H3[j][k];

								nz[matrixRow]++;

								//===========================
								if(nz[matrixRow]==model.Ss.row[matrixRow].nzLength-1){
									model.Ss.row[matrixRow].extend(ext);
								}
								//===========================
							}
						}
						else{


							model.Ss.row[matrixRow].addToNz(model.H3[j][k],m);
						}
					}


				}			
			}
		}



		for(int i=0;i<model.numberOfUnknownEdges;i++){
			model.Ss.row[i].sortAndTrim(nz[i]);
		}


		if(model.eddyTimeIntegMode==0 || model.eddyTimeIntegMode==1){

			if(!model.AC)
				model.Ss.times(1.0/model.dt);


		}



	}


	public void setA_Phi_couplingMat(Model model){


		model.Ps=getPs(model);


		model.Qs=getQs(model);



	}

	public void setRHS(Model model){


		if(model.phiCoils!=null){	
			double[] losses=new double[1];
			setRHS(model,losses);

			return;	
		}

		setRHS(model,true);
		
		if(model.hasBunif){
			
			TimeFunction tf=model.timeFunctions[model.unifBTimeId];
			double time=model.getCurrentTime();
	
			double factor=tf.getValue(time);
		
			model.RHS=model.RHS.sub(model.RHS_boundary.times(factor));
			model.scaleKnownEdgeAL(factor);

		}

	}

	public void setRHS(Model model,boolean includeM){		

		model.RHS=new Vect(model.numberOfUnknowns);

		int matrixRow=0, rowEdgeNumb;

		double[] Cj=new double[model.nElEdge];

		boolean hasM, hasJ;

		Vect J=null;


		for(int ir=1;ir<=model.numberOfRegions;ir++){


			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				hasJ=model.element[i].hasJ();

				hasM=model.element[i].hasM();

				if(!hasJ && !hasM) continue;


				if(hasJ){
					J=model.element[i].getJ();

				}


				this.calc.He(model,0,0,i,false,false,hasJ,hasM);


				if(hasJ ){
					for(int j=0;j<model.nElEdge;j++){

						if(model.dim==2)
							Cj[j]=J.el[2]*model.Cj2d[j];
						else{
							Cj[j]=J.dot(model.Cj[j]);

						}


					}

				}


				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;


					//===========  right-hand side
					if( hasM && includeM  ){					
						model.RHS.el[matrixRow]+=model.C[j];	
					}

					if(hasJ ){
						model.RHS.el[matrixRow]+=Cj[j];	

					}


				}
			}
		}


	}
	
	public void setRHS(Model model,double[] totalLoss){		

		if(model.seperateCoil){
			setRHSSep(model);
			return;
		}
		model.RHS=new Vect(model.numberOfUnknowns);

		double[] loss=new double[1];

		int matrixRow=0, rowEdgeNumb;

		Vect elemV=new Vect(model.nElEdge);

		
		double[] coilLosses=new double[model.phiCoils.length];

		for(int ir=1;ir<=model.numberOfRegions;ir++){


			int coilIndex=model.coilInicesByRegion[ir];
			if(coilIndex<0) continue;

			double conductivity=model.phiCoils[coilIndex].getConductivity();
			double current=model.phiCoils[coilIndex].getCurrent();

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

	
				elemV=calc.elemVectPhiCoil(model,i,loss, 1.);

				elemV.timesVoid(current);

				if(conductivity>0)
				coilLosses[coilIndex]+=loss[0]*current*current/conductivity;

				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<model.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnown ) continue;


					matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;


					//===========  right-hand side

					model.RHS.el[matrixRow]+=elemV.el[j];	


				}
			}
		}

		//model.RHS.show();
		totalLoss[0]=0;
		for(int ic=0;ic<model.phiCoils.length;ic++){
			
			double current=model.phiCoils[ic].getCurrent();
			double res=coilLosses[ic]/(current*current);
			model.phiCoils[ic].setResistance(res);
			util.pr("Coil Resistance= "+res);
			totalLoss[0]+=coilLosses[ic];

		}
	}


	private  SpMat getQs(Model model){

		double eps=1e-12;
		Mat He;
		int m,nodeNumber,columnIndex=0,matrixRow;
		int[] nz=new int[model.numberOfVarNodes];

		SpMat Qs=new SpMat(model.numberOfVarNodes, model.numberOfVarNodes,model.nNodNod);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				if(!model.element[i].isConductor()) continue;

				He=this.calc.Qe(model,i);


				int[] vertNumb=model.element[i].getVertNumb();

				for(int j=0;j<model.nElVert;j++){



					//	if(!model.node[vertNumb[j]].isPhiVar() || model.node[vertNumb[j]].isPhiKnown() ) continue;

					matrixRow=model.nodeVarIndex[vertNumb[j]]-1;
					if(matrixRow<0) continue;
					for(int k=0;k<model.nElVert;k++){

						nodeNumber=vertNumb[k];
						columnIndex=model.nodeVarIndex[nodeNumber]-1;			
						if(columnIndex<0) continue;
						if(columnIndex>matrixRow) continue;


						m=util.search(Qs.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	

							if(abs(He.el[j][k])>eps ){
								Qs.row[matrixRow].index[nz[matrixRow]]=columnIndex;
								Qs.row[matrixRow].el[nz[matrixRow]++]=He.el[j][k];
							}

						}

						else{

							Qs.row[matrixRow].addToNz(He.el[j][k],m);
						}

					}			
				}
			}
		}

		Qs.sortAndTrim(nz);

		return Qs;

	}


	private  SpMat getPs(Model model){

		double ePs=1e-12;
		double[][] He;

		int ext=6;

		int m,edgeNumber,nodeNumber,columnIndex,matrixRow;
		int[] nz=new int[model.numberOfVarNodes];

		SpMat Ps=new SpMat(model.numberOfVarNodes, model.numberOfUnknownEdges,model.nNodEd);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				if(!model.element[i].isConductor()) continue;

				He=this.calc.Pe(model,i);


				int[] vertNumb=model.element[i].getVertNumb();
				int[] edgeNumb=model.element[i].getEdgeNumb();


				for(int j=0;j<model.nElVert;j++){
					nodeNumber=vertNumb[j];
					if(!model.node[nodeNumber].isPhiVar() || model.node[nodeNumber].isPhiKnown()  ) continue;
					matrixRow=model.nodeVarIndex[nodeNumber]-1;
					if(matrixRow<0) continue;
					for(int k=0;k<model.nElEdge;k++){		
						edgeNumber=edgeNumb[k];
						if(model.edge[edgeNumber].edgeKnown) continue;

						columnIndex=model.edgeUnknownIndex[edgeNumber]-1;
						if(columnIndex<0) continue;
						m=util.search(Ps.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	
							if(abs(He[j][k])>ePs ){
								Ps.row[matrixRow].index[nz[matrixRow]]=columnIndex;
								Ps.row[matrixRow].el[nz[matrixRow]++]=He[j][k];
								//===========================
								if(nz[matrixRow]==Ps.row[matrixRow].nzLength-1){
									Ps.row[matrixRow].extend(ext);
								}
								//===========================
							}

						}

						else{

							Ps.row[matrixRow].addToNz(He[j][k],m);
						}

					}			
				}
			}
		}
		Ps.sortAndTrim(nz);




		return Ps;
	}


	private double hatFunc(double tt,int px,double W,double period){
		double result=0;


		tt=(tt-px*W)%period;

		if(tt<0) tt+=period;


		if(tt<W || period-tt<W){
			if(tt<W && period-tt>W)
				result=(W-tt)/W;
			else if(tt>W && period-tt<W)
				result=1-(period-tt)/W;
		}

		return result;

	}

	public void setMagBCUniform(Model model){

		if(model.hasBunif){
	
		
			double Bx=model.unifB.el[0];
			double By=model.unifB.el[1];
			double Bz=0;
			if(model.dim==3)
				Bz=model.unifB.el[2];
			double Ax,Ay,Az;
			double x,y,z;
			Vect A=new Vect(3);
			for(int i=1;i<=model.numberOfEdges;i++){

				if(!model.edge[i].edgeKnown) continue;
		
				Vect edgeVect=model.edge[i].node[1].getCoord().sub(model.edge[i].node[0].getCoord());

				Vect center=model.edge[i].node[1].getCoord().add(model.edge[i].node[0].getCoord()).times(.5);
				x=center.el[0];
				y=center.el[1];
		
				Az=y*Bx-x*By;
				if(model.dim==3){
					z=center.el[2];
					Ax=x*By;
					Ay=y*Bz;
					A=new Vect(Ax,Ay,Az);
				}else{
					A=new Vect(0,0,Az);
				}
			
				double a=edgeVect.dot(A);
				
				model.edge[i].setKnownA(a);

			}
		}


	}
	
 public void setRHSSep(Model model){		
		model.RHS=new Vect(model.numberOfUnknowns);
		
		model.network.solveStatic(model);

		
		for(int ic=0;ic<model.phiCoils.length;ic++){

			if(model.phiCoils[ic].load!=null) {
				
				double factor=model.phiCoils[ic].getCurrent();
				
				model.RHS=	model.RHS.add(model.phiCoils[ic].load.times(factor));
			}

		}

	}

	public void annexNetworkCoupling(Model model){	
		
		Network network=model.network;
		
		int nCurrents=network.no_unknown_currents;
		int row=model.numberOfUnknowns-nCurrents;

		Elem elem;
		
		SpVect[] couplings=new SpVect[nCurrents];
		
		for (int j = 0; j < network.indep_elems.length; j++)
		{
			elem =network.indep_elems[j];

			if (elem.type == ElemType.CPS) continue;
			

			int jx=elem.unknown_seq_no;


			Vect v=new Vect(model.numberOfUnknowns);
			for (int k = 0; k < network.numElements; k++){
				double f=network.tiesetMat.el[j][k];
				if(f!=0){
					Elem elem2=network.elems[k];
					if(elem2.type==ElemType.FEM){
						Vect load=model.phiCoils[elem2.fem_index].load;
						v=v.add(load.times(f));
					}
				
				}
			}
			
			for (int k = 0; k <nCurrents;k++)
				v.el[v.length-1+k]=network.PRPt.el[jx][k];
			
			couplings[j]=new SpVect(v.times(1./model.dt));
			
			model.Hs.row[row+jx]=couplings[j].deepCopy();
			
		}


		

	}

}
