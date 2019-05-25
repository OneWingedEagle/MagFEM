package femSolver;

import static java.lang.Math.abs;

import java.util.Arrays;

import fem.Calculator;
import fem.Model;
import fem.Network;
import fem.Network.ElemType;
import fem.PhiCoil;
import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;


public class StaticElectricSolver{
	int stepNumb;

	public SpMat phi_matrix;
	public SpMat t_matrix;
	public SpMat conductiveMat;
	public Vect RHS;
	private Calculator calc;
	private int[] phiVarIndex;
	
	public int numberOfUnknownPhis,numberOfUnknowns,nCurrents;

	
	public StaticElectricSolver(){

	}

	public void setMatrix(Model model){

		this.calc=new Calculator(model);

		setPhiMat(model);
		setTmat(model);
		annexTmat();

		
	}

	public Vect solve(Model model ){

		Vect x=new Vect(numberOfUnknowns);

		if(numberOfUnknowns==0 || RHS.norm()<1e-12) return x;

		SpMat L=new SpMat();


		model.solver.terminate(false);

	//	RHS.show();
		SpMat Ks=conductiveMat.deepCopy();

			//Ks.diagSym().show();

		Vect Ci=Ks.scale(RHS);
	//	util.pr(Ks.maxAbs());
	////util.pr("RHS --------------------------> "+RHS.norm());


		L=Ks.ichol();
		
		double errMax=1e-11;
		
		

	//	io.Console.redirectOutput(model.main.gui.iccgArea);

		util.pr("\n"+this.getClass().getName()+" starts...\n");

		x=model.solver.ICCG(Ks,L, RHS,errMax,model.iterMax);
		
		util.pr("\n"+this.getClass().getName()+" ends.\n");


	//	io.Console.redirectOutput(model.main.gui.dataArea);
		//util.pr(model.main.gui.iccgArea.getText());

		x.timesVoid(Ci);

		//x.show();
		//util.pr("x --------------------------> "+x.norm());
		

		return x;

	}

	public  void setPhiMat(Model model){


		int[] coilIndices=new int[model.numberOfRegions+1];

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			coilIndices[ir]=-1;

		for(int ic=0;ic<model.phiCoils.length;ic++){
			int nr=model.phiCoils[ic].regNo;
			coilIndices[nr]=ic;

		}
		
		double eps=1e-12;

		int ext=10;
		Mat Ke;
		int m,nodeNumber,columnIndex=0,matrixRow;
		int[] nz=new int[numberOfUnknownPhis];

		phi_matrix=new SpMat(numberOfUnknownPhis, numberOfUnknownPhis,model.nNodNod);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			if(coilIndices[ir]<0) continue;
		//	if(!model.region[ir].isConductor) continue;

			double conductivity=1;//model.phiCoils[coilIndices[ir]].conductivity;
			if(coilIndices[ir]>0){
				double conductivity1=model.phiCoils[0].conductivity;
				double conductivity2=model.phiCoils[coilIndices[ir]].conductivity;
			if(conductivity1>0) conductivity=conductivity2/conductivity1;
			}
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				//if(!model.element[i].isConductor()) continue;

				Ke=this.calc.elemPhiMat(model,i);
		
				Ke=Ke.times(conductivity);

				int[] vertNumb=model.element[i].getVertNumb();

				for(int j=0;j<model.nElVert;j++){

					if(!model.node[vertNumb[j]].isPhiVar() || model.node[vertNumb[j]].isPhiKnown() ) continue;

					matrixRow=phiVarIndex[vertNumb[j]];
					for(int k=0;k<model.nElVert;k++){
						nodeNumber=vertNumb[k];
						columnIndex=phiVarIndex[nodeNumber];								

						if(columnIndex==-1 || columnIndex>matrixRow) continue;	
				
						m=util.search(phi_matrix.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	

							if(abs(Ke.el[j][k])>eps ){
		
								//===========================
								if(nz[matrixRow]==phi_matrix.row[matrixRow].nzLength-1){
									phi_matrix.row[matrixRow].extend(ext);
								}
								//===========================

								phi_matrix.row[matrixRow].index[nz[matrixRow]]=columnIndex;
								phi_matrix.row[matrixRow].el[nz[matrixRow]++]=Ke.el[j][k];
							}

						}

						else{

							phi_matrix.row[matrixRow].addToNz(Ke.el[j][k],m);

						}

					}			
				}
			}
		}

		phi_matrix.sortAndTrim(nz);
		//phi_matrix.show();
	}
	
	public  void setTmat(Model model){

		Network network=model.network;
		
			t_matrix=new SpMat(nCurrents,numberOfUnknowns);
	
			for(int n=0;n<network.indep_elems.length;n++){
				int i=network.indep_elems[n].unknown_seq_no;
				if(i==-1) continue;
			t_matrix.row[i]=new SpVect(numberOfUnknowns,model.phiCoils.length+i+1);

			int jx=0;
			int ic=0;
			for(int j=0;j<network.numElements;j++){
				if(network.elems[j].type==ElemType.FEM){
				if(network.tiesetMat.el[n][j]!=0){
				t_matrix.row[i].index[jx]=numberOfUnknownPhis-model.phiCoils.length+ic;
				PhiCoil coil=model.phiCoils[network.elems[j].fem_index];
	
				double turns=coil.getNumTurns();

				t_matrix.row[i].el[jx]=-network.tiesetMat.el[n][j]*turns;
				jx++;
				}
				ic++;			
				}
				
			}
			for(int j=0;j<=i;j++){
		//	if(p>=i || p==-1) continue;
			t_matrix.row[i].index[model.phiCoils.length+j]=numberOfUnknownPhis+j;

			t_matrix.row[i].el[model.phiCoils.length+j]=-network.PRPt.el[n][j];
			}
			}
			//t_matrix.matForm().show();;

			//t_matrix.shownz();;
	}
	

	private void annexTmat(){

		conductiveMat=new SpMat(numberOfUnknowns);
		
		for(int i=0;i<numberOfUnknownPhis;i++){
			conductiveMat.row[i]=new SpVect(numberOfUnknowns,phi_matrix.row[i].nzLength);

			conductiveMat.row[i].index=phi_matrix.row[i].index;
			conductiveMat.row[i].el=phi_matrix.row[i].el;
		}

		for(int i=0;i<t_matrix.nRow;i++){
			int matrixRow=i+numberOfUnknownPhis;

			conductiveMat.row[matrixRow]=new SpVect(numberOfUnknowns,t_matrix.row[i].nzLength);
		for(int j=0;j<t_matrix.row[i].nzLength;j++){
			
			conductiveMat.row[matrixRow].index=t_matrix.row[i].index;
			
			conductiveMat.row[matrixRow].el=t_matrix.row[i].el;
		
		}
	}
		//t_matrix.matForm().show();;
		//conductiveMat.size();

	}


	public void setRHS0(Model model){		


		RHS=new Vect(numberOfUnknowns);
		Network network=model.network;

		boolean byCPS=true;
		int rowIndex=-1;
		for(int j=0;j<network.indep_elems.length;j++){
			if(network.indep_elems[j].type==ElemType.VPS){
				 rowIndex=network.indep_elems[j].unknown_seq_no+numberOfUnknownPhis;

				byCPS=false;
				 break;
			}
		}
			
		if(!byCPS &&rowIndex>=0){

			double vps_volatge=1;
			RHS.el[rowIndex]=-vps_volatge;
		}
		
		if(byCPS){
			
			Vect indepCurrents=new Vect(network.indep_elems.length);
	
			for (int j = 0; j<network.indep_elems.length; ++j){
				if(network.indep_elems[j].type==ElemType.CPS){
					int time_id=network.indep_elems[j].time_id;

					double value=0;
					if(model.timeFunctions!=null && time_id>0)
					  value=model.timeFunctions[time_id].getValue(0);
					indepCurrents.el[j]=value;
				}
			}
			
			Vect rhsCurrents=network.PRPt.mul(indepCurrents);

			
			for (int j = 0; j<network.indep_elems.length; ++j){
				if(network.indep_elems[j].type!=ElemType.CPS){
					
					 rowIndex=this.numberOfUnknownPhis+network.indep_elems[j].unknown_seq_no;


					RHS.el[rowIndex]+=rhsCurrents.el[j];
				}
			}
		

			
			Vect allCurrents=network.tiesetMat.transp().mul(indepCurrents);
			
			for (int j = 0; j<network.numElements; ++j)
			{
				if(network.elems[j].type!=ElemType.CPS) {
					network.elems[j].I=allCurrents.el[j];
				}
			}
			

			for(int j=0;j<network.elems.length;j++){
				if(network.elems[j].type==ElemType.FEM){
			
					PhiCoil coil=model.phiCoils[network.elems[j].fem_index];
										
					rowIndex=this.phiVarIndex[coil.infaceNodes[0]];
					double current=network.elems[j].I;
					//double turns=coil.getNumTurns();
		
					RHS.el[rowIndex]+=current;
					

				}
			}
		}
	
	//	RHS.show();

	}



	public  void setBoundaryCondition(Model model){




		int[] coilIndices=new int[model.numberOfRegions+1];

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			coilIndices[ir]=-1;

		for(int ic=0;ic<model.phiCoils.length;ic++){
			int nr=model.phiCoils[ic].regNo;
			coilIndices[nr]=ic;

		}

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			//if(!model.region[ir].isConductor) continue;

			if(coilIndices[ir]<0) continue;
	

			PhiCoil coil=model.phiCoils[coilIndices[ir]];

			int[] infaceNodes1=new int[model.numberOfNodes];
			boolean[] nc=new boolean[model.numberOfNodes];

			int nx=0;
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

			//	if(!model.element[i].isConductor()) continue;

			
				if(model.element[i].getRegion()!=coil.regNo) continue;

		
				int[] vertNumb=model.element[i].getVertNumb();

				for(int k=0;k<model.nElVert;k++){

					int nodeNumber=vertNumb[k];

					model.node[nodeNumber].setPhiVar(true);

					Vect coord=model.node[nodeNumber].getCoord();

					boolean onFace1=withinBox(coil.faceBox[0],coord,coil.faceCoordType[0]);

					if(onFace1){
						
						model.node[nodeNumber].setPhiKnown(true);
						model.node[nodeNumber].setPhi(1.0);
						if(nc[nodeNumber]==false){
							infaceNodes1[nx++]=nodeNumber;
							nc[nodeNumber]=true;
							
						}

					}

					else if(model.dim==3){ 
						boolean onFace2=withinBox(coil.faceBox[1],coord,coil.faceCoordType[1]);
				
						if(onFace2){

						if(nc[nodeNumber]==false){
						model.node[nodeNumber].setPhiKnown(true);
						model.node[nodeNumber].setPhi(1e-10);
					}
					}

					}

				}

			}

			coil.infaceNodes=Arrays.copyOf(infaceNodes1, nx);
				//util.hshow(coil.infaceNodes);
				//util.pr(coil.infaceNodes.length);
		}




		setPhiIndices(model);
		
		setCurrentIndices(model);
	}



	public void setPhiIndices(Model model){


		int[] coilIndices=new int[model.numberOfRegions+1];

		for(int ir=1;ir<=model.numberOfRegions;ir++)
			coilIndices[ir]=-1;

		for(int ic=0;ic<model.phiCoils.length;ic++){
			int nr=model.phiCoils[ic].regNo;
			coilIndices[nr]=ic;
		}

		int ix=0;
		phiVarIndex=new int[model.numberOfNodes+1];
		
		for(int i=0;i<=model.numberOfNodes;i++)
			phiVarIndex[i]=-1;


		for(int ir=1;ir<=model.numberOfRegions;ir++){


		//	if(!model.region[ir].isConductor) continue;

			if(coilIndices[ir]<0) continue;


			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				//if(!model.element[i].isConductor()) continue;

				int[] vertNumb=model.element[i].getVertNumb();


				for(int k=0;k<model.nElVert;k++){
					int nodeNumber=vertNumb[k];
					if(model.node[nodeNumber].isPhiVar() && !model.node[nodeNumber].isPhiKnown())
					{
				
						if(phiVarIndex[nodeNumber]==-1)
							phiVarIndex[nodeNumber]=ix++;


					}
				}
			}
		}

			for(int ic=0;ic<model.phiCoils.length;ic++){
				PhiCoil coil=model.phiCoils[ic];	

				phiVarIndex[coil.infaceNodes[0]]=ix++;
				model.node[coil.infaceNodes[0]].setPhiKnown(false);
				for(int i=1;i<coil.infaceNodes.length;i++){
					model.node[coil.infaceNodes[i]].setPhiKnown(false);
					phiVarIndex[coil.infaceNodes[i]]=phiVarIndex[coil.infaceNodes[0]];
				}

			
			}
			numberOfUnknownPhis=ix;
				

		
	//	util.pr("------->>>     "+numberOfUnknownPhis);
	}

	public void setCurrentIndices(Model model){
		
		nCurrents=model.network.no_unknown_currents;
		numberOfUnknowns=numberOfUnknownPhis;
		for(int i=0;i<nCurrents;i++)
			numberOfUnknowns++;
		
	}

	

	public void setSolution(Model model, Vect x){
		
		for(int i=1;i<=model.numberOfNodes;i++){
			int index=phiVarIndex[i];	
		
			if(index>=0){
				model.node[i].setPhi(x.el[index]);	
			}
		}
		
		setNetworkSol(model,x);


	}
	
	
	private boolean withinBox(double[] box, Vect coord1,int coordType){
		
		Vect coord=coord1.deepCopy();
		int dim=coord1.length;

		if(coordType==1){
			Vect v2=coord1.v2();
			double r=v2.norm();
			double tt=util.getAng(v2)*180/Math.PI;
			if(dim==3) coord=new Vect(r,tt,coord1.el[2]);
			else coord=new Vect(r,tt);

		}
			
		boolean result=false;

		if(coord.el[0]>=box[0] &&coord.el[0]<=box[1]&&
				coord.el[1]>=box[2] &&coord.el[1]<=box[3] &&
				(dim==2 || coord.el[2]>=box[4] &&coord.el[2]<=box[5])){
					result=true;
				}

				
				return result;
	}
	
	
	
	public  void setNetworkSol(Model model, Vect x){

		
		Network network=model.network;
		Vect indepCurrents=new Vect(network.indep_elems.length);
		int ix=x.length-nCurrents;
		for (int j = 0; j<network.indep_elems.length; ++j)
		{
			if(network.indep_elems[j].type!=ElemType.CPS) {
				network.indep_elems[j].I=x.el[ix+network.indep_elems[j].unknown_seq_no];
			indepCurrents.el[j]=network.indep_elems[j].I;
			}else{
				int time_id=network.indep_elems[j].time_id;

				double value=0;
				if(model.timeFunctions!=null && time_id>0)
				  value=model.timeFunctions[time_id].getValue(0);
				indepCurrents.el[j]=value;
				network.elems[j].I=value;
			}
		}
		
		Vect allCurrents=network.tiesetMat.transp().mul(indepCurrents);

		for (int j = 0; j<network.numElements; ++j)
		{
			if(network.elems[j].type!=ElemType.CPS) {
				network.elems[j].I=allCurrents.el[j];
				util.pr(network.elems[j].id+" "+network.elems[j].I);
			}
		}
		
	}


}





