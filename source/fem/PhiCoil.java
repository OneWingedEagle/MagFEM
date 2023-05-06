package fem;
import static java.lang.Math.*;

import java.util.Arrays;

import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;



public class PhiCoil {

	public int regNo,index;
	public double conductivity;
	private double numTurns;
	public int constPotIndex;
	public double volatge,current,resistance;
	
	public double[][] faceBox;
	public int[] faceCoordType;
	
	private Calculator calc;
	private int[] phiVarIndex;
	public int[] infaceNodes;
	private int numberOfUnknownPhis;
	
	private SpMat matrix;
	public Vect load;

	
	public PhiCoil(int nr)
	{
		regNo=nr;
		
		faceBox=new double[2][6];
		faceCoordType=new int[2];
	}
	
	public void setMatrix(Model model){
		
		this.calc=new Calculator(model);
		
		setPhiMat(model);
		
	}

	
	public void setSigma(double sig){
		conductivity=sig;
		current=1.;
		
	}
	public void setNumTurns(double nt){
		numTurns=nt;


	}
	
	public void setResistance(double res){
		resistance=res;


	}
	public double getResistance(){
		return resistance;


	}

	public int getRegNo(){return regNo;}
	
	public double getNumTurns(){return numTurns;}
	
	public double  getConductivity(){return conductivity;}
	
	double  getCurrent(){return current;}
	

	public  void setPhiMat(Model model){


		double eps=1e-12;

		int ext=10;
		Mat Ke;
		int m,nodeNumber,columnIndex=0,matrixRow;
		int[] nz=new int[numberOfUnknownPhis];

		matrix=new SpMat(numberOfUnknownPhis, numberOfUnknownPhis,model.nNodNod);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			if(ir!=this.regNo) continue;

						
			double conductivity=1.;


			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

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
				
						m=util.search(matrix.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	

							if(abs(Ke.el[j][k])>eps ){
		
								//===========================
								if(nz[matrixRow]==matrix.row[matrixRow].nzLength-1){
									matrix.row[matrixRow].extend(ext);
								}
								//===========================

								matrix.row[matrixRow].index[nz[matrixRow]]=columnIndex;
								matrix.row[matrixRow].el[nz[matrixRow]++]=Ke.el[j][k];
							}

						}

						else{

							matrix.row[matrixRow].addToNz(Ke.el[j][k],m);

						}

					}			
				}
			}
		}

		matrix.sortAndTrim(nz);
		//matrix.show();
	}


	public  void setBoundaryCondition(Model model){



		for(int ir=1;ir<=model.numberOfRegions;ir++){


			if(ir!=this.regNo) continue;
	


			int[] infaceNodes1=new int[model.numberOfNodes];
			boolean[] nc=new boolean[model.numberOfNodes];

			int nx=0;
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
	
				int[] vertNumb=model.element[i].getVertNumb();
		
				for(int k=0;k<model.nElVert;k++){

					int nodeNumber=vertNumb[k];

					model.node[nodeNumber].setPhiVar(true);

					Vect coord=model.node[nodeNumber].getCoord();

					boolean onFace1=withinBox(this.faceBox[0],coord,this.faceCoordType[0]);

					if(onFace1){
						
						model.node[nodeNumber].setPhiKnown(true);
						model.node[nodeNumber].setPhi(1.0);
						if(nc[nodeNumber]==false){
							infaceNodes1[nx++]=nodeNumber;
							nc[nodeNumber]=true;
							
						}

					}

					else if(model.dim==3){ 
						boolean onFace2=withinBox(this.faceBox[1],coord,this.faceCoordType[1]);
				
						if(onFace2){

						if(nc[nodeNumber]==false){
						model.node[nodeNumber].setPhiKnown(true);
						model.node[nodeNumber].setPhi(0);
					}
					}

					}

				}

			}

			infaceNodes=Arrays.copyOf(infaceNodes1, nx);

		}




		setPhiIndices(model);
		
	}



public void setPhiIndices(Model model){


	int ix=0;
	phiVarIndex=new int[model.numberOfNodes+1];
	
	for(int i=0;i<=model.numberOfNodes;i++)
		phiVarIndex[i]=-1;


	for(int ir=1;ir<=model.numberOfRegions;ir++){


		if(ir!=this.regNo) continue;


		for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


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



			phiVarIndex[this.infaceNodes[0]]=ix++;
			model.node[this.infaceNodes[0]].setPhiKnown(false);
			for(int i=1;i<this.infaceNodes.length;i++){
				model.node[this.infaceNodes[i]].setPhiKnown(false);
				phiVarIndex[this.infaceNodes[i]]=phiVarIndex[this.infaceNodes[0]];
			}

		
		
			numberOfUnknownPhis=ix;
			

	
//	util.pr("------->>>     "+numberOfUnknownPhis);
}


public Vect solve(Model model ){

Vect x=new Vect(numberOfUnknownPhis);

if(numberOfUnknownPhis==0) return x;

SpMat L=new SpMat();


model.solver.terminate(false);

	SpMat Ks=matrix.deepCopy();


	Vect b=new Vect(numberOfUnknownPhis);

	b.el[b.length-1]=current;

	Vect Ci=Ks.scale(b);

	L=Ks.ichol();

	x=model.solver.ICCG(Ks,L, b,model.errCGmax*1e-3,model.iterMax);

	x.timesVoid(Ci);
	
	for(int i=1;i<=model.numberOfNodes;i++){
		int index=phiVarIndex[i];	
	
		if(index>=0){
			model.node[i].setPhi(x.el[index]);	
		}
	}
	
	return x;

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


public void makeLoadVector(Model model){	
	
	this.setBoundaryCondition(model);
	
	this.setMatrix(model);
	
	this.solve(model);

	load=new Vect(model.numberOfUnknowns);

	double[] loss=new double[1];

	int matrixRow=0, rowEdgeNumb;

	Vect elemV=new Vect(model.nElEdge);

	
	double coilLoss=0;

	for(int ir=1;ir<=model.numberOfRegions;ir++){


		if(ir!=this.regNo) continue;


		for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


			elemV=calc.elemVectPhiCoil(model,i,loss, 1.);

			elemV.timesVoid(1);

			if(conductivity>0)
				coilLoss+=loss[0];

			int[] edgeNumb=model.element[i].getEdgeNumb();

			for(int j=0;j<model.nElEdge;j++){
				rowEdgeNumb=edgeNumb[j];

				if(model.edge[rowEdgeNumb].edgeKnown ) continue;


				matrixRow=model.edgeUnknownIndex[rowEdgeNumb]-1;


				//===========  right-hand side

				load.el[matrixRow]+=elemV.el[j];	


			}
		}
	}
	
	if(conductivity>0){
		coilLoss*=numTurns*numTurns/conductivity;
	}

	load.timesVoid(numTurns);

		double res=coilLoss;
		this.setResistance(res);
		
		util.pr("Coil Resistance= "+res);

	
}

}