package fem;
import static java.lang.Math.*;

import java.util.Arrays;

import math.Mat;
import math.SpMat;
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
	
	private SpMat phiMat;
	private Vect RHS;
	boolean byCircuit;
	
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
		byCircuit=true;
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

		double eps=1e-8;
		
		int ext=10;
		Mat Ke;
		int m,nodeNumber,columnIndex=0,matrixRow;
		int[] nz=new int[numberOfUnknownPhis];

		phiMat=new SpMat(numberOfUnknownPhis, numberOfUnknownPhis,model.nNodNod);

		for(int ir=1;ir<=model.numberOfRegions;ir++){

			if(ir!=regNo) continue;
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				if(!model.element[i].isConductor()) continue;

				Ke=this.calc.elemPhiMat(model,i);
			
				int[] vertNumb=model.element[i].getVertNumb();

				for(int j=0;j<model.nElVert;j++){

					if(!model.node[vertNumb[j]].isPhiVar() || model.node[vertNumb[j]].isPhiKnown() ) continue;

					matrixRow=phiVarIndex[vertNumb[j]]-1;
					for(int k=0;k<model.nElVert;k++){
						nodeNumber=vertNumb[k];
						columnIndex=phiVarIndex[nodeNumber]-1;								

						if(columnIndex==-1 || columnIndex>matrixRow) continue;	

						
						m=util.search(phiMat.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	
							
							if(abs(Ke.el[j][k])>eps ){
								
								//===========================
								if(nz[matrixRow]==phiMat.row[matrixRow].nzLength-1){
									phiMat.row[matrixRow].extend(ext);
								}
								//===========================
								
								phiMat.row[matrixRow].index[nz[matrixRow]]=columnIndex;
								phiMat.row[matrixRow].el[nz[matrixRow]++]=Ke.el[j][k];
							}

						}

						else{

							phiMat.row[matrixRow].addToNz(Ke.el[j][k],m);
						}

					}			
				}
			}
		}
		if(byCircuit){
		matrixRow=numberOfUnknownPhis-1;
		nz[matrixRow]=2;
		//phiMat.row[matrixRow].hshow();
		double Rext=1e-10;
		phiMat.row[matrixRow].index[0]=matrixRow-1;
		phiMat.row[matrixRow].el[0]=-1;
		phiMat.row[matrixRow].index[1]=matrixRow;
		phiMat.row[matrixRow].el[1]=-Rext;
		}
		
		phiMat.sortAndTrim(nz);

		util.pr("PHIMAT");
		phiMat.shownz();
		//phiMat.show();
	}


	public void setRHS(Model model,double factor){		

		if(byCircuit){
		int indx=numberOfUnknownPhis-1;
		phiMat.row[indx].el[1]=-1e12; // openning port
		}
		
		RHS=new Vect(numberOfUnknownPhis);

		Vect elemRHS=new Vect(model.nElVert);
		
		int matrixRow=0;
				
		boolean isConductive;
		
	//	if(!byCircuit)
		for(int ir=1;ir<=model.numberOfRegions;ir++){

			if(ir!=regNo) continue;
			
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				isConductive=model.element[i].isConductor();
			

				if(!isConductive) continue;
	
				int[] vertNumb=model.element[i].getVertNumb();

				
				boolean elemHasPhiVar=false;
				for(int k=0;k<model.nElVert;k++){
					int nodeNumber=vertNumb[k];
					
					if(model.node[nodeNumber].isPhiVar() && !model.node[nodeNumber].isPhiKnown() ){
						elemHasPhiVar=true;
						break;
					}
				}
					
				if(!elemHasPhiVar) continue;
					
				elemRHS=this.calc.elemPhiVect(model,i);
			
				elemRHS=elemRHS.times(factor);
			
				
					for(int k=0;k<model.nElVert;k++){
						int nodeNumber=vertNumb[k];
						
						if(model.node[nodeNumber].isPhiVar() && !model.node[nodeNumber].isPhiKnown() ){
							matrixRow=phiVarIndex[nodeNumber]-1;
							RHS.el[matrixRow]+=-elemRHS.el[k];
						}
					}

				

			}


	}
		

	}
	
	public void setRHS(Model model){		

		Mat Ke;
		
		RHS=new Vect(numberOfUnknownPhis);

		Vect elemPhiVec=new Vect(model.nElVert);
		
		Vect elemRHS=new Vect(model.nElVert);
		
		int matrixRow=0, columnIndex;
				
		boolean isConductive;
		
		if(!byCircuit)
		for(int ir=1;ir<=model.numberOfRegions;ir++){

			if(ir!=regNo) continue;

			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

				
				int[] vertNumb=model.element[i].getVertNumb();

				
				boolean elemHasKnownPhi=false;
				for(int k=0;k<model.nElVert;k++){
					int nodeNumber=vertNumb[k];
					if(model.node[nodeNumber].isPhiVar() && model.node[nodeNumber].isPhiKnown() ){
						elemHasKnownPhi=true;
						break;
					}
				}
					
				if(!elemHasKnownPhi) continue;
					
				Ke=this.calc.elemPhiMat(model,i);

				elemPhiVec.timesVoid(0);

	
					for(int k=0;k<model.nElVert;k++){
						int nodeNumber=vertNumb[k];
						if(model.node[nodeNumber].isPhiVar() && model.node[nodeNumber].isPhiKnown() ){
							elemPhiVec.el[k]=model.node[nodeNumber].getPhi();
						}
					}

					elemRHS=Ke.mul(elemPhiVec);
				
					for(int k=0;k<model.nElVert;k++){
						int nodeNumber=vertNumb[k];
						
						if(model.node[nodeNumber].isPhiVar() && !model.node[nodeNumber].isPhiKnown() ){
							
							matrixRow=phiVarIndex[nodeNumber]-1;
							RHS.el[matrixRow]+=-elemRHS.el[k];
						}
					}

				

			}


	}


		if(byCircuit){
			RHS.el[RHS.length-1]=-current;
		}
	}


public  void setBoundaryCondition(Model model){

	int[] infaceNodes1=new int[model.numberOfNodes];
	boolean[] nc=new boolean[model.numberOfNodes];
	
	
	int nx=0;

	for(int ir=1;ir<=model.numberOfRegions;ir++){


		if(ir!=regNo) continue;
		
		for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

			if(!model.element[i].isConductor()) continue;
			
			int[] vertNumb=model.element[i].getVertNumb();

		for(int k=0;k<model.nElVert;k++){
			
				int nodeNumber=vertNumb[k];
				
				model.node[nodeNumber].setPhiVar(true);
	
				Vect coord=model.node[nodeNumber].getCoord();
			

				if(coord.el[0]>=faceBox[0][0] &&coord.el[0]<=faceBox[0][1]&&
					coord.el[1]>=faceBox[0][2] &&coord.el[1]<=faceBox[0][3] &&
						(model.dim==2 || coord.el[2]>=faceBox[0][4] &&coord.el[2]<=faceBox[0][5])){
					model.node[nodeNumber].setPhiKnown(true);
					model.node[nodeNumber].setPhi(1.);
				
					if(nc[nodeNumber]==false){
					infaceNodes1[nx++]=nodeNumber;
					nc[nodeNumber]=true;
					}
					
				}

				else if(coord.el[0]>=faceBox[1][0] &&coord.el[0]<=faceBox[1][1] &&
							coord.el[1]>=faceBox[1][2] &&coord.el[1]<=faceBox[1][3] &&
								(model.dim==2 || coord.el[2]>=faceBox[1][4] &&coord.el[2]<=faceBox[1][5])){
			//		else if(model.node[nodeNumber].getCoord(0)<.001){
					model.node[nodeNumber].setPhiKnown(true);
					model.node[nodeNumber].setPhi(1e-10);
			
								
		}

		}
		
	}
		
	}
	
	infaceNodes=Arrays.copyOf(infaceNodes1, nx);
	
	//util.show(infaceNodes);
	
	setPhiIndices(model);
}


public void setPhiIndices(Model model){

	int ix=1;
	phiVarIndex=new int[model.numberOfNodes+1];
	
	for(int ir=1;ir<=model.numberOfRegions;ir++){


		if(ir!=regNo) continue;
		
		for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){

			if(!model.element[i].isConductor()) continue;

			
			int[] vertNumb=model.element[i].getVertNumb();

			
		for(int k=0;k<model.nElVert;k++){
				int nodeNumber=vertNumb[k];
		if(model.node[nodeNumber].isPhiVar() && !model.node[nodeNumber].isPhiKnown())
		{

			if(phiVarIndex[nodeNumber]==0)
			phiVarIndex[nodeNumber]=ix++;


		}
		}
		}
	}

	if(byCircuit){
	phiVarIndex[infaceNodes[0]]=ix++;
	model.node[infaceNodes[0]].setPhiKnown(false);
	for(int i=1;i<infaceNodes.length;i++){
		model.node[infaceNodes[i]].setPhiKnown(false);
		phiVarIndex[infaceNodes[i]]=phiVarIndex[infaceNodes[0]];
	}
	
	ix++; // for current
	}


	numberOfUnknownPhis=ix-1;

	//util.pr("------->>>     "+numberOfUnknownPhis);
}


public Vect solve(Model model ){

Vect x=new Vect(numberOfUnknownPhis);

if(numberOfUnknownPhis==0) return x;

SpMat L=new SpMat();


model.solver.terminate(false);

//RHS.show();
SpMat Ks=phiMat.deepCopy();



	Vect Ci=Ks.scale(RHS);


	L=Ks.ichol(1.);


	x=model.solver.ICCG(Ks,L, RHS,model.errCGmax*1e-4,model.iterMax);


	x.timesVoid(Ci);
	
	for(int i=1;i<=model.numberOfNodes;i++){
		int index=phiVarIndex[i]-1;	

		if(index>=0){
			model.node[i].setPhi(x.el[index]);	
		}
	}

	return x;

}
}