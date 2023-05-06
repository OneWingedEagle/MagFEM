package fem;

import static java.lang.Math.abs;

import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;

/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 20, 2012.
 */
public class TMatrix {

	private int nEdEd,nNodEd,nEdNod,
	nNodNod,numberOfRegions;
	private int nElVert,nElEdge,dim;
	private Calculator calc;

	public  TMatrix(){}

	public   TMatrix(Model model){


		this.nElVert=model.nElVert;
		this.nElEdge=model.nElEdge;
		this.dim=model.dim;

		this.nEdEd=model.nEdEd;
		this.nNodEd=model.nNodEd;
		this.nNodNod=model.nNodNod;
		this.nEdNod=model.nEdNod;
		this.numberOfRegions=model.numberOfRegions;

		this.calc=new Calculator(model);
	}


	public void setTMat(Model model){
		setTMatEdge(model);

	}

	public void setTMatEdge(Model model){	



		
		
		double eps=1e-6,cPB=model.cpb; 
		double[][] H=new double[this.nElEdge][this.nElEdge];
		int m,columnEdgeNumb,columnIndex,matrixRow=0, rowEdgeNumb,ext=6;


		int[] nz=new int[model.numberOfUnknownT];
		model.Cs=new SpMat(model.numberOfUnknownT);
		
		
		for(int i=0;i<model.numberOfUnknownT;i++){
			
			model.Cs.row[i]=new SpVect(model.numberOfUnknownT,this.nEdEd);
		}

		model.bT=new Vect(model.numberOfUnknownT);

		for(int ir=1;ir<=this.numberOfRegions;ir++){
			if(!model.region[ir].hasJ) continue;
		
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){


				H=this.calc.Te(model,i);

				int[] edgeNumb=model.element[i].getEdgeNumb();

				for(int j=0;j<this.nElEdge;j++){
					rowEdgeNumb=edgeNumb[j];

					if(model.edge[rowEdgeNumb].edgeKnownT ) continue;


					matrixRow=model.T_unknownIndex[rowEdgeNumb]-1;

					
					for(int k=0;k<this.nElEdge;k++){

						columnEdgeNumb=edgeNumb[k];

						//===== Periodic BC ================

						if(model.edge[columnEdgeNumb].aPBC ||model.edge[rowEdgeNumb].aPBC){

							cPB=model.edge[columnEdgeNumb].getnPBC()*model.edge[rowEdgeNumb].getnPBC();
							
							H[j][k]*=cPB;
						}	



						//===========  right-hand side 1

						if(model.edge[columnEdgeNumb].edgeKnownT  ){
								model.bT.el[matrixRow]-=H[j][k]*model.edge[columnEdgeNumb].T;
					
								continue;
						}


		
						//=======================

			
						columnIndex=model.T_unknownIndex[columnEdgeNumb]-1;
						if(columnIndex>matrixRow) continue;
						m=util.search(model.Cs.row[matrixRow].index,nz[matrixRow]-1,columnIndex);
						if(m<0)
						{	

							if(abs(H[j][k])>eps  ){	

								model.Cs.row[matrixRow].index[nz[matrixRow]]=columnIndex;

								model.Cs.row[matrixRow].el[nz[matrixRow]]=H[j][k];
								
								

								nz[matrixRow]++;
								
								//===========================
								if(nz[matrixRow]==model.Cs.row[matrixRow].nzLength-1){
									model.Cs.row[matrixRow].extend(ext);
								}
								//===========================
							}
						}
						else{

							model.Cs.row[matrixRow].addToNz(H[j][k],m);
							
						}


					}			
				}
			}
		}



		for(int i=0;i<model.numberOfUnknownT;i++){
			model.Cs.row[i].sortAndTrim(nz[i]);
		}

	


	}





}
