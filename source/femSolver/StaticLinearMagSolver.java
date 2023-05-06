package femSolver;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.log10;

import fem.Model;
import fem.TimeFunction;
import io.Writer;
import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;


public class StaticLinearMagSolver{
	int stepNumb;
	boolean usePrev=true;

	public StaticLinearMagSolver(){	}

	public Vect solve(Model model,int step,Vect x_init){

		
		this.stepNumb=step;

		if(!usePrev || x_init.length==0){
			x_init=new Vect(model.numberOfUnknowns);

		}

		if(model.numberOfUnknowns==0)
			return new Vect(model.numberOfUnknowns);
	
	SpMat L=new SpMat();

	Vect x=new Vect(model.numberOfUnknowns);

	model.solver.terminate(false);

	

if(step==0){
	model.setMagMat();

}
	
model.magMat.setRHS(model);


	
	//	util.pr("|RHS|="+model.RHS.norm());

	SpMat  Ks=model.Hs.deepCopy();
	
	/*//model.RHS.show();
	//Ks.shownz();
	//Ks.plot();
	Mat M1=Ks.matForm();
	for(int i=0;i<M1.nRow;i++)
		for(int j=i+1;j<M1.nCol;j++)
			M1.el[i][j]=M1.el[j][i];

	Mat M=new Mat(model.numberOfEdges,model.numberOfEdges);
	M.eye();
	
	for(int i=0;i<model.numberOfEdges;i++){
		for(int j=0;j<=i;j++){
	int rind=model.edgeUnknownIndex[i+1]-1;
	int cind=model.edgeUnknownIndex[j+1]-1;
	if(rind>=0 && cind>=0)
		M.el[i][j]=M1.el[rind][cind];
		}
	}
	
	SpMat Ms=new SpMat(M);
//	Ms.plot();
	
	for(int i=1;i<=model.numberOfEdges;i++){
	
		if(model.edge[i].edgeKnown) util.pr(i);
	}
	
	Vect vs=new Vect(model.numberOfEdges);
	for(int i=0;i<model.numberOfEdges;i++){

	int rind=model.edgeUnknownIndex[i+1]-1;
	if(rind>=0)
	vs.el[i]=model.RHS.el[rind];
	}*/
//Ms.shownz();
	//M.show();
	
//	vs.show();
	
	Vect Ci=Ks.scale(model.RHS);

	//util.pr("|RHS|="+model.RHS.norm());
	x_init.timesVoid(Ci.inv());
		
	L=Ks.ichol();



	if(model.RHS.abs().max()>1e-8){
			if(!usePrev || model.xp==null){
				x=model.solver.ICCG(Ks,L, model.RHS,model.errCGmax,model.iterMax,x_init);
				//x=model.solver.CG(Ks, model.RHS,model.errCGmax,model.iterMax,x_init);
			}
			else{
				x=model.solver.ICCG(Ks,L, model.RHS,model.errCGmax,model.iterMax,model.xp);
				//x=model.solver.err0ICCG(Ks,L, model.RHS,1e-2*model.errCGmax,model.iterMax,model.xp);	
			
			}
	}else{
		
		x=new Vect(model.numberOfUnknowns);
		
		model.solver.totalIter++;
		model.solver.errs.add(0.);
		model.solver.totalIter++;
		model.solver.errs.add(log10(model.errCGmax));
		model.solver.errs.add(0.);
		if(model.hasBunif) model.scaleKnownEdgeAL(0);
	}
		


		model.xp=x.deepCopy();


		x.timesVoid(Ci);
		
		//util.pr("|x|="+x.norm());

boolean unif=false;

if(unif){
x.zero();
Vect u=new Vect(0,0,1);

double Bx=1;
double By=0;
double Bz=0;
if(model.dim==3)
	Bz=model.unifB.el[2];
double Ax,Ay,Az;
double x1,y,z;
Vect A=new Vect(3);
int[] edgeDirs=new int[1+model.numberOfEdges];

for(int i=1;i<=model.numberOfElements;i++){
	boolean[] edgeDir=model.element[i].getEdgeReverse();
	int[] ne=model.element[i].getEdgeNumb();
	for(int j=0;j<model.nElEdge;j++)
		if(edgeDir[j])
		edgeDirs[ne[j]]=-1;
		else	edgeDirs[ne[j]]=1;
}

for(int i=1;i<=model.numberOfEdges;i++){


	Vect edgeVect=model.edge[i].node[1].getCoord().sub(model.edge[i].node[0].getCoord());

	Vect center=model.edge[i].node[1].getCoord().add(model.edge[i].node[0].getCoord()).times(.5);
	
	x1=center.el[0];
	y=center.el[1];

	Az=y*Bx-x1*By;
	if(model.dim==3){
		z=center.el[2];
		Ax=x1*By;
		Ay=y*Bz;
		A=new Vect(Ax,Ay,Az);
	}else{
		A=new Vect(0,0,Az);
	}


	double a=edgeVect.dot(A);
	

	if(model.edgeUnknownIndex[i]>0)
		x.el[model.edgeUnknownIndex[i]-1]=a;
		else{
		model.edge[i].setA(a);
		}




}
	

	}
	

/*for(int i=1;i<=model.numberOfEdges;i++){

//util.pr("Edge "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")= "+model.edge[i].getA());
util.pr("Edge "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")");

}*/

boolean edgeNodes=false;
if(edgeNodes)
for(int i=1;i<=model.numberOfElements;i++){

	int[] ne=model.element[i].getEdgeNumb();
	//util.pr("Edge "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")= "+model.edge[i].getA());
	util.pr("element "+i);
	util.pr(" edges:");
	for(int j=0;j<model.nElEdge;j++)
		util.ph(ne[j]+", ");
	util.pr("");
	
	int[] nn=model.element[i].getVertNumb();
	//util.pr("Edge "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")= "+model.edge[i].getA());
	util.pr(" nodes:");
	for(int j=0;j<model.nElVert;j++)
		util.ph(nn[j]+",");
	util.pr("");
	
	
	//util.pr("element "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")");


	}
	

	model.setSolution(x);	

	
		System.out.println("Bmax ( linear analysis): "+model.Bmax);
		


		return x;



}




}
