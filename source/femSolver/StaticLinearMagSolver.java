package femSolver;

import static java.lang.Math.PI;
import static java.lang.Math.abs;

import fem.Model;
import io.Writer;
import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;


public class StaticLinearMagSolver{
	int stepNumb;
	boolean usePrev=false;

	public StaticLinearMagSolver(){	}

	public Vect solve(Model model,int step,Vect x_init){
		
		this.stepNumb=step;
		
		if(x_init.length==0)
			x_init=new Vect(model.numberOfUnknowns);


	
	SpMat L=new SpMat();

	Vect x=new Vect(model.numberOfUnknowns);

	model.solver.terminate(false);


	model.magMat.setRHS(model);


if(step==0){
	model.setMagMat();

}
	
	//=== known values go to right hand side 

		model.RHS=model.RHS.sub(model.HkAk);
	//	model.RHS.show();

		//util.pr("|RHS|="+model.RHS.norm());

	SpMat  Ks=model.Hs.deepCopy();
	
	//model.RHS.show();
	//Ks.shownz();
	//Ks.plot();

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
		}

		else
			x=new Vect(x_init.length);

		model.xp=x.deepCopy();


		x.timesVoid(Ci);

boolean unif=false;
		
if(unif){
		Vect u=new Vect(0,0,1);
	for(int i=1;i<=model.numberOfEdges;i++){

Vect v=model.edge[i].node[1].getCoord().sub(model.edge[i].node[0].getCoord());
//v.hshow();
Vect xx=model.edge[i].node[1].getCoord().add(model.edge[i].node[0].getCoord());
double a=v.dot(u.times(xx.el[0]/2));
if(model.edgeUnknownIndex[i]>0)
x.el[model.edgeUnknownIndex[i]-1]=a;
else
model.edge[i].setA(a);
	}

//util.pr("Edge "+i+" ("+model.edge[i].node[0].id+" --->"+model.edge[i].node[1].id+")= "+model.edge[i].getA());

	}
	
	

	model.setSolution(x);	
	
	//model.setB();
/*	for(int i=1;i<=model.numberOfElements;i++)
		if(i>1)
		for(int j=0;j<model.nElEdge;j++)
			model.edge[model.element[i].getEdgeNumb(j)].setA(0);*/
	
		System.out.println("Bmax ( linear analysis): "+model.Bmax);
		


		return x;



}




}
