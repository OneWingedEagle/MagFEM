package femSolver;


import static java.lang.Math.PI;

import fem.Model;
import math.Complex;
import math.SpMat;
import math.SpMatComp;
import math.SpVectComp;
import math.Vect;
import math.VectComp;
import math.util;

public class ACMagSolver {


	public ACMagSolver(){	}

	public VectComp solve(Model model, int step ){


		//SpMat L=new SpMat();

	//	Vect x=new Vect(model.numberOfUnknowns);

		model.solver.terminate(false);


		model.magMat.setRHS(model);



		model.setMagMat();

		
		util.pr("frequency ="+model.freq);
		double  w=2.*Math.PI*model.freq;
		double  wr=-1./w;
		SpMatComp Ks=new SpMatComp(model.Hs,model.Ss.timesNew(w));

		if(model.analysisMode==2){
	
			for(int i=0;i<model.numberOfVarNodes;i++){
			
				SpVectComp spv=new SpVectComp(model.Ps.row[i]);
				spv=spv.augh(new SpVectComp(model.Qs.row[i]).times(new Complex(0,wr)));
				Ks.row[i+model.numberOfUnknownEdges]=spv.deepCopy();
			}



		}
		//model.Qs.shownz();
	//	model.RHS.show();

		//Ks.diagSym().show();
		Ks.setSymHerm(1); // symmetric but not Hermitian

	//	Ks.shownz();
		VectComp  b=new VectComp(model.RHS);
	

		int m=b.length;
		model.Ci=Ks.scale(b);



		SpMatComp Ls=Ks.ichol(1.2);
		
		Ls.setSymHerm(1); // symmetric but not Hermitian


		VectComp xc;

		if(b.norm()>1e-8){
			xc=model.solver.COICCG(Ks,Ls,b,model.errCGmax,model.iterMax,new VectComp(m),1,false);
			//xc=model.solver.COCG(Ks,b,model.errCGmax,model.iterMax,new VectComp(m),1,false);
		}
		else
			xc=new VectComp(m);

		xc.timesVoid(model.Ci);	
	

		return xc;



	}


}
