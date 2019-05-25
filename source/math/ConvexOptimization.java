package math;

import static java.lang.Math.pow;
public class ConvexOptimization {
	
	Mat Af0,Af,Ah;
	Vect bf0,bf,bh, xopt;
	double cf0;
	double minF0;
	
	// nx: number of variables, nineq: number of inequality constraints,  neq: number of equality constraints
	int nx,nineq, neq;

	
	

	ConvexOptimization(){	}
	
public static void main(String[] args) throws Exception{
	
	ConvexOptimization cv=new ConvexOptimization();

	//  f0=1/2(x^T)(Af0)(x)+(bf0^T)(x)+cf0
	
	//  inequality constraints  (Af)x+bf<0 , equal to:
										//Af(0,0)*x0+ Af(0,1)*x1+...+ Af(0,nx)*xnx+bf(0)<0
										//Af(1,0)*x0+ Af(1,1)*x1+...+ Af(1,nx)*xnx+bf(1)<0
										//.
										//.
										//Af(nineq,0)*x0+ Af(nineq,1)*x1+...+ Af(nineq,nx)*xnx+bf(nineq)<0

	
	// quality constraints  (Ah)(x)+bh=0, equal to:
										//Ah(0,0)*x0+ Ah(0,1)*x1+...+ Ah(0,nx)*xnx+bh(0)=0
										//Ah(1,0)*x0+ Ah(1,1)*x1+...+ Ah(1,nx)*xnx+bh(1)=0
										//.
										//.
										//Ah(neq,0)*x0+ Ah(neq,1)*x1+...+ Ah(neq,nx)*xnx+bj(neq)=0
	
	
	
	double[][] aa={{1,0},{0,1}};
	cv.Af0=new Mat(aa);
	
	double[] cc={-1,-1};
	cv.bf0=new Vect(cc);
	
	cv.cf0=0;
	
	double[][] bb={{-1,0},{1,0},{0,-1},{0,1}};
	cv.Af=new Mat(bb);
	
	double[] dd={6,-10,6,-10};
	cv.bf=new Vect(dd);

	
	Vect x=cv.convexMin();
	
	x.hshow();


	
}


	public Vect convexMin(){
		
		 nx=bf0.length;
		 nineq=0;
		if(bf!=null)
			nineq=bf.length;
		
		 neq=0;
		if(bh!=null)
			neq=bh.length;
		
		
		Vect x=new Vect(nx);
		Vect lam=new Vect(nineq);
		Vect v=new Vect(neq);
		
		double stp0=.5;
		double m=1;
		Vect xp;
		double eps=1e-6;
		
		x=new Vect(0,0);
		
		for(int i=1;i<1000;i++){
			xp=x.deepCopy();
			 x=minPointLagr(lam,v);

			 double er=x.sub(xp).norm();
			if(i>1 && er<eps) break;
			 Vect[] s=gradDual(x);
			
			double stp=stp0;//(m+stp0)/(i+m);
			
			if(nineq>0)
			 lam=lam.add(s[0].times(stp));
			 //s[0].hshow();
			for(int k=0;k<nineq;k++)
		
			if(lam.el[k]<0) lam.el[k]=0;
			
			if(neq>0)
			v=v.add(s[1].times(stp));
		
			 double mf=f0( x,Af0,bf0,cf0);
		//	util.pr(i+" "+mf+"  "+er+"  "+lam.el[0]);

			 x.hshow();
		}
		
		return x;

		
	}

	public double f0(Vect x,Mat Af0,Vect bf0, double cf0){
		
		double f=0;
		
		
		f=x.dot(Af0.mul(x))/2+bf0.dot(x)+cf0;
		
		return f;
		
	}
	
	public Vect fi(Vect x){
		
		
		Vect f=Af.mul(x).add(bf);

		return f;
		
	}
	
	public Vect hi(Vect x){
		
				
		Vect f=Ah.mul(x).add(bh);

		return f;
		
	}



public Vect minPointLagr(Vect lam,Vect v){
	
	MatSolver ms=new MatSolver();

	Vect r=new Vect(bf0.length);
	if(nineq>0)
		r=r.add(Af.transp().mul(lam));
	if(neq>0)
		r=r.add(Ah.transp().mul(lam));
	Vect xmin=ms.gaussel(Af0, r.times(-1));
				
	return xmin;
}

public Vect[] gradDual(Vect x){
	
	Vect[] g=new Vect[2];
	if(nineq>0)
	 g[0]=fi(x);
	 if(neq>0)
	 g[1]=hi(x);
	return g;
}
}
