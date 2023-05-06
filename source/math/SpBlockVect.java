package math;

import java.util.Arrays;


public class SpBlockVect {
	
	public Mat[] el;
	public int [] index;
	public int nzLength;
	public int length;

	
	public SpBlockVect(){}	
	public SpBlockVect(int I){
		length=I;
	}
	
	public SpBlockVect(int I,int L){
		length=I;
		nzLength=L;
		el=new Mat[L];
		index=new int[L];
	}
	
	public SpBlockVect(int I,int L,int m, int n){
		System.out.println(I+" "+L+" "+m+" "+n);

		length=I;
		nzLength=L;
		el=new Mat[L];
		for(int i=0;i<L;i++)
			el[i]=new Mat(m,n);
		index=new int[L];
	}
	

	
	public void addToNz(Mat A,int nzr){
		int dim=A.nCol;
		for(int i=0;i<dim;i++)
			for(int j=0;j<dim;j++)
				el[nzr].el[i][j]+=A.el[i][j];
	}
	
	public int getNzLength(){
		return nzLength;
	
	}
	
	public int getLength(){
		return length;
	
	}
	
	
	public SpBlockVect times(double a){
		
		SpBlockVect b=this;
		for(int i=0;i<b.nzLength;i++){
		b.el[i]=el[i].times(a);
	
		}
		
		return b;
	}
	
	public void trim(int nzLnew){
		
		SpBlockVect b=new SpBlockVect(length,nzLnew);
		for(int i=0;i<b.nzLength;i++){
		b.el[i]=el[i];
		b.index[i]=index[i];
		}

		el=null;
		index=null;
		el=b.el;
		index=b.index;
		nzLength=nzLnew;

	}
	
	public void trim(){
		
		int nzL=0;
		for(int i=0;i<nzLength;i++)
			if(el[i].norm()!=0) nzL++;
		
		SpBlockVect b=new SpBlockVect(length,nzL);
		
		for(int i=0;i<b.nzLength;i++){
		b.el[i]=el[i];
		b.index[i]=index[i];
		}

		el=null;
		index=null;
		el=b.el;
		index=b.index;
		nzLength=nzL;

		
	}
	
	public  void sortAndTrim(int L){
		SpBlockVect sorted=new SpBlockVect(length,L);
		
		sorted.index=Arrays.copyOf(index, L);
		Arrays.sort(sorted.index);
		int m;
		for(int i=0;i<L;i++){
			m=Arrays.binarySearch(sorted.index,index[i]);
			sorted.el[m]=el[i];
		}
		el=sorted.el;
		index=sorted.index;
		nzLength=L;
			
	}
	
	
	public void addSame(SpBlockVect v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		for(int i=0;i<nzLength;i++)
			el[i]=el[i].add(v.el[i]);

	}
	
	
	public void addSmaller(SpBlockVect v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		for(int i=0;i<nzLength;i++)
			for(int j=0;j<v.nzLength;j++)
				if(index[i]==v.index[j])
					el[i]=el[i].add(v.el[i]);

	}
	
	
	public SpBlockVect deepCopy(){
		
		int I1=getLength();
		int L1=getNzLength();
		SpBlockVect w=new SpBlockVect(I1,L1);
		for(int i=0;i<nzLength;i++){
			w.el[i]=el[i];
			w.index[i]=index[i];
		}
		w.length=length;
		w.nzLength=nzLength;
		
	return w;
	}

	
	public double norm(){
		double s=0;
		for(int i=0;i<nzLength;i++)
		s+=el[i].norm();
		
		return Math.sqrt(s);
	}

	
	public void showr(){
		
		util.show(index);
	}

	public Mat  getEl(int i){
		Mat eli=new Mat();
		for(int j=0;j<nzLength;j++){
			if(index[j]==i){
				eli=el[j];
				break;
			}
		}
		return eli;
	}
	
	public void extend(int L){
		SpBlockVect vs=this.deepCopy();
		nzLength+=L;
		el=new Mat[nzLength];
		this.index=new int[nzLength];
		for(int i=0;i<vs.nzLength;i++){
		this.el[i]=vs.el[i];
		this.index[i]=vs.index[i];
		}
	}
	
}
