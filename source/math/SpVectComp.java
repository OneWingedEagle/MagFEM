package math;

import java.util.Arrays;


public class SpVectComp {
	
	public Complex[] el;
	public int [] index;
	public int nzLength;
	public int length;

	
	public SpVectComp(){}	
	
	public SpVectComp(int I){
		length=I;
	}
	
	
	public SpVectComp(Complex[] el){
		
		int L=el.length;
		length=L;
		int nz=0;
		for(int i=0;i<L;i++){
			if(!el[i].equals(new Complex(0,0))) nz++;
		}
		nzLength=nz;
		this.el=new Complex[nz];
		index=new int[nz];
		nz=0;
		for(int i=0;i<L;i++){
			if(!el[i].equals(new Complex(0,0))) {
				this.el[nz]=el[i].deepCopy();
				index[nz++]=i;
			}
		}
	}
	
	public SpVectComp(double[] el1){
		
		int L=el1.length;
		length=L;
		int nz=0;
		for(int i=0;i<L;i++){
			if(el1[i]!=0) nz++;
		}
		nzLength=nz;
		this.el=new Complex[nz];
		index=new int[nz];
		nz=0;
		for(int i=0;i<L;i++){
			if(el1[i]!=0) {
				this.el[nz]=new Complex(el1[i],0);
				index[nz++]=i;
			}
		}
	}
	
	public SpVectComp(SpVect sv,Complex a){
		
		int L=sv.length;
		length=L;
	
		nzLength=sv.nzLength;
		this.el=new Complex[nzLength];
		index=new int[nzLength];
		for(int i=0;i<nzLength;i++){
		{
				this.el[i]=new Complex(sv.el[i],0).times(a);
				index[i]=sv.index[i];
			}
		}
	}
	
	
	public SpVectComp(SpVect sv){
		
		int L=sv.length;
		length=L;
	
		nzLength=sv.nzLength;
		this.el=new Complex[nzLength];
		index=new int[nzLength];
		for(int i=0;i<nzLength;i++){
		{
				this.el[i]=new Complex(sv.el[i],0);
				index[i]=sv.index[i];
			}
		}
	}
	
	
	public SpVectComp(double[] el1,Complex a){
		
		int L=el.length;
		length=L;
		int nz=0;
		for(int i=0;i<L;i++){
			if(el1[i]!=0) nz++;
		}
		nzLength=nz;
		this.el=new Complex[nz];
		index=new int[nz];
		nz=0;
		for(int i=0;i<L;i++){
			if(el1[i]!=0) {
				this.el[nz]=new Complex(el1[i],0).times(a);
				index[nz++]=i;
			}
		}
	}
	
	public SpVectComp(int I,int L){
		length=I;
		nzLength=L;
		el=new Complex[L];
		index=new int[L];
	}
	
	
	public Complex dot(SpVectComp v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		
		Complex dot=new Complex(0,0);
		Complex vel;
		for(int i=0;i<nzLength;i++){
			vel=v.getEl(index[i]);
			if(!vel.equals(new Complex(0,0)))
				dot=dot.add(el[i].times(vel));
	
		}
		
		return dot;
	}
	
	

	public Complex dot(SpVect v,int n){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		int pass=0;
		Complex dot=new Complex(0,0);
		for(int i=0;i<nzLength;i++){
			for(int j=pass+1;j<v.nzLength;j++)
				if(v.index[j]==index[i] && index[i]<=n)
				{
					pass=j;
					dot=dot.add(el[i].times(v.el[j]));
				break;
				}
				
			}
		
		return dot;
	}
	
	
	
	public Complex dot(Vect v,int n){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		int j;
		Complex dot=new Complex(0,0);		
		
		for(int i=0;i<nzLength;i++){
			j=index[i];
			if(j>n) break;
			dot=dot.add(el[i].times(v.el[j]));
				}
		
		return dot;
	}
	
	public Complex dot(VectComp v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		int j;
		Complex dot=new Complex(0,0);		
		
		for(int i=0;i<nzLength;i++){
			j=index[i];

			dot=dot.add(el[i].times(v.el[j]));
				}
		
		return dot;
	}
	
	
	public Complex dotUp(Vect v,int n){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		int j;
		Complex dot=new Complex(0,0);		
		for(int i=0;i<nzLength;i++){
			j=index[i];
			if(j<n) continue;
			dot=dot.add(el[i].times(v.el[j]));

				}
		
		return dot;
	}
	

	
	
	public Complex dot(int n,Complex[] D,int[] index){
	
		Complex dot=new Complex();		
		int i=0;
			while(i<nzLength && index[i]<=n){
		
			dot=dot.add(D[index[i]].times(el[i].times(el[i].conj())));
			i++;
			}

		
		return dot;
		
		
	}
	
	public Complex dotSym(int n,Complex[] D,int[] index){
		
		Complex dot=new Complex();		
		int i=0;
			while(i<nzLength && index[i]<=n){
		
			dot=dot.add(D[index[i]].times(el[i].times(el[i])));
			i++;
			}

		
		return dot;
		
		
	}
	
	
	
	public Complex dot(SpVectComp v,int n,Complex[] D,int[] index1, int[] index2){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		
		int indi;
		Complex dot=new Complex(0,0);		
		int i=0,j=0;
		while(i<nzLength && index1[i]<=n){
			indi=index1[i];
			j=0;
			while(j<v.nzLength && index2[j]<=n){
				if(index2[j]==indi){
					dot=dot.add(el[i].times(v.el[j].conj()).times(D[index[i]]));
				}
				j++;
			}
		i++;
		}
		
		return dot;
	}

	public Complex dotSym(SpVectComp v,int n,Complex[] D,int[] index1, int[] index2){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		
		int indi;
		Complex dot=new Complex(0,0);		
		int i=0,j=0;
		while(i<nzLength && index1[i]<=n){
			indi=index1[i];
			j=0;
			while(j<v.nzLength && index2[j]<=n){
				if(index2[j]==indi){
					dot=dot.add(el[i].times(v.el[j]).times(D[index[i]]));
				}
				j++;
			}
		i++;
		}
		
		return dot;
	}

	
	public void setEl(Complex a, int cl,int nzr){
		el[nzr]=a.deepCopy();
		index[nzr]=cl;
	
	}

	
	public void setNzEl(Complex a,int nzr){
		el[nzr]=a.deepCopy();
	}
	
	public void addToNz(Complex a,int nzr){
		el[nzr]=el[nzr].add(a);
	}
	
	public int getNzLength(){
		return nzLength;
	
	}
	
	public int getLength(){
		return length;
	
	}
	
	
	public SpVectComp times(double a){
		
		SpVectComp b=this.deepCopy();
		for(int i=0;i<b.nzLength;i++){
		b.el[i]=el[i].times(a);
	
		}
		
		return b;
	}
	
	public SpVectComp times(Complex a){
		
		SpVectComp b=this.deepCopy();
		for(int i=0;i<b.nzLength;i++){
		b.el[i]=el[i].times(a);
	
		}
		
		return b;
	}
	
	public void trim(int nzLnew){
		
		SpVectComp b=new SpVectComp(length,nzLnew);
		for(int i=0;i<b.nzLength;i++){
		b.el[i]=el[i].deepCopy();
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
			if(!el[i].equals(new Complex(0,0))) nzL++;
		
		SpVectComp b=new SpVectComp(length,nzL);
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
		SpVectComp sorted=new SpVectComp(length,L);
		
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

	
	public void addSame(SpVectComp v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		for(int i=0;i<nzLength;i++)
			el[i]=el[i].add(v.el[i]);

	}
	
	
	public void addSmaller(SpVectComp v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		for(int i=0;i<nzLength;i++)
			for(int j=0;j<v.nzLength;j++)
				if(index[i]==v.index[j])
					el[i]=el[i].add(v.el[i]);
	}
	

	public SpVectComp augh(SpVectComp v){
		int I1=getLength();
		int L1=getNzLength();
		

		SpVectComp w=new SpVectComp(I1+v.getLength(),L1+v.getNzLength());
		for(int i=0;i<nzLength;i++){
			w.el[i]=el[i].deepCopy();
			w.index[i]=index[i];
		}
		
		for(int i=0;i<v.nzLength;i++){
			w.el[i+L1]=v.el[i];
			w.index[i+L1]=v.index[i]+I1;
		}
		return w;
	}
	
public SpVectComp deepCopy(){
		
		int I1=getLength();
		int L1=getNzLength();
		SpVectComp w=new SpVectComp(I1,L1);
		for(int i=0;i<nzLength;i++){
			w.el[i]=el[i].deepCopy();
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
	
	public void Ddotv(Vect D){
		if(length!=D.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		for(int i=0;i<nzLength;i++)
			el[i]=el[i].times(D.el[index[i]]);
	}
	
	public void show(){
		
		for(int i=0;i<nzLength;i++)
			System.out.format("(%4d)   %12.4f  %12.4f\n",index[i], el[i].re,  el[i].im);
	}
	
	public void show(Complex[] el){
		
		for(int i=0;i<el.length;i++)
		//	System.out.format("%12.4f +j %12.4f\n", getEl(i).re,getEl(i).im);
			getEl(i).show();
	//	System.out.format("%12.4f\n", getEl(i).re);
	}
	
	public void hshow(){
		
		for(int i=0;i<length;i++)
			System.out.format("%12.4f  %12.4f", getEl(i).re, getEl(i).im);
		System.out.println();
	}
	
	public void shownz(){
		
		for(int i=0;i<nzLength;i++)
			System.out.format("(%4d)   %12.4f  %12.4f\n",index[i], el[i].re,  el[i].im);
	}
	
public void showr(){
		
		util.show(index);
	}

	public Complex getEl(int i){
		Complex eli=new Complex();
		for(int j=0;j<nzLength;j++){
			if(index[j]==i){
				eli=el[j].deepCopy();
				eli.show();
				break;
			}
		}
		return eli;
	}
	
	public void extend(int L){
		SpVectComp vs=this.deepCopy();
		nzLength+=L;
		el=new Complex[nzLength];
		index=new int[nzLength];
		for(int i=0;i<vs.nzLength;i++){
		this.el[i]=vs.el[i].deepCopy();
		this.index[i]=vs.index[i];
		}
	}
	
}
