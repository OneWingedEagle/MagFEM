package math;

import java.util.Arrays;


import static java.lang.Math.*;

public class SpVect {
	
	public double[] el;
	public int [] index;
	public int nzLength;
	public int length;

	
	public SpVect(){}	
	public SpVect(int I){
		length=I;
	}
	
	public SpVect(Vect  v){
		int L=v.el.length;
		length=L;
		int nz=0;
		for(int i=0;i<L;i++){
			if(v.el[i]!=0) nz++;
		}
		nzLength=nz;
		this.el=new double[nz];
		index=new int[nz];
		nz=0;
		for(int i=0;i<L;i++){
			if(v.el[i]!=0) {
				this.el[nz]=v.el[i];
				index[nz++]=i;
			}
		}
		
	}
	

	
	public SpVect(double el[]){
		
		int L=el.length;
		length=L;
		int nz=0;
		for(int i=0;i<L;i++){
			if(el[i]!=0) nz++;
		}
		nzLength=nz;
		this.el=new double[nz];
		index=new int[nz];
		nz=0;
		for(int i=0;i<L;i++){
			if(el[i]!=0) {
				this.el[nz]=el[i];
				index[nz++]=i;
			}
		}
	}
	
	public SpVect(int I,int L){
		length=I;
		nzLength=L;
		el=new double[L];
		index=new int[L];
	}
	
	public SpVect(int I,int L,int flag){
		length=I;
		nzLength=L;
		if(flag==0)
		el=new double[L];
		else
		if(flag==1)
		index=new int[L];
	}

	
	public Vect add(Vect a){
		
		Vect b=a;
		for(int i=0;i<length;i++)
		b.el[index[i]]+=el[i];
		
		return b;
	}
	
	public Vect vectForm(){
	
		Vect b=new Vect(length);
		for(int i=0;i<this.nzLength;i++){
		b.el[index[i]]=el[i];
		}
		
		return b;
	}
	
	public void setEl(double a, int cl,int nzr){
		el[nzr]=a;
		index[nzr]=cl;
	
	}

	
	public void setNzEl(double a,int nzr){
		el[nzr]=a;
	}
	
	public void addToNz(double a,int nzr){
		el[nzr]+=a;
	}
	
	public int getNzLength(){
		return nzLength;
	
	}
	
	public int getLength(){
		return length;
	
	}
	
	
	public SpVect times(double a){
		
		SpVect b=this.deepCopy();
		for(int i=0;i<b.nzLength;i++){
		b.el[i]=el[i]*a;
	
		}
		
		return b;
	}
	
	public void trim(int nzLnew){
		
		SpVect b=new SpVect(length,nzLnew);
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
			if(el[i]!=0) nzL++;
		SpVect b=new SpVect(length,nzL);
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
		SpVect sorted=new SpVect(length,L);
		
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
	
	public double dot(SpVect v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		
		double dot=0;
		double vel;
		for(int i=0;i<nzLength;i++){
			vel=v.getEl(index[i]);
			if(vel!=0)
				dot+=el[i]*vel;
	
		}
		
		return dot;
	}
	
	

	public double dot(SpVect v,int n){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		int pass=0;
		double dot=0;
		for(int i=0;i<nzLength;i++){
			for(int j=pass+1;j<v.nzLength;j++)
				if(v.index[j]==index[i] && index[i]<=n)
				{
					pass=j;
				dot+=el[i]*v.el[j];
				break;
				}
				
			}
		
		return dot;
	}
	
	public double dot(Vect v,int n){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		int j;
		double dot=0;
		for(int i=0;i<nzLength;i++){
			j=index[i];
			if(j>n) break;
				dot+=el[i]*v.el[j];
				}
		
		return dot;
	}
	
	
	public double dotUp(Vect v,int n){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		int j;
		double dot=0;
		for(int i=0;i<nzLength;i++){
			j=index[i];
			if(j<n) continue;
				dot+=el[i]*v.el[j];
				}
		
		return dot;
	}
	
	
	public double dot(int n,double[] D,int[] index){
		
		double dot=0;
		int i=0;
			while(i<nzLength && index[i]<=n){
			dot+=pow(el[i],2)*D[index[i]];
			i++;
			}

		
		return dot;
	}
	
	
	public double dot(SpVect v,int n,double[] D,int[] index1, int[] index2){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		
		int indi;
		double dot=0;
		int i=0,j=0;
		while(i<nzLength && index1[i]<=n){
			indi=index1[i];
			j=0;
			while(j<v.nzLength && index2[j]<=n){
				if(index2[j]==indi){
			
					dot+=el[i]*v.el[j]*D[indi];
				}
				j++;
			}
		i++;
		}
		
		return dot;
	}
	
	public void addSame(SpVect v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

		for(int i=0;i<nzLength;i++)
			el[i]+=v.el[i];

	}
	
	
	public SpVect addGeneral(SpVect v){
		
		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		
	 SpVect sum= new SpVect(this.vectForm().add(v.vectForm()));
	 
	 return sum;
	 
	}
	
	public SpVect subGeneral(SpVect v){
		
		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");

	 SpVect sum= new SpVect(this.vectForm().sub(v.vectForm()));
	 
	 return sum;
	 
	}
	
	public void addSmaller(SpVect v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		for(int i=0;i<nzLength;i++)
			for(int j=0;j<v.nzLength;j++)
				if(index[i]==v.index[j])
					el[i]+=v.el[j];
	}
	


	public SpVect augh(SpVect v){
		int I1=getLength();
		int L1=getNzLength();

		SpVect w=new SpVect(I1+v.getLength(),L1+v.getNzLength());
		

		for(int i=0;i<nzLength;i++){
			w.el[i]=el[i];
			w.index[i]=index[i];
		}
		
		for(int i=0;i<v.nzLength;i++){
			w.el[i+L1]=v.el[i];
			w.index[i+L1]=v.index[i]+I1;
		}
		return w;
	}
	
public SpVect deepCopy(){
		
		int I1=getLength();
		int L1=getNzLength();
		SpVect w=new SpVect(I1,L1);
		for(int i=0;i<nzLength;i++){
			w.el[i]=el[i];
			w.index[i]=index[i];
		}
		w.length=length;
		w.nzLength=nzLength;
		
	return w;
	}

	
	public double dot(Vect v){

		if(length!=v.length) throw new IllegalArgumentException("Vectrs have different lengths.");
		
		double dot=0;
		for(int i=0;i<nzLength;i++){
			if(v.el[index[i]]!=0)
				dot+=el[i]*v.el[index[i]];
	
		}
		
		return dot;
	}
	
	public double norm(){
		double s=0;
		for(int i=0;i<nzLength;i++)
		s+=el[i]*el[i];
		
		return Math.sqrt(s);
	}
	
	public void Ddotv(Vect D){
		if(length!=D.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		for(int i=0;i<nzLength;i++)
			el[i]=el[i]*D.el[index[i]];
	}
	
	public void show(){
		
		for(int i=0;i<length;i++)
			System.out.format("%12.4f\n", getEl(i));
	}
	
	public void hshow(){
		
		for(int i=0;i<length;i++)
			System.out.format("%12.4f\t", getEl(i));
		System.out.println();
	}
	
	public void shownz(){
		
		util.hshow(el);
	}
	
	public void shownzA(){
		System.out.format(" Sparse Vector ********  dimention: %d  size: %d \n",length, nzLength);
	
				for(int j=0;j<nzLength;j++)
				System.out.format(" ( %d )  %25.12e\n",index[j], el[j]);
	}
	
	
public void showr(){
		
		util.hshow(index);
	}

	public double getEl(int i){
		double eli=0;
		for(int j=0;j<nzLength;j++){
			if(index[j]==i){
				eli=el[j];
				break;
			}
		}
		return eli;
	}
	
	public void extend(int L){
		SpVect vs=this.deepCopy();
		nzLength+=L;
		el=new double[nzLength];
		index=new int[nzLength];
		for(int i=0;i<vs.nzLength;i++){
		this.el[i]=vs.el[i];
		this.index[i]=vs.index[i];
		}
	}
	
}
