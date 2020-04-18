package math;


import static java.lang.Math.*;

import java.util.Arrays;

public class Vect{
	public Vect(){};
	public double[] el;
	public int length;	
	

	public static void main(String[] args){

		int m=1,K=100;
		double[][] lim= new double[m][2];
		
		lim[0][0]=0;
		lim[0][1]=5;
/*		lim[1][0]=0;
		lim[1][1]=8;*/
/*		lim[2][0]=0;
		lim[2][1]=8;
		lim[3][0]=0;
		lim[3][1]=5;*/
/*		lim[4][0]=0;
		lim[4][1]=8;
		lim[5][0]=0;
		lim[5][1]=6;
		lim[6][0]=0;
		lim[6][1]=10;*/

		
		double[] dd=new double[m];
		for(int h=0;h<m;h++)
			dd[h]=(lim[h][1]-lim[h][0])/K;
		
		Vect x=new Vect(m);
		int[] c=new int[m];
		double f=0;
		
		for(c[0]=0;c[0]<K;c[0]++)
			//for(c[1]=0;c[1]<K;c[1]++)
				/*for(c[2]=0;c[2]<K;c[2]++)
					for(c[3]=0;c[3]<K;c[3]++)*/
		/*				for(c[4]=0;c[4]<K;c[4]++)
							for(c[5]=0;c[5]<K;c[5]++)
							for(c[6]=0;c[6]<K;c[6]++)*/
							{
									for(int h=0;h<m;h++)
									x.el[h]=c[h]*dd[h];
									double p=1;
									for(int h=0;h<m;h++)
										p*=dd[h];
									
								f+=ff(x)*p;
							}
		
util.pr("integral1= "+f);


f=0;
int L=1000;
for(int k=0;k<L;k++){
for(int h=0;h<m;h++)
	x.el[h]=lim[h][0]+random()*(lim[h][1]-lim[h][0]);
/*if(k==0)
	f+=ff(x);
else
f+=(f/(k)+ff(x))/2;*/
f+=ff(x);
}

double p=1;
for(int h=0;h<m;h++)
	p*=(lim[h][1]-lim[h][0]);
util.pr("integral2= "+f*p/L);
		int N=50;
Vect v1=new Vect(N);
	v1.rand();
Vect v2=v1.deepCopy();
//v1.bubble();
v2.quickSort();
System.out.println("done");

Mat A=new Mat(5,5);
A.symrand(5);
A.ddom(1,1);

for(int i=0;i<A.nRow;i++)
	for(int j=0;j<A.nRow;j++)
		if(A.el[i][j]<.5) A.el[i][j]=0;
		else  A.el[i][j]=floor(5* A.el[i][j]);
A.show();




SpMat As=new SpMat(A);

SpMat Bs=new SpMat(A.nRow,A.nRow);
for(int i=0;i<A.nRow;i++){
	int nz=0;
	for(int j=0;j<As.row[i].nzLength;j++){
		nz++;
		if(As.row[i].index[j]==i) break;
	}
	 Bs.row[i]=new SpVect(A.nRow,nz);
	 
	 for(int j=0;j<nz;j++){
		 Bs.row[i].index[j]= As.row[i].index[j];
		 Bs.row[i].el[j]= As.row[i].el[j];
	 }
} 
/*Bs.show();
Bs.shownz();*/

SpMat Cs=Bs.reflect(8);
Cs.show();

	}
	
	

	public static double ff(Vect x){
		double f=0;
		
		//f=3*x.el[0]*pow(x.el[2],2)*x.el[3]*tan(x.el[6])+pow(x.el[1],1.5)*log(Math.abs(x.el[4]))/10+x.el[1]*x.el[2]*x.el[5];
		f=x.dot(x)*x.dot(x);/*+x.outer(x).mul(x).norm2()*/;
		f=exp(-2*x.el[0]);
		
		return f;
	}
	
	public static double ff2(Vect x){
		double f=0;
		
	for(int i=0;i<x.length;i++)
		f+=exp(-Math.abs(x.el[i]));
		return f;
	}
	
	public Vect(int I){
		this.length=I;
		this.el=new double[I];
	}

	public Vect(double[] array){
		this.length=array.length;
		this.el=Arrays.copyOf(array, length);
	}

	public Vect(double x,double y,double z){
		length=3;
		el=new double[length];
		el[0]=x; el[1]=y; el[2]=z;
	}
	
	public Vect(double x,double y){
		length=2;
		el=new double[length];
		el[0]=x; el[1]=y;;
	}
	
	
	
	public void set(double[] u){
	
		this.length=u.length;		
		this.el=Arrays.copyOf(u,this.length);
	}
	
	public Vect fill(SpVect vs){
		Vect v=new Vect(length);
		for(int i=0;i<vs.nzLength;i++)
			v.el[vs.index[i]]=vs.el[i];
		return v;
		}
	
	public Vect deepCopy(){
		Vect w=new Vect(length);
		w.el=Arrays.copyOf(el, length);
		return w;
		}

	public void set(int[] u){
		for(int i=0;i<u.length;i++)
		this.el[i]=(double)u[i];
		this.length=u.length;
	}
	
	public Vect add(Vect v){
	
		if(this.length!=v.length) throw new NullPointerException("vectrs have different lengths");
		Vect w=new Vect(v.length);
		for(int i=0;i<v.length;i++)
			w.el[i]=this.el[i]+v.el[i];
		return w;
	}
	
	public Vect sub(Vect v){
		
		if(this.length!=v.length) throw new NullPointerException("vectrs have different lengths");
		Vect w=new Vect(v.length);
		for(int i=0;i<v.length;i++)
			w.el[i]=this.el[i]-v.el[i];
		return w;
	}
	
	public Vect aug(Vect v){
		
		int L=this.length+v.length;
		Vect w=new Vect(L);
		for(int i=0;i<this.length;i++)
			w.el[i]=this.el[i];
    	for(int i=0;i<v.length;i++)
			w.el[i+this.length]=v.el[i];
    	return w;
	}
	
	
	public void rand(){
	for(int i=0;i<length;i++)
			el[i]=random();

	}
	
	public Vect ones(int I){
		
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=1;
		return v;
	}

	public Vect flr(){
		int I=length;
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=floor(el[i]);
		return v;
	}
	
public Vect rand(int I,double a, double b){
		
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=a+(b-a)*random();
		return v;
	}

public Vect mul(Vect u){
	if(length!=u.length) throw new IllegalArgumentException("vectrs have different lengths");
	Vect v=new Vect(length);
	for(int i=0;i<v.length;i++)
		v.el[i]=el[i]*u.el[i];
	return v;
}

public Mat outer(Vect u){
	if(length!=u.length) throw new IllegalArgumentException("vectrs have different lengths");
Mat M=new Mat(length,length);
for(int i=0;i<length;i++)
	for(int j=0;j<length;j++)
		M.el[i][j]=this.el[i]*u.el[j];
	return M;
}
public Mat outer(){
	Vect u=this.deepCopy();
	return outer(u);
}

public Mat mul(Mat A){
	if(length!=A.nRow) throw new IllegalArgumentException("Array dimensions must agree.");
	Mat M=new Mat(A.size());
	for(int i=0;i<A.nRow;i++)
		for(int j=0;j<A.nCol;j++)
		M.el[i][j]=A.el[i][j]*this.el[i];
	return M;
}

	public void rand(double a, double b){

		for(int i=0;i<length;i++)
			el[i]=a+(b-a)*random();
	}
	
public Vect rand(int I,int a, int b){
		
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=a+(b-a)*random();
		return v;
	}


	public Vect linspace(double a, double b,int I){
		
		Vect v=new Vect(I);
		double d=(b-a)/(I-1);
		v.el[0]=a;
		for(int i=1;i<I-1;i++)
			v.el[i]=a+i*d;
		v.el[I-1]=b;
		return v;
	}
	
	public Vect linspace(double a, double b,double d){
		double r=(b-a)/d;
		int N=(int)r-1;
		
		if(N<2) N=2;
		return linspace(a,b,N);
	}
	
	public Vect sqspace(double a, double b,int I){
		
		Vect v=new Vect(I);
		v.el[0]=0;
		for(int i=1;i<I;i++)
			v.el[i]=v.el[i-1]+1+0.2*Math.abs(i-(double)I/2);
		v=v.times((b-a)/v.el[I-1]);
		for(int i=0;i<I;i++)
			v.el[i]+=a;
		return v;
	}
	
public Vect randspace(double a, double b,int I,double r){
		
		Vect v=new Vect();
		double d=(b-a)/(I-1);
		v=v.rand(I,-r*d,r*d);	
		v=v.add(v.linspace(a,b,I));
		v.el[0]=a;
		v.el[I-1]=b;

		return v;
	}
	
	public Vect times(double a){
		
		Vect v=new Vect(this.length);
		for(int i=0;i<this.length;i++)
			v.el[i]=a*this.el[i];
		return v;
	}
	
	public void timesVoid(double a){
		
		for(int i=0;i<this.length;i++)
			this.el[i]=a*this.el[i];
	}
	
public Vect add(double a){
		
		Vect v=new Vect(length);
		for(int i=0;i<length;i++)
			v.el[i]=a+el[i];
		return v;
	}

	public Vect times(int a){
		
		Vect v=new Vect(this.length);
		for(int i=0;i<this.length;i++)
			v.el[i]=a*this.el[i];
		return v;
	}
	
	public Vect abs(){
		
		Vect v=new Vect(length);
		for(int i=0;i<length;i++)
			v.el[i]=Math.abs(el[i]);
		return v;
	}
	
	public Vect sqrt(){
		
		Vect v=new Vect(this.length);
		for(int i=0;i<length;i++){
			if(this.el[i]<0) throw new IllegalArgumentException("Square root error: Entry"+Integer.toString(i) +" is less than zero.");
			v.el[i]=Math.sqrt(this.el[i]);
		}
		return v;
	}
	
	public double sum(){
		 double sum=0;
		for(int i=0;i<length;i++)
			sum+=el[i];
			
		return sum;
	}
	
	public double max(){
		 double max=el[0];
		for(int i=1;i<length;i++)
			if(el[i]>max)
				max=el[i];
			
		return max;
	}
	
	
	
	public double min(){
		 double min=el[0];
		for(int i=1;i<length;i++)
			if(el[i]<min)
				min=el[i];
			
		return min;
	}
	
	public double dot(Vect u){
	
		if(this.length!=u.length) throw new IllegalArgumentException("vectrs have different lengths");
		double s=0;
		for(int i=0;i<u.length;i++)
			s=s+this.el[i]*u.el[i];
		return s;
	}
	
	public double norm2(){
		
		double s=0;
		for(int i=0;i<this.length;i++)
			s=s+this.el[i]*this.el[i];
		return s;
	}
	
	public Vect inv(){
		Vect v=new Vect(this.length);
		for(int i=0;i<length;i++){
			if(this.el[i]==0) throw new IllegalArgumentException("Divided by zerp error: Entry "+Integer.toString(i) +" is zero.");
			v.el[i]=1/el[i];
		}
		return v;
	}
	public void timesVoid(Vect D){
		if(length!=D.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		for(int i=0;i<length;i++)
			el[i]=el[i]*D.el[i];
	}
	
	public void extend(Vect D){
		Vect temp=this.deepCopy();
		
		this.length+=D.length;

		el=new double[this.length];
		
		for(int i=0;i<temp.length;i++)
			el[i]=temp.el[i];
		
		int ix=0;
		for(int i=0;i<D.length;i++)
			el[ix++]=temp.el[i];
	}

	public Vect times(Vect D){
		Vect v=new Vect(this.length);
		if(length!=D.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		for(int i=0;i<length;i++)
			v.el[i]=el[i]*D.el[i];
		
		return v;
	}

	public Vect v3(){
		int L=this.length;
		if(L>2) return this;
		Vect v=new Vect(3);
		for(int i=0;i<length;i++)
			v.el[i]=el[i];
		
		return v;
	}
	public void zero(){
		for(int i=0;i<length;i++)
			this.el[i]=0;
	}
	public void zero(int start,int end){
		for(int i=start;i<=end;i++)
			this.el[i]=0;
	}
	public Vect v2(){
		return new Vect(el[0],el[1]);
	}
	
	public Vect div(Vect D){
		Vect v=new Vect(this.length);
		if(length!=D.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		for(int i=0;i<length;i++){
			if(D.el[i]==0){throw new IllegalArgumentException("Divided by zerp error: Entry "+Integer.toString(i) +" is zero.");}
			v.el[i]=el[i]/D.el[i];
		}
		return v;
	}
	
	public Vect div0(Vect D){
		Vect v=new Vect(this.length);
		if(length!=D.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		for(int i=0;i<length;i++){
			if(D.el[i]==0) v.el[i]=0;
			else
			v.el[i]=el[i]/D.el[i];
		}
		return v;
	}
	
	public Vect cross(Vect u){
		
		if(length>3) throw new IllegalArgumentException("Cross product is not defined for dimentions higher than 3");
		if(length!=u.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		if(length==2)
			return new Vect(0,0,el[0]*u.el[1]-el[1]*u.el[0]);
		
		Vect w=new Vect(3);
			w.el[0]=el[1]*u.el[2]-el[2]*u.el[1];
			w.el[1]=el[2]*u.el[0]-el[0]*u.el[2];
			w.el[2]=el[0]*u.el[1]-el[1]*u.el[0];
			return w;
		}

	
	public double norm(){
	double s=0;
		for(int i=0;i<this.length;i++)
			s=s+this.el[i]*this.el[i];
		return Math.sqrt(s);
	}
	
	public void normalize(){
		double nr=this.norm();
		
		if(nr==0) 	this.timesVoid(0);
		double normRev=1/nr;
		this.timesVoid(normRev);
	}
	
	public Vect normalized(){
		double nr=this.norm();
		
		if(nr==0) return new Vect(this.length);
		double normRev=1/nr;
		return this.times(normRev);
	}
	
	public void show(){
		int I;
		I=this.length;
		for(int i=0;i<I;i++)
				System.out.format("%10.8e\n",this.el[i]);
		System.out.println();

	}
	public void length(){

		System.out.println(length);

	}
	public void hshow(){
		int I;
		I=this.length;
		double am=this.abs().max();
		int nf=(int)log10(am);
		if(nf<5)
		for(int i=0;i<I;i++)
				System.out.format("%15.8f\t",this.el[i]);
		else if(nf<10)
			for(int i=0;i<I;i++)
				System.out.format("%20.8f\t",this.el[i]);
		else
			for(int i=0;i<I;i++)
				System.out.format("%25.8f\t",this.el[i]);
		System.out.println();

	}

	
	public int[] bubble(){
		if(length<1) throw new IllegalArgumentException("Null vector");
		int[] indice=new int[length];
		for(int i=0;i<length;i++)
			indice[i]=i;
		
		if(length==1) return indice;
		    int n = el.length;
		    for (int pass=0; pass < n-1; pass++) { 
		        for (int i=0; i < n-1-pass; i++) {
		            if (el[i] > el[i+1]) {
		                double temp = el[i];  el[i] = el[i+1];  el[i+1] = temp;
		               int tempi=indice[i]; indice[i]=indice[i+1]; indice[i+1]=tempi;
		            }
		        }
		    }
	
		
		return indice;
		}
	
	
	int partition(double arr[], int left, int right)

	{
	      int i = left, j = right;
	      double tmp;
	      double pivot = arr[(left + right) / 2];
	      while (i <= j) {
	            while (arr[i] < pivot)
	                  i++;
	            while (arr[j] > pivot)
	                  j--;
	            if (i <= j) {
	                  tmp = arr[i];
	                  arr[i] = arr[j];
	                  arr[j] = tmp;
	                  i++;
	                  j--;
	            }
	      };

	      return i;
	}

	public void quickSort(double arr[], int left, int right) {
	      int index = partition(arr, left, right);
	      if (left < index - 1)
	            quickSort(arr, left, index - 1);
	      if (index < right)
	            quickSort(arr, index, right);
	}
	
	public void quickSort() {
	      int index = partition(this.el, 0, this.length-1);
	      int left=0;
	      if (left < index - 1)
	            quickSort(this.el, left, index - 1);
	      if (index < this.length-1)
	            quickSort(this.el, index, this.length-1);

	}

	
	private double trunc(double a,int n){
		return floor(a*pow(10,n))/pow(10,n);
	}
	public void trunc(int n){
		for(int i=0;i<length;i++)
			el[i]=trunc(el[i],n);
	}
	
}
