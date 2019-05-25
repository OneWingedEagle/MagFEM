package math;


import static java.lang.Math.*;

import java.util.Arrays;

public class VectComp implements Cloneable{
	public VectComp(){};
	public Complex[] el;
	public int length;	
	

	public static void main(String[] args){}
	
	

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
	
	public VectComp(int I){
		this.length=I;
		this.el=new Complex[I];
		for(int i=0;i<I;i++)
			el[i]=new Complex();
	}

	public VectComp(Complex[] array){
		this.length=array.length;
		this.el=Arrays.copyOf(array, length);
	}

	public VectComp(Vect v){
		this.length=v.length;
		this.el=new Complex[v.length];
		for(int i=0;i<v.length;i++)
		this.el[i]=new Complex(v.el[i], 0);
	}

	
	
	
	public VectComp deepCopy(){
		VectComp w=new VectComp(length);
		w.el=Arrays.copyOf(el, length);
		return w;
		}

	
	public VectComp add(VectComp v){
	
		if(this.length!=v.length) throw new NullPointerException("vectrs have different lengths");
		VectComp w=new VectComp(v.length);
		for(int i=0;i<v.length;i++)
			w.el[i]=this.el[i].add(v.el[i]);
		return w;
	}
	
	public VectComp sub(VectComp v){
		
		if(this.length!=v.length) throw new NullPointerException("vectrs have different lengths");
		VectComp w=new VectComp(v.length);
		for(int i=0;i<v.length;i++)
			w.el[i]=this.el[i].sub(v.el[i]);
		return w;
	}
	
	public VectComp aug(VectComp v){
		
		int L=this.length+v.length;
		VectComp w=new VectComp(L);
		for(int i=0;i<this.length;i++)
			w.el[i]=this.el[i];
    	for(int i=0;i<v.length;i++)
			w.el[i+this.length]=v.el[i];
    	return w;
	}
	
	


public VectComp mul(VectComp u){
	if(length!=u.length) throw new IllegalArgumentException("vectrs have different lengths");
	VectComp v=new VectComp(length);
	for(int i=0;i<v.length;i++)
		v.el[i]=el[i].times(u.el[i]);
	return v;
}



	

	public VectComp times(double a){
		
		VectComp v=new VectComp(this.length);
		for(int i=0;i<this.length;i++)
			v.el[i]=this.el[i].times(a);
		return v;
	}
	
	public VectComp times(Complex a){
		
		VectComp v=new VectComp(this.length);
		for(int i=0;i<this.length;i++)
			v.el[i]=this.el[i].times(a);
		return v;
	}
	
	public void timesVoid(double a){
		
		for(int i=0;i<this.length;i++)
			this.el[i]=this.el[i].times(a);
	}
	

	
	public Complex dot(VectComp u){
	
		if(this.length!=u.length) throw new IllegalArgumentException("vectrs have different lengths");
		Complex s=new Complex(0,0);
		for(int i=0;i<u.length;i++)
			s=s.add(this.el[i].times(u.el[i]));
		return s;
	}
	
	public double norm2(){
		
		double s=0;
		for(int i=0;i<this.length;i++)
			s=s+this.el[i].times(this.el[i].conj()).re;
		return s;
	}
	
	public Vect abs(){
		
		 Vect vn=new Vect(this.length);
		for(int i=0;i<this.length;i++)
			vn.el[i]=this.el[i].norm();
		
		return vn;
	}
	
	public void timesVoid(Vect D){
		if(length!=D.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		for(int i=0;i<length;i++)
			el[i]=el[i].times(D.el[i]);
	}

	public VectComp times(Vect D){
		VectComp v=new VectComp(this.length);
		if(length!=D.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		for(int i=0;i<length;i++)
			v.el[i]=el[i].times(D.el[i]);
		
		return v;
	}
	
	public VectComp conj(){
		VectComp v=new VectComp(this.length);
		
		for(int i=0;i<length;i++)
			v.el[i]=el[i].conj();
		
		return v;
	}


	public VectComp div(Vect D){
		VectComp v=new VectComp(this.length);
		if(length!=D.length) throw new IllegalArgumentException("vectrs have different lengths");
		
		for(int i=0;i<length;i++){
			if(D.el[i]==0){throw new IllegalArgumentException("Divided by zerp error: Entry "+Integer.toString(i) +" is zero.");}
			v.el[i]=el[i].times(1.0/D.el[i]);
		}
		return v;
	}
	


	
	public double norm(){
	double s=0;
		for(int i=0;i<this.length;i++)
			s=s+this.el[i].norm2();
		return Math.sqrt(s);
	}
	

	
	public void show(){
		int I;
		I=this.length;
		for(int i=0;i<I;i++)
			this.el[i].show();
		System.out.println();

	}
	public void length(){

		System.out.println(length);

	}

}
