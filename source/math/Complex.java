package math;


import static java.lang.Math.*;

import java.text.DecimalFormat;
import java.util.Arrays;

public class Complex{
	
	public double re,im;
	public Complex(){
			};
	
	public Complex(double a, double b){
	
	re=a; im=b;
	
	}
	
	public Complex(double[] a){
		
		re=a[0]; im=a[1];
		
		}
	    
	public double  norm(){
		
		Vect v=new Vect(re,im) ;
		double r=v.norm();
		return r;
		}
	
	public double  norm2(){
		
		Vect v=new Vect(re,im) ;
		double r=v.dot(v);
		return r;
		}

	public double  ang(){
		
		Vect v=new Vect(re,im) ;
		double ang=util.getAng(v);		
		return ang;
		}
	
	public Complex  add(Complex c){
		
		Complex z=new Complex(c.re+this.re,c.im+this.im);	
		return z;
		}

	public Complex  sub(Complex c){
		
		Complex z=new Complex(this.re-c.re,this.im-c.im);	
		return z;
		}
	
public Complex  conj(){
		
		Complex z=new Complex(this.re,-this.im);	
		return z;
		}

public Complex  inv(){
	
	Complex z=this.conj().times(1.0/this.norm2());	
	return z;
	}

public Complex  pow(double p){
	
	double zn=Math.pow(this.norm(),p);
	double phi=this.ang()*p;
	Complex z=new Complex(zn*Math.cos(phi), zn*Math.sin(phi));
	return z;
	}

	public Complex  times(double c){
		
		Complex z=new Complex(c*this.re,c*this.im);	
		return z;
		}
	
	public Complex  deepCopy(){
		
		Complex z=new Complex(this.re,this.im);	
		return z;
		}
	
	public Complex  times(Complex c){
		double a=this.re*c.re-this.im*c.im;
		double b=this.re*c.im+this.im*c.re;
		Complex z=new Complex(a,b);	
		return z;
		}
	
	public void  show(){
		
		DecimalFormat df=new DecimalFormat("0.000E00");
		String s="+";
		if(this.im<0) s="-";
		util.pr(df.format(this.re)+" "+s+"j"+df.format(abs(this.im)));
		}
	
	public static void main(String[] args){
		
		Complex c=new Complex(1,1);
		
		double r=c.norm();
		util.pr(r);
	}
	
	
}
