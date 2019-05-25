package math;


import static java.lang.Math.*;

import java.util.Arrays;

public class IntVect implements Cloneable{
	public IntVect(){};
	public int[] el;
	public int length;	
		
	

	public IntVect(int I){
		this.length=I;
		this.el=new int[I];
	}

	public IntVect(int[] array){
		this.length=array.length;
		this.el=Arrays.copyOf(array, length);
	}

	public IntVect(int x,int y,int z){
		length=3;
		el=new int[length];
		el[0]=x; el[1]=y; el[2]=z;
	}
	
	public IntVect(int x,int y){
		length=2;
		el=new int[length];
		el[0]=x; el[1]=y;;
	}
	


	
	public IntVect deepCopy(){
		IntVect w=new IntVect(length);
		w.el=Arrays.copyOf(el, length);
		return w;
		}

	public void set(int[] u){
		for(int i=0;i<u.length;i++)
		this.el[i]=u[i];
		this.length=u.length;
	}
	
	public void extend(int L){
		int[] a=Arrays.copyOf(el,length);
		length+=L;
		el=new int[length];
		for(int i=0;i<a.length;i++)
		this.el[i]=a[i];
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

	
	public Vect ones(int I){
		
		Vect v=new Vect(I);
		for(int i=0;i<I;i++)
			v.el[i]=1;
		return v;
	}

	public void hshow(){

util.hshow(this.el);

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
		                int temp = el[i];  el[i] = el[i+1];  el[i+1] = temp;
		               int tempi=indice[i]; indice[i]=indice[i+1]; indice[i+1]=tempi;
		            }
		        }
		    }
	
		
		return indice;
		}
	
	
	int partition(int arr[], int left, int right)

	{
	      int i = left, j = right;
	      int tmp;
	      int pivot = arr[(left + right) / 2];
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

	public void quickSort(int arr[], int left, int right) {
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

}
