package math;
import java.util.Arrays;

import math.SpVect;
import static java.lang.Math.*;
import math.Vect;

public class SpMatAsym  {

	public SpVect[] row;
	public int nRow;
	public SpMatAsym(){}

	public SpMatAsym(int nRow){
		this.nRow=nRow;
		row=new SpVect[nRow];
	}

	public SpMatAsym(int nRow,int I){
		this.nRow=nRow;
		row=new SpVect[nRow];
		for(int i=0;i<nRow;i++)
			row[i]=new SpVect(I);
	}

	public SpMatAsym(int nRow,int I,int L){
		this.nRow=nRow;
		row=new SpVect[nRow];
		for(int i=0;i<nRow;i++)
			row[i]=new SpVect(I,L);
	}

	public SpMatAsym (Mat M){
		nRow=M.nRow;
		row=new SpVect[nRow];
		for(int i=0;i<nRow;i++){
			row[i]=new SpVect(M.el[i]);
		}
	}

	public SpMatAsym deepCopy(){
		SpMatAsym M=new SpMatAsym(nRow);
		for(int i=0;i<nRow;i++){
			M.row[i]=row[i].deepCopy();
		}

		return M;
	}
	

	
	public void sortAndTrim(int[]L){

		for(int i=0;i<nRow;i++)
			row[i].sortAndTrim(L[i]);
	}

	public void trim(int[] L){

		for(int i=0;i<nRow;i++)
			row[i].trim(L[i]);
	}

	public void trim(){

		for(int i=0;i<nRow;i++)
			row[i].trim();
	}
	

	public SpMat iLDU(Vect D,double domfact){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		SpMat L=new SpMat(nRow);
		double s;

		int j=0,kd=0;
		s=0;

		for(int i=0;i<nRow;i++){
			kd=row[i].nzLength;
			L.row[i]=new SpVect(nRow,kd,0);

			j=0;

			for(int k=0;k<kd;k++){
				j=row[i].index[k];
				s=L.row[j].dot(j-1,D.el,row[j].index);
				D.el[j]=domfact*row[i].el[kd-1]-s;

				s=L.row[i].dot(L.row[j],j-1,D.el,row[i].index,row[j].index);
				L.row[i].el[k]=(row[i].el[k]-s)/D.el[j];
			}

		}

	
		return L;


	}


	public Mat matForm(){
		Mat M=new Mat(nRow,getnCol());
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				M.el[i][row[i].index[k]]=row[i].el[k];
				//M.el[row[i].index[k]][i]=row[i].el[k];
			}

		return M;
	}

	
	

	public SpMat transpose(int cLmax){

		SpMat T=new SpMat(getnCol(),nRow,cLmax);
		int[] nz=new int[nRow];
		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].nzLength;j++){
				int cl=row[i].index[j];
				T.row[cl].el[nz[cl]]=row[i].el[j];
				T.row[cl].index[nz[cl]++]=i;
			}

		T.sortAndTrim(nz);
		return T;
	}


	public int getnRow(){

		return nRow;
	}

	public int getnCol(){

		return row[0].getLength();
	}

	public int getNzWidth(){

		return row[0].getNzLength();
	}

	public void size(){

		System.out.format("%10d%10d%10d\n",getnRow(),getnCol(),getNzWidth());
	}


	public Vect mul(Vect u){
		Vect v=new Vect(getnRow());

		int j;
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				j=row[i].index[k];
				v.el[i]+=this.row[i].el[k]*u.el[j];
			}
		return v;
	}
	



	public void addSmaller(SpMat M){

		for(int i=0;i<nRow;i++)
			row[i].addSmaller(M.row[i]);
	}
	
	public SpMatAsym addNew(SpMat M){
		SpMatAsym B=this.deepCopy();
		for(int i=0;i<nRow;i++)
			B.row[i].addSmaller(M.row[i]);
		
		return B;
	}
	
	public SpMatAsym addAugNew(SpMat M){
		
		int I1=this.nRow;
		SpMatAsym B=new SpMatAsym(I1,I1);
		
		for(int i=0;i<nRow;i++){
			
		int L1=this.row[i].getNzLength();
		int L2=M.row[i].getNzLength();
		
		SpVect w=new SpVect(I1,L1+L2);
		
		for(int j=0;j<this.row[i].nzLength;j++){
			w.el[j]=this.row[i].el[j];
			w.index[j]=this.row[i].index[j];
		}
		
		for(int j=0;j<M.row[i].nzLength;j++){
			w.el[j+L1]=M.row[i].el[j];
			w.index[j+L1]=M.row[i].index[j]+row[i].index[L1-1];
		}
			B.row[i]=w.deepCopy();
		}
		return B;
	}





	public Vect diag(){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		Vect D=new Vect(nRow);
		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].getNzLength();j++)
				if(row[i].index[j]==i){
					D.el[i]=row[i].el[j];
					break;
				}
		return D;
	}
	
	public SpMatAsym offDiag(){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");
		SpMatAsym Ls=this.deepCopy();
		
		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].getNzLength();j++){
				if(row[i].index[j]==i)
					Ls.row[i].el[j]=0;
			}
			

		return Ls;
	}


	public double max(){
		double elMax=0;
		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].getNzLength();j++)
				if(row[i].el[j]>elMax) elMax=row[i].el[j];

		return elMax;
	}
	
	public double maxAbs(){
		double elMax=0;
		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].getNzLength();j++)
				if(abs(row[i].el[j])>elMax) elMax=row[i].el[j];

		return elMax;
	}


	

	public Vect scale(Vect b){

		int I=getnRow();
		int J=getnCol();
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		Vect invDS=this.diag().abs().sqrt().inv();

		DMD(invDS);

		b.timesVoid(invDS);

		return invDS;
	}
	public Vect scale(){

		int I=getnRow();
		int J=getnCol();
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		Vect invDS=this.diag().abs().abs().sqrt().inv();

		DMD(invDS);


		return invDS;
	}





	public void DMD(Vect D){}


	public void times(double a){
		for(int i=0;i<nRow;i++)
			row[i]=row[i].times(a);
	}

	public SpMatAsym timesNew(double a){
		SpMatAsym B=this.deepCopy();
		B.times(a);
		return B;
	}


	public void show(){
		Mat M=matForm();
		M.show();
	}

	public void shownz(){

		for(int i=0;i<nRow;i++)
			row[i].shownz();
	}

	public void showcl(){

		for(int i=0;i<nRow;i++)
			row[i].showr();
	}




	public int search(int[] A,int ic,int a){
		int m=-1;
		for(int i=0;i<=ic;i++){
			if(A[i]==a){
				m=i;
				break;
			}
		}
		return m;
	}

	public int search(int[] A,int a){
		int m=-1;
		for(int i=0;i<A.length;i++){
			if(A[i]==a){
				m=i;
				break;
			}
		}
		return m;
	}

}
