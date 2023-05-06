package math;
import java.util.Arrays;

import math.SpVect;
import static java.lang.Math.*;
import math.Vect;
import fem.Model;

public class SpBlockMat  {

	public SpBlockVect[] row;
	private int nRow;
	public SpBlockMat(){}

	public SpBlockMat(int nRow){
		this.nRow=nRow;
		row=new SpBlockVect[nRow];
	}

	public SpBlockMat(int nRow,int I){
		this.nRow=nRow;
		row=new SpBlockVect[nRow];
		for(int i=0;i<nRow;i++)
			row[i]=new SpBlockVect(I);
	}

	public SpBlockMat(int nRow,int I,int L){
		this.nRow=nRow;
		row=new SpBlockVect[nRow];
		for(int i=0;i<nRow;i++)
			row[i]=new SpBlockVect(I,L);
	}

	public SpBlockMat(int nRow,int I,int L,int m, int n){
		this.nRow=nRow;
		row=new SpBlockVect[nRow];
		for(int i=0;i<nRow;i++)
			row[i]=new SpBlockVect(I,L,m,n);
	}

	public SpBlockMat deepCopy(){
		SpBlockMat M=new SpBlockMat(nRow);
		for(int i=0;i<nRow;i++)
			M.row[i]=row[i].deepCopy();

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



	public Mat matForm(){
		Mat M=new Mat(3*getnRow(),3*getnCol());
		for(int i=0;i<nRow;i++)
			for(int id=0;id<3;id++)
				for(int k=0;k<row[i].nzLength;k++)
					for(int jd=0;jd<=((k==row[i].nzLength-1)?id:2);jd++){
						M.el[3*i+id][3*row[i].index[k]+jd]=row[i].el[k].el[id][jd];
					}

		return M;
	}

	/**
	 * TODO Converts the sparse block matrix into a real block matrix.
	 *
	 * @return
	 */
	public SpMat spMatForm(){
		SpMat M=new SpMat(3*getnRow());
		int rowNumb,nzColumn;
		for(int i=0;i<getnRow();i++){
			for(int id=0;id<3;id++){
				rowNumb=3*i+id;
				M.row[rowNumb]=new SpVect(3*getnCol(),3*(row[i].nzLength-1)+1+id);
				for(int k=0;k<row[i].nzLength;k++)
					for(int jd=0;jd<=((k==row[i].nzLength-1)?id:2);jd++){
						nzColumn=3*k+jd;
						M.row[rowNumb].el[nzColumn]=row[i].el[k].el[id][jd];
						M.row[rowNumb].index[nzColumn]=3*row[i].index[k]+jd;


					}

			}
		}

		return M;
	}

	public SpMat spMatForm(Model model,int[][] ind){
		
		int dim=model.dim;
		int nComp=model.numberOfUnknownUcomp;

		int[] nz=new int[ind.length];
		for(int i=0;i<getnRow();i++)
			for(int k=0;k<row[i].nzLength;k++)
				for(int p=0;p<dim;p++)
				if(ind[row[i].index[k]][p]>0)
					nz[i]++;

		SpMat M=new SpMat(nComp);
		
		int rowIndex=0,nzColumn=0;
		for(int i=0;i<getnRow();i++){

			for(int id=0;id<dim;id++){

				if(ind[i][id]==0) continue;
				rowIndex=ind[i][id]-1;
				nzColumn=0;
				int c=0;

				for(int q=id+1;q<dim;q++) if(ind[i][q]>0) c++;
		
				M.row[rowIndex]=new SpVect(nComp,nz[i]-c);
				for(int k=0;k<row[i].nzLength;k++){

					for(int jd=0;jd<=((k==row[i].nzLength-1)?id:dim-1);jd++){
						if(ind[row[i].index[k]][jd]==0) 
										continue;
						M.row[rowIndex].el[nzColumn]=row[i].el[k].el[id][jd];
						M.row[rowIndex].index[nzColumn]=ind[row[i].index[k]][jd]-1;
						nzColumn++;
											}
					}

			}
		}

		return M;
	}



	public void reflect(){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");
		int[] nz=new int[nRow];

		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].nzLength;j++){
				if(row[i].index[j]>=i) break;
				nz[i]++;
			}

		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].nzLength;j++){
				int cl=row[i].index[j];
				if(cl>=i)
					break;
				row[cl].el[nz[cl]]=row[i].el[j];
				row[cl].index[nz[cl]]=i;
				nz[cl]++;

			}

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


	public void add(SpBlockMat M){

		for(int i=0;i<nRow;i++)
			row[i].addSmaller(M.row[i]);
	}




	public void times(double a){
		for(int i=0;i<nRow;i++)
			row[i]=row[i].times(a);
	}

	public void show(){
		SpMat M=spMatForm();
		M.show();
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
