package math;
import java.util.Arrays;


import static java.lang.Math.*;

public class SpMatComp  {

	public SpVectComp[] row;
	public int nRow;
	public boolean hermit=true, symm=false;
	public SpMatComp(){}

	public SpMatComp(int nRow){
		this.nRow=nRow;
		row=new SpVectComp[nRow];
	}

	public SpMatComp(int nRow,int I){
		this.nRow=nRow;
		row=new SpVectComp[nRow];
		for(int i=0;i<nRow;i++)
			row[i]=new SpVectComp(I);
	}

	public SpMatComp(int nRow,int I,int L){
		this.nRow=nRow;
		row=new SpVectComp[nRow];
		for(int i=0;i<nRow;i++)
			row[i]=new SpVectComp(I,L);
	}

	public SpMatComp (Mat Mr, Mat Mm, double factor){
		
		
		nRow=Mr.nRow;
		row=new SpVectComp[nRow];
		
		int[] nz=new int[nRow];
		for(int i=0;i<nRow;i++){
			for(int j=0;j<=i;j++)
				if(Mr.el[i][j]!=0 || Mm.el[i][j]!=0) nz[i]++;
		}
		
		for(int i=0;i<nRow;i++){
			row[i]=new SpVectComp(nRow,nz[i]);

			int nz2=0;

			for(int j=0;j<=i;j++)
				if(Mr.el[i][j]!=0 || Mm.el[i][j]!=0) {
					row[i].el[nz2]=new Complex(Mr.el[i][j],factor*Mm.el[i][j]);
					row[i].index[nz2]=j;
					nz2++;
				}
		}
		
	}
	
	public SpMatComp (Mat M){
		
		
		nRow=M.nRow;
		row=new SpVectComp[nRow];
		
		for(int i=0;i<nRow;i++){
			row[i]=new SpVectComp(M.el[i]);

		}
	
	}
	
	public SpMatComp(SpMat Mr, SpMat Mm){

		nRow=Mr.nRow;
		row=new SpVectComp[nRow];
		int nCol=Mr.getnCol();
		for(int i=0;i<Mr.nRow;i++){
			if(Mr.row[i]==null){
				row[i]=new SpVectComp(nCol);
				continue;
			}


			row[i]=new SpVectComp(nCol,Mr.row[i].nzLength);
			
			for(int j=0;j<row[i].nzLength;j++){
			
					double ur=Mr.row[i].el[j];
					double um=0;
					for(int k=0;k<Mm.row[i].nzLength;k++){
						if(Mm.row[i].index[k]==Mr.row[i].index[j]){
							um=Mm.row[i].el[k];
							break;
						}
					}
		
					row[i].el[j]=new Complex(ur,um);
					row[i].index[j]=Mr.row[i].index[j];
		}
			}
		
		
	}
	public SpMatComp(SpMat Mr){

		nRow=Mr.nRow;
		row=new SpVectComp[nRow];
		int nCol=Mr.getnCol();
		for(int i=0;i<Mr.nRow;i++){
			if(Mr.row[i]==null){
				row[i]=new SpVectComp(nCol);
				continue;
			}
			
			row[i]=new SpVectComp(nCol,Mr.row[i].nzLength);
			
			for(int j=0;j<row[i].nzLength;j++){
			
					double ur=Mr.row[i].el[j];
	
					row[i].el[j]=new Complex(ur,0);
					row[i].index[j]=Mr.row[i].index[j];
		}
			}
		
		
	}

	public SpMatComp deepCopy(){
		SpMatComp M=new SpMatComp(nRow);
		for(int i=0;i<nRow;i++){
			if(this.row[i]!=null)
			M.row[i]=row[i].deepCopy();
		}

		return M;
	}
	
	public SpMatComp eye(int L){
		SpMatComp M=new SpMatComp(L,L,1);
		for(int i=0;i<L;i++){
			M.row[i].el[0]=new Complex(1,0);
			M.row[i].index[0]=i;
		}

		return M;
	}
	
	
	public void setSymHerm(int  b){

		if(b==0){
		this.hermit=false;
		this.symm=false;
		}
		else if(b==1){
			this.hermit=false;
			this.symm=true;
			}
		else if(b==2){
			this.hermit=true;
			this.symm=false;
			}
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

	public SpMatComp ichol(){
		return ichol(1.2);

	}
	
/*	public SpMatComp chol(){
	
		Mat M=this.matForm(true).chol();
		return new SpMatComp(M);

	}*/

	
	public SpMatComp ichol(double domfact){
	
		if(this.symm)
			return  icholSymComp(domfact);
		
		util.pr("Computing Hermitian Preconditioner..");

		int[] nz1=new int[nRow];
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				int j=this.row[i].index[k];
				nz1[j]++;			
			}

		SpMatComp upInd=new SpMatComp(nRow);
		for(int i=0;i<nRow;i++){
			upInd.row[i]=new SpVectComp(nRow,nz1[i]);
		}
		
		int[] nz2=new int[nRow];
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				int j=this.row[i].index[k];
				upInd.row[j].index[nz2[j]++]=i;
			}

	
		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		SpMatComp L=new SpMatComp(nRow);
		Complex[] D=new Complex[nRow];
		int[] nz3=new int[nRow];

		int kd=0;
		Complex s=new Complex();

		for(int i=0;i<nRow;i++){
			kd=this.row[i].nzLength;	
			L.row[i]=new SpVectComp(nRow,kd);
			L.row[i].index=this.row[i].index.clone();
		
		}
			
		for(int j=0;j<nRow;j++){
			kd=row[j].nzLength;

			s=L.row[j].dot(j-1,D,row[j].index);
			D[j]=row[j].el[kd-1].times(domfact).sub(s);
			//D[j].show();
			kd=nz1[j];

			for(int k=0;k<kd;k++){
				int i=upInd.row[j].index[k];	
				if(i==j){
					L.row[i].el[nz3[i]]=new Complex(1,0);
					continue;
				}
				s=L.row[i].dot(L.row[j],j-1,D,L.row[i].index,L.row[j].index);
				L.row[i].el[nz3[i]]=this.row[i].el[nz3[i]].sub(s).times(D[j].inv());
				
				nz3[i]++;

			}

		}
	
		int j=0;
		for(int i=0;i<nRow;i++){
			for(int k=0;k<row[i].nzLength ;k++){
				j=row[i].index[k];
				if(j==i){
			
					L.row[i].el[k]=D[i].pow(0.5).times(domfact);
				}
				else{
					
					L.row[i].el[k]=L.row[i].el[k].times(D[j].pow(0.5));
				}
			}
	
		}


		return L;

	}
	
	public SpMatComp icholSymComp(double domfact){

		
		util.pr("Computing Symmetric Complex Preconditioner..");

		int[] nz1=new int[nRow];
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				int j=this.row[i].index[k];
				nz1[j]++;			
			}

		SpMatComp upInd=new SpMatComp(nRow);
		for(int i=0;i<nRow;i++){
			upInd.row[i]=new SpVectComp(nRow,nz1[i]);
		}
		
		int[] nz2=new int[nRow];
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				int j=this.row[i].index[k];
				upInd.row[j].index[nz2[j]++]=i;
			}

	
		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		SpMatComp L=new SpMatComp(nRow);
		Complex[] D=new Complex[nRow];
		int[] nz3=new int[nRow];

		int kd=0;
		Complex s=new Complex();

		for(int i=0;i<nRow;i++){
			kd=this.row[i].nzLength;	
			L.row[i]=new SpVectComp(nRow,kd);
			L.row[i].index=this.row[i].index.clone();
		
		}
			
		for(int j=0;j<nRow;j++){
			kd=row[j].nzLength;

			s=L.row[j].dotSym(j-1,D,row[j].index);
			Complex tt=row[j].el[kd-1].times(domfact);

			D[j]=tt.sub(s);

			kd=nz1[j];

			for(int k=0;k<kd;k++){
				int i=upInd.row[j].index[k];	
				if(i==j){
					L.row[i].el[nz3[i]]=new Complex(1,0);
					continue;
				}
				s=L.row[i].dotSym(L.row[j],j-1,D,L.row[i].index,L.row[j].index);
				L.row[i].el[nz3[i]]=this.row[i].el[nz3[i]].sub(s).times(D[j].inv());
				nz3[i]++;

			}
			
					

		}
		
		
		
		int j=0;
		for(int i=0;i<nRow;i++){
			for(int k=0;k<row[i].nzLength ;k++){
				j=row[i].index[k];
				if(j==i){
					L.row[i].el[k]=D[i].pow(0.5).times(domfact);

				}
				else{
					
					L.row[i].el[k]=L.row[i].el[k].times(D[j].pow(0.5));
				}
			}
	
		}



		return L;

	}




	public Mat[] matForm(){
		
		Mat[] M=new Mat[2];
		
		M[0]=new Mat(nRow,getnCol());
		M[1]=new Mat(nRow,getnCol());
		
		for(int i=0;i<nRow;i++){

			for(int k=0;k<row[i].nzLength;k++){
				int j=row[i].index[k];
				M[0].el[i][j]=row[i].el[k].re;
				M[1].el[i][j]=row[i].el[k].im;
				if(i!=j){
				if(symm){
					M[0].el[j][i]=row[i].el[k].re;
					M[1].el[j][i]=row[i].el[k].im;

				}
				else if(this.hermit){
					M[0].el[j][i]=row[i].el[k].re;
					M[1].el[j][i]=-row[i].el[k].im;

				}
				}
			}
			}

		return M;
	}

	public Mat matFormLow(){
		Mat M=new Mat(nRow,getnCol());
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				M.el[i][2*row[i].index[k]]=row[i].el[k].re;
				M.el[i][2*row[i].index[k]+1]=row[i].el[k].im;
			}

		return M;
	}
	
	public Mat matFormAssym(){
		Mat M=new Mat(nRow,getnCol());
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				M.el[i][2*row[i].index[k]]=row[i].el[k].re;
				M.el[i][2*row[i].index[k]+1]=row[i].el[k].im;
			}

		return M;
	}
	
	


	public SpMatComp transpose(int cLmax){

		SpMatComp T=new SpMatComp(getnCol(),nRow,cLmax);
		int[] nz=new int[nRow];
		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].nzLength;j++){
				int cl=row[i].index[j];
				T.row[cl].el[nz[cl]]=row[i].el[j].deepCopy();
				T.row[cl].index[nz[cl]++]=i;
			}

		T.sortAndTrim(nz);
		return T;
	}

	public SpMatComp reflect(int clmax){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");
		int[] nz=new int[nRow];

		SpMatComp A=this.transpose(clmax);
		

		for(int i=0;i<nRow;i++){
			SpVectComp rx=new SpVectComp(nRow,clmax);

			for(int j=0;j<row[i].nzLength;j++){
				int cl=row[i].index[j];
				rx.el[nz[i]]=row[i].el[j];
				rx.index[nz[i]]=cl;
				nz[i]++;
				
			}
			for(int j=1;j<A.row[i].nzLength;j++){
				int cl=A.row[i].index[j];
				if(cl>i){
				rx.el[nz[i]]=A.row[i].el[j];
				rx.index[nz[i]]=cl;
				nz[i]++;
				}
			}

			A.row[i]=rx.deepCopy();
			
			}
		A.sortAndTrim(nz);

		return A;
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


	public VectComp mulLower(SpVectComp u){
		VectComp v=new VectComp(getnRow());
		for(int i=0;i<v.length;i++)
			v.el[i]=row[i].dot(u);
		return v;
	}

	public VectComp mulLower(VectComp u){
		VectComp v=new VectComp(getnRow());
		for(int i=0;i<v.length;i++)
			v.el[i]=row[i].dot(u);
		return v;
	}

	public VectComp mul(VectComp u){

		VectComp v=new VectComp(getnRow());
		VectComp w=new VectComp(getnRow());
		int j;
		for(int i=0;i<nRow;i++){
			
			if(row[i]==null) continue;
			v.el[i]=new Complex();
			w.el[i]=new Complex();
			for(int k=0;k<row[i].nzLength;k++){
				j=row[i].index[k];
				if(j>i) break;
				v.el[i]=v.el[i].add(row[i].el[k].times(u.el[j]));
				if(j!=i){
					if(this.symm)
					w.el[j]=w.el[j].add(row[i].el[k].times(u.el[i]));
					else if(this.hermit)
						w.el[j]=w.el[j].add(row[i].el[k].conj().times(u.el[i]));
				}
			}
		}
		return v.add(w);
	}
	
	public VectComp smul(VectComp u){
		

		VectComp v=new VectComp(getnRow());
		VectComp w=new VectComp(getnRow());
		int j;
		for(int i=0;i<nRow;i++){
			if(row[i]==null) continue;
			for(int k=0;k<row[i].nzLength;k++){
				j=row[i].index[k];
				if(j>i) break;
				v.el[i]=v.el[i].add(row[i].el[k].times(u.el[j]));
				if(j!=i)
					w.el[j]=w.el[j].add(row[i].el[k].times(u.el[i]));
			}
		}
		return v.add(w);
	}



	public void addSmaller(SpMatComp M){

		for(int i=0;i<nRow;i++)
			row[i].addSmaller(M.row[i]);
	}
	
	public SpMatComp addNew(SpMatComp M){
		SpMatComp B=this.deepCopy();
		for(int i=0;i<nRow;i++)
			if(M.row[i]!=null && B.row[i]!=null)
			B.row[i].addSmaller(M.row[i]);
		
		return B;
	}
	
	public void add(VectComp D){
		if(nRow!=D.length) throw new IllegalArgumentException("Dimensions do not agree.");
		for(int i=0;i<nRow;i++)
			row[i].el[row[i].nzLength-1]=row[i].el[row[i].nzLength-1].add(D.el[i]);
	}


	public void augh(SpMatComp N){
		for(int i=0;i<nRow;i++)
			row[i]=row[i].augh(N.row[i]);
	}

	public SpMatComp augv(SpMatComp N){
		SpMatComp MN=new SpMatComp(nRow+N.getnRow());
		for(int i=0;i<nRow;i++)
			MN.row[i]=row[i];
		for(int i=0;i<N.getnRow();i++)
			MN.row[i+nRow]=N.row[i];
		return MN;
	}

	public VectComp diag(){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		VectComp D=new VectComp(nRow);
		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].getNzLength();j++)
				if(row[i].index[j]==i){
					D.el[i]=row[i].el[j];
					break;
				}
		return D;
	}


	public VectComp diagSym(){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		VectComp D=new VectComp(nRow);
		for(int i=0;i<nRow;i++)		{
			D.el[i]=row[i].el[row[i].nzLength-1];
		}

		return D;
	}
	
	public void addToDiag(VectComp D){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");
		if(nRow!=D.length) throw new IllegalArgumentException("Arrays dimensions do not square.");

		for(int i=0;i<nRow;i++)		
			this.row[i].el[row[i].nzLength-1]=this.row[i].el[row[i].nzLength-1].add(D.el[i]);

	}
	
	public SpMatComp offDiagSym(){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");
		SpMatComp Ls=new SpMatComp(nRow);

		for(int i=0;i<nRow;i++)	{
			int nnz=this.row[i].nzLength-1;
			Ls.row[i]=new SpVectComp(this.nRow,nnz);
			Ls.row[i].el=Arrays.copyOf(this.row[i].el,nnz);
			Ls.row[i].index=Arrays.copyOf(this.row[i].index,nnz);
		}
			

		return Ls;
	}



	public void times(double a){
		for(int i=0;i<nRow;i++)
			row[i]=row[i].times(a);
	}
	
	public void times(Complex a){
		for(int i=0;i<nRow;i++)
			row[i]=row[i].times(a);
	}

	public SpMatComp timesNew(double a){
		SpMatComp B=this.deepCopy();
		B.times(a);
		return B;
	}
	
	public SpMatComp timesNew(Complex a){
		SpMatComp B=this.deepCopy();
		B.times(a);
		return B;
	}


	public void show(){
		Mat M[]=matForm();


		String s="";
		
		for (int i=0;i<M[0].nRow;i++){
			for (int j=0;j<M[0].nCol;j++){
				if(M[1].el[i][j]<0) s="-";
				else
					s="+";
				System.out.print(M[0].el[i][j]+s+"j"+abs(M[1].el[i][j])+"\t");
			}
			System.out.println();
		}
	}
	
	public void plot(){


		Mat[] M=matForm();

		M[0].plot();
		
	}
	

	public void showLower(){
		Mat M=matFormLow();
		M.show();
	}

	public void shownz(){

		for(int i=0;i<nRow;i++)
		for(int j=0;j<row[i].nzLength;j++)
			System.out.format("(% d % d)   %25.12e  %25.12e\n",i,row[i].index[j], row[i].el[j].re,  row[i].el[j].im);
	
	}

	public void showcl(){

		for(int i=0;i<nRow;i++)
			row[i].showr();
	}


	public Vect scale(VectComp b){

		int I=getnRow();
		int J=getnCol();
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		Vect invDS=this.diagSym().abs().sqrt().inv();
		
		DMD(invDS);

		b.timesVoid(invDS);

		return invDS;
	}
	public Vect scale(){

		int I=getnRow();
		int J=getnCol();
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		Vect invDS=this.diagSym().abs().sqrt().inv();
		DMD(invDS);


		return invDS;
	}





	public void DMD(Vect D){

		int I=getnRow();

		for(int i=0;i<I;i++)
			for(int j=0;j<row[i].nzLength;j++){
				row[i].el[j]=row[i].el[j].times(D.el[row[i].index[j]]);
			}

		for(int i=0;i<I;i++)
			for(int j=0;j<row[i].nzLength;j++)
				row[i].el[j]=row[i].el[j].times(D.el[i]);

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
