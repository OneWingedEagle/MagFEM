package math;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import math.SpVect;
import static java.lang.Math.*;
import math.Vect;

public class SpMat  {

	public SpVect[] row;
	public int nRow;
	public boolean lower=true;
	public SpMatSolver solver=new SpMatSolver();
	public SpMat(){}

	public SpMat(int nRow){
		this.nRow=nRow;
		row=new SpVect[nRow];
	}

	public SpMat(int nRow,int I){
		this.nRow=nRow;
		row=new SpVect[nRow];
		for(int i=0;i<nRow;i++)
			row[i]=new SpVect(I);
	}

	public SpMat(int nRow,int I,int L){
		this.nRow=nRow;
		row=new SpVect[nRow];
		for(int i=0;i<nRow;i++)
			row[i]=new SpVect(I,L);
	}

	public SpMat (Mat M,boolean sym){
		
		
		nRow=M.nRow;
		row=new SpVect[nRow];
		
		int[] nz=new int[nRow];
		for(int i=0;i<nRow;i++){
			for(int j=0;j<=i;j++)
				if(M.el[i][j]!=0) nz[i]++;
		}
		for(int i=0;i<nRow;i++){
			row[i]=new SpVect(nRow,nz[i]);

			int nz2=0;
			for(int j=0;j<=i;j++)
				if(M.el[i][j]!=0) {
					row[i].el[nz2]=	M.el[i][j];
					row[i].index[nz2]=i;
					nz2++;
				}
		}
		
	}
	
	public SpMat (Mat M){
		
		
		nRow=M.nRow;
		row=new SpVect[nRow];
		
		for(int i=0;i<nRow;i++){
			row[i]=new SpVect(M.el[i]);

		}
	
	}

	public SpMat deepCopy(){
		SpMat M=new SpMat(nRow);
		for(int i=0;i<nRow;i++){
			if(this.row[i]!=null)
			M.row[i]=row[i].deepCopy();
		}

		return M;
	}
	
	public SpMat eye(int L){
		SpMat M=new SpMat(L,L,1);
		for(int i=0;i<L;i++){
			M.row[i].el[0]=1;
			M.row[i].index[0]=i;
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

	public SpMat ichol(){
		return ichol(1.2);

	}
	
	public SpMat chol(){
	
		Mat M=this.matForm(true).chol();
		return new SpMat(M);

	}

	public SpMat icholw(double domfact){
		util.pr("Weak Precond.");
		
		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		SpMat T=new SpMat(nRow);
		double[] D=new double[nRow];
		double s;

		int j=0,kd=0;
		s=0;

		for(int i=0;i<nRow;i++){
			kd=row[i].nzLength;
			T.row[i]=new SpVect(nRow,kd,0);

			j=0;

			for(int k=0;k<kd;k++){
				j=row[i].index[k];
				s=T.row[j].dot(j-1,D,row[j].index);
				D[j]=abs(domfact*row[i].el[kd-1]-s);

				s=T.row[i].dot(T.row[j],j-1,D,row[i].index,row[j].index);
				T.row[i].el[k]=(row[i].el[k]-s)/D[j];
			}

		}

		for(int i=0;i<nRow;i++){
			j=0;
			for(int k=0;(k<row[i].nzLength && j<=i);k++){
				j=row[i].index[k];
				if(i==j){
					T.row[i].el[k]=sqrt(D[i])/domfact;
				}
				else{
					T.row[i].el[k]=T.row[i].el[k]*sqrt(D[i]);
				}
			}
		}


		return T;

	}
	
	public SpMat ichol(double domfact){
		util.pr("Computing Preconditioner..");

		int[] nz1=new int[nRow];
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				int j=this.row[i].index[k];
				nz1[j]++;			
			}

		SpMat upInd=new SpMat(nRow);
		for(int i=0;i<nRow;i++){
			upInd.row[i]=new SpVect(nRow,nz1[i],1);
		}
		
		int[] nz2=new int[nRow];
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				int j=this.row[i].index[k];
				upInd.row[j].index[nz2[j]++]=i;
			}

	
		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		SpMat L=new SpMat(nRow);
		double[] D=new double[nRow];
		int[] nz3=new int[nRow];
		
		double s;

		int kd=0;
		s=0;

		for(int i=0;i<nRow;i++){
			kd=this.row[i].nzLength;	
			L.row[i]=new SpVect(nRow,kd,0);
			L.row[i].index=this.row[i].index.clone();
		
		}
			
		for(int j=0;j<nRow;j++){
			kd=row[j].nzLength;

			s=L.row[j].dot(j-1,D,row[j].index);
			D[j]=abs(domfact*row[j].el[kd-1]-s);

			kd=nz1[j];

			for(int k=0;k<kd;k++){
				int i=upInd.row[j].index[k];	
				//=========
				
				if(i==j){
					L.row[i].el[nz3[i]]=1;
					continue;
				}
				
				//=======
				s=L.row[i].dot(L.row[j],j-1,D,L.row[i].index,L.row[j].index);
				L.row[i].el[nz3[i]]=(this.row[i].el[nz3[i]]-s)/D[j];
				nz3[i]++;

			}
			
					

		}
		
		int j=0;
		for(int i=0;i<nRow;i++){
			for(int k=0;k<row[i].nzLength ;k++){
				j=row[i].index[k];
				if(j==i){
					L.row[i].el[k]=sqrt(D[i])/domfact;
				}
				else{
					L.row[i].el[k]*=sqrt(D[j]);
				}
			}
	
		}

		return L;

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


	public int numbEigLessThan(double mu){
		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		SpMat A=this.deepCopy();
		A.lamShift(mu);
		A.scale();
		Vect D=new Vect(A.nRow);

		A.iLDU(D,1);
		
		int nneg=0;
			for(int i=0;i<A.nRow;i++)
				if(D.el[i]<0) nneg++;
	
		return nneg;
	}
	
	public int numbEigLessThan(SpMat Ms,double mu){
		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		SpMat A=this.addGeneral(Ms.timesNew(-mu));

		//A.scale();
		A.times(1e-9);

		Vect D=new Vect(A.nRow);

		A.iLDU(D,1);

		int nneg=0;
			for(int i=0;i<A.nRow;i++)
				if(D.el[i]<0) nneg++;
	
		return nneg;
	}
	
	
	public Mat matForm(){
		
		return matForm(false);
		
	}

	public Mat matForm(boolean sym){
		Mat M=new Mat(nRow,getnCol());
		
		for(int i=0;i<nRow;i++){

			for(int k=0;k<row[i].nzLength;k++){
				M.el[i][row[i].index[k]]=row[i].el[k];

				if(sym)
				M.el[row[i].index[k]][i]=row[i].el[k];
			}
			}

		return M;
	}

	public Mat matFormLow(){
		Mat M=new Mat(nRow,getnCol());
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				M.el[i][row[i].index[k]]=row[i].el[k];
			}

		return M;
	}
	
	public Mat matFormAssym(){
		Mat M=new Mat(nRow,getnCol());
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				M.el[i][row[i].index[k]]=row[i].el[k];
			}

		return M;
	}

	public SpMat band(){
//this.show();
		HashSet R = new HashSet();
		List<Integer> Q=new ArrayList<Integer>();

	int ix=0;
	
	SpMat band=new SpMat(this.nRow);

	int N=this.nRow;
	boolean[] counted=new boolean[N];

/*	int df0=N;
	int m=0;
	for(int i=0;i<this.nRow;i++)
		if(this.row[i].nzLength<df0)
			m=i;


		int n0=m;
//	util.pr(n0);
		R.add(n0);

		int[] arr=this.row[n0].index;
		util.hshow(arr);
	Vect df=new Vect(arr.length);
		for(int i=0;i<arr.length;i++)
			df.el[i]=this.row[arr[i]].nzLength;
		int[] ind=df.bubble();
		
		for(int i=0;i<arr.length;i++)
			ind[i]=i;
		
		for(int i=0;i<arr.length;i++)
	//		Q.add(arr[i]);
		Q.add(arr[ind[i]]);*/
		


		while(R.size()<N){
			
			if(Q.size()==0){
				int df0=N;
				int m=0;
				for(int i=0;i<this.nRow;i++){
					if(!counted[i]){
					if(this.row[i].nzLength<df0)
						m=i;
					}
				}
				m=0;


					int n0=m;
					R.add(n0);

					int[] arr=this.row[n0].index;
					Vect df=new Vect(arr.length);
					for(int i=0;i<arr.length;i++)
						df.el[i]=this.row[arr[i]].nzLength;
					int[] ind=df.bubble();
					
		/*			for(int i=0;i<arr.length;i++)
						ind[i]=i;*/
					
					for(int i=0;i<arr.length;i++){
						Integer I=new Integer(arr[ind[i]]);
						if(!counted[arr[ind[i]]] && !Q.contains(I)){
							Q.add(arr[ind[i]]);
							//counted[arr[ind[i]]]=true;
						}
					}
			}

			//if(Q.size()==0) break;
/*			for(int i=0;i<Q.size();i++){
				int nx=Q.get(i);
				System.out.print(" "+nx);
			}*/
	
		util.pr("");
		Integer n=Q.get(0);
		Q.remove(n);

	
			
			if(!counted[n]){
			R.add(n);
			counted[n]=true;
			int[] arr=this.row[n].index;
			Vect df=new Vect(arr.length);
			for(int i=0;i<arr.length;i++)
				df.el[i]=this.row[arr[i]].nzLength;
			int[] ind=df.bubble();
			
			for(int i=0;i<arr.length;i++){
				Integer I=new Integer(arr[ind[i]]);
				if(!counted[arr[ind[i]]]/* && !Q.contains(I)*/)
					Q.add(arr[ind[i]]);
			}

			}
			

	
			//util.pr(R.size()+" sR "+nnz);
			
			}
		
	
	
		ArrayList<Integer> rr = new ArrayList<Integer>(R);
		for(int i=0;i<R.size();i++){
			int nx=rr.get(N-1-i);
			System.out.println(" "+nx);
		}
		
		
		for(int i=0;i<1*this.nRow;i++)
		{
			int rx=rr.get(i);
			band.row[rx]=this.row[i].deepCopy();
			for(int j=0;j<band.row[rx].nzLength;j++)
				band.row[rx].index[j]=rr.get(this.row[i].index[j]);
		//	band.row[rx].el=this.row[i].el;
		}
	return band;

	}
	
	
	public int[] RCMorder(){
		return RCMorder(this);
	}
	
	public int[] RCMorder(SpMat K){
		
		List<Integer> R = new ArrayList<Integer>();
		List<Integer> Q=new ArrayList<Integer>();

	
	
	int N=this.nRow;
	boolean[] counted=new boolean[N];


		while(R.size()<K.nRow){
		
			if(Q.size()==0){
				int df0=N;
				int m=0;
				for(int i=0;i<K.nRow;i++){
					if(!counted[i]){
					if(K.row[i].nzLength<df0)
						m=i;
					}
				}
					int n0=m;

					int[] arr=K.row[n0].index;
					Vect df=new Vect(arr.length);
					for(int i=0;i<arr.length;i++)
						df.el[i]=K.row[arr[i]].nzLength;
					int[] ind=df.bubble();
					

					
					for(int i=0;i<arr.length;i++){
						Integer I=new Integer(arr[ind[i]]);
						if(!counted[arr[ind[i]]] && !Q.contains(I)){
							Q.add(arr[ind[i]]);
				
						}
					}
			}


		Integer n=Q.get(0);
		Q.remove(n);


	if(!counted[n]){
			R.add(n);
			counted[n]=true;
			int[] arr=K.row[n].index;
			Vect df=new Vect(arr.length);
			for(int i=0;i<arr.length;i++)
				df.el[i]=K.row[arr[i]].nzLength;
			int[] ind=df.bubble();

			
			for(int i=0;i<arr.length;i++){
				Integer I=new Integer(arr[ind[i]]);
				if(!counted[arr[ind[i]]] && !Q.contains(I))
					Q.add(arr[ind[i]]);
			}

			}

			
			}
		
	
	
		int[] rr=new int[R.size()];
		for(int i=0;i<R.size();i++){
			rr[i]=R.get(N-1-i);
			
			
		
		}


	return rr;			

	}
	
	
	
	public SpMat RCM(){
		

		SpMat K=new SpMat();
		
		if(this.lower)
			K=this.reflect(200);
		else K=this;

		int[] rr=this.RCMorder(K);
		
		int[] map=new int[rr.length];
		for(int i=0;i<rr.length;i++){
			map[rr[i]]=i;
			
		
		}
	
		SpMat Krcm=new SpMat(this.nRow);
		Krcm.lower=this.lower;

		for(int i=0;i<this.nRow;i++){
			int nn=rr[i];
			Krcm.row[i]=this.row[nn].deepCopy();
			for(int j=0;j<Krcm.row[i].nzLength;j++)
				Krcm.row[i].index[j]=map[this.row[nn].index[j]];
	
			
		}

		return Krcm;
		
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

	public SpMat reflect(int clmax){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");
		int[] nz=new int[nRow];

		SpMat A=this.transpose(clmax);
		

		for(int i=0;i<nRow;i++){
			SpVect rx=new SpVect(nRow,clmax);

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


	public Vect mulLower(SpVect u){
		Vect v=new Vect(getnRow());
		for(int i=0;i<v.length;i++)
			v.el[i]=row[i].dot(u);
		return v;
	}

	public Vect mulLower(Vect u){
		Vect v=new Vect(getnRow());
		for(int i=0;i<v.length;i++)
			v.el[i]=row[i].dot(u);
		return v;
	}

	public Vect smul(Vect u){
		

		Vect v=new Vect(getnRow());
		Vect w=new Vect(getnRow());
		int j;
		for(int i=0;i<nRow;i++){
			if(row[i]==null) continue;
			for(int k=0;k<row[i].nzLength;k++){
				j=row[i].index[k];
				if(j>i) break;
				v.el[i]+=row[i].el[k]*u.el[j];
				if(j!=i)
					w.el[j]+=row[i].el[k]*u.el[i];
			}
		}
		return v.add(w);
	}
	
	public Vect amul(Vect u){

		Vect v=new Vect(getnRow());
		int j;
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				j=row[i].index[k];
				v.el[i]+=row[i].el[k]*u.el[j];
				
			}
		return v;
	}
	
	
	public Vect mulUpper(Vect u){
		Vect w=new Vect(getnRow());
		int j;
		for(int i=0;i<nRow;i++)
			for(int k=0;k<row[i].nzLength;k++){
				j=row[i].index[k];
				if(j!=i)
					w.el[j]+=row[i].el[k]*u.el[i];
			}
		return w;
	}
	

	public Mat smul( Mat X){
		return smul(X,new boolean[X.nCol]);
	}

	public Mat smul(Mat X,boolean[] lock){
		Mat Y=new Mat(X.size());
		Vect v=new Vect(X.nRow);
		for(int i=0;i<X.nCol;i++){
			if(lock[i]) continue;
			for(int j=0;j<X.nRow;j++)
				v.el[j]=X.el[j][i];
			v=smul(v);
			for(int j=0;j<X.nRow;j++)
				Y.el[j][i]=v.el[j];
		}

		return Y;
	}


	public Vect mulLow(Vect u){
		Vect v=new Vect(getnRow());
		for(int i=0;i<v.length;i++){
			v.el[i]=row[i].dot(u,i);
		}
		return v;
	}

	public Vect mulUp(Vect u){
		Vect v=new Vect(getnRow());
		for(int i=0;i<v.length;i++){
			v.el[i]=row[i].dotUp(u,i);
		}
		return v;
	}

	public void addSmaller(SpMat M){

		for(int i=0;i<nRow;i++)
			row[i].addSmaller(M.row[i]);
	}
	
	public SpMat addGeneral(SpMat M){
		SpMat B=this.deepCopy();
		for(int i=0;i<nRow;i++)
			if(M.row[i]!=null && B.row[i]!=null)
			B.row[i].addGeneral(M.row[i]);
		
		return B;
	}
	
	public SpMat addSmallerNew(SpMat M){
		SpMat B=this.deepCopy();
		for(int i=0;i<nRow;i++)
			if(M.row[i]!=null && B.row[i]!=null)
			B.row[i].addSmaller(M.row[i]);
		
		return B;
	}
	
	public void add(Vect D){
		if(nRow!=D.length) throw new IllegalArgumentException("Dimensions do not agree.");
		for(int i=0;i<nRow;i++)
			row[i].el[row[i].nzLength-1]+=D.el[i];
	}

	public void lamShift(double mu){
		for(int i=0;i<nRow;i++)
			row[i].el[row[i].nzLength-1]-=mu;
	}

	public void augh(SpMat N){
		for(int i=0;i<nRow;i++)
			row[i]=row[i].augh(N.row[i]);
	}

	public SpMat augv(SpMat N){
		SpMat MN=new SpMat(nRow+N.getnRow());
		for(int i=0;i<nRow;i++)
			MN.row[i]=row[i];
		for(int i=0;i<N.getnRow();i++)
			MN.row[i+nRow]=N.row[i];
		return MN;
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

	public Vect diagSym(){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");

		Vect D=new Vect(nRow);
		for(int i=0;i<nRow;i++)		{
			D.el[i]=row[i].el[row[i].nzLength-1];
		}

		return D;
	}
	
	public void addToDiag(Vect D){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");
		if(nRow!=D.length) throw new IllegalArgumentException("Arrays dimensions do not square.");

		for(int i=0;i<nRow;i++)		
			this.row[i].el[row[i].nzLength-1]+=D.el[i];

	}
	
	public SpMat offDiagSym(){

		if(nRow!=getnCol()) throw new IllegalArgumentException("Matrix is not square.");
		SpMat Ls=new SpMat(nRow);

		for(int i=0;i<nRow;i++)	{
			int nnz=this.row[i].nzLength-1;
			Ls.row[i]=new SpVect(this.nRow,nnz);
			Ls.row[i].el=Arrays.copyOf(this.row[i].el,nnz);
			Ls.row[i].index=Arrays.copyOf(this.row[i].index,nnz);
		}
			

		return Ls;
	}

	public Vect scale(Vect b, Vect x_init){
		int I=getnRow();
		int J=getnCol();
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");

		Vect vec=this.diagSym().abs().sqrt();
		Vect invDS=new Vect(I);
		for(int i=0;i<I;i++)
			if(vec.el[i]!=0) invDS.el[i]=1./vec.el[i];
			else  invDS.el[i]=1.;		
		DMD(invDS);

		b.timesVoid(invDS);
		x_init.timesVoid(invDS);

		return invDS;
	}

	public Vect scale(Vect b){

		int I=getnRow();
		int J=getnCol();
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");

		Vect vec=this.diagSym().abs().sqrt();
		Vect invDS=new Vect(I);
		for(int i=0;i<I;i++)
			if(vec.el[i]!=0) invDS.el[i]=1./vec.el[i];
			else  invDS.el[i]=1.;		
		DMD(invDS);

		b.timesVoid(invDS);

		return invDS;
	}
	
	public Vect scale(){

		int I=getnRow();
		int J=getnCol();
		if(I!=J) throw new IllegalArgumentException("Matrix is not square");
		
		Vect vec=this.diagSym().abs().sqrt();
		Vect invDS=new Vect(I);
		for(int i=0;i<I;i++)
			if(vec.el[i]!=0) invDS.el[i]=1./vec.el[i];
			else  invDS.el[i]=1.;		

		DMD(invDS);


		return invDS;
	}





	public void DMD(Vect D){

		int I=getnRow();

		for(int i=0;i<I;i++)
			for(int j=0;j<row[i].nzLength;j++){
				row[i].el[j]*=D.el[row[i].index[j]];
			}

		for(int i=0;i<I;i++)
			for(int j=0;j<row[i].nzLength;j++)
				row[i].el[j]*=D.el[i];

	}


	public void times(double a){
		for(int i=0;i<nRow;i++)
			if(row[i]!=null)
			row[i]=row[i].times(a);
	}

	public SpMat timesNew(double a){
		SpMat B=this.deepCopy();
		B.times(a);
		return B;
	}


	public void show(String format){
		Mat M=matForm();
		M.show(format);
	}
	public void show(){
		Mat M=matForm();
		M.show();
	}

	
	public void plot(){
		int N=this.nRow;
		//if(N>10000) N=6000;
		if(this.nRow<=N)
			util.plot(this);
		else{
			util.pr("Matrix too large to plot. It is plotted partially.");
			
			SpMat Ms=new SpMat(N,N);
			for(int i=0;i<Ms.nRow;i++){
				Ms.row[i]=this.row[i].deepCopy();
				
				}
			
			util.plot(Ms);
		}

		
	}
	

	public void showLower(){
		Mat M=matFormLow();
		M.show();
	}

	public void shownz(){

		for(int i=0;i<nRow;i++)
			for(int j=0;j<row[i].nzLength;j++)
				System.out.format(" ( %d  %d )  %25.12e\n",i,row[i].index[j], row[i].el[j]);
	}

	public int numNonzeros(){
		int nz=0;
		for(int i=0;i<nRow;i++)
			if(row[i]!=null) nz+=row[i].nzLength;
	
	return nz;
	}


	public void shownzA(){

		for(int i=0;i<nRow;i++)
				for(int j=0;j<row[i].nzLength;j++)
					if(row[i].el[j]!=0)
				System.out.format(" ( %d  %d )  %25.12e\n",i,row[i].index[j], row[i].el[j]);
	}
	
	public void showcl(){

		for(int i=0;i<nRow;i++)
			row[i].showr();
	}


	public SpMat QR(SpMat R) {
		SpMat Q=deepCopy();
		double c;
		Q.transpose(100);
		Vect v=new Vect(),e=new Vect(nRow),vr=new Vect(nRow);
		v=Q.row[0].vectForm();
		c=v.norm();
		v=v.times(1/c);
		//R.el[0][0]=c;
		Q.row[0]=new SpVect(v);


		for(int i=1;i<nRow;i++){
			v=Q.row[i].vectForm();
			vr=v;
			for(int j=0;j<i;j++){
				e=Q.row[j].vectForm();
				c=v.dot(e);
				//R.el[j][i]=c;
				v=v.sub(e.times(c));
			}

			v.normalize();
			//R.el[i][i]=vr.dot(v);;	
			Q.row[i]=new SpVect(v);

		}

		Q.transpose(100);

		return Q;

	}


	public Vect eigPow(double errc) {

		Vect v=new Vect(nRow);
		v.rand();
		v.normalize();
		Vect vp;
		double err=1;
		int k=0;
		while(err>errc){
			k++;
			if(k>1000) {  break;}				
			vp=v;
			v=smul(v);
			v.normalize();
			err=v.sub(vp).max();

		}
		System.out.println(" iterantions: "+k+" error: "+  err);

		return v;

	}

	public double eigMax(double errc) {



		Vect v=new Vect(nRow);
		v.rand();
		v.normalize();
		double lam=1,lamr;
		double err=1;
		int k=0;
		while(err>errc){
			k++;
			if(k>1000) {  break;}				
			v=smul(v);
			v.normalize();
			
			lamr=lam;
			lam=this.smul(v).norm();
			if(lam!=0)
			err=abs(lam-lamr)/abs(lam);
			else if(lam!=0)
			err=abs(lam-lamr)/abs(lamr);
			else
				err=abs(lam-lamr)/1e-20;

		}
		System.out.println(" iterantions: "+k+" error: "+  err);

		return smul(v).norm();
	}

	public double eigMin(double errc, SpMatSolver solver) {


		Vect v=new Vect(nRow);
		Vect xini=new Vect(nRow);
		v.rand();
		v.normalize();
		Vect b;
		double lam=1,lamr;
		double err=1;
		int k=0;
		SpMat A=this.deepCopy();
		Vect Ci=A.scale(xini);
		SpMat L=A.ichol();
		//SpMat L=A.iLDU(new Vect(nRow),1.1);
		xini=new Vect(nRow);
		while(	!solver.terminate &&  err>errc){
			k++;
			if(k>1000) {  break;}			
			b=v.times(Ci);
			v=solver.silentICCG(A,L, b,1e-6,1000,xini);
			v.timesVoid(Ci);
			v.normalize();
			lamr=lam;
			lam=this.smul(v).norm();
			if(lam!=0)
			err=abs(lam-lamr)/abs(lam);
			else if(lam!=0)
			err=abs(lam-lamr)/abs(lamr);
			else
				err=abs(lam-lamr)/1e-20;
				
			System.out.println(" iterantions: "+k+" error: "+  err);
		}

		System.out.println(" iterantions: "+k+" error: "+  err);
		return smul(v).norm();

	}

	public double eigMin(SpMat Ms, double errc, SpMatSolver solver) {

		Vect v=new Vect(nRow);
		Vect xini=new Vect(nRow);
		v.rand();
		v.normalize();
		Vect vp,b;

		double err=1;
		int k=0;
		SpMat A=this.deepCopy();
		Vect Ci=A.scale(xini);
		SpMat L=A.ichol();
		xini=new Vect(nRow);
		while(	!solver.terminate &&  err>errc){
			k++;
			if(k>1000) {  break;}			
			vp=v;
			b=Ms.smul(v);		

			b.timesVoid(Ci);
			v=solver.silentICCG(A,L, b,1e-6,1000,xini);
		
			v.timesVoid(Ci);
			v.normalize();
			err=v.sub(vp).abs().max();
			System.out.println(" iterantions: "+k+" error: "+  err);
		}

		System.out.println(" iterantions: "+k+" error: "+  err);
		return smul(v).norm()/Ms.smul(v).norm();

	}

	public Mat eigVectLanc(int m){
		Mat Q=new Mat(nRow,m);
		Mat V=new Mat(nRow,m);
		Mat T=tridiagLanc(m,V);
		Q=V.mul(T.eigVect());

		return Q;
	}


	public Vect eigLanc(int m){
		Mat T=tridiagLanc(m);
		Vect lam= T.eigVal();
		lam.quickSort();
		return lam;
	}
	
	public Vect eigLancLower(int m){
		Mat T=tridiagLancLower(m);
		Vect lam1= T.eigVal();
		lam1=lam1.inv();
		lam1.quickSort();
		int nx=1;
		double[] temp=new double[m];
		temp[0]=lam1.el[0];
		for(int i=1;i<m;i++)
			if(abs(abs(lam1.el[i])/abs(lam1.el[i-1])-1)>1e-2){			
				temp[nx]=lam1.el[i];
				nx++;
			}
		Vect lam=new Vect(nx);	
		for(int i=0;i<lam.length;i++)
			lam.el[i]=temp[i];
		
		return lam;
	}


	public Mat tridiagLanc(int m,Mat V) {


		Vect v=new Vect(nRow);
		Vect r=new Vect(nRow);
		Vect vp=new Vect(nRow);
		double[] a=new double[m];
		double[] b=new double[m+1];
		r.rand();
		b[0]=r.norm();
		for(int i=0;i<m;i++){
			vp=v;
			v=r.times(1.0/b[i]);
			r=smul(v).sub(vp.times(b[i]));
			a[i]=r.dot(v);
			r=r.sub(v.times(a[i]));
			b[i+1]=r.norm();

			for(int j=0;j<nRow;j++)
				V.el[j][i]=v.el[j];

		}

		Mat P=new Mat(m,m);
		P.el[0][0]=a[0];
		P.el[0][1]=b[1];
		P.el[m-1][m-2]=b[m-1];
		P.el[m-1][m-1]=a[m-1];
		for(int i=1;i<m-1;i++){
			P.el[i][i-1]=b[i];
			P.el[i][i]=a[i];
			P.el[i][i+1]=b[i+1];
		}

		return P;
	}

	public Mat tridiagLanc(int m) {


		Vect v=new Vect(nRow);
		Vect r=new Vect(nRow);
		Vect vp=new Vect(nRow);
		double[] a=new double[m];
		double[] b=new double[m+1];
		r.rand();
		b[0]=r.norm();
		for(int i=0;i<m;i++){
			vp=v;
			v=r.times(1.0/b[i]);
			r=smul(v).sub(vp.times(b[i]));
			a[i]=r.dot(v);
			r=r.sub(v.times(a[i]));
			b[i+1]=r.norm();


		}

		Mat P=new Mat(m,m);
		P.el[0][0]=a[0];
		P.el[0][1]=b[1];
		P.el[m-1][m-2]=b[m-1];
		P.el[m-1][m-1]=a[m-1];
		for(int i=1;i<m-1;i++){
			P.el[i][i-1]=b[i];
			P.el[i][i]=a[i];
			P.el[i][i+1]=b[i+1];
		}

		return P;
	}

	public Mat tridiagLancLower(int m) {

		SpMat Ks=this.deepCopy();
		Vect Ci=Ks.scale();
		SpMat L=Ks.ichol();

		Vect v=new Vect(nRow);
		Vect r=new Vect(nRow);
		Vect vp=new Vect(nRow);
		double[] a=new double[m];
		double[] b=new double[m+1];
		r.rand();
		b[0]=r.norm();
		for(int i=0;i<m;i++){
			vp=v.deepCopy();
			v=r.times(1.0/b[i]);
			r=v.deepCopy();
			r.times(Ci);
			r=solver.silentICCG(this, L, r, 1e-8, 1000, new Vect(nRow));
			r.times(Ci);
		
			r=r.sub(vp.times(b[i]));
			a[i]=r.dot(v);
			r=r.sub(v.times(a[i]));
			b[i+1]=r.norm();
			
			if(i%10==0)
				System.out.println("Lanczos method: " +i+" out of "+m);


		}

		Mat P=new Mat(m,m);
		P.el[0][0]=a[0];
		P.el[0][1]=b[1];
		P.el[m-1][m-2]=b[m-1];
		P.el[m-1][m-1]=a[m-1];
		for(int i=1;i<m-1;i++){
			P.el[i][i-1]=b[i];
			P.el[i][i]=a[i];
			P.el[i][i+1]=b[i+1];
		}

		return P;
	}

	public Vect eigSubspace(int p,double errMax){
		
		return eigSubspace(p,errMax,solver);
		
	}

	public Vect eigSubspace(int p,double errMax, SpMatSolver solver){


		int n=this.nRow;
		int q=min(min(p+8,2*p),n);

		Mat X=new Mat(n,q);
		

		Mat X_hat=new Mat(n,q);
		Vect x=new Vect(n);
		Vect errv=new Vect(q);
		Vect lamqr=new Vect(q);
		Vect lamq=new Vect(q);
		lamq=lamq.ones(q);
		
		int[] index=new int[q];
		for(int i=0;i<q;i++)
			index[i]=i;
		
		X.setCol(new Vect().ones(n), 0);
		for(int i=1;i<q-1;i++)
			X.el[i][i]=1;
		x.rand();
		X.setCol(x, q-1);
		
		boolean[] locked=new boolean[q];
	
		Mat Q=new Mat(q,q);
		Mat I=new Mat(q,q);
		I.eye();
		Mat K_hat=new Mat();
		Mat M_hat=new Mat();
		int iter;
		int nconv=0;
		int[] fallen=new int[q];
		SpMat Ks=this.deepCopy();
		Vect Ci=Ks.scale();
		SpMat L=Ks.ichol();
		//SpMat L=Ks.iLDU(new Vect(nRow),1.2);
		
		for(iter=0;iter<200;iter++){
			for(int j=0;j<q;j++){				
				x=X.getColVect(index[j]);	
				if(!locked[j])
				{
				x.timesVoid(Ci);
				x=solver.silentICCG(Ks,L,x,1e-6,1000,new Vect(nRow));	

				x.timesVoid(Ci);
				}
				else{
					if(lamq.el[j]!=0)
					x=x.times(1/lamq.el[j]);
				}

				if( solver.terminate) break;	
			
				X_hat.setCol(x,index[j]);		
			}
			K_hat=X_hat.transp().mul(X);
			M_hat=X_hat.transp().mul(X_hat);	
			K_hat.eigGenSym(M_hat,Q);
			X=X_hat.mul(Q);

				lamqr=lamq.deepCopy();
				lamq=this.smul(X).normCol();
				index=lamq.bubble();

				errv=lamq.sub(lamqr).abs().div0(lamq);


				for( int i=0;i<p;i++)
						if(!locked[i] && errv.el[i]<errMax){
						fallen[i]++;
						if(fallen[i]>0){
							locked[i]=true;
							nconv++;
						}

						}
					
				//errv.hshow();
				System.out.println("Subspace Iteration: " +iter+"   | error: "+errv.el[p-1]+   "  |  Numver of converged eigenvalues: "+nconv);
				
			if( nconv==p || solver.terminate) {  break;}
			
		}
	

		Vect lam=new Vect(p);		
		for(int i=0;i<p;i++)
			lam.el[i]=lamq.el[i];

	
		
		lam.quickSort();
		
		return lam;	

	}
	
	
	
	public Vect eigSubspace(SpMat Ms,int p,double errMax,SpMatSolver solver){
		return eigSubspace(Ms,new Mat(Ms.nRow,p), p, errMax,solver);
	}

	public Vect eigSubspace(SpMat Ms,Mat P,int p,double errMax, SpMatSolver solver){

		int n=this.nRow;
		int q=min(min(p+8,2*p),n);

		Mat X=new Mat(n,q);
		

		Mat X_hat=new Mat(n,q);
		Vect x=new Vect(n);
		Vect errv=new Vect(q);
		Vect lamqr=new Vect(q);
		Vect lamq=new Vect(q);
		lamq=lamq.ones(q);
		
		int[] index=new int[q];
		for(int i=0;i<q;i++)
			index[i]=i;
		
		X.setCol(new Vect().ones(n), 0);
		for(int i=1;i<q-1;i++)
			X.el[i][i]=1;
		x.rand();
		X.setCol(x, q-1);
		
		boolean[] locked=new boolean[q];
	
		Mat Q=new Mat(q,q);
		Mat I=new Mat(q,q);
		I.eye();
		Mat K_hat=new Mat();
		Mat M_hat=new Mat();
		int iter;
		int nconv=0;
		SpMat Ks=this.deepCopy();
		Vect Ci=Ks.scale();
		SpMat L=Ks.ichol();
		Mat X2;
		for(iter=0;iter<200;iter++){
			X2=Ms.smul(X);
			for(int j=0;j<q;j++){				
				
				if(!locked[j])
				{
				x=X2.getColVect(index[j]);	
				x.timesVoid(Ci);
				x=solver.silentICCG(Ks,L,x,1e-6,1000,new Vect(nRow));	
				x.timesVoid(Ci);
				}

				else{
					x=X.getColVect(index[j]);	
					if(lamq.el[j]!=0)
					x=x.times(1/lamq.el[j]);
				}

				if( solver.terminate) break;	
			
				X_hat.setCol(x,index[j]);		
			}
			
			K_hat=X_hat.transp().mul(X2);
			M_hat=X_hat.transp().mul(Ms.smul(X_hat));	
			K_hat.eigGenSym(M_hat,Q);
			X=X_hat.mul(Q);

				lamqr=lamq.deepCopy();
				lamq=this.smul(X).normCol().div0(Ms.smul(X).normCol());
				index=lamq.bubble();
		
				
				errv=lamq.sub(lamqr).abs().div0(lamq);


				for( int i=0;i<p;i++)
					if(!locked[i] && errv.el[i]<errMax){
					
						locked[i]=true;
						nconv++;

					}
						
					
				//errv.hshow();
				System.out.println("Subspace Iteration: " +iter+"   | error: "+errv.el[p-1]+   "  |  Numver of converged eigenvalues: "+nconv);
				
			if( nconv==p || solver.terminate) {  break;}
		}

		Vect lam=new Vect(p);		
		for(int i=0;i<p;i++)
			lam.el[i]=lamq.el[i];

		System.out.println("Subspace Iteration: " +iter+"   error: "+errv.el[p-1]);
		
		lam.quickSort();
		
		return lam;	
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
	
	
	public Vect solveICCG(Vect b){

		Vect x=new Vect(this.nRow);
	

		 SpMatSolver solver=new SpMatSolver();

		 solver.terminate(false);
		
		SpMat  Ks=this.deepCopy();


		Vect Ci=Ks.scale(b);

			
		SpMat L=Ks.ichol();

			if(b.abs().max()>1e-8){
				x=solver.ICCG(Ks,L, b,1e-6,2000,x);
				
			//	x=solver.CG(Ks,b,1e-6,2000,x);
					
			}
		x.timesVoid(Ci);
			
		return x;
	}

}
