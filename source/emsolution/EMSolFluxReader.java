package emsolution;


import io.Loader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import math.Mat;
import math.Vect;
import math.util;
import fem.Model;



public class EMSolFluxReader {

	String regex="[:; ,\\t]+";


	public static void main(String[] args){

		EMSolFluxReader x=new EMSolFluxReader();
		
		//x.fetchB();
		//x.fetchJQelem(6513);
		//x.fetchEQAll();
	//	x.fetchMagAll();
	//	x.fetchMagAllEls();
		
		//x.getFluxAtlasOrig(3,41514,1);
		
		
		String modelFile=System.getProperty("user.dir")+"\\EMSol\\Hexa.txt";
		String mapFile=System.getProperty("user.dir")+"\\EMSol\\map.txt";
		
		x.getFluxAtlas(modelFile,mapFile,3,1);
		
	
	}
	
	public void fetchB(){
		
		int elNumb=26942;

		
			int nSteps=20;
	
			String fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\classic\\magnetic";

			Mat BB=extractFlux( fluxFile,3,nSteps,  elNumb);

			BB.show();
			
			util.plot(BB.getColVect(0));

		
		
	}
	
	
	public void fetchBH(){
				
		int elNumb=390;
			int nSteps=21;
			
		String bhfolder="C:\\Works\\Problems and Models\\Large-scale Test";
		Mat BH=	getBHcurve( bhfolder,3,nSteps, elNumb,0);
//		BH.show(); 
		Mat HB=new Mat(BH.size());
		HB.setCol(BH.getColVect(1), 0);
		HB.setCol(BH.getColVect(0), 1);
		HB.show();
		util.plot(BH);
////		util.plot(BH.getColVect(0));
		util.plot(BH.getColVect(1));
	}
	
	public void fetchJQelem(int elNumb){
		
		double[] time=new double[10000];
		int[] nStepsr=new int[1];
		//int elNumb=6513;
	

		
			String fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\classic\\current";
		
			
			int nD=1;
			Mat[] JJD=new Mat[nD+1];
			double[] dE=new double[nD+1];
			double[] dQ=new double[nD+1];
			int kd=1;
			JJD[0]=extractJQ( fluxFile,3, elNumb,nStepsr, time);
			Vect v=new Vect().ones(JJD[0].nCol);
			double Eclassic=JJD[0].mul(v).norm();


			int nSteps=nStepsr[0];
			util.pr(nSteps);
			Vect timex=new Vect(nSteps);
			for(int i=0;i<nSteps;i++)
				timex.el[i]=time[i];
			
			
			for(int i=1;i<=nD;i++){
				kd*=2;
	
			fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\domain"+kd+"\\current";
			
			
			JJD[i]=extractJQ( fluxFile,3, elNumb,nStepsr, time);
			
		
			dE[i]=Math.abs(JJD[i].mul(v).norm()-Eclassic)/Eclassic*100;
			
			
			util.pr("Error (%): "+dE[i]);
			}
			
		
		
		
	}
	
	public void fetchEQAll(){
		
	

			double[] time=new double[10000];
			int[] nStepsr=new int[1];
		
			String fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\classic\\current";

			
			int nD=1;
			double[][] JJD=new double[nD+1][2];
			double[] dE=new double[nD+1];
			double[] dQ=new double[nD+1];
			int kd=1;
			JJD[0]=extractEQ( fluxFile,3,nStepsr,time);
			
			int nSteps=nStepsr[0];
			util.pr(nSteps);
			Vect timex=new Vect(nSteps);
			for(int i=0;i<nSteps;i++)
				timex.el[i]=time[i];
			


			for(int i=1;i<=nD;i++){
				kd*=2;
	
			fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\domain"+kd+"\\current";
			
			
			JJD[i]=extractEQ(fluxFile,3,nStepsr,time);
			
		
			dE[i]=Math.abs(JJD[i][0]-JJD[0][0])/JJD[0][0]*100;
			
			dQ[i]=Math.abs(JJD[i][1]-JJD[0][1])/JJD[0][1]*100;
			
			
			util.pr(" dEJ Error (%): "+dE[i]);
			util.pr(" dQ Error (%): "+dQ[i]);
			}
		
		
	}
	
	public void fetchMagAll(){
		
		

		double[] time=new double[10000];
		int[] nStepsr=new int[1];
		

		String fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\classic\\magnetic";
		
		
		int nD=1;
		double[][] JJD=new double[nD+1][2];
		double[] dE=new double[nD+1];
		double[] dQ=new double[nD+1];
		int kd=1;
		JJD[0]=extractMag( fluxFile,3,nStepsr,time);
		
		int nSteps=nStepsr[0];

		Vect timex=new Vect(nSteps);
		for(int i=0;i<nSteps;i++)
			timex.el[i]=time[i];
		
		util.pr(" EMagClassic : "+JJD[0][0]);
		util.pr(" EMagClassicBm2 : "+JJD[0][1]);

		for(int i=1;i<=nD;i++){
			kd*=2;

		fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\domain"+kd+"\\magnetic";
		
		
		JJD[i]=extractMag(fluxFile,3,nStepsr,time);
		
		util.pr(" EMagDomain"+kd+" : "+JJD[i][0]);
		util.pr(" EMagDomain"+kd+"Bm2 : "+JJD[i][1]);
	
		dE[i]=Math.abs(JJD[i][0]-JJD[0][0])/JJD[0][0]*100;
		
		dQ[i]=Math.abs(JJD[i][1]-JJD[0][1])/JJD[0][1]*100;
		
		
		util.pr(" dEm Error (%): "+dE[i]);
		util.pr(" XX Error (%): "+dQ[i]);
		}
	
	
}
	
	public void fetchMagAllEls(){
		
		int dim=3;
		

		double[] time=new double[10000];
		int[] nStepsr=new int[1];
		

		String fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\classic\\magnetic";
		
		Mat[] mag0=extractMagEls( fluxFile,dim,time);
		

			
		
		fluxFile="C:\\Works\\Problems and Models\\Large-scale Test\\domain2\\magnetic";
		Mat[] mag1=extractMagEls( fluxFile,dim,time);

		
/*		for(int i=0;i<mag0[0].nRow;i++){
			if(Math.abs(mag1[0].el[i][0]-mag0[0].el[i][0])>1e-3){
			util.hshow(mag0[0].el[i]);
			util.hshow(mag1[0].el[i]);
			util.pr("");
			}
			
		}*/
		


		
		//util.plot(mag1[0].getColVect(0).sub(mag0[0].getColVect(0)));

	//	mag0[0].show();
		
		Mat diff=mag1[0].sub(mag0[0]);
		
		util.pr(diff.getColVect(0).abs().sum());
		
		Mat diffExpanded=new Mat(35000,2);

		double den=0;
		double en=0;
		Vect B1=new Vect(3);
		Vect B2=new Vect(3);
		Vect dB1=new Vect(3);
		for(int i=0;i<diff.nRow;i++){
			B1=new Vect(mag0[0].el[i][1],mag0[0].el[i][2],mag0[0].el[i][3]);
			B2=new Vect(mag1[0].el[i][1],mag1[0].el[i][2],mag1[0].el[i][3]);
			dB1=new Vect(diff.el[i][1],diff.el[i][2],diff.el[i][3]);
			//dB1=new Vect(mag1[0].el[i][1],mag0[0].el[i][2],mag0[0].el[i][3]);
			
			diffExpanded.el[(int)mag0[0].el[i][0]][0]=B1.dot(B1);
			diffExpanded.el[(int)mag1[0].el[i][0]][1]=B2.dot(B2);
		//	diffExpanded.el[(int)mag0[0].el[i][0]][0]=mag0[0].el[i][4];
		//	diffExpanded.el[(int)mag1[0].el[i][0]][1]=mag1[0].el[i][4];
			
		/*	if(diff.el[i][0]<.1){
			en+=B1.dot(B1);
			
			den+=dB1.dot(dB1);
			}
			else{
				
				util.hshow(mag0[0].el[i]);
				util.hshow(mag1[0].el[i]);
				util.pr("");
			}*/
		}
		
		util.plot(diffExpanded.getColVect(0));
		util.plot(diffExpanded.getColVect(1));
		//util.plot(diffExpanded.getColVect(1).sub(diffExpanded.getColVect(0)));
		for(int i=0;i<diffExpanded.nRow;i++){
			den+=Math.abs(diffExpanded.el[i][1]-diffExpanded.el[i][0]);
		}
			
		
		util.pr(en+" <<<<<<<<<<<<<<< "+den);
		
		util.pr(" <<<<end %  "+den/en*100);
		
	/*	String fout="C:\\Works\\Problems and Models\\Large-scale Test\\classic\\magneticDiff";

		try{
			DecimalFormat formatter= new DecimalFormat("0.000000000E000");
			DecimalFormat intFormatter= new DecimalFormat("00000000");
			
			PrintWriter pwBun = new PrintWriter(new BufferedWriter(new FileWriter(fout)));	
			
			pwBun.println("STEP    (2I5,E18.0)");
			pwBun.println("    1    1  1.-003");
			pwBun.println("EVAL   1(I8,6E17.0)");
			for(int i=0;i<diff.nRow;i++){
				
				pwBun.print((int)mag0[0].el[i][0]+"\t");
				for(int k=1;k<diff.nCol;k++)
				pwBun.print(formatter.format(diff.el[i][k])+"\t");
				
				pwBun.println();
				
			}
			pwBun.println("\t"+(-1));
			
			pwBun.close();


	System.out.println("Flux was written to "+fout);

		}
	catch(Exception e){System.err.println("error");	e.printStackTrace(); }
		*/
		

}
	

	

	public Mat extractFlux(String bbf,int dim, int nSteps, int nelem){

		Vect[] B=new Vect[nSteps];
		
		for(int i=0;i<nSteps;i++)
		B[i]=new Vect(dim);
		
		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
	
			for(int i=0;i<nSteps;i++){
						
			while(!line.startsWith("STEP")){
			line=br.readLine();
					}
		
			String ss="";
			while(!ss.equals(Integer.toString(nelem))){
				line=br.readLine();
				sp=line.split(regex);

				ss=sp[1];
				
						}
			sp=line.split(regex);
			for(int j=0;j<dim;j++)
				B[i].el[j]=Double.parseDouble(sp[2+j]);
		
			}

				
		
	br.close();
	fr.close();
	
	Mat result=new Mat(nSteps,dim);
	for(int i=0;i<nSteps;i++)
		result.el[i]=B[i].el;
	
return result;
	
		}


		catch(Exception e){
			System.err.println("error");	e.printStackTrace(); 
			return null;
		}
		

	
	}
	
	
	public Mat getBHcurve(String bhfolder,int dim, int nSteps, int nelem,double angdeg){

		
		String fileB=bhfolder+"\\magnetic";
		String fileH=bhfolder+"\\magnetization";
		
		Mat BH1=new Mat(nSteps,2);
		Mat BH=null;
		Vect B=new Vect(dim);
		Vect H=new Vect(dim);
		Vect er=new Vect(dim);
		er.el[0]=Math.cos(angdeg*Math.PI/180);
		er.el[1]=Math.sin(angdeg*Math.PI/180);
		
			String regex="[ ,\\t]+";
		try{

			
			
			File f=new File(fileB);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
	
			int ix=0;
			
			for(int i=0;i<nSteps;i++){
					
			while(!line.startsWith("STEP") ){
				line=br.readLine();
	

				if(line==null) break;
			
					}
			if(line==null) break;

			String ss="";
			while(!ss.equals(Integer.toString(nelem))){
	
				line=br.readLine();

				sp=line.split(regex);

				ss=sp[1];
				
						}
			sp=line.split(regex);
			for(int j=0;j<dim;j++)
				B.el[j]=Double.parseDouble(sp[2+j]);
		
	ix++;

			BH1.el[i][1]=B.dot(er);
			//BH.el[i][1]=B.el[1];
		
		
			}
			
			 f=new File(fileH);
			 fr=new FileReader(f);
			 br = new BufferedReader(fr);
			 line="";
			 sp=new String[15];
	
				for(int i=0;i<nSteps;i++){
					
					while(!line.startsWith("STEP") ){
						line=br.readLine();
				

						if(line==null) break;
					
							}
					if(line==null) break;
				
					String ss="";
					while(!ss.equals(Integer.toString(nelem))){

				line=br.readLine();
				sp=line.split(regex);

				ss=sp[1];
				
						}
			sp=line.split(regex);
			for(int j=0;j<dim;j++)
				H.el[j]=Double.parseDouble(sp[2+j]);
		
	

			BH1.el[i][0]=H.dot(er);
			//BH.el[i][0]=H.el[0];
		
		
			}
	
		
				BH=new Mat(ix,2);
				for(int i=0;i<ix;i++)
					BH.el[i]=BH1.el[i];
			
			//util.plot(BH.getColVect(1));

		//	BH.getColVect(1).show();
			
			br.close();
			fr.close();

	
	
		}


		catch(Exception e){System.err.println("error");	e.printStackTrace(); }

return BH;
	
	
	}
	
	public Mat extractJQ(String bbf,int dim,int nelem,int[] nStepsr,double[] time){

		int nSteps=10000;
		Vect[] B=new Vect[nSteps];
		
		double[] Q=new double[nSteps];
		
		for(int i=0;i<nSteps;i++)
		B[i]=new Vect(dim);
		
		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
		
			
			int stepCount=0;

		
		
			while(line!=null){
				
				while(!line.startsWith("STEP")){
				line=br.readLine();
				if(line==null) break;
				}
		
				line=br.readLine();
				if(line==null) break;
				
				sp=line.split(regex);

				time[stepCount]=Double.parseDouble(sp[3]);
			
	
					while(!line.startsWith("EVAL")){
					line=br.readLine();
			
							}
					
				
		

			int ne=0;
			boolean elFound=false;
			boolean isInt=true;
			while(!elFound){
				
				line=br.readLine();
				if(line==null) break;
				
				sp=line.split(regex);
				
			try{
			ne=Integer.parseInt(sp[1]);
			}
			catch(Exception e){
			isInt=false;
			}
			
			if(isInt && ne==nelem)
				elFound=true;
			else continue;

			
		
			for(int j=0;j<dim;j++)
				B[stepCount].el[j]=Double.parseDouble(sp[2+j]);
			
			Q[stepCount]=Double.parseDouble(sp[2+dim]);

		
			stepCount++;
			}
			}

		
	br.close();
	fr.close();
	
	nStepsr[0]=stepCount;
	
	Mat result=new Mat(nStepsr[0],dim+1);
	for(int i=0;i<nStepsr[0];i++){
		for(int j=0;j<dim;j++)
		result.el[i][j]=B[i].el[j];
		
		result.el[i][dim]=Q[i];
	}
	
return result;
	
		}


		catch(Exception e){
			System.err.println("error");	e.printStackTrace(); 
			return null;
		}
		

	
	}
	
	public double[] extractEQ(String bbf,int dim,int[] nStepsr,double[] time){

	
		
		double[] EQ=new double[2];

		
		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
		
			
			int stepCount=0;
			int nelCount=0;

		
		
			while(line!=null){
				
				while(!line.startsWith("STEP")){
				line=br.readLine();
				if(line==null) break;
	
				}
		
				line=br.readLine();
				if(line==null) break;
				
				sp=line.split(regex);

				time[stepCount]=Double.parseDouble(sp[3]);
				stepCount++;
	
					while(!line.startsWith("EVAL")){
					line=br.readLine();
			
							}
					
					if(stepCount!=1) continue;
		

			int ne=0;
			boolean elNum=true;
			
			nelCount=0;
			
			while(elNum){
				
				line=br.readLine();
				if(line==null) break;
				
				sp=line.split(regex);
				
			try{
			ne=Integer.parseInt(sp[1]);
			}
			catch(Exception e){
				elNum=false;
			}
			
			if(ne<1)
				elNum=false;
			


		
			if(!elNum) continue;
			
			
			nelCount++;

			
		//	if(Math.abs(ne)>00000) elNum=false;
			
	
			Vect B=new Vect(dim);
			for(int j=0;j<dim;j++)
				B.el[j]=Double.parseDouble(sp[2+j]);
			
			EQ[0]+=B.norm2();
		//	util.pr(sp[2+dim]);
			double Q=Double.parseDouble(sp[2+dim]);

			EQ[1]+=Q;
			
			
			line=br.readLine();
			
			sp=line.split(regex);
				
			
			}
			
			}

			nStepsr[0]=stepCount;
		
	br.close();
	fr.close();


		
return EQ;
	
		}


		catch(Exception e){
			System.err.println("error");	e.printStackTrace(); 
			return null;
		}
		

	
	}
	
	public double[] extractMag(String bbf,int dim,int[] nStepsr,double[] time){

	
		
		double[] EQ=new double[2];

		
		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
		
			
			int stepCount=0;
			int nelCount=0;

		
		
			while(line!=null){
				
				while(!line.startsWith("STEP")){
				line=br.readLine();
				if(line==null) break;
	
				}
		
				line=br.readLine();
				if(line==null) break;
				
				sp=line.split(regex);

				time[stepCount]=Double.parseDouble(sp[3]);
	
					while(!line.startsWith("EVAL")){
					line=br.readLine();
			
							}
					
					//if(stepCount!=0) continue;

		

			int ne=0;
			boolean elNum=true;
			
			nelCount=0;
			
			while(elNum){
				
				line=br.readLine();
				if(line==null) break;
				
				sp=line.split(regex);
				
			try{
			ne=Integer.parseInt(sp[1]);
			}
			catch(Exception e){
				elNum=false;
			}
			
			if(ne<1)
				elNum=false;
			


		
			if(!elNum) continue;
			
			
			nelCount++;
			
		
		//	if(Math.abs(ne)>00000) elNum=false;
			
	
			Vect B=new Vect(dim);
			for(int j=0;j<dim;j++)
				B.el[j]=Double.parseDouble(sp[2+j]);
			
			EQ[0]+=B.norm2();
		//	util.pr(sp[2+dim]);
			double Q=Double.parseDouble(sp[2+dim]);

			EQ[1]+=Q*Q;

			
			line=br.readLine();
			
			sp=line.split(regex);
				
			
			}
			
			//util.pr(nelCount);
			//util.pr(EQ[1]);
			
			
			stepCount++;
			}

			nStepsr[0]=stepCount;
		
	br.close();
	fr.close();


		
return EQ;
	
		}


		catch(Exception e){
			System.err.println("error");	e.printStackTrace(); 
			return null;
		}
		

	
	}
	
	
public Mat[] extractMagEls(String bbf,int dim,double[] time){

	
		
		double[] EQ=new double[2];
		Mat[] results1=new Mat[1000];

		
		String regex="[ ,\\t]+";
		try{

			File f=new File(bbf);
			FileReader fr=new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line="";
			String[] sp=new String[15];
		
			
			int stepCount=0;
			int nelCount=0;

		
		
			while(line!=null){
				
				while(!line.startsWith("STEP")){
				line=br.readLine();
				if(line==null) break;
	
				}
		
				line=br.readLine();
				if(line==null) break;
				
				sp=line.split(regex);

				time[stepCount]=Double.parseDouble(sp[3]);
				
	
					while(!line.startsWith("EVAL")){
					line=br.readLine();
			
							}
					
					

					if(stepCount>0) break;

			int ne=0;
			boolean elNum=true;
			
			nelCount=0;
			
			int NeMax=50000;
			
			Mat data=new Mat(NeMax,dim+2);
			
			while(elNum){
				
				line=br.readLine();

				if(line==null) break;
				
				sp=line.split(regex);
				
			try{
			ne=Integer.parseInt(sp[1]);
			}
			catch(Exception e){
				elNum=false;
			}
			
			if(ne<1)
				elNum=false;
			


		
			if(!elNum) continue;
	
			
		
		//	if(Math.abs(ne)>00000) elNum=false;
			data.el[nelCount][0]=ne;
	
			Vect B=new Vect(dim);
			for(int j=0;j<dim;j++){
				B.el[j]=Double.parseDouble(sp[2+j]);
				data.el[nelCount][j+1]=B.el[j];
			}
			
			EQ[0]+=B.norm2();
		//	util.pr(sp[2+dim]);
			double Q=Double.parseDouble(sp[2+dim]);
			data.el[nelCount][dim+1]=Q;

			EQ[1]+=Q;

			
			line=br.readLine();
			
			sp=line.split(regex);
			
			
			
			nelCount++;
			
			}



			Vect numbs=data.getColVect(0);

			int[] index=numbs.bubble();
		
			results1[stepCount]=new  Mat(nelCount,dim+2);
			for(int i=0;i<nelCount;i++){
				results1[stepCount].el[nelCount-1-i]=data.el[index[NeMax-1-i]];
				//results1[stepCount].el[i]=data.el[i];
			}
			
	
			//numbs.hshow();
			//util.pr(numbs.sum());

			
			stepCount++;
		
			}

				
	br.close();
	fr.close();
	
	Mat[] results=new Mat[stepCount];
	for(int i=0;i<stepCount;i++)
		results[i]=results1[i];


		
return results;
	
		}


		catch(Exception e){
			System.err.println("error");	e.printStackTrace(); 
			return null;
		}
		

	
	}

public void getFluxAtlasOrig(String bbf,int dim, int nEmax, int nStepMax){




	Model model;
	if(dim==3)
		model=new Model(1,nEmax,1,"hexahedron");
	else
		model=new Model(1,nEmax,1,"quadrangle");
	
	
	String regex="[ ,\\t]+";
	try{

		File f=new File(bbf);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		String line="";
		String[] sp=new String[15];
	
		
		int stepCount=0;
		int nelCount=0;

	
	
		while(line!=null){
			
			while(!line.startsWith("STEP")){
			line=br.readLine();
			if(line==null) break;

			}
	
			line=br.readLine();
			if(line==null) break;
			
			sp=line.split(regex);

	

				while(!line.startsWith("EVAL")){
				line=br.readLine();
						}
				
				

				//if(stepCount>0) break;

		int ne=0;
		boolean elNum=true;
		
		nelCount=0;
		
		
	
		while(elNum){
			
			line=br.readLine();

			if(line==null) break;
			
			sp=line.split(regex);
			
		try{
		ne=Integer.parseInt(sp[1]);
		}
		catch(Exception e){
			elNum=false;
		}
		
		if(ne<1)
			elNum=false;
		


	
		if(!elNum) continue;

		
	
	//	if(Math.abs(ne)>00000) elNum=false;

		Vect B=new Vect(dim);
		for(int j=0;j<dim;j++){
			B.el[j]=Double.parseDouble(sp[2+j]);
		}
		
		model.element[ne].setB(B);
		
	

		
		line=br.readLine();
		
		sp=line.split(regex);
		
		
		
		nelCount++;
		
		
		}


		String fout=System.getProperty("user.dir")+"\\EMSol\\fluxes\\flux"+stepCount+".txt";
		
		model.writeB(fout);


		stepCount++;
		
		if(stepCount==nStepMax) break;
	
		}

			
br.close();
fr.close();


	

	}


	catch(Exception e){
		System.err.println("error");	e.printStackTrace(); 
	}

}

public void getFluxAtlas(String modelFile,String mapFile,int dim, int nStepMax){


	String bbf=util.getFile();
	
		
		Model model=new Model(modelFile);
		
		
		double[][] map=new Loader().loadArrays(model.numberOfElements,2,mapFile);
		
		int maxel=100000;
		int[] map2=new int[maxel];
		
		for(int i=0;i<map.length;i++)
			map2[(int)map[i][0]]=(int)map[i][1];
		util.show(map2);
			
/*	Model model;
	if(dim==3)
		model=new Model(1,nEmax,1,"hexahedron");
	else
		model=new Model(1,nEmax,1,"quadrangle");
	*/


	
	String regex="[ ,\\t]+";
	try{

		File f=new File(bbf);
		FileReader fr=new FileReader(f);
		BufferedReader br = new BufferedReader(fr);
		String line="";
		String[] sp=new String[15];
	
		
		int stepCount=0;
		int nelCount=0;

	
	
		while(line!=null){
			
			while(!line.startsWith("STEP")){
			line=br.readLine();
			if(line==null) break;

			}
	
			line=br.readLine();
			if(line==null) break;
			
			sp=line.split(regex);

	

				while(!line.startsWith("EVAL")){
				line=br.readLine();
						}
				
				

				//if(stepCount>0) break;

		int ne=0;
		boolean elNum=true;
		
		nelCount=0;
		
		
	
		while(elNum){
			
			line=br.readLine();

			if(line==null) break;
			
			sp=line.split(regex);
			
		try{
		ne=Integer.parseInt(sp[1]);
		}
		catch(Exception e){
			elNum=false;
		}
		
		if(ne<1)
			elNum=false;
		


	
		if(!elNum) continue;

		
	
	//	if(Math.abs(ne)>00000) elNum=false;

		Vect B=new Vect(dim);
		for(int j=0;j<dim;j++){
			B.el[j]=Double.parseDouble(sp[2+j]);
		}
		
		
		nelCount++;

		int nemap=map2[ne];
		if(nemap<=model.numberOfElements)
		model.element[nemap].setB(B);
		
	

		
		line=br.readLine();
		
		sp=line.split(regex);
		
		
		

		
		
		}


		String fout=System.getProperty("user.dir")+"\\EMSol\\fluxes\\flux"+stepCount+".txt";
		
		model.writeB(fout);


		stepCount++;
	
		if(stepCount==nStepMax) break;
		
		}

			
br.close();
fr.close();


	

	}


	catch(Exception e){
		System.err.println("error");	e.printStackTrace(); 
	}
	


}
}
