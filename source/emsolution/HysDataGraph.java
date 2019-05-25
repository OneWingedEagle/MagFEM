package emsolution;

import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.cos;
import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Scanner;

import math.Mat;
import math.Vect;
import math.util;
import fem.Model;



public class HysDataGraph {

	Mat[][] BH;
	Mat rotLoss;
	Mat[][] BHAni;
	int nSet;
	double[] Bs, Hs,angle;
	int[] nInitial,nMajor,nSymLoop,nDescending,nAscending,nTotCurves,nAni,Lani;
	String regex="[:; ,\\t]+";


	public static void main(String[] args){

		HysDataGraph pg=new HysDataGraph();

		pg.loadHysData();
	//	pg.combineXYHysData();
	//	pg.loadAltHysData();
	//	pg.loadDoshishaAltHysData();
		//pg.loadRotHysData();
	//	pg.loadAngData();
		//pg.loadAngSymData();
	//	pg.loadHysDataOld();
	}

	public boolean loadHysData(){

		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_B入力対称ループhts_data\\hys_data";

		String file="C:\\Works\\PlayModel\\hys_data";
		//String file="C:\\Works\\HVID\\hys_data";
	//	String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		


		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			
			line=br.readLine();
			
			sp=line.split(regex);
			if(sp.length<3) this.nSet=1;
			else this.nSet=Integer.parseInt(sp[2]);
			
			this.Bs=new double[this.nSet];
			this.Hs=new double[this.nSet];
			this.nInitial=new int[this.nSet];
			this.nMajor=new int[this.nSet];
			this.nSymLoop=new int[this.nSet];
			this.nDescending=new int[this.nSet];
			this.nAscending=new int[this.nSet];
			this.nTotCurves=new int[this.nSet];
			this.nAni=new int[this.nSet];
			this.Lani=new int[this.nSet];
			this.angle=new double[this.nSet];
			

			BH=new Mat[nSet][];
			BHAni=new Mat[nSet][];
	

				line=br.readLine();
				if(line.startsWith("*"))
					line=br.readLine();
				

				sp=line.split(regex);	

				this.Bs[0]=Double.parseDouble(sp[0]);
				for(int ia=0;ia<this.nSet;ia++){
				this.Hs[ia]=Double.parseDouble(sp[ia+1]);
				}

				line=br.readLine();
				if(line.startsWith("*"))
					line=br.readLine();
				
				sp=line.split(regex);	
				this.nInitial[0]=Integer.parseInt(sp[0]);
				this.nMajor[0]=Integer.parseInt(sp[1]);
				this.nSymLoop[0]=Integer.parseInt(sp[2]);
				this.nDescending[0]=Integer.parseInt(sp[3]);
				this.nAscending[0]=Integer.parseInt(sp[4]);

				this.nTotCurves[0]=this.nInitial[0]+this.nMajor[0]+this.nSymLoop[0]+this.nDescending[0]+this.nAscending[0];
				

				for(int ia=0;ia<nSet;ia++)
				BH[ia]=new Mat[this.nTotCurves[0]];

				line=br.readLine();
				if(line.startsWith("*"))
					line=br.readLine();
				
				sp=line.split(regex);	
				for(int ia=0;ia<nSet;ia++)
					this.angle[ia]=Double.parseDouble(sp[ia]);

				
				for(int i=0;i<nTotCurves[0];i++){
				line=br.readLine();
				if(line.startsWith("*"))
					line=br.readLine();
				int L=Integer.parseInt(line);
				
				for(int ia=0;ia<nSet;ia++)
				BH[ia][i]=new Mat(L,2);

				for( int p=0;p<L;p++){
					line=br.readLine();
					sp=line.split(regex);	
					double[] hp=getCSV(line);
					for(int ia=0;ia<nSet;ia++){
						BH[ia][i].el[p][0]=hp[ia];
						BH[ia][i].el[p][1]=hp[nSet];
					}
					

				}

				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				int nLoss=Integer.parseInt(line);
				rotLoss=new Mat(nLoss,2);
				line=br.readLine();

				for( int i=0;i<nLoss;i++){
					line=br.readLine();
					double[] BL=getCSV(line);
					rotLoss.el[i][0]=BL[0];
					rotLoss.el[i][1]=BL[1];
				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				double[] dd=this.getCSV(line);
				
				int ia=0;
				Lani[ia]=(int)dd[0];
				nAni[ia]=(int)dd[1];
			
				line=br.readLine();

				BHAni[ia]=new Mat[nAni[ia]];

				for( int i=0;i<nAni[ia];i++)
					BHAni[ia][i]=new Mat(Lani[ia],3);

				for( int i=0;i<Lani[ia];i++){
					line=br.readLine();
			
					double[] BL=getCSV(line);
					double B=BL[0];
		
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][0]=B;
						BHAni[ia][j].el[i][1]=BL[j+1];
					}
				}

				line=br.readLine();
				for( int i=0;i<Lani[ia];i++){
					line=br.readLine();
					double[] BL=getCSV(line);
			
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][2]=BL[j+1];
					}
				}



			
			//util.plotBunch(BH[0]);
		//	util.plotBunch(BH[0],6,11);
			
		//createAngleDepData();

			
			boolean ani=false;
			
			if(ani){
				
				Mat[] BHaniT=new Mat[BHAni[0].length];
				
	
				for(int i=0;i<BHAni[0].length;i++){
					
					BHaniT[i]=new Mat(BHAni[0][i].nRow,2);
					
					double ang=i*Math.PI/18;
					
					Vect er=new Vect(Math.cos(ang),Math.sin(ang));

				for(int j=0;j<BHAni[0][i].nRow;j++){
					double Ht=new Vect(BHAni[0][i].el[j][1], BHAni[0][i].el[j][2]).dot(er);
					BHaniT[i].el[j][0]=Ht;
					BHaniT[i].el[j][1]=BHAni[0][i].el[j][0];
					
				}
			
				}
				


			}
			
		
			String filex="C:\\Works\\PlayModel\\hys_dataNewFormat";
		//	this.writeHystDataColumns(BH, filex);
			
			int nSetOut=7;
			Mat[][] BHs=new Mat[nSetOut][BH[0].length];
			
			for(int i=0;i<nSetOut;i++)
				for(int j=0;j<BH[i].length;j++)
					BHs[i][j]=BH[i][j].deepCopy();
			
			PlayModel2D pm=new PlayModel2D();
			
			pm.writeHystData(BHs, filex);
			
	/*	
		//	util.plotBunch(BHs[2]);
	
			boolean distill=false;
			if(distill){
			Mat[][] BHdist=distHysData(BH);
			
			
			for(int i=0;i<BHdist.length;i++)
				for(int j=0;j<BHdist[i].length;j++)
					for(int k=0;k<BHdist[i][j].nRow;k++)
						BHdist[i][j].el[k][1]*=1+Math.abs(k-BHdist[i][j].nRow/2)*.005;
						
			//util.plotBunch(BHdist[0]);
			
		boolean	average=false;
	
			if(average){

			Mat[] BHsAv=new Mat[BH[0].length];
			
			Mat[] BHsAni=new Mat[0];
	
			//Mat[] BHsrAv=new Mat[BH[0].length];
			for(int i=0;i<BHsAv.length;i++){
				Mat M=BH[0][i].deepCopy();
				for(int j=1;j<nSet;j++)
					M=M.add(BH[j][i]);
		
				BHsAv[i]=M.times(1.0/nSet);

				//BHsrAv[i]=pm.getHBij(BHsAv[i],0);
			}
			
			String fileav="C:\\Works\\HVID\\hys_dataHAvdist";
			this.writeHystDataAv(BHsAv, BHsAni, fileav);

				util.plotBunch(BHsAv);
			}

			
			//Mat[] hysDataAv=new Mat[BH[0].length];
			
		
			if(BHdist.length==1){
				String filed="C:\\Works\\HVID\\hys_dataHAvdist";
			this.writeHystDataAv(BHdist[0], BHAni[0], filed);
			}
			else{
				String filed="C:\\Works\\HVID\\hys_dataHdist";
			new PlayModel2D().writeHystData(BHdist, filed);
			}
			}*/
		
			return true;

		}

		catch(IOException e){System.err.println("Error in loading BH data file.");
		return false;
		}

	}	
	
	
	public boolean combineXYHysData(){

		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_B入力対称ループhts_data\\hys_data";

		String file="C:\\Works\\PlayModelAni\\hys_datax";
		//String file="C:\\Works\\HVID\\hys_data";
	//	String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		

		Mat[][][] BHxy=new Mat[2][][];
		
		for(int k=0;k<2;k++){
			
			if(k==1)
			file="C:\\Works\\PlayModelAni\\hys_datay";

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			
			line=br.readLine();
			
			sp=line.split(regex);
			if(sp.length<3) this.nSet=1;
			else this.nSet=Integer.parseInt(sp[2]);
			
			this.Bs=new double[this.nSet];
			this.Hs=new double[this.nSet];
			this.nInitial=new int[this.nSet];
			this.nMajor=new int[this.nSet];
			this.nSymLoop=new int[this.nSet];
			this.nDescending=new int[this.nSet];
			this.nAscending=new int[this.nSet];
			this.nTotCurves=new int[this.nSet];
			this.nAni=new int[this.nSet];
			this.Lani=new int[this.nSet];
			this.angle=new double[this.nSet];
			
		

			BH=new Mat[nSet][];
			BHAni=new Mat[nSet][];
	

				line=br.readLine();
				if(line.startsWith("*"))
					line=br.readLine();
				

				sp=line.split(regex);	

				this.Bs[0]=Double.parseDouble(sp[0]);
				for(int ia=0;ia<this.nSet;ia++){
				this.Hs[ia]=Double.parseDouble(sp[ia+1]);
				}

				line=br.readLine();
				if(line.startsWith("*"))
					line=br.readLine();
				
				sp=line.split(regex);	
				this.nInitial[0]=Integer.parseInt(sp[0]);
				this.nMajor[0]=Integer.parseInt(sp[1]);
				this.nSymLoop[0]=Integer.parseInt(sp[2]);
				this.nDescending[0]=Integer.parseInt(sp[3]);
				this.nAscending[0]=Integer.parseInt(sp[4]);

				this.nTotCurves[0]=this.nInitial[0]+this.nMajor[0]+this.nSymLoop[0]+this.nDescending[0]+this.nAscending[0];
				
				
				if(k==0){
					BHxy=new Mat[2][nSet][this.nTotCurves[0]];
				}

				for(int ia=0;ia<nSet;ia++)
				BH[ia]=new Mat[this.nTotCurves[0]];

				line=br.readLine();
				if(line.startsWith("*"))
					line=br.readLine();
				
				sp=line.split(regex);	
				for(int ia=0;ia<nSet;ia++)
					this.angle[ia]=Double.parseDouble(sp[ia]);

				
				for(int i=0;i<nTotCurves[0];i++){
				line=br.readLine();
				if(line.startsWith("*"))
					line=br.readLine();
				int L=Integer.parseInt(line);
				
				for(int ia=0;ia<nSet;ia++)
				BH[ia][i]=new Mat(L,2);

				for( int p=0;p<L;p++){
					line=br.readLine();
					sp=line.split(regex);	
					double[] hp=getCSV(line);
					for(int ia=0;ia<nSet;ia++){
						BH[ia][i].el[p][0]=hp[ia];
						BH[ia][i].el[p][1]=hp[nSet];
					}
					

				}

				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				int nLoss=Integer.parseInt(line);
				rotLoss=new Mat(nLoss,2);
				line=br.readLine();

				for( int i=0;i<nLoss;i++){
					line=br.readLine();
					double[] BL=getCSV(line);
					rotLoss.el[i][0]=BL[0];
					rotLoss.el[i][1]=BL[1];
				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				double[] dd=this.getCSV(line);
				
				int ia=0;
				Lani[ia]=(int)dd[0];
				nAni[ia]=(int)dd[1];
			
				line=br.readLine();

				BHAni[ia]=new Mat[nAni[ia]];

				for( int i=0;i<nAni[ia];i++)
					BHAni[ia][i]=new Mat(Lani[ia],3);

				for( int i=0;i<Lani[ia];i++){
					line=br.readLine();
			
					double[] BL=getCSV(line);
					double B=BL[0];
		
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][0]=B;
						BHAni[ia][j].el[i][1]=BL[j+1];
					}
				}

				line=br.readLine();
				for( int i=0;i<Lani[ia];i++){
					line=br.readLine();
					double[] BL=getCSV(line);
			
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][2]=BL[j+1];
					}
				}


			boolean ani=false;
			
			if(ani){
				
				Mat[] BHaniT=new Mat[BHAni[0].length];
				
	
				for(int i=0;i<BHAni[0].length;i++){
					
					BHaniT[i]=new Mat(BHAni[0][i].nRow,2);
					
					double ang=i*Math.PI/18;
					
					Vect er=new Vect(Math.cos(ang),Math.sin(ang));

				for(int j=0;j<BHAni[0][i].nRow;j++){
					double Ht=new Vect(BHAni[0][i].el[j][1], BHAni[0][i].el[j][2]).dot(er);
					BHaniT[i].el[j][0]=Ht;
					BHaniT[i].el[j][1]=BHAni[0][i].el[j][0];
					
				}
			
				}
				


			}
		}

		catch(IOException e){System.err.println("Error in loading BH data file.");
		return false;
		}
		
	
		
		for(int ia=0;ia<this.nSet;ia++)
			for(int i=0;i<BH[0].length;i++)
				BHxy[k][ia][i]=BH[ia][i].deepCopy();
		
		}
		
		file="C:\\Works\\PlayModelAni\\hys_dataxy";
		
		util.plotBunch(BHxy[0][4]);	
		new PlayModel2D().writeHystDataXY(BHxy[0], BHxy[1], file, this.nInitial[0]);
		
		return true;


	}	
	
	
	public void loadRotHysData(){


		String file="C:\\Works\\PlayModel\\hysRotation";

		
	 HystDataLoader loader=new HystDataLoader();
	 
		

	int nTot=17;


	Mat BHij=new Mat(loader.loadArrays(1024,4*nTot,file));

		
		Mat[] BB=new Mat[nTot];
		Mat[] HH=new Mat[nTot];
		
		for(int i=0;i<BB.length;i++){
			BB[i]=new Mat(BHij.nRow,2);
			HH[i]=new Mat(BHij.nRow,2);
			
			BB[i].setCol(BHij.getColVect(4*i), 0);
			BB[i].setCol(BHij.getColVect(4*i+1), 1);
			
			HH[i].setCol(BHij.getColVect(4*i+2), 0);
			HH[i].setCol(BHij.getColVect(4*i+3), 1);
		}
	
		
	//util.plotBunch(HH,9);
	util.plot(HH[4]);
	HH[16].show();
		
		//util.plot(BB[12]);

		


	}
	
	public void loadAltHysData(){

		int nTot=18;
		int nset=7;
	PlayModel2D pm=new PlayModel2D();		
	Mat[][] BHsDns=new Mat[nset][nTot];
	 
	


	
	for(int ia=0;ia<nset;ia++){
		int deg=ia*15;

		String file="C:\\Works\\PlayModel\\KitaoData\\symmetricData\\hysAlt"+deg;

		
	 HystDataLoader loader=new HystDataLoader();



	Mat BHij=new Mat(loader.loadArrays(1024,4*nTot,file));


		Mat[] BB=new Mat[nTot];
		Mat[] HH=new Mat[nTot];
		
		for(int i=0;i<nTot;i++){
	
			BB[nTot-1-i]=new Mat(BHij.nRow,2);
			HH[nTot-1-i]=new Mat(BHij.nRow,2);
			
			BB[nTot-1-i].setCol(BHij.getColVect(4*i), 0);
			BB[nTot-1-i].setCol(BHij.getColVect(4*i+1), 1);
			
			HH[nTot-1-i].setCol(BHij.getColVect(4*i+2), 0);
			HH[nTot-1-i].setCol(BHij.getColVect(4*i+3), 1);
		}

		Vect er=new Vect(Math.cos(deg*Math.PI/180),Math.sin(deg*Math.PI/180));
		Vect en=new Vect(-Math.sin(deg*Math.PI/180),Math.cos(deg*Math.PI/180));
		 for(int i=0;i<nTot;i++)
		 {
			 
			 int skip=8;
			
			 int L=BB[i].nRow/skip;
			Mat  BH1=new Mat(L,2);
		

			 for(int j=0;j<L;j++){
				 BH1.el[j][1]=new Vect(BB[i].el[j*skip]).dot(er);
				// BH1.el[j][0]=new Vect(HH[i].el[j*skip]).dot(er);
				//BH1.el[j][0]=new Vect(HH[i].el[j*skip]).dot(er);
				// BH1.el[j][0]=new Vect(HH[i].el[j*skip]).dot(er);
				 BH1.el[j][0]=HH[i].el[j*skip][1];
				// BH1.el[j][0]=HH[i].el[j*skip][1];
			 }

			 
			 int n1=0;
			while( BH1.el[n1+1][1]>BH1.el[n1][1]){n1++;}
	
			int n2=n1+1;
			while( BH1.el[n2+1][1]<=BH1.el[n2][1]){n2++;}
util.pr(n1+" "+n2);
			Mat  BH2=new Mat(n2-n1+1,2);
			 
			 for(int j=0;j<BH2.nRow;j++){
				 BH2.el[j][0]=BH1.el[j+n1][0];
				 BH2.el[j][1]=BH1.el[j+n1][1];
			 }
			 
			 int L2=BH2.nRow;
			 
			 int Lp=2*(nTot-i)+1;
			 Vect B= new Vect().linspace(BH2.el[0][1],BH2.el[L2-1][1],Lp);
		
			
			 Mat  BH=new Mat(Lp,2);
			 for(int j=0;j<Lp;j++){
				 BH.el[j][1]=B.el[j];
				 BH.el[j][0]=pm.getH(BH2, B.el[j]);
				// BH.el[j][0]=HH[i].el[j*skip][0];
			 }
			 
			 BHsDns[ia][i]=BH.deepCopy();
		 }
		 
		
		 }

	
	//util.plotBunch(HH,9);
//	util.plot(HH[4]);
	//HH[16].show();
		
		//util.plot(BB[12]);
	Mat[][] BHs=new Mat[nset][nTot];
	Mat[][] BHs2=new Mat[nset][nTot];
	 for(int ia=0;ia<nset;ia++){
 
	 for(int i=0;i<nTot;i++)
	 {
			/* Vect B=new Vect().linspace(BHsDns[ia][i].el[0][1], BHsDns[ia][i].el[BHsDns[ia][i].nRow-1][1], L);

		 for(int j=0;j<L;j++){
			 BHs[ia][i].el[j][1]=B.el[j];
			 BHs[ia][i].el[j][0]=pm.getH(BHsDns[ia][i], B.el[j]);
		 }*/
		 
		 BHs[ia][i]= BHsDns[ia][i].deepCopy();
		 BHs2[ia][i]= BHsDns[ia][nTot-1-i].deepCopy();

	 }
	 
	 }
	
	// BHs[0][0].show();
	
	// util.plotBunch(BHs[0],4);
	// util.plotBunch(BHs[3]);
	 
	 util.plotBunch(BHs2[3],10);
	 
	//String fileout="C:\\Works\\HVID\\KitaoData\\symmetricData\\hys_data"+deg;
	String fileout="C:\\Works\\PlayModel\\hys_dataAllDensey";
	pm.writeHystData(BHs, fileout,0);

	}
	
	
	
	public void loadDoshishaAltHysData(){
		
		int deg=0;
		
		 HystDataLoader loader=new HystDataLoader();
		
		 
		 
		 int nset=7;
	
		 
		 Mat[][] hys=new Mat[nset][];
	
		 
		 for(int ia=0;ia<nset;ia++){
			 
			 deg=ia*15;
		 String file="C:\\Works\\PlayModel\\KitaoData\\symmetricData\\sym"+deg+".dat";


		 Mat[][] syms=loader.loadDataSym(file);
		 
		 hys[ia]=new Mat[syms[1].length+1];
				 
		 Mat init=new Mat(syms[1].length+1,2);
		 for(int i=1;i<init.nRow;i++){
			 init.el[i][0]=syms[1][i-1].el[0][0];
			 init.el[i][1]=syms[1][i-1].el[0][1];
		 }
		 
		 hys[ia][0]=init;
		 
		 for(int i=1;i<hys[ia].length;i++)
		 hys[ia][hys[ia].length-i]=syms[1][i-1];
		 
		 }
		 
		util.plotBunch(hys[0]);
		
		PlayModel2D pm=new PlayModel2D();

		
		Mat[][] BHs=new Mat[nset][];

		 int L=50;
		 
		 for(int ia=0;ia<nset;ia++){
			 int nc=hys[ia].length;
			 BHs[ia]=new Mat[nc];
			 
		 for(int i=0;i<nc;i++)
		 {
/*			 BHs[ia][i]=new Mat(L,2);
			 Vect B=new Vect().linspace(hys[ia][i].el[0][1], hys[ia][i].el[hys[ia][i].nRow-1][1], L);
	
			 for(int j=0;j<L;j++){
				 BHs[ia][i].el[j][1]=B.el[j];
				 BHs[ia][i].el[j][0]=pm.getH(hys[ia][i], B.el[j]);
			 }*/
			 
			 BHs[ia][i]= hys[ia][i].deepCopy();

		 }
		 }
		
		
		 util.plotBunch(BHs[0]);
		// util.plotBunch(BHs[3]);
		 
		// util.plotBunch(BHs[0]);
		 
		//String fileout="C:\\Works\\HVID\\KitaoData\\symmetricData\\hys_data"+deg;
		String fileout="C:\\Works\\PlayModel\\KitaoData\\symmetricData\\hys_dataAll";
		pm.writeHystData(BHs, fileout);
		
}
	
	
	
	public void loadAngData(){

		int deg=60;
		
		Vect er=new Vect(Math.cos(deg*180.0/Math.PI),Math.sin(deg*180/Math.PI));

		//String file="C:\\Works\\HVID\\KitaoData\\deg45";

		String file="C:\\Works\\HVID\\KitaoData\\deg"+deg;


	 HystDataLoader loader=new HystDataLoader();
	 

	Mat BHij=new Mat(loader.loadArrays(1024,72,file));
	
		int nc=72/4;
		
		double err=1e-4;
		
		Mat[] BB=new Mat[nc];
		Mat[] HH=new Mat[nc];
		
		Mat[] BH=new Mat[nc];
		Mat[] BHs=new Mat[nc];
		
		for(int i=0;i<BB.length;i++){
			BB[i]=new Mat(BHij.nRow,2);
			HH[i]=new Mat(BHij.nRow,2);
			
			BB[i].setCol(BHij.getColVect(4*i+0), 0);
			BB[i].setCol(BHij.getColVect(4*i+1), 1);
			
			HH[i].setCol(BHij.getColVect(4*i+2), 0);
			HH[i].setCol(BHij.getColVect(4*i+3), 1);
			
			int L=BB[i].nRow;
			
			BH[i]=new Mat(L,2);
			
			for(int j=0;j<L;j++){
				
				BH[i].el[j][0]=new Vect(HH[i].el[j]).dot(er);
				BH[i].el[j][1]=new Vect(BB[i].el[j]).dot(er);
			}
			

			double Bmax=BH[i].getColVect(1).max();
			double Bmin=BH[i].getColVect(1).min();
			
			int n1=0;
			int n2=0;
			int jx=0;

			
			while(BH[i].el[jx][1]<-err+Bmax /*|| BH[i].el[jx+1][0]>=BH[i].el[jx][0]*/){jx++;}
			n1=jx;

			while(BH[i].el[jx][1]>err+Bmin /*|| BH[i].el[jx+1][0]<=BH[i].el[jx][0]*/){jx++;}
			
			n2=jx;
			
			int Ls=n2-n1;
			
			BHs[i]=new Mat(Ls,2);
			
			for(int j=0;j<Ls;j++){
				
				BHs[i].el[j]=BH[i].el[j+n1];
			}
			
			
		}
	

		
	//util.plotBunch(HH,1);
	util.plotBunch(BH,10);
/*	BHs[0].show();	
	BH[0].show();*/
		//util.plot(HH[9]);

		


	}
	

	
	
	Mat[][] distHysData(Mat[][] BH){
		
		Mat[][] distBH=new Mat[BH.length][];
		
		PlayModel2D pm=new PlayModel2D();
		for(int k=0;k<distBH.length;k++){
			
			distBH[k]=new Mat[BH[k].length];
		
			
		for(int i=0;i<distBH[k].length;i++){
			int Lx=50;
			if (i==0) Lx=Lx/2;
			distBH[k][i]=new Mat(Lx,2);
			Vect B=new Vect();
			if(i==0){
				int jx=0;
				while(jx<BH[k][i].nRow-1 &&BH[k][i].el[jx+1][0]>=BH[k][i].el[jx][0]){jx++;};
		
				B=new Vect().linspace(0, BH[k][i].el[jx][1], Lx);
				
				if(B.el[Lx-1]<BH[k][i].el[BH[k][i].nRow-1][1]) B.el[Lx-1]=BH[k][i].el[BH[k][i].nRow-1][1];
			}
			else if(i==1){
				
				int jx1=0;
				while(jx1<BH[k][i].nRow-1&& BH[k][i].el[jx1+1][0]>=BH[k][i].el[jx1][0]){jx1++;};
				
				int r=BH[k][i].nRow;
				int jx2=r-1;
				while(jx2>0 &&BH[k][i].el[jx2-1][0]<=BH[k][i].el[jx2][0]){jx2--;};

				B=new Vect().linspace(BH[k][i].el[jx1][1], BH[k][i].el[jx2][1], Lx);
				
				if(B.el[Lx-1]>-BH[k][0].el[BH[k][0].nRow-1][1]) B.el[Lx-1]=-BH[k][0].el[BH[k][0].nRow-1][1];
				if(B.el[0]<BH[k][0].el[BH[k][0].nRow-1][1]) B.el[0]=BH[k][0].el[BH[k][0].nRow-1][1];
			}
			else{
				B=new Vect().linspace(BH[k][i].el[0][1], BH[k][i].el[BH[k][i].nRow-1][1], Lx);
			}

			for(int j=0;j<Lx;j++)
			{
				double H=pm.getH(BH[k][i], B.el[j]);
				distBH[k][i].el[j][0]=H;
			
				distBH[k][i].el[j][1]=B.el[j];
			}
			
			if(i==0) distBH[k][i].el[0][0]=0;
			
		}
		}
		

		return distBH;

		
	}
	
	public void writeHystDataAv(Mat[] BHs,Mat[] BHani, String file){
		
		int nSet=1;
		
		double Bseff= BHs[0].el[BHs[0].nRow-1][1];
		
		int nAni=BHani.length;
		int Lani=0;
		if(nAni>0)
		 Lani=BHani[0].nRow;
		int nInit=1;
		int nMajor=1;
		int nTot=BHs.length;
		int nSymLoops=nTot-2;
		int nDescending=0;
		int nAscending=0;
		
		double Hseff=BHs[0].el[BHs[0].nRow-1][0];

		DecimalFormat dfB=new DecimalFormat("#.00");
		DecimalFormat dfH=new DecimalFormat("#.0");


			try{
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));		


				pw.println(1+"\t"+1+"\t"+nSet+"\t"+0);
				pw.println("*Bs*Hs*");
				pw.println(Bseff+"\t"+Hseff);

				pw.println("* 初磁化曲線数 * メジャーループ数 * 対称ループ数 * 下降曲線数 * 上昇曲線数 *");
				pw.println(nInit+"\t"+nMajor+"\t"+nSymLoops+"\t"+nDescending+"\t"+nAscending);

				for(int i=0;i<nTot;i++){
					pw.println("*xxx");
					pw.println(BHs[i].nRow);
					for(int j=0;j<BHs[i].nRow;j++)
						pw.println(BHs[i].el[j][0]+"\t"+BHs[i].el[j][1]);
				}

				pw.println("* ----- 回転ヒステリシス損");
				pw.println("* B数 *");
				pw.println("0");
				pw.println("* B * 損失");
				pw.println("* ----- 異方性");
				pw.println("* B数 * 角度数 *");
				pw.println(Lani+"\t"+nAni);
				pw.println("* B * H ･････ *　磁化容易軸");
				
				for(int i=0;i<Lani;i++){
					pw.print(dfB.format(BHani[0].el[i][0])+"\t");
					for(int j=0;j<nAni;j++){
						pw.print(dfH.format(BHani[j].el[i][1])+"\t");
					}
					pw.println();
				}
				pw.println("* B * H ･････ *　磁化困難軸");
				for(int i=0;i<Lani;i++){
					pw.print(dfB.format(BHani[0].el[i][0])+"\t");
					for(int j=0;j<nAni;j++){
						pw.print(dfH.format(BHani[j].el[i][1])+"\t");
					}
					pw.println();
				}
				//BHani[0].show();

				util.pr("Simulated angle-dependent hysteresis data was written to "+file+".");

				pw.close();
			}
			catch(IOException e){}
			

	}
	
	public boolean loadHysDataOld(){

		/*String file=util.getFile();
	if(file==null || file.equals("") )return false;*/
		//	String file="C:\\Works\\HVID\\folder1\\data\\A_B入力対称ループhts_data\\hys_data";

		String file="C:\\Works\\HVID\\hys_dataHX";
		//String file="C:\\Works\\HVID\\hys_data";
	//	String file=System.getProperty("user.dir") + "\\hys_dataH.txt";
		


		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			
			line=br.readLine();
			
			sp=line.split(regex);
			if(sp.length<3) this.nSet=1;
			else this.nSet=Integer.parseInt(sp[2]);
			
			this.Bs=new double[this.nSet];
			this.Hs=new double[this.nSet];
			this.nInitial=new int[this.nSet];
			this.nMajor=new int[this.nSet];
			this.nSymLoop=new int[this.nSet];
			this.nDescending=new int[this.nSet];
			this.nAscending=new int[this.nSet];
			this.nTotCurves=new int[this.nSet];
			this.nAni=new int[this.nSet];
			this.Lani=new int[this.nSet];

			BH=new Mat[nSet][];
			BHAni=new Mat[nSet][];
	
		
		
			fr=new FileReader(file);
			br = new BufferedReader(fr);

			for(int ia=0;ia<this.nSet;ia++){

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				sp=line.split(regex);	

				this.Bs[ia]=Double.parseDouble(sp[0]);
				this.Hs[ia]=Double.parseDouble(sp[1]);
				line=br.readLine();
				line=br.readLine();
				sp=line.split(regex);	
				this.nInitial[ia]=Integer.parseInt(sp[0]);
				this.nMajor[ia]=Integer.parseInt(sp[1]);
				this.nSymLoop[ia]=Integer.parseInt(sp[2]);
				this.nDescending[ia]=Integer.parseInt(sp[3]);
				this.nAscending[ia]=Integer.parseInt(sp[4]);

				this.nTotCurves[ia]=this.nInitial[ia]+this.nMajor[ia]+this.nSymLoop[ia]+this.nDescending[ia]+this.nAscending[ia];
				BH[ia]=new Mat[this.nTotCurves[ia]];

				line=br.readLine();
				line=br.readLine();
				int L=Integer.parseInt(line);

				BH[ia][0]=new Mat(L,2);

				for( int p=0;p<L;p++){
					line=br.readLine();
					sp=line.split(regex);	
					double[] bh=getCSV(line);
					BH[ia][0].el[p][0]=bh[0];
					BH[ia][0].el[p][1]=bh[1];

				}



				line=br.readLine();
				int L1=0;
				for( int ip=1;ip<this.nTotCurves[ia];ip++){
					line=br.readLine();
					if(line.startsWith("*")) {line=br.readLine();};
					sp=line.split(regex);	
					//if(sp.length==1)
					L1=Integer.parseInt(sp[0]);

					BH[ia][ip]=new Mat(L1,2);

					for( int i=0;i<L1;i++){
						line=br.readLine();

						double[] bh=getCSV(line);

						BH[ia][ip].el[i][0]=bh[0];
						BH[ia][ip].el[i][1]=bh[1];
					}


				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();

				int nLoss=Integer.parseInt(line);
				rotLoss=new Mat(nLoss,2);
				line=br.readLine();

				for( int i=0;i<nLoss;i++){
					line=br.readLine();
					double[] BL=getCSV(line);
					rotLoss.el[i][0]=BL[0];
					rotLoss.el[i][1]=BL[1];
				}

				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				double[] dd=this.getCSV(line);
				
				Lani[ia]=(int)dd[0];
				nAni[ia]=(int)dd[1];
			
				line=br.readLine();

				BHAni[ia]=new Mat[nAni[ia]];

				for( int i=0;i<nAni[ia];i++)
					BHAni[ia][i]=new Mat(Lani[ia],3);

				for( int i=0;i<Lani[ia];i++){
					line=br.readLine();
			
					double[] BL=getCSV(line);
					double B=BL[0];
		
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][0]=B;
						BHAni[ia][j].el[i][1]=BL[j+1];
					}
				}

				line=br.readLine();
				for( int i=0;i<Lani[ia];i++){
					line=br.readLine();
					double[] BL=getCSV(line);
			
					for(int j=0;j<nAni[ia];j++){
						BHAni[ia][j].el[i][2]=BL[j+1];
					}
				}


			}
			
			//util.plotBunch(BH[0]);
			//util.plotBunch(BH[9],3);
			
		//createAngleDepData();

			
			boolean ani=false;
			
			if(ani){
				
				Mat[] BHaniT=new Mat[BHAni[0].length];
				
	
				for(int i=0;i<BHAni[0].length;i++){
					
					BHaniT[i]=new Mat(BHAni[0][i].nRow,2);
					
					double ang=i*Math.PI/18;
					
					Vect er=new Vect(Math.cos(ang),Math.sin(ang));

				for(int j=0;j<BHAni[0][i].nRow;j++){
					double Ht=new Vect(BHAni[0][i].el[j][1], BHAni[0][i].el[j][2]).dot(er);
					BHaniT[i].el[j][0]=Ht;
					BHaniT[i].el[j][1]=BHAni[0][i].el[j][0];
					
				}
			
				}
				

			
			
		//	util.plotBunch(BHaniT,0);

			}
			
		
			String filex="C:\\Works\\PlayModel\\hys_dataNewFormat";
		//	this.writeHystDataColumns(BH, filex);
			
			
			int nSetOut=10;
			Mat[][] BHs=new Mat[nSetOut][BH[0].length];
			
			for(int i=0;i<nSetOut;i++)
				for(int j=0;j<BH[i].length;j++)
					BHs[i][j]=BH[i][j].deepCopy();
			
			PlayModel2D pm=new PlayModel2D();
			
			pm.writeHystData(BHs, filex);
			
			
			util.plotBunch(BHs[0]);	
			util.plotBunch(BHs[4]);
		//	util.plotBunch(BHs[2]);
	/*
			boolean distill=false;
			if(distill){
			Mat[][] BHdist=distHysData(BH);
			
			
			for(int i=0;i<BHdist.length;i++)
				for(int j=0;j<BHdist[i].length;j++)
					for(int k=0;k<BHdist[i][j].nRow;k++)
						BHdist[i][j].el[k][1]*=1+Math.abs(k-BHdist[i][j].nRow/2)*.005;
						
			//util.plotBunch(BHdist[0]);
			
		boolean	average=false;
	
			if(average){

			Mat[] BHsAv=new Mat[BH[0].length];
			
			Mat[] BHsAni=new Mat[0];
	
			//Mat[] BHsrAv=new Mat[BH[0].length];
			for(int i=0;i<BHsAv.length;i++){
				Mat M=BH[0][i].deepCopy();
				for(int j=1;j<nSet;j++)
					M=M.add(BH[j][i]);
		
				BHsAv[i]=M.times(1.0/nSet);

				//BHsrAv[i]=pm.getHBij(BHsAv[i],0);
			}
			
			String fileav="C:\\Works\\HVID\\hys_dataHAvdist";
			this.writeHystDataAv(BHsAv, BHsAni, fileav);

				util.plotBunch(BHsAv);
			}

			
			//Mat[] hysDataAv=new Mat[BH[0].length];
			
		
			if(BHdist.length==1){
				String filed="C:\\Works\\HVID\\hys_dataHAvdist";
			this.writeHystDataAv(BHdist[0], BHAni[0], filed);
			}
			else{
				String filed="C:\\Works\\HVID\\hys_dataHdist";
			this.writeHystData(BHdist, filed);
			}
			}*/
		
			return true;

		}

		catch(IOException e){System.err.println("Error in loading BH data file.");
		return false;
		}

	}	
	


	private double[] getCSV(String line){

		String[] sp=line.split(regex);	

		int p0=0;
		if(sp[0].equals(""))
		{
			p0=1;
		}
		int L=sp.length-p0;

		double[] v=new double[L];

		for( int p=0;p<L;p++){

			v[p]=Double.parseDouble(sp[p+p0]);
		}

		return v;
	}


}
