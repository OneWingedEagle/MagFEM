package io;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Scanner;

import javax.swing.JOptionPane;

import materialData.BHCurve;
import materialData.CurrentWaveForm;
import math.Complex;
import math.Mat;
import math.SpMat;
import math.SpVect;
import math.Vect;
import math.util;
import fem.*;
import fem.Network.ElemType;


/**
 * TODO Put here a description of what this class does.
 *
 * @author Hassan.
 *         Created Aug 15, 2012.
 */
public class Loader {

	private String regex="[:; ,=\\t]+";
	public String regex2="[\\[\\]\\s: )(,=\\t]+";

	private String logFile;
	

	public void loadMesh(Model model, String bunFilePath){

		model.meshFilePath=bunFilePath;
		
		
		try{
			FileReader fr=new FileReader(bunFilePath);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;

			String elType=br.readLine();
			model.setElType(elType);

			br.readLine();	
			
			line=br.readLine();
			sp=line.split(regex);
			if(!sp[0].equals(""))		
				model.numberOfNodes=Integer.parseInt(sp[0]);
			else
				model.numberOfNodes=Integer.parseInt(sp[1]);


			br.readLine();		
			line=br.readLine();

			sp=line.split(regex);
			if(!sp[0].equals(""))		
				model.numberOfElements=Integer.parseInt(sp[0]);
			else
				model.numberOfElements=Integer.parseInt(sp[1]);
		
			model.element=new Element[model.numberOfElements+1];
			for(int i=1;i<=model.numberOfElements;i++){
				model.element[i]=new Element(elType);
			}
			
		
			model.node=new Node[model.numberOfNodes+1];

			for(int i=1;i<=model.numberOfNodes;i++)
				model.node[i]=new Node(i, model.dim);

			br.readLine();		
			line=br.readLine();
			sp=line.split(regex);

			int nRegs=0;
			if(!sp[0].equals(""))		
				nRegs=Integer.parseInt(sp[0]);
			else
				nRegs=Integer.parseInt(sp[1]);
			
			if(model.numberOfRegions<nRegs) model.numberOfRegions=nRegs;

			model.region=new Region[model.numberOfRegions+1];
			for(int i=1;i<=model.numberOfRegions;i++)
				model.region[i]=new Region(model.dim);

			br.readLine();

			line=br.readLine();


			model.scaleFactor=Double.parseDouble(line);

			double factor=1.0/model.scaleFactor;

			for(int i=1;i<=model.numberOfElements;i++){
				line=br.readLine();
				sp=line.split(regex);
				int k=0;
				for(int j=0;j<sp.length;j++){
					if(!sp[j].equals(""))		
						model.element[i].setVertNumb(k++, Integer.parseInt(sp[j]));		
				}
			}

			
			Vect z=new Vect(model.dim);
			

			
			for(int i=1;i<=model.numberOfNodes;i++){
				line=br.readLine();
				sp=line.split(regex);
				int k=0;
				
				for(int j=0;j<model.dim;j++)
					if(!sp[j].equals(""))
						z.el[k++]=Double.parseDouble(sp[j])*factor;
				
	
				model.node[i].setCoord(z);

	
				}
			
			
			model.setBounds();
			
		
			
				for(int i=1;i<=nRegs;i++){
					
					line=br.readLine();
					if(line==null) line="1,0,x"+i;
					sp=line.split(regex);
					String[] str=new String[10];
					int k=0;
					for(int j=0;j<sp.length;j++){
						if(!sp[j].equals(""))
							str[k++]=sp[j];
						
					}
					model.region[i].setFirstEl(Integer.parseInt(str[0]));	
					
					model.region[i].setLastEl(Integer.parseInt(str[1]));
					model.region[i].setName(str[2]);
					model.region[i].setMaterial(str[2]);

			}
				
				
			for(int i=nRegs+1;i<=model.numberOfRegions;i++){
	
				model.region[i].setFirstEl(1);	
					
					model.region[i].setLastEl(0);
					model.region[i].setName("xxx");
					model.region[i].setMaterial("xmat");

			}
				
			
				//==============
			
				//=========
				
			System.out.println();
			System.out.println("Loading mesh file completed.");

			br.close();
			fr.close();

			for(int ir=1;ir<=nRegs;ir++)
				for( int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++)
					model.element[i].setRegion(ir);
			
			model.setMaxDim();
			
			model.setFemCalc();

		//	util.pr(model.maxDim);			
				
		}
		catch(IOException e){
			e.printStackTrace();//System.err.println("Error in loading model file.");
		}


	}	

	public void loadData(Model model,String dataFilePath){
		
		int dim0=model.dim;
		
		//logFile=model.logFilePath;


		try{
			BufferedReader br = new BufferedReader(new FileReader(dataFilePath));
			String line;
			line=getNextDataLine(br,"// DATA TYPE (0: Magnetic)");

			int dataType =getIntData(line);
			model.dataType=dataType;
			
			line=getNextDataLine(br,"// DIMENSION (2: 2D, 3: 3D, 4: Axisymmetric 2D)");

			int dim =getIntData(line);
			if(dim==4){
				dim=2;
				model.axiSym=true;
			}

			line=getNextDataLine(br,"// COORDINATE (0: Cartesian, 1: Cylindrical)");

			int coordCode =getIntData(line);
			model.coordCode=coordCode;
			
			if(dim!=dim0){
				System.err.println("Mesh and Data do not match in dimension: "+dim0+" and "+dim);
			}

				setDataMag( model,br);
		
				br.close();
		}

		catch(IOException e){
			e.printStackTrace();
		}

		System.out.println();
		System.out.println("Loading data file completed.");
		
		util.setLogFile(null);

	}
	
	public void setDataMag(Model model,BufferedReader br){
		String line;
		String s;
		int dim=model.dim;

		try {

	
			line=getNextDataLine(br,"// ANALYSIS MODE (0: Magnetostatic, 1:  A-method,  2: A-fi-method ");
			

			int am =getIntData(line);
			model.analysisMode=am;

			line=getNextDataLine(br,"// NONLINEAR (0: Linear , 1: Nonliear ");
			
			boolean nonlin=getBooleanData(line);;
			
			model.setNonLin(nonlin);


			line=getNextDataLine(br,"// AC (1: AC , 0: Time domain) // FREQ (if AC =1) ");
			
			
			String[] sp0=line.split(this.regex);	

			int ib0=0;
			if(sp0[0].equals("")) ib0=1;
			model.AC=getBooleanData(sp0[ib0++]);
			
			if (ib0<sp0.length){
				
			double f0=Double.parseDouble(sp0[ib0++]);
			
			model.setFreq(f0);
			}

			line=getNextDataLine(br,"// NUMBER OF REGIONS ");

			int nRegions =getIntData(line);		
			if(nRegions!=model.numberOfRegions){
				System.out.println("Mesh and Data do not match in the number of regions: "+model.numberOfRegions+" and "+nRegions);
			}
		
			for(int ir=1;ir<=model.numberOfRegions;ir++){
		
			line=getNextDataLine(br,"// *REGION_ID * BH_ID * MU * SIGMA , MAGNETIZATION [Mx, My, Mz] ");

			readAndSetRegMagPropery(model,ir,line);
			}
			

			model.BCtype=new int[model.nBoundary];
			for(int j=0;j<model.nBoundary;j++)
				model.BCtype[j]=-1;
			model.PBCpair=new int[model.nBoundary];

			Vect B;
			int[] bcData=new int[2];
			
			for(int j=0;j<model.nBoundary;j++){

				line=getNextDataLine(br,"// * BOUNDRAY CONDITION (D: Drichlet, N: Neumann) *");

				if(model.BCtype[j]>-1) continue;
				bcData=getBCdata(line);
				model.BCtype[j]=bcData[0];
				
				model.PBCpair[j]=bcData[1];
				if(model.BCtype[j]>1){
				
					model.BCtype[model.PBCpair[j]]=bcData[0];
					model.PBCpair[model.PBCpair[j]]=j;

				 model.hasPBC=true;

				}

								
			}
			
			line=getNextDataLine(br,"// *NUMBER OF GIVEN CURRENT DENSITY * ");

			int numbRegsWithJ=Integer.parseInt(line);
			for(int j=0;j<numbRegsWithJ;j++){
				if(model.axiSym)
					line=getNextDataLine(br,"// * COIL_ID * 0. *  0. * Jy *");
				else
					line=getNextDataLine(br,"// * COIL_ID * Jx *  Jy * Jz *");
		
				String[] sp=line.split(this.regex);	

				int ib=0;
				if(sp[0].equals("")) ib=1;
				int nr=Integer.parseInt(sp[ib++]);
				Vect J=new Vect(3);
				
				for(int k=0;k<3;k++)
					J.el[k]=Double.parseDouble(sp[ib++]);
				
				model.region[nr].setJ(J);

				
			}

			line=getNextDataLine(br,"//NUMBER OF COILS ");
			
			int numCoils=Integer.parseInt(line);
			if(numCoils<0){
				numCoils=-numCoils;
				model.seperateCoil=false;
			}
			
			if(numCoils>0){
				model.hasJ =true;

				model.phiCoils=new PhiCoil[numCoils];
				model.coilInicesByRegion=new int[model.numberOfRegions+1];
				for(int i=1;i<=model.numberOfRegions;i++)
					model.coilInicesByRegion[i]=-1;

				
				for(int j=0;j<numCoils;j++){
					line=getNextDataLine(br,"// * REGION_ID *  TURNS * SIGMA *");
					String[] sp=line.split(this.regex);	

					int ib=0;
					if(sp[0].equals("")) ib=1;
					int nr=Integer.parseInt(sp[ib++]);
					model.phiCoils[j]=new PhiCoil(nr);
					model.phiCoils[j].index=j;
					model.coilInicesByRegion[nr]=j;
					
					double turns=Double.parseDouble(sp[ib++]);
					model.phiCoils[j].setNumTurns(turns);
					double sigma=0;
					
					if(ib<=sp.length){
					 sigma=Double.parseDouble(sp[ib++]);
					}
						
					model.phiCoils[j].setSigma(sigma);
					model.region[nr].setSigma(new Vect(3));
					
					double[][] boxdata=new double[2][6];
					int[] boxCoordType=new int[2];
					
					for(int k=0;k<2;k++){
					if(k==0)
						line=getNextDataLine(br,"//BOX OF COIL INPUT FACE NODES ");
					else
						line=getNextDataLine(br,"//BOX OF COIL OUTPUT FACE NODES ");

					ib=0;
					if(sp[0].equals("")) ib=1;
					sp=line.split(this.regex);	

					String ss=sp[ib++];
					if(ss.equals("x") || ss.equals("r"))
					{
						boxdata[k][0]=Double.parseDouble(sp[ib++]);
						boxdata[k][1]=Double.parseDouble(sp[ib++]);
						if ( ss.equals("r")) boxCoordType[k]=1;
						
					}
					ss=sp[ib++];
					if(ss.equals("y") || ss.equals("t"))
					{
						boxdata[k][2]=Double.parseDouble(sp[ib++]);
						boxdata[k][3]=Double.parseDouble(sp[ib++]);
						
					}
					ss=sp[ib++];
					if(ss.equals("z"))
					{
						boxdata[k][4]=Double.parseDouble(sp[ib++]);
						boxdata[k][5]=Double.parseDouble(sp[ib++]);
					}
					model.phiCoils[j].faceBox[k]=boxdata[k];
					model.phiCoils[j].faceCoordType[k]=boxCoordType[k];
					}
					
				}	
	
			}

			line=getNextDataLine(br,"// UNIFORM FIELD TIME_ID ");

			model.unifBTimeId=getIntData(line);

			if(model.unifBTimeId>0){
				model.hasBunif=true;

				line=getNextDataLine(br,"// Bx  By 0 ");

				double[] array=getCSV(line);
				
				model.unifB=new Vect(array);
				if(model.unifB.length==3 && model.unifB.el[2]!=0){
					System.out.println("!!!!!!!!!!! Uniform field in Z direction not available. !!!!!!!!!!");
					System.out.println("!!!!!!!!!!! Calculation abandoned. !!!!!!!!!!");
					wait(10000*10000);

				}
			}
			line=getNextDataLine(br,"//NETWORK (CIRCUT) ");

			if(line.equals("NETWORK")){
				Network network=new Network();
				network.read(this, br);

				model.network=network;
				
			for(int j=0;j<network.numElements;j++){

					if(network.elems[j].type==ElemType.FEM){

						network.elems[j].fem_index=network.elems[j].fem_id-1;
					
						
					}
				}
			}
			
			line=getNextDataLine(br,"// NUM TIME FUNCTIONS");

			int numTimeFuncs=getIntData(line);	
			if(numTimeFuncs>0){
				
			model.timeFunctions=new TimeFunction[numTimeFuncs+1];

				
			for(int j=0;j<numTimeFuncs;j++){

				line=getNextDataLine(br,"// TIME ID // TYPE");

				String[] sp=line.split(regex);
				int ib=0;
				if(sp[ib].equals("")) ib++;
				int id=Integer.parseInt(sp[ib++]);
				int type=0;
				if(ib<sp.length)
				 type=Integer.parseInt(sp[ib++]);
				if(type==0){

				line=getNextDataLine(br,"// AMPLITUDE // PERIOD // PHASE");
			
				
				sp=line.split(regex);
				ib=0;
				if(sp[ib].equals("")) ib++;
				
				double amp= Double.parseDouble(sp[ib++]);
				
				double per= Double.parseDouble(sp[ib++]);
				
				double phase=0;
				if(ib<sp.length)
				 phase= Double.parseDouble(sp[ib++]);
				model.timeFunctions[id]=new TimeFunction(id,amp,per,phase);

				}else if(type==1){
					line=getNextDataLine(br,"// C_DC * C_RAMP * PERIOD * C_COS * C_SIN * T_EXP * //(f(t)= C_DC+exp(T_EXP*t)(C_COS(2*PI*t/PERIOD)+C_SIN(2*PI*t/PERIOD)");
					
					sp=line.split(regex);
					ib=0;
					if(sp[ib].equals("")) ib++;

					double a0= Double.parseDouble(sp[ib++]);
					
					double a1=0;
					if(ib<sp.length)
						a1= Double.parseDouble(sp[ib++]);
					
					double period=1;
					if(ib<sp.length)
						period= Double.parseDouble(sp[ib++]);
					double a2=0;
					if(ib<sp.length)
						a2= Double.parseDouble(sp[ib++]);
					
					double a3=0;
					if(ib<sp.length)
						a3= Double.parseDouble(sp[ib++]);
					double a4=0;
					if(ib<sp.length)
					a4= Double.parseDouble(sp[ib++]);
					model.timeFunctions[id]=new TimeFunction(id,a0,a1,period,a2,a3,a4);
				}
				
			
			}
			}


			line=getNextDataLine(br,"//DELTA_TIME");
			{
				String sp[]=line.split(regex);
				int ib=0;
				if(sp[ib].equals("")) ib++;
				model.dt=Double.parseDouble(sp[ib++]);	
				if(ib<sp.length)
					model.initialTime=Double.parseDouble(sp[ib++]);	
				else
					model.initialTime=0;

			}

			line=getNextDataLine(br,"//STEP_BEGIN  *  STEP_END * INTERVAL");
			if(line!=null){
				int nSteps=1;
				String sp[]=line.split(regex);
				int L=sp.length;
				
				int n1=0,n2=0, d=1;
		
					n1=Integer.parseInt(sp[L-3]);
					n2=Integer.parseInt(sp[L-2]);
					d=Integer.parseInt(sp[L-1]);
					if(d!=0)
					 nSteps=(n2-n1)/d+1;
			
			
			
					model.setnTsteps(nSteps);
					model.nBegin=n1;
					model.nEnd=n2;
					model.nInc=d;



				}


			line=getNextDataLine(br,"//SAVE_FLUX * SAVE_CURRENT");

			String sp[]=line.split(regex);
			int ib=0;
			if(sp[ib].equals("")) ib++;
	
			 model.saveFlux=false;
			 model.saveJe=false;
			 if(Integer.parseInt(sp[ib++])==1)  model.saveFlux=true;
			 if(ib<sp.length)
			 if(Integer.parseInt(sp[ib++])==1)  model.saveJe=true;


			line=getNextDataLine(br,"//NUMBER OF BH_DATA");

			int numBHdata=getIntData(line);	
			

			for(int j=0;j<numBHdata;j++){

				line=getNextDataLine(br,"// * H * B* ");

				sp=line.split(regex);
				ib=0;
				if(sp[ib].equals("")) ib++;
				int bhID=Integer.parseInt(sp[ib]);
				for(int ir=1;ir<=model.numberOfRegions;ir++)
					if(model.region[ir].BHnumber==bhID)
						model.region[ir].setMaterial(sp[ib+1]);
				
				
			}
				


	if(model.axiSym) model.height=2*Math.PI;
		
			model.setHasJ();

			model.setHasM();

			model.setNonLinToElememts();
			
			model.setEdge();


			model.setElementsParam();

			
			model.setBounds();

		
		if(model.coordCode==1) {
			
			model.cpb=1;

			for(int j=0;j<model.nBoundary;j++){

				if(model.BCtype[j]==3) model.cpb=-1;
			}

		}
		
		//=====================


		model.setNodeOnBound();
		
		if(model.hasPBC) model.mapPBC();
				
			for(int ir=1;ir<=model.numberOfRegions;ir++){
				if(model.region[ir].circuit ) {
					model.circuit=true;
					break;
				}
				
			}

				
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		model.network.constructNetwork();

		
		//model.setForceCalc();

}

	
		
	public void loadBH(BHCurve BH,String mateName) throws Exception{
		double[][] BH1=new double[200][2];

		String file = System.getProperty("user.dir") + "\\BH\\"+mateName+".txt";

		try{
			Scanner scr=new Scanner(new FileReader(file));
			while(scr.hasNext()){

				while(!scr.next().equals("begin")){}
				int j=0;
				String s=scr.next();
				while(!s.equals("end")){
					BH1[j][0]=Double.parseDouble(s);
					s=scr.next();
					BH1[j++][1]=Double.parseDouble(s);
					s=scr.next();
				}
				BH.length=j;
			}

			BH.BH=new double[BH.length][2];
			for(int k=0;k<BH.length;k++)
				BH.BH[k]=BH1[k];

			scr.close();
		}
		catch(IOException fnf){
			throw new Exception(fnf);
		}
		BH.setGradBH();

	}

	private int[] getBCdata(String line){
	
		
		int[] bctp=new int[2];
		String[] sp=line.split(regex);	
		int pair=-1,bct=0;
	
		int ib=0;
		if(sp[0].equals("")) ib++;
		if(line.startsWith("N")){
				bct=0;
				
		}
		else if(line.startsWith("D")){
			bct=1;
			
	}
		else if(line.startsWith("X")){
			bct=-1;
			
	}
		else if(line.startsWith("PS")){
			bct=2;
			pair=Integer.parseInt(sp[ib+1])-1;
			
	}

		else if(line.startsWith("PA")){

				bct=3;
				pair=Integer.parseInt(sp[ib+1])-1;
			}
		else{

			bct=0;
		}
		
		bctp[0]=bct;
		bctp[1]=pair;
		
		return bctp;
	}


	private double getScalarData(String line){
		String[] sp=line.split(regex);	
		return Double.parseDouble(sp[sp.length-1]);
	}

	private int getIntData(String line){
		String[] sp=line.split(regex);	
		return Integer.parseInt(sp[sp.length-1]);
	}

	private boolean getBooleanData(String line){
		boolean b=false;
		String[] sp=line.split(regex);	
		
		if(sp[sp.length-1].startsWith("t") || sp[sp.length-1].equals("1"))	
			b=true;
		else 	if(sp[sp.length-1].startsWith("f") || sp[sp.length-1].equals("0"))	
			b=false;
		else {
	
				System.out.println("!!!!!!!!!!! Bad input! Job sopps!!!!!!!!!!");
	
				wait(10000*10000);

			}		
		return b;

	}

	private void readAndSetRegMagPropery(Model model,int i,String line){
		int is=0;
		String[] sp=line.split(this.regex);	

		int ib=0;
		if(sp[0].equals("")) ib=1;
		int ir=Integer.parseInt(sp[ib++]);
	
		int BH_id=Integer.parseInt(sp[ib++]);
	
		model.region[ir].BHnumber=BH_id;
		
		Vect v=new Vect(model.dim);
		
		double mu=Double.parseDouble(sp[ib++]);
		
		for(int k=0;k<v.length;k++){
			v.el[k]=mu;
		}
	
	//	if(model.region[ir].BHnumber==0)
		model.region[ir].setMur(v);
		
		double sigma=Double.parseDouble(sp[ib++]);
		Vect v3=new Vect(3);
		for(int k=0;k<v3.length;k++){
			v3.el[k]=sigma;
		}

		model.region[ir].setSigma(v3);

		if(ib<sp.length){
			Vect M=new Vect(model.dim);
			for(int k=0;k<v.length;k++){
				M.el[k]=Double.parseDouble(sp[ib++]);;
			}
			model.region[ir].setM(M);
		}
		
		//=========== 
		if(model.region[ir].stranded) {


			model.stranded=true;
			
			if(model.dim==2)
			model.region[ir].windingSurf=model.getRegionArea(ir);
			else
			{
				model.region[ir].nloop=93;
			//	model.region[ir].windingSurf=model.getRegionAreaXY(ir);
				model.region[ir].windingSurf=model.getRegionVolume(ir)/.05;
			}
				
			/*
			model.region[ir].terminalVoltage0=Double.parseDouble(sp[is++]);
			model.region[ir].setFreq(Double.parseDouble(sp[is++]));

			model.region[ir].phase0=Double.parseDouble(sp[is++])*Math.PI/180;
			model.region[ir].setWireRes(Double.parseDouble(sp[is++]));

			if(is<sp.length)
				model.region[ir].nloop=Double.parseDouble(sp[is++]);
				else
					model.region[ir].nloop=100;

			
			if(is<sp.length)
				model.region[ir].circuit=sp[is++].startsWith("t");
			
			if(model.region[ir].circuit){

				if(is<sp.length)
				model.region[ir].curMap1=Integer.parseInt(sp[is++]);
			
			if(is<sp.length)
				model.region[ir].currCoef1=Double.parseDouble(sp[is++]);
			else
				model.region[ir].currCoef1=1;
			}*/
			
			model.region[ir].NtS=model.region[ir].nloop/model.region[ir].windingSurf;
		//	model.region[ir].NtS=462357.4142194;
 
			

		}

		
	}
	


public double[] loadArray(){
	String file=util.getFile();
	if(file==null || file.equals("") )  throw new NullPointerException("file not found.");
	return loadArray(file);
}

public double[] loadArray(String arrayPath){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;

		int N=100000;
		
		double[] x1=new double[N];
		
		int i=0;
		line=br.readLine();
		while(line!=null){
			if(i>N) break;
			x1[i++]=Double.parseDouble(line);
			line=br.readLine();
			
		}

		double[] x=Arrays.copyOf(x1, i);
		
	br.close();
	fr.close();
			return x;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}	


public Mat loadMatSymm(int n,String arrayPath){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;
	
		Mat A=new Mat(n,n);
		
		for(int i=0;i<n;i++){
			line=br.readLine();
			if(line==null) continue;
			double[] x=getCSV(line);
			for(int j=0;j<=i;j++){
				A.el[i][j]=x[j];	
				if(i!=j)
					A.el[j][i]=x[j];
			}
	
		}

		br.close();
		fr.close();
		
		return A;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}



public double[][] loadArrays(int n, int m,String arrayPath){
return loadArrays( n,  m, arrayPath,  0);
}

public double[][] loadArrays(int n, int m,String arrayPath, int skip){

	try{
		FileReader fr=new FileReader(arrayPath);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		for(int i=0;i<skip;i++){
			line=br.readLine();
		}
		
		double[][] A=new double[n][m];
		
		for(int i=0;i<n;i++){
			line=br.readLine();
			if(line==null) continue;
			double[] x=getCSV(line);
			for(int j=0;j<m;j++)
				A[i][j]=x[j];
	
			
			
		}

		br.close();
		fr.close();
			return A;
			
	}
	catch(IOException e){
		e.printStackTrace();//System.err.println("Error in loading model file.");
	}


	return null;
}


public double[] getCSV(String line){
	
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

public int[] getCSInt(String line){
	String[] sp=line.split(regex);	
	int L=sp.length;
	int[] v=new int[L];
	for( int p=0;p<L;p++)
				v[p]=Integer.parseInt(sp[p]);

	return v;
}

public String getNextDataLine(BufferedReader br) throws IOException{
	String line="";
	while(true){
		line=br.readLine();
		if(line==null) break;
		if(!line.startsWith("/")) break;
	}
return line;
}


public String getNextDataLine(BufferedReader br,String title) throws IOException{
	String line="";
	//try{
		//PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(logFile,true)));
	util.pr(title);
	//pw.println(title);
	
	while(true){
		line=br.readLine();
		if(line==null) break;
		if(!line.startsWith("/")) break;
	}
	util.pr(line);
	

	//	pw.println(line);
	//	pw.close();
	//} catch(IOException e){System.err.println("IOException: " + e.getMessage());}

return line;
}
public double outputLoss(Model model,String file,int step,double phase_time,boolean append){
	

	double totalLoss=0;
	try{
		PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(file,append)));
		pw.println("Step: "+step+"\t phase/time: "+phase_time);
		pw.println();
		util.pr("step "+step);
		pw.println(" ****** ============= ******** ==================== ******");
			pw.println("Joule Losses [W]");
			pw.println();
			util.pr("Joule Losses [W]");
			for(int ir=1;ir<=model.numberOfRegions;ir++){
				if((model.region[ir].isConductor && model.analysisMode>0) || model.coilInicesByRegion[ir]>=0){
			double  loss=0;
			
			loss=model.obtainLoss(ir);
			
			totalLoss+=loss;
			pw.format("%5d %12.5e\n",ir,loss);
			pw.println();
			System.out.format("%5d %25.12e\n",ir,loss);

				}
			}
			pw.format("total %12.5e",totalLoss);
			pw.println();
	
		pw.close();
	} catch(IOException e){System.err.println("IOException: " + e.getMessage());}

	
	//System.out.println(" Temperature was written to "+file);

	
	return totalLoss;

}


public void wait(int ms){
	try {
		Thread.sleep(ms);
	} catch (InterruptedException e) {
		e.printStackTrace();
	}
}



}
