package emsolution;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import io.Loader;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Frame;
import java.awt.Toolkit;


import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.Files;
import java.text.DecimalFormat;

import javax.swing.JFrame;
import javax.swing.JPanel;

import materialData.CurrentWaveForm;
import math.Complex;
import math.Mat;
import math.Vect;
import math.util;

import org.math.plot.Plot2DPanel;




public class LoopRunRelease {

	DecimalFormat df=new DecimalFormat("0.00E00");
	DataExtractor dex=new DataExtractor();

	String regex="[:; . ,\\t]+";
	String regex_no_point="[:;  ,\\t]+";

	boolean[] success;
	
	boolean inputsReady=false;
	
	boolean twoMeshes=false;

	String[][] inputMesh=new String[2][2];
	static int geom=1;


	

	String caseFolder,execDir;


	public  LoopRunRelease(){

	};


	public static void main(String[] args) throws IOException{

		//Loader loader=new Loader();
		
		LoopRunRelease x=new LoopRunRelease();
		x.sweepComplexPerm();
		//.x.sweepForFineMesh();
		//x.sweepACwithRelease();

	}

	public void sweepComplexPerm() throws IOException{
		
	//	LoopRunRelease x=new LoopRunRelease();

		Loader loader=new Loader();

		//x.execDir="C:/Works/2017 works/Homogenization/one-trun/actual-100coils/self_current/axisym/homogenized";
				//execDir="C:/Works/2017 works/Homogenization/LitzCell/cell-36wires-2018/homogenizedModel/axisym";
				//execDir="C:/Works/2017 works/Homogenization/LitsCell/cell-36wires-2018/homogenizedModel/axisym/FEM3x3/homogenized";
				//execDir="C:/Works/2017 works/Homogenization/LitsCell/cell-36wires-2018/homogenizedModel/axisym/singleWire/homogenized";
				
				//execDir="C:/Works/2018 Works/homogenization-by-CLN/Litz Wire 6x6/homogenizedModel/axisym";
				//execDir="C:/Works/2018 Works/homogenization-by-CLN/5x10LitzWire/experimentalCoilCalculations/homogenizedModel/ModifieldCoilDimension";
			//	execDir="C:/Works/2018 Works/homogenization-by-CLN/5x10LitzWire/experimentalCoilCalculations/homogenizedModel/RougherMesh/veryRough";
				execDir="C:/Works/2018 Works/homogenization-by-CLN/by_Biot_Savart/IsolatedCylinder/2D/d=1mm/Homogenization of Single Wire/FEM-CLN/Homogenized model/rout=.5mm";
				execDir="C:/Works/2018 Works/homogenization-by-CLN/by_Biot_Savart/IsolatedCylinder-office/2D/d=1mm/Homogenization of Single Wire/FEM-CLN/Homogenized model/3x3-wires-close/Homogenized";

				String folder="C:/Users/Hassan Ebrahimi/Desktop";
				
				String mur_vs_freq=execDir+"/mur_vs_freq";
			
				Mat mur_c_all=new Mat(loader.loadArrays(17, 3, mur_vs_freq,2));
				//Mat mur_c_all=new Mat(loader.loadArrays(25, 3, mur_vs_freq));
				
				
				CurrentWaveForm mur_real_all= new CurrentWaveForm(mur_c_all.getColVect(0),mur_c_all.getColVect(1));
			
				CurrentWaveForm mur_imag_all= new CurrentWaveForm(mur_c_all.getColVect(0),mur_c_all.getColVect(2));
				
				mur_real_all.setPeriodic(false);
				mur_imag_all.setPeriodic(false);

			//	DecimalFormat df=new DecimalFormat("0x0E00");

				File fRef=new File(execDir+"/inputRef");
				File fInput1=new File(execDir+"/input1");
				File fInput=new File(execDir+"/input");
				
				boolean yasuda=false;
				
				int nf=17;
				
				if(yasuda)
					nf=10;
				
				double fStart=1e2;
				double fEend=1e6;
				double factor = 1;
				if (nf > 1)
					factor = Math.pow(fEend/ fStart, 1. / (nf - 1));
				
				Mat mur_c=new Mat(nf,3);
				
				Vect freqs=new Vect(nf);
				double freq=fStart;
				for(int i=0;i<nf;i++){
					freqs.el[i]=freq;
					freq*=factor;
				}
				
				//===== 
				
				if(yasuda){
				freqs= new Vect(nf);
				
				int ix=0;
				freqs.el[ix++]=1e0;
				freqs.el[ix++]=1e3;
				freqs.el[ix++]=1e4;
				freqs.el[ix++]=5e4;
				freqs.el[ix++]=8.5e4;
				freqs.el[ix++]=1e5;
				freqs.el[ix++]=2e5;
				freqs.el[ix++]=3e5;
				freqs.el[ix++]=4e5;
				freqs.el[ix++]=5e5;
				}
				//======
				
				boolean showMurOnly =false;
				//freqs.hshow();

				for(int k=0;k<freqs.length;k++){
			
			
				 freq=freqs.el[k];
				mur_c.el[k][0]=freq;
				mur_c.el[k][1]=mur_real_all.getI(freq);
				mur_c.el[k][2]=mur_imag_all.getI(freq);
			
				if(!showMurOnly){
				setFrequency(fRef,fInput1,freq);

				setComplexPermeability(fInput1,fInput,freq,mur_c.el[k][1],-mur_c.el[k][2]);

					runRelease(folder,"");
					
					
				//	try{
					util.wait(500);
						File fout=new File(execDir+"/output");
						
						double[] impedance=getImpedance(fout.getPath());
						double[] freq_impedance=new double[3];
						freq_impedance[0]=freq;
						freq_impedance[1]=impedance[0];
						freq_impedance[2]=impedance[1];
				
						util.hshow(freq_impedance);
						
					

				}
				else
				util.hshow(mur_c.el[k]);
						//String st=x.df.format(freq);
						//st.SetharAt(1,'p');
						//File fout_x=new File(x.execDir+"/output"+x.df.format(freq));
					//	util.copyFile(fout, fout_x);

					//} catch (IOException e1) {
						// TODO Auto-generated catch block
					//	e1.printStackTrace();
					//}
				}
					
	}
	
	public void sweepACwithRelease() throws IOException{
		
				double volume=1.139*1.139*18*1e-9;
				double Bave=1;
			//	execDir="C:/Works/2017 works/Homogenization/LitsCell/cell-36wires-2018/homogenizedModel/axisym/FEM";
				//execDir="C:/Works/2017 works/Homogenization/LitsCell/cell-36wires-2018/homogenizedModel/axisym/FEM5wiresInRow";
				//execDir="C:/Works/2018 Works/homogenization-by-CLN/one-turn/diam0.8/AC";
			//	execDir="C:/Works/2018 Works/homogenization-by-CLN/5x10LitzWire/NewMeshWithQuads/3D-AC_SLIDE";
				execDir="C:/Works/2018 Works/homogenization-by-CLN/by_Biot_Savart/IsolatedCylinder-home/2D/d=1mm/Homogenization of Single Wire/FEM-CLN/Homogenized model/rout=.52mm";

				String folder="C:/Users/Hassan Ebrahimi/Desktop";
				
				



				File fRef=new File(execDir+"/inputRef");
				File fInput1=new File(execDir+"/input1");
				File fInput=new File(execDir+"/input");
				
				int nf=5;
				double fStart=1e2;
				double fEend=1e6;
				double factor = 1;
				if (nf > 1)
					factor = Math.pow(fEend/ fStart, 1. / (nf - 1));
				
				Mat mur_c=new Mat(nf,3);
				
				Vect freqs=new Vect(nf);
				double freq=fStart;
				for(int i=0;i<nf;i++){
					freqs.el[i]=freq;
					freq*=factor;
				}


				for(int k=0;k<freqs.length;k++){
			
			
				 freq=freqs.el[k];
				mur_c.el[k][0]=freq;
	
			
				setFrequency(fRef,fInput1,freq);
				setTimeEvolutionData(fInput1,fInput,1.,1./freq,0);

				runRelease(folder,"");
					
					
				//	try{
					util.wait(500);
						File fout=new File(execDir+"/output");
						
						double W=getW_AC_Release(fout.getPath());
						double Q=getQ_AC_Release(fout.getPath());
						
					//	util.pr("W: "+W+"    Q: "+Q);
						//Complex mur=Bave*volume
						double omega=2*PI*freq;
								Complex denum=new Complex(W,Q/omega);
								Complex num=new Complex(Bave*volume/2/(PI*4e-7),0);
								Complex mur=num.times(denum.inv());
					//	COMPLEX mur = Bav2Vol / P_jw / PMU0 / 2.;

					//	double[] mur=getComplexPermeabilityRelease(fout.getPath());
						double[] freq_mur=new double[3];
						freq_mur[0]=freq;
						freq_mur[1]=mur.re;
						freq_mur[2]=mur.im;
				
						
						freq_mur[1]=Q;
						freq_mur[2]=W;
						util.hshow(freq_mur);
						
					

				}
			
				
					
	}

	public  boolean runRelease(String sourceFolder, String destFolder){


		boolean success=false;



		ProcessBuilder builder = new ProcessBuilder();

		//builder.command("EMSolBatchCversion_x64.exe","input");
		
		
		builder.command(sourceFolder+"/EMSolutionComplex.exe","-b","-m","-f","input");

		builder.directory(new File(execDir));

		


		Process pp = null;
		try {
			pp = builder.start();
		} catch (IOException e) {
			e.printStackTrace();
		}
		// Process has started here
		try {
			pp.waitFor();
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		try {
			// i = 143
			int i = pp.exitValue();
		} catch( IllegalThreadStateException e){
			e.printStackTrace();
		}

		pp.destroy();


		return success;


	}						

	//************************************************



	public   String getFile(int mode,String message){
		String filePath="";
		FileDialog fd;
		Frame f=new Frame();
		if(mode==0)
			fd= new FileDialog(f,message,FileDialog.LOAD);
		else
			fd= new FileDialog(f,message,FileDialog.SAVE);
		fd.setVisible(true);
		fd.toFront();
		String Folder=fd.getDirectory();
		String File = fd.getFile();
		if(Folder!=null && File!=null)
		{

			filePath=Folder+"\\"+File;

		}
		f.dispose();
		fd.dispose();

		return filePath;
	}



private void setComplexPermeability(File source, File dest,double freq, double mur_real,double mur_complex ){

	try{
		FileReader fr=new FileReader(source);
		BufferedReader br = new BufferedReader(fr);

		FileWriter fr2=new FileWriter(dest);
		PrintWriter pw = new PrintWriter(fr2);

		String line="";
		String line1="";

		int ix=0;
		boolean manip=false;
		while(line!=null){
			line=br.readLine();

			if(line==null) break;

			line1=util.dropLeadingSpaces(line);

			if(manip && !line1.startsWith("*")){

				pw.println("* "+df.format(freq));

				line=mur_real+"\t"+mur_complex;

				manip=false;

			}
			pw.println(line);


			if(line.startsWith("* MU_Re *")){
				manip=true;
			}

		}

		fr2.close();

		pw.close();

		br.close();
		fr.close();


	}

	catch(IOException e2){e2.printStackTrace();}
}





private void setFrequency(File source, File dest,double freq){

	try{
		FileReader fr=new FileReader(source);
		BufferedReader br = new BufferedReader(fr);

		FileWriter fr2=new FileWriter(dest);
		PrintWriter pw = new PrintWriter(fr2);

		String line="";
		String line1="";

		int ix=0;
		boolean manip=false;
		while(line!=null){
			line=br.readLine();

			if(line==null) break;

			line1=util.dropLeadingSpaces(line);

			if(manip && !line1.startsWith("*")){

				line=df.format(freq);



				manip=false;

			}

			pw.println(line);


			if(line.startsWith("* FREQUENCY *")){
				manip=true;
			}

		}

		fr2.close();

		pw.close();

		br.close();
		fr.close();


	}

	catch(IOException e2){e2.printStackTrace();}
}

private void setTimeEvolutionData(File source, File dest, double amp, double period,double phase){

	try{
		FileReader fr=new FileReader(source);
		BufferedReader br = new BufferedReader(fr);

		FileWriter fr2=new FileWriter(dest);
		PrintWriter pw = new PrintWriter(fr2);

		String line="";
		String line1="";

		int ps_ind=0;
		int ix=0;
		boolean wasSet=false;
		boolean manip=false;
		while(line!=null){
			line=br.readLine();

			if(line==null) break;

			line1=util.dropLeadingSpaces(line);

			if(manip && !line1.startsWith("*")){
				//line=br.readLine();
				//util.pr(line);
				//pw.println(line);
				line=df.format(amp)+"   "+df.format(period)+"   "+df.format(phase);
				wasSet=true;
				//break;
				manip=false;

			}

			pw.println(line);


			if(!wasSet && line.startsWith("* TIME_ID * OPTION *")){
				line=br.readLine();
				//util.pr(line);
				pw.println(line);
					manip=true;
			}

		}

		fr2.close();

		pw.close();

		br.close();
		fr.close();


	}

	catch(IOException e2){e2.printStackTrace();}
}



public double[] getImpedance(String file){



	double[] impedance=new double[2];

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("***  Power Sources   ********")){}
			if(line==null) {break;}
			line=br.readLine();
			line=br.readLine();
			sp=line.split(regex_no_point);
	
			int ib=0;
			if(sp[0]=="") ib++;
			double Zimag=Double.parseDouble(sp[3+ib]);
			
			impedance[1]=Zimag;
	
			
			while((line=br.readLine())!=null && !line.startsWith("***  Power Sources   ********")){}
			
			if(line==null) {break;}
			line=br.readLine();
			line=br.readLine();
			sp=line.split(regex_no_point);
			 ib=0;
			if(sp[0]=="") ib++;
			double Zreal=Double.parseDouble(sp[3+ib]);
			
			impedance[0]=Zreal;
	
			break;
		}

		br.close();
		fr.close();



	}

	catch(IOException e){System.err.println("Error in loading output file.");

	}
	return impedance;


}


public double[] getImpedanceNetwork(String file){



	double[] impedance=new double[2];

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("***  Network elements ")){}
			if(line==null) {break;}
			line=br.readLine();
			line=br.readLine();
			sp=line.split(regex_no_point);
	
			int ib=0;
			if(sp[0]=="") ib++;
			double Zimag=Double.parseDouble(sp[3+ib]);
			
			impedance[1]=Zimag;
	
			
			while((line=br.readLine())!=null && !line.startsWith("***  Network elements ")){}
			
			if(line==null) {break;}
			line=br.readLine();
			line=br.readLine();
			sp=line.split(regex_no_point);
			 ib=0;
			if(sp[0]=="") ib++;
			double Zreal=Double.parseDouble(sp[3+ib]);
			
			impedance[0]=Zreal;
	
			break;
		}

		br.close();
		fr.close();



	}

	catch(IOException e){System.err.println("Error in loading output file.");

	}
	return impedance;


}

public double getQ_AC_Release(String file){



	double Q=0;

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** Total Joule")){}
			if(line==null) {break;}
			
			while((line=br.readLine())!=null && !line.startsWith("        Total")){}
			//line=br.readLine();
			//line=br.readLine();
			sp=line.split(regex_no_point);
	
			int ib=0;
			if(sp[0]=="") ib++;
			double Qimag=Double.parseDouble(sp[2+ib]);
			
			Q+=.5*Qimag;
	
			
			while((line=br.readLine())!=null && !line.startsWith("*** Total Joule")){}
			if(line==null) {break;}
			
			while((line=br.readLine())!=null && !line.startsWith("        Total")){}
			//line=br.readLine();
			//line=br.readLine();
			sp=line.split(regex_no_point);
	
			 ib=0;
			if(sp[0]=="") ib++;
			double Qreal=Double.parseDouble(sp[2+ib]);
			
			Q+=.5*Qreal;
				
			break;
		}

		br.close();
		fr.close();



	}

	catch(IOException e){System.err.println("Error in loading output file.");

	}
	return Q;


}

public double getW_AC_Release(String file){



	double W=0;

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*** Total magnetic en")){}
			if(line==null) {break;}
			
			while((line=br.readLine())!=null && !line.startsWith("        Total")){}
			//line=br.readLine();
			//line=br.readLine();
			sp=line.split(regex_no_point);
	
			int ib=0;
			if(sp[0]=="") ib++;

			double Wimag=Double.parseDouble(sp[2+ib]);
			W+=Wimag;
	
			
			while((line=br.readLine())!=null && !line.startsWith("*** Total magnetic en")){}
			if(line==null) {break;}
			
			while((line=br.readLine())!=null && !line.startsWith("        Total")){}
			//line=br.readLine();
			//line=br.readLine();
			sp=line.split(regex_no_point);
	
			 ib=0;
			if(sp[0]=="") ib++;
			double Wreal=Double.parseDouble(sp[2+ib]);

			W+=Wreal;
				
			break;
		}

		br.close();
		fr.close();



	}

	catch(IOException e){System.err.println("Error in loading output file.");

	}
	return W;


}

}

