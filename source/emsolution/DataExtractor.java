package emsolution;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import math.Vect;


public class DataExtractor {
	String regex="[:; ,\\t]+";

	public String[][] loadDataString(String file, String key,int nDataMax,int nTmax,PrintWriter stderr){
		
		return  loadDataString( file,  key,1, nDataMax, nTmax, stderr);
	}
	
	
	public String[][] loadDataString(String file, String key, int nLinesAfter,int nDataMax,int nTmax,PrintWriter stderr){

	
		
		String[][] data1=new String[nDataMax][nTmax];

		try{
			FileReader fr=new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line;
			String s;
			String[] sp;
			int ix=0;
			int sourceIndex=0;
			while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith(key)){}
		
			if(line==null) {break;}
			
			for(int i=0;i<nLinesAfter;i++){
	
				line=br.readLine();
				}
			if(line==null) {break;}
			

			
			sourceIndex=0;


			data1[sourceIndex++][ix]=line;
			line=br.readLine();

			while(line!=null && !line.startsWith("*")){
				sp=line.split(regex);
				if(sp.length>1 && nLinesAfter==4){
					line="";
					for(int k=1;k<sp.length;k++)
						line=line+sp[k]+"  ";

				}
				if(sp.length<2) break;
				data1[sourceIndex][ix]=line;

			sourceIndex++;
		
			line=br.readLine();
			
			}
			
			ix++;
			}
			
	br.close();
	fr.close();

	if(ix==0) return null;

	String[][] data=new String[sourceIndex][ix];

	for(int i=0;i<sourceIndex;i++)
		for(int j=0;j<ix;j++){
		data[i][j]=data1[i][j];


	}


	return data;
		}

		catch(IOException e){stderr.println(e.getMessage());
		return null;
		}

		
	}	

public Vect loadTimesSteps(int[] stepNumb,String file,PrintWriter stderr){


	Vect time1=new Vect(stepNumb.length);

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;
		int ix=0;
		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*                            Step No.   ")){}
			if(line==null) {break;}

			sp=line.split(regex);

			stepNumb[ix]=Integer.parseInt(sp[3]);


			time1.el[ix]=Double.parseDouble(sp[5]);

			ix++;


		}

		br.close();
		fr.close();

		if(ix==0) return null;

		Vect time=new Vect(ix);

		for(int i=0;i<ix;i++){
			time.el[i]=time1.el[i];
		}

		return time;
	}

	catch(IOException e){stderr.println(e.getMessage());
	return null;
	}



}



public int[] getIterNumb(String file,PrintWriter stderr){



	int[] iterNumbs=new int[2];

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("**** Total number of ICCG")){}
			if(line==null) {break;}

			//line=br.readLine();

			//util.pr(line);
			//	line=br.readLine();
			sp=line.split(regex);

			iterNumbs[0]=Integer.parseInt(sp[6]);

			while((line=br.readLine())!=null && !line.startsWith("**** Total number of nonlinear")){}
			if(line==null) {break;}

			sp=line.split(regex);

			iterNumbs[1]=Integer.parseInt(sp[6]);
			//line=br.readLine();
			//sp=line.split(regex);

		}

		br.close();
		fr.close();




	}

	catch(IOException e){stderr.println(e.getMessage());
	}

	return iterNumbs;
}

public int[] getNodeElNumbs(String file,PrintWriter stderr){



	int[] neNumbs=new int[2];

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("                         total number of nodes")){}
			if(line==null) {break;}

			//line=br.readLine();

			//util.pr(line);
			//	line=br.readLine();
			sp=line.split(regex);

			neNumbs[0]=Integer.parseInt(sp[6]);

			while((line=br.readLine())!=null && !line.startsWith("                         total number of volume elements")){}
			if(line==null) {break;}

			sp=line.split(regex);

			neNumbs[1]=Integer.parseInt(sp[7]);
			//line=br.readLine();
			//sp=line.split(regex);

		}

		br.close();
		fr.close();




	}

	catch(IOException e){stderr.println(e.getMessage());
	}

	return neNumbs;
}	

public String[] getComputationTimesAndDate(String file,PrintWriter stderr){



	String[] tt=new String[11];

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*                     CPU Times")){}
			if(line==null) {break;}

			tt[0]=line;
			for(int i=1;i<7;i++){
				line=br.readLine();
				tt[i]=line;
			}

			for(int i=0;i<3;i++){
				line=br.readLine();
			}

			for(int i=8;i<11;i++){
				line=br.readLine();
				tt[i]=line;
				line=br.readLine();
			}

		}

		br.close();
		fr.close();



	}

	catch(IOException e){stderr.println(e.getMessage());
	}
	return tt;


}

public int[] getNodeElNumbs(String file){



	int[] neNumbs=new int[2];

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("                         total number of nodes")){}
			if(line==null) {break;}

			//line=br.readLine();

			//util.pr(line);
			//	line=br.readLine();
			sp=line.split(regex);

			neNumbs[0]=Integer.parseInt(sp[6]);

			while((line=br.readLine())!=null && !line.startsWith("                         total number of volume elements")){}
			if(line==null) {break;}

			sp=line.split(regex);

			neNumbs[1]=Integer.parseInt(sp[7]);
			//line=br.readLine();
			//sp=line.split(regex);

		}

		br.close();
		fr.close();




	}

	catch(IOException e){System.err.println("Error in loading output file.");
	}

	return neNumbs;
}	

public String[] getComputationTimesAndDate(String file){



	String[] tt=new String[12];

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("*                     CPU Times")){}
			if(line==null) {break;}

			tt[0]=line;
			for(int i=1;i<7;i++){
				line=br.readLine();
				tt[i]=line;
			}

			for(int i=0;i<3;i++){
				line=br.readLine();
			}

			for(int i=8;i<11;i++){
				line=br.readLine();
				tt[i]=line;
				line=br.readLine();
			}
			
			line=br.readLine();
			tt[11]=line;

		}

		br.close();
		fr.close();



	}

	catch(IOException e){System.err.println("Error in loading output file.");

	}
	return tt;


}	

public int getRefine(String file){
	int refine=0;


	int L=file.length();

	int ix=0;
	while(file.charAt(L-1-ix)!='\\'){
		ix++;
		if(ix>L-2) break;
	}


	String checkfile="";
	for(int i=0;i<L-ix;i++)
		checkfile=checkfile+file.charAt(i);
	checkfile=checkfile+"check";

	try{
		FileReader fr=new FileReader(checkfile);
		BufferedReader br = new BufferedReader(fr);
		String line;


		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("* INPUT_MESH_FILE *  UNIT  * NO_MESH *")){}
			if(line==null) {break;}
			line=br.readLine();
			String[] sp2=line.split(regex);

			int ib=0;
			if(sp2[0].equals("")) ib=1;

			if(sp2.length-ib>5) refine=Integer.parseInt(sp2[ib+5]);


		}

		br.close();
		fr.close();


	}
	catch(IOException e){};

	return refine;

}

public boolean isRelease(String file){

	boolean rel=false;


	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);

		String line;
		for(int i=0;i<100;i++){
			line=br.readLine();
			if(line==null) {break;}
			if(line.startsWith("      No. of Parallel")){
				rel=true; 
				br.close();
				fr.close();
				break;
				
			}
		}
	

		br.close();
		fr.close();

	}

	catch(IOException e){System.err.println("Error in loading output file.");
	}


			return rel;

}

public int[] getIterNumb(String file){



	int[] iterNumbs=new int[2];

	try{
		FileReader fr=new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		String s;
		String[] sp;

		while((line=br.readLine())!=null){
			while((line=br.readLine())!=null && !line.startsWith("**** Total number of ICCG")){}
			if(line==null) {break;}

			//line=br.readLine();

			//util.pr(line);
			//	line=br.readLine();
			sp=line.split(regex);

			iterNumbs[0]=Integer.parseInt(sp[6]);

			while((line=br.readLine())!=null && !line.startsWith("**** Total number of nonlinear")){}
			if(line==null) {break;}

			sp=line.split(regex);

			iterNumbs[1]=Integer.parseInt(sp[6]);
			//line=br.readLine();
			//sp=line.split(regex);

		}

		br.close();
		fr.close();




	}

	catch(IOException e){System.err.println("Error in loading output file.");
	}

	return iterNumbs;
}


}
