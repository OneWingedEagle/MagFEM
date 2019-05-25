
package emsolution;

import java.io.*;
import java.util.*;

import math.util;


public class ironloss3dani{
	
	public static void main(String args[]) throws IOException{
		
		
		int fluxFormat=1; // 0 for java program results, 1 for fortran 
		
		int m;
		int npoint=0,nelem=0,nedge=0,nblock=0;
		int nstep=0,nperiod;
		double factor=0.;
		int nod[][];
		double xyz[][];
		int nbl[][];
		String name[];
		double flux[][];
		double eddy[][][];
		double AbsB[],BR[],BT[],Angle[];
		double zl[],eddyl[];
		String meshname,fluxname,eddyname,lossname,lossfilename,totallossname;
		String line,s;
		double totalloss=0.0;
		BufferedReader myReader=new BufferedReader(new InputStreamReader(System.in),1);
/*
		System.out.println("Input Mesh File");
		meshname=myReader.readLine();
		System.out.println("Input Flux File");
		fluxname=myReader.readLine();
		System.out.println("Output totalLoss");
		lossname=myReader.readLine();
		System.out.println("Output totaltotalLoss");
		totallossname=myReader.readLine();
		System.out.println("Output Loss File");
		lossfilename=myReader.readLine();
		*/
	
		meshname="mesh";
		if(fluxFormat==0)
		fluxname="fluxJava";
		else
		fluxname="fluxFort";
		lossname="loss";
		System.out.println("Output totaltotalLoss");
		totallossname="totLoss";
		lossfilename="lossout";
		
		nstep=1;
		nperiod=1;
		System.out.println(meshname+".txtファイルを開いています！！");

		BufferedReader br=new BufferedReader(new FileReader(meshname+".txt"));
		for(int icount=0;icount<2;icount=icount){
			line = br.readLine();
			if(line.equals("mesh")){
				icount=icount+1;
				line = br.readLine();
				StringTokenizer tk=new StringTokenizer(line," ");
				m = tk.countTokens();
				for(int j=1; j<=m; j++){
					try{
						s = tk.nextToken();
						if(j==1) {npoint	=Integer.parseInt(s);}    //
						if(j==2) {nelem		=Integer.parseInt(s);}    //
						if(j==3) {nedge		=Integer.parseInt(s);}    //座標の単位設定
						if(j==4) {factor	=Double.parseDouble(s);}  //
					}catch(NoSuchElementException e0){
						System.out.println("No Token");
					}
				}
			}
			if(line.equals("nbl")){
				icount=icount+1;
				line = br.readLine();
				nblock=Integer.parseInt(line);
			}
		}
		br.close();
		flux = new double[nelem+1][3+1]; //得られた磁束密度
		xyz = new double[npoint+1][3+1];					//xyz座標
		nod = new int[nelem+1][8+1];
		nbl = new int[nblock+1][3+1];
		name = new String[nblock+1];
		
		BufferedReader br1=new BufferedReader(new FileReader(meshname+".txt"));
		for(int icount=0;icount<2;icount=icount){
			line = br1.readLine();
			if(line.equals("mesh")){
				icount=icount+1;
				line = br1.readLine();
				for(int i=1; i<=npoint; i++){
					line = br1.readLine();
					StringTokenizer tm=new StringTokenizer(line," ");
					m = tm.countTokens();
					for(int j=1; j<=m; j++){
						try{
							s = tm.nextToken();
							xyz[i][j]=Double.parseDouble(s)/factor;
						}catch(NoSuchElementException e1){
							System.out.println("No Token");
						}
					}
				}
			
				for(int i=1; i<=nelem; i++){
					line = br1.readLine();
					StringTokenizer tm=new StringTokenizer(line," ");
					m = tm.countTokens();
					for(int j=1; j<=m; j++){
						try{
							s = tm.nextToken();
							nod[i][j]=Integer.parseInt(s);
						}catch(NoSuchElementException e1){
							System.out.println("No Token");
						}
					}
				}
			}
			if(line.equals("nbl")){
				icount=icount+1;
				line = br1.readLine();
				for(int ib=1;ib<=nblock;ib=ib+1){
					line = br1.readLine();
					StringTokenizer tk=new StringTokenizer(line," ");
					m = tk.countTokens();
					for(int j=1; j<=m; j++){
						try{
							s = tk.nextToken();
							if(j==1) {nbl[ib][1]=Integer.parseInt(s);}    //
							if(j==2) {nbl[ib][2]=Integer.parseInt(s);}    //
							if(j==3) {nbl[ib][3]=Integer.parseInt(s);}    //
							if(j==4) {name[ib]=s;}
						}catch(NoSuchElementException e0){
							System.out.println("No Token");
						}
					}
				}
			}
		}

		int nstart=nstep-nperiod+1;
		
		BufferedReader br2=null;
		
		if(fluxFormat==0)
		 br2=new BufferedReader(new FileReader(fluxname+".txt"));
		else if(fluxFormat==1)
		 br2=new BufferedReader(new FileReader(fluxname+".dat"));
		
		if(fluxFormat==1){
			br2.readLine();
			br2.readLine();
		}
		
		line = br2.readLine();
		for(int is=1;is<nstart;is++){
			for(int ie=1;ie<=nelem;ie++){
				line = br2.readLine();
			}
		}
		for(int ie=1;ie<=nelem;ie++){
			line = br2.readLine();
			StringTokenizer tk=null;
		
			if(fluxFormat==0)
			tk=new StringTokenizer(line," ");
			else if(fluxFormat==1)
			tk=new StringTokenizer(line,",");
			
			m = tk.countTokens();
			for(int j=1; j<=m; j++){
				try{
					s = tk.nextToken();
					if(j==1) {flux[ie][1]=Double.parseDouble(s);}    //
					if(j==2) {flux[ie][2]=Double.parseDouble(s);}    //
					if(j==3) {flux[ie][3]=Double.parseDouble(s);}    //
				}catch(NoSuchElementException e0){
					System.out.println("No Token");
				}
			} 
		}
	for(int i=0;i<20;i++)
		util.hshow(flux[i]);
		
		System.out.println("nelem="+nelem);
		System.out.println("nod[nelem][8]="+nod[nelem][8]);
		System.out.println("npoint="+npoint);
		System.out.println("xyz[npoint][3]="+xyz[npoint][3]);
		System.out.println("nblock="+nblock);
		System.out.println("nbl[nblock][2]="+nbl[nblock][2]);
		System.out.println("name[1]="+name[1]);
		System.out.println("flux[nelem][3]="+flux[nelem][3]);
		System.out.println("factor="+factor);

		BR   = new double[nelem+1]; //圧延方向の磁束密度
		BT   = new double[nelem+1]; //直角方向の磁束密度
		AbsB = new double[nelem+1]; //磁束密度の合成値
		Angle= new double[nelem+1]; // 圧延方向との角度
		String s1="50H230", s2="50H400";
		int score1=1,ncore1=5160;
		int score2=5161,ncore2=5590;
		
		
		
		//圧延方向と直角方向の定義,圧延方向との角度の計算
		for(int ib=1;ib<=nblock;ib=ib+1){
			if(name[ib].equals(s1)){
				for(int ie=nbl[ib][1];ie<=nbl[ib][2];ie++){
				}
			}
			//磁束密度の合成値の計算square root calculation
			for(int ie=1;ie<=ncore2;ie++){
				AbsB[ie]=Math.sqrt(flux[ie][1]*flux[ie][1]+flux[ie][2]*flux[ie][2]+flux[ie][3]*flux[ie][3]);
			}
		}
			// loss data
			
			double ready1[],Wiron1[],W1[],B1[],Iloss1[];
			ready1 = new double[5+1];
			ready1[1]=0.;
			ready1[2]=Math.PI/8.;
			ready1[3]=Math.PI/4.;
			ready1[4]=Math.PI*3./8.;
			ready1[5]=Math.PI/2.;
			B1 = new double[41+1];
			Iloss1 = new double[41+1];

			B1[1]=0.00; B1[2]=0.05; B1[3]=0.10; B1[4]=0.15; B1[5]=0.20;
			B1[6]=0.25; B1[7]=0.30; B1[8]=0.35; B1[9]=0.40; B1[10]=0.45;
			B1[11]=0.50;B1[12]=0.55;B1[13]=0.60;B1[14]=0.65;B1[15]=0.70;
			B1[16]=0.75;B1[17]=0.80;B1[18]=0.85;B1[19]=0.90;B1[20]=0.95;
			B1[21]=1.00;B1[22]=1.05;B1[23]=1.10;B1[24]=1.15;B1[25]=1.20;
			B1[26]=1.25;B1[27]=1.30;B1[28]=1.35;B1[29]=1.40;B1[30]=1.45;
			B1[31]=1.50;B1[32]=1.55;B1[33]=1.60;B1[34]=1.65;B1[35]=1.70;
			B1[36]=1.75;B1[37]=1.80;B1[38]=1.85;B1[39]=1.90;B1[40]=1.95;
			B1[41]=2.0;

			double[] Iloss2= //方向性電磁鋼帯のP88の鉄損データと比較せよ。[W/kg]
				{0,0,0,0,0,0,0.13,0.17,0.21,0.25,0.31,0.33,0.41,0.49,0.53,0.61,0.69,0.78,0.82,0.90,1.00,1.10,1.15,1.27,1.38,
					1.50,1.60,1.77,1.92,2.10,2.25,2.38,2.50,2.62,2.80,2.95,3.00,3.07,3.12,3.18,3.20};


			for(int j=1;j<=41;j++){
				Iloss1[j]=Iloss2[j-1];
			}

	// interpolation
			Wiron1 = new double[nelem+1];
			W1     = new double[nelem+1];

			for(int ie=score1;ie<=ncore1;ie++){
				//	for(int i=1;i<=5;i++){
				for(int j=2;j<=41;j++){
					if((AbsB[ie]-1.E-8)>=B1[j]&&(AbsB[ie]-1.E-8)<=B1[j+1]){
						Wiron1[ie]=Iloss1[j+1]-(Iloss1[j+1]-Iloss1[j])*(B1[j+1]-AbsB[ie])/(B1[j+1]-B1[j]);
			
					}
					if((AbsB[ie]-1.E-8)<B1[2]){
						Wiron1[ie]=Iloss1[2]*AbsB[ie]/B1[2];
					}
				}
			}

				//			if(name[ib].equals(s2)){
				//				for(int ie=nbl[ib][1];ie<=nbl[ib][2];ie++){
				//					BR[ie]=Math.abs(flux[ie][2]);
				//					BT[ie]=Math.sqrt(flux[ie][3]*flux[ie][3]+flux[ie][1]*flux[ie][1]);
				//					BT[ie]=Math.abs(flux[ie][1]);
				//					Angle[ie]=Math.atan(BT[ie]/BR[ie]);
				
				

				// loss data
				
				double ready[],Wiron2[],W2[],B2[],Iloss3[];
				//		ready = new double[5+1];
				//		ready[1]=0.;
				//		ready[2]=Math.PI/8.;
				//		ready[3]=Math.PI/4.;
				//		ready[4]=Math.PI*3./8.;
				//		ready[5]=Math.PI/2.;
				B2 = new double[41+1];
				Iloss3 = new double[41+1];
				B2[1]=0.00; B2[2]=0.05; B2[3]=0.10; B2[4]=0.15; B2[5]=0.20;
				B2[6]=0.25; B2[7]=0.30; B2[8]=0.35; B2[9]=0.40; B2[10]=0.45;
				B2[11]=0.50;B2[12]=0.55;B2[13]=0.60;B2[14]=0.65;B2[15]=0.70;
				B2[16]=0.75;B2[17]=0.80;B2[18]=0.85;B2[19]=0.90;B2[20]=0.95;
				B2[21]=1.00;B2[22]=1.05;B2[23]=1.10;B2[24]=1.15;B2[25]=1.20;
				B2[26]=1.25;B2[27]=1.30;B2[28]=1.35;B2[29]=1.40;B2[30]=1.45;
				B2[31]=1.50;B2[32]=1.55;B2[33]=1.60;B2[34]=1.65;B2[35]=1.70;
				B2[36]=1.75;B2[37]=1.80;B2[38]=1.85;B2[39]=1.90;B2[40]=1.95;
				B2[41]=2.0;

				double[] Iloss4=    //方向性電磁鋼帯のP88の鉄損データと比較せよ。[W/kg]
					{0,0,0,0,0.14,0.195,0.25,0.32,0.38,0.46,0.53,0.62,0.70,0.80,0.88,1.00,1.18,1.20,1.30,1.41,1.55,1.68,1.80,1.97,2.10,
						2.27,2.45,2.65,2.90,3.38,3.60,3.80,3.95,4.01,4.10,4.20,4.23,4.31,4.34,4.38,4.39};


				for(int j=1;j<=41;j++){
					Iloss3[j]=Iloss4[j-1];
				}

				// interpolation
				Wiron2 = new double[nelem+1];
				W2     = new double[nelem+1];

				for(int ie=score2;ie<=ncore2;ie++){
					//	for(int i=1;i<=5;i++){
						
					for(int j=2;j<41;j++){
						if((AbsB[ie]-1.E-8)>=B2[j]&&(AbsB[ie]-1.E-8)<=B2[j+1]){
							Wiron2[ie]=Iloss3[j+1]-(Iloss3[j+1]-Iloss3[j])*(B2[j+1]-AbsB[ie])/(B2[j+1]-B2[j]);
						}
					}
				}
				//		System.out.println("BR[9]="+BR[9]);
				//		System.out.println("BT[9]="+BT[9]);
				//		System.out.println("Angle[9]="+Angle[9]);
				//		System.out.println("BR[2288]="+BR[2288]);
				//		System.out.println("BT[2288]="+BT[2288]);
				//		System.out.println("Angle[2288]="+Angle[2288]);

				//		System.out.println("AbsB[2288]="+AbsB[2288]);		
				// 圧延方向との角度別の異方性鉄損データ：磁束密度B[35],鉄損Iloss[5][35]
				// 0℃圧延方向、22.5℃、45℃、67.5℃、90℃直角方向
			
				double Wiron[][],W[],B[],Iloss[][];
				ready = new double[5+1];
				ready[1]=0.;
				ready[2]=Math.PI/8.;
				ready[3]=Math.PI/4.;
				ready[4]=Math.PI*3./8.;
				ready[5]=Math.PI/2.;
				B = new double[40+1];
				Iloss = new double[5+1][40+1];
				B[1]=0.0; B[2]=0.05; B[3]=0.1; B[4]=0.15; B[5]=0.20; B[6]=0.25;
				B[7]=0.30; B[8]=0.35; B[9]=0.40; B[10]=0.45; B[11]=0.50; B[12]=0.55;
				B[13]=0.60;	B[14]=0.65;	B[15]=0.70;	B[16]=0.75;	B[17]=0.80;	B[18]=0.85;
				B[19]=0.90;	B[20]=0.95;	B[21]=1.00;	B[22]=1.05;	B[23]=1.10;	B[24]=1.15;
				B[25]=1.20;	B[26]=1.25;	B[27]=1.30;	B[28]=1.35;	B[29]=1.40;	B[30]=1.45;
				B[31]=1.50;	B[32]=1.55;	B[33]=1.60;	B[34]=1.65;	B[35]=1.70;
				B[36]=1.75;	B[37]=1.8;	B[38]=1.85;	B[39]=1.9;	B[40]=2.2;

				double[][] Iloss22={    //方向性電磁鋼帯のP88の鉄損データと比較せよ。[W/kg]
						{0.,0.009,0.018,0.027,0.036,0.045,0.054,0.063,0.072,0.081,0.09,
							0.108,0.128,0.15,0.171,0.19,0.22,0.245,0.275,0.32,0.335,0.362,
							0.42,0.44,0.48,0.517,0.56,0.59,0.645,0.675,0.73,0.77,0.83,0.88,0.98,
							1.05,1.2,1.35,1.6,3.1},

							{0.,0.0222985,0.044597,0.0668955,0.089194,0.1114925,0.133791,0.1560895,
								0.178388,0.2006865,0.222985075,0.26758209,0.317134328,0.371641791,
								0.423671642,0.470746269,0.545074627,0.607014925,0.681343284,0.792835821,
								0.83,0.89,0.97,1.06,1.16,1.27,1.4,1.54,1.72,1.92,2.17,2.3,2.43,2.56,2.68,
								2.8,2.92,3.04,3.16,3.88},

								{0.,0.0416418,0.0832836,0.1249254,0.1665672,0.208209,0.2498508,0.2914926,
									0.3331344,0.3747762,0.41641791,0.499701493,0.592238806,0.694029851,0.79119403,
									0.879104478,1.017910448,1.13358209,1.27238806,1.480597015,1.55,1.7,1.85,2,2.2,
									2.4,2.6,2.84,3.02,3.2,3.35,3.45,3.62,3.8,3.86,4.16,4.34,4.52,4.7,5.78},

									{0.,0.0505075,0.101015,0.1515225,0.20203,0.2525375,0.303045,0.3535525,0.40406,
										0.4545675,0.505074627,0.606089552,0.718328358,0.841791045,0.959641791,
										1.066268657,1.234626866,1.374925373,1.543283582,1.795820896,1.88,2.04,2.21,
										2.4,2.6,2.79,2.99,3.19,3.35,3.44,3.6,3.66,3.88,4.02,4.2,4.38,4.56,4.74,4.92,5.42},

										{0.,0.0454,0.0908,0.1362,0.1816,0.227,0.2724,0.3178,0.3632,0.4086,0.454029851,
											0.544835821,0.645731343,0.756716418,0.862656716,0.958507463,1.109850746,1.235970149,
											1.387313433,1.614328358,1.69,1.85,1.99,2.12,2.2,2.3,2.5,2.85,3.1,3.42,3.7,3.92,4.1,4.22,4.35,
											4.48,4.61,4.74,4.87,5.52}
				};
				for(int i=1;i<=5;i++){
					for(int j=1;j<=40;j++){
						Iloss[i][j]=Iloss22[i-1][j-1];
					}
				}
				System.out.println("ready[2]="+ready[2]);		

				System.out.println("Iloss[5][9]="+Iloss[5][9]);		
				System.out.println("Iloss2[4][8]="+Iloss22[4][8]);			
	
/*
				
				//角度別の異方性鉄損データのBによる線形補間
				Wiron = new double[5+1][nelem+1];
				W     = new double[nelem+1];
				for(int ie=1;ie<=ncore;ie++){
					for(int i=1;i<=5;i++){
						for(int j=2;j<=39;j++){
							if((AbsB[ie]-1.E-8)>=B[j]&&(AbsB[ie]-1.E-8)<=B[j+1]){
								Wiron[i][ie]=Iloss[i][j+1]-(Iloss[i][j+1]-Iloss[i][j])*(B[j+1]-AbsB[ie])/(B[j+1]-B[j]);
							}
							if((AbsB[ie]-1.E-8)<B[2]){
								Wiron[i][ie]=Iloss[i][2]*AbsB[ie]/B[2];
							}
						}
					}
				}
				
				
				
				//		System.out.println("Wiron[1][2288]="+Wiron[1][2288]);	
				//		System.out.println("Wiron[2][2288]="+Wiron[2][2288]);	
				//		System.out.println("Wiron[3][2288]="+Wiron[3][2288]);	
				//		System.out.println("Wiron[4][2288]="+Wiron[4][2288]);	
				//		System.out.println("Wiron[5][2288]="+Wiron[5][2288]);	
				//角度間の異方性鉄損データの角度による補間
		
				for(int ie=1;ie<=ncore;ie++){
					for(int i=1;i<=4;i++){
						if((Angle[ie]-1.E-8)>=ready[i]&&(Angle[ie]-1.E-8)<=ready[i+1]){
							W[ie]=((ready[i+1]-Angle[ie])*Wiron[i][ie]+(Angle[ie]-ready[i])*Wiron[i+1][ie])/(ready[i+1]-ready[i]); //内分
						}
					}
				}

				// 0℃圧延方向、22.5℃、45℃、67.5℃、90℃直角方向

				//		System.out.println("ready[2]="+ready[2]);		

				//		System.out.println("Iloss[5][9]="+Iloss[5][9]);		
				//		System.out.println("Iloss2[4][8]="+Iloss2[4][8]);			





				//		System.out.println("ready[2]="+ready[2]);		

				//		System.out.println("Iloss[5][9]="+Iloss[5][9]);		
				//		System.out.println("Iloss2[4][8]="+Iloss2[4][8]);			

*/


				//		System.out.println("W[2288]="+W[2288]);		
				double V[],Losspervol[],Lossinvol[],V1[],Losspervol1[],Lossinvol1[],V2[],Losspervol2[],Lossinvol2[];
				V          =new double[nelem+1];
				Losspervol =new double[nelem+1];
				Lossinvol  =new double[nelem+1];
				V1         =new double[nelem+1];
				Losspervol1=new double[nelem+1];
				Lossinvol1 =new double[nelem+1];
				V2         =new double[nelem+1];
				Losspervol2=new double[nelem+1];
				Lossinvol2 =new double[nelem+1];
				double Vtotal=0.;
				double hen1,hen2,hen3;
				for(int ie=1;ie<=ncore2;ie++){
					V2[ie]=Math.abs((xyz[nod[ie][5]][1]-xyz[nod[ie][6]][1])*(xyz[nod[ie][8]][2]-xyz[nod[ie][5]][2])*(xyz[nod[ie][5]][3]-xyz[nod[ie][1]][3]));
					hen1=Math.abs(xyz[nod[ie][5]][1]-xyz[nod[ie][6]][1]);
					hen2=Math.abs(xyz[nod[ie][8]][2]-xyz[nod[ie][5]][2]);
					hen3=Math.abs(xyz[nod[ie][5]][3]-xyz[nod[ie][1]][3]);
					V1[ie]=hen1*hen2*hen3;
					Vtotal=Vtotal+V1[ie];   //コアの体積
					Losspervol1[ie]=Wiron1[ie]*7650.;    //磁気特性の密度[kg/dm^3]（1000*kg/m^3）と先程の損失[W/kg]の積
					Lossinvol1[ie]=Wiron1[ie]*7650.*V1[ie];     //[W/m^3]*[m^3]=[W]
					Losspervol2[ie]=Wiron2[ie]*7650.;    //磁気特性の密度[kg/dm^3]（1000*kg/m^3）と先程の損失[W/kg]の積
					Lossinvol2[ie]=Wiron2[ie]*7650.*V2[ie];     //[W/m^3]*[m^3]=[W]	
					Lossinvol[ie]=Lossinvol1[ie]+Lossinvol2[ie];
					//		    System.out.println("Lossinvol"+Lossinvol[ie]);

					totalloss=totalloss+Lossinvol[ie];	//[W]		
					//		System.out.println("totalloss="+totalloss);

				}	
				totalloss=totalloss/Vtotal;    //[W]/[m^3]

				System.out.println("Vtotal="+Vtotal);	
				System.out.println("totalloss2="+totalloss);	
				// bmaxを探す
				double bmax[],bb1=0.,bb2=0.;
				bmax = new double[2+1];
				for(int ie=1;ie<=ncore2;ie++){
					bb1=AbsB[ie];
					bb2=Math.abs(flux[ie][2]);
					if(bb1>=bmax[1]) {
						bmax[1]=bb1;
					}
					if(bb2>=bmax[2]) {
						bmax[2]=bb2;
					}
				}
				System.out.println("bmax[1]="+bmax[1]);		
				System.out.println("bmax[2]="+bmax[2]);		

				System.out.println(lossname+".txtファイルを書き込みます！！");
				PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(lossname+".txt")));
				for(int ie=1;ie<=ncore2;ie++){
					pw.println("損失＝"+ie+","+Lossinvol[ie]);
				}
				pw.println("体積＝"+Vtotal);
				pw.close();

				System.out.println(lossname+".txtファイルを書き込みます！！");
				PrintWriter pw2 = new PrintWriter(new BufferedWriter(new FileWriter(totallossname+".txt")));			
				pw2.println("損失＝"+totalloss);
				pw2.println("体積＝"+Vtotal);
				pw2.close();

				System.out.println(lossfilename+".txtファイルを書き込みます！！");

				PrintWriter pw1 = new PrintWriter(new BufferedWriter(new FileWriter(lossfilename+".txt")));	
				pw1.println("flux");
				pw1.println("3");
				pw1.println(nelem);
				for(int ie=1;ie<=nelem;ie++){		
					pw1.println(0.+","+0.+","+(Wiron1[ie]+Wiron2[ie]));
				}
				pw1.close();
			}

		
	}
