package fem;
import java.io.File;

import java.text.DecimalFormat;

import femSolver.StaticElectricSolver;
import main.Main;
import math.Vect;
import math.util;

public class RunMag {

	private DecimalFormat formatter=new DecimalFormat("0.00");

	private Vect xp;
	int[] nn=new int[1];

	public static void main(String[] args){
		new Main();
	}

	public void runMag(Model model, Main main){
		

		
			String folder=model.resultFolder;


			int nTsteps=model.nTsteps;
		
			Vect T=new Vect(nTsteps);
			int ix=0;
			
			int nBegin=model.nBegin;
			int nEnd=model.nEnd;
			

			int inc=model.nInc;


			
			String dispFolder=model.resultFolder;
			String fluxFolder=model.resultFolder;

			if(model.saveForce){
	/*		
			
					dispFolder = folder+"\\forces";
				File dfolder = new File(dispFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();
*/
			}

			if(model.saveFlux){
/*					fluxFolder = folder+"\\fluxes";
				
				File dfolder = new File(fluxFolder);
				if(dfolder.exists())
					util.deleteDir(dfolder);
				dfolder.mkdir();
*/
			}

			
			if(nTsteps>1)
				model.writeFiles=false;
			else
				model.writeFiles=true;

			boolean writeFiles=model.writeFiles;
			
			
			Vect x=new Vect();

			main.gui.lbX[0].setText("rot ang. ");
			main.gui.lbX[1].setText("Trq. ");

			
			model.solveCoils();

	
				for(int i=nBegin;i<=nEnd;i+=inc){
					
				double t0=0;
		
				
				main.gui.tfX[0].setText((i)+"/"+nEnd);
				
				model.setJ0(t0+i*model.dt);	
				if(model.phiCoils!=null)
				model.phiCoils[0].current=Math.cos(2*model.freq*i*model.dt);

					if(!model.loadFlux){
						
						
							model.setMagBC();
							
						
							if(i==nBegin){
								model.writer.reportData(model);
							}

							
							if(!model.nonLin /*|| i==nBegin*//* ||model.Q!=null*/){
								
							
									if(!model.loadPrevMag){
										model.saveAp();		
										x=model.solveMagLin(i-nBegin,x);	
									

									}
									else
									{
							
									}

								}
								else
									x=this.xp;
						
								
							if(model.nonLin ){

								if(x==null) x=new Vect(model.numberOfUnknowns);
									if(model.analysisMode>0 && i!=nBegin )	{	
								
										model.saveAp();				
									}
						
									x=model.solveNonLinear(x,true,i-nBegin);
								}

						if(model.analysisMode>0)
							model.setJe();
				
						if(x!=null)
							this.xp=x.deepCopy();

							if(model.saveFlux){
								String fluxFile = fluxFolder+"\\flux"+i+".txt";
								if(i==nBegin)
									model.writeMesh(fluxFolder+"\\bun"+i+".txt");
			
								model.writeB(fluxFile);
								
		
							}


					//	model.resetReluctForce();
					//	model.setReluctForce();
						//
						//model.setMSForce();
						
						
				}


					
			
					/*if(model.dim==2)
						model.setTorque(0,model.rm,1);
						else
							model.setTorque(model.r1,1,1);

					
					
					util.pr("torque >>>>>>>"+model.TrqZ);*/
					
				//	T.el[ix++]=model.TrqZ;
					
	
					
					
								//************************************************
				//	writeFiles=true;

					if(model.phiCoils==null)
						model.writeJ0(folder+"\\J"+i+".txt");
					
					double loss=0;
										
					if(model.saveFlux){
						model.setSolution(x);	

						String fluxFile = folder+"\\flux"+i+".txt";
						model.writeB(fluxFile);
					}
					if(model.saveJe){
						String JeFile =  folder+"\\Je"+i+".txt";

			
						model.setJe();	
	
						model.writeJe(JeFile);
						for(int ir=1;ir<=model.numberOfRegions;ir++){
							if(model.region[ir].isConductor){
						 loss=model.obtainLoss(ir);
						util.pr("Joule Loss of region "+ir+"=  "+loss);
							}
						}
					}
					
					T.el[ix++]=loss;
				
					
					if(model.saveForce){
					//	String dispFile =dispFolder+"\\force"+mangs+".txt";
						String forceFile =dispFolder+"\\force"+i+".txt";
						model.writeNodalField(forceFile,model.forceCalcMode);
					}


					main.gui.tfX[1].setText(this.formatter.format(model.TrqZ));

					if(model.solver.terminate){

						break;
					}

			if(model.solver.terminate) break;


			}
			//	util.pr(model.main.gui.iccgArea.getText());

				//util.plot(T);
			//	T.show();
				
				Vect errs=new Vect(model.solver.totalIter);
				
				for(int i=0;i<errs.length;i++)
					errs.el[i]=model.solver.errs.get(i);
				
				util.plot("error",errs.el,"ICCG Convergence");

			//	util.plot(errs);

				boolean writeInit=false;
				
				if(writeInit){
			
			
				String initfile = folder+"\\initxMag.txt";
				double[][] arr=new double[model.numberOfUnknowns+1+1][2];
				
				Vect vk=model.getUnknownA();

				Vect vp=model.getUnknownAp();

				for(int j=0;j<vp.length;j++){
					arr[j][0]=vp.el[j];
					arr[j][1]=vk.el[j];
				}
				

				for(int i=0;i<model.numberOfUnknownCurrents;i++){

					int nr=model.unCurRegNumb[i];
					
					arr[i+vp.length][0]=model.region[nr].currentp;
					arr[i+vp.length][1]=model.region[nr].current;


				}

				arr[arr.length-2][0]=model.vNeutral;
				arr[arr.length-2][1]=0;
				
				arr[arr.length-1][0]=0;
				arr[arr.length-1][1]=0;



				model.writer.writeArray(arr,initfile);
			
				}




	}

}
