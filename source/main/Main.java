package main;
import math.*;
import static java.lang.Math.*;
import io.Console;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Scanner;

import javax.swing.JOptionPane;

import java.io.FileDescriptor;
import java.io.FileOutputStream;

import components.GUI;
import fem.Model;
import fem.RunCLN;
import fem.RunMag;
import fem.RunMagAC;

	public class Main implements ActionListener{

	public GUI gui;

	private Model model;
	private Thread thread;
	private int iterMax=1000;	
	private double  errMax;
	private String path = System.getProperty("user.dir");
	public  boolean console=true,dated=false;

	public Main()
	{		
		this.model=new Model();
		this.model.main=this;
		this.gui=new GUI(this.path);
		this.gui.Run.addActionListener(this);
		this.gui.bTerminate.addActionListener(this);
		this.gui.setVisible(true);
		//this.gui.Run.doClick();
	}	
	
	public static void main(String[] args){
	
		new Main();
	}

		
	public void runGUI(){
		
		if(console){
			Console.redirectOutput(this.gui.dataArea);
		}

		
		this.model.meshFilePath=this.gui.tfMeshFile.getText();
		this.model.dataFilePath=this.gui.tfDataFile.getText();
	
		this.model.resultFolder = new File(model.meshFilePath).getParentFile().getAbsolutePath();
		
		DateFormat dateFormat = new SimpleDateFormat("MM.dd.HH.mm.ss");
		
		Date date = new Date();
		
		String suff=dateFormat.format(date);

		if(dated) this.model.resultFolder=this.model.resultFolder+suff;

		
		model.logFilePath=model.resultFolder+ "\\log.txt";

		util.setLogFile(model.logFilePath);
		
		if(util.getLogFile()!=null){
		try{
			DateFormat dateFormat1 = new SimpleDateFormat("YY/MM/dd  HH:mm:ss");

			PrintWriter pw=new PrintWriter(new BufferedWriter(new FileWriter(util.getLogFile())));
			pw.println("***  log file *** started at "+ dateFormat1.format(date)+" ****");
			pw.println("================================================================");
			pw.close();
		} catch(IOException e){System.err.println("IOException: " + e.getMessage());}
		}

		this.errMax=Double.parseDouble(this.gui.tfErrorMax.getText());
		this.iterMax=Integer.parseInt(this.gui.tfIterMax.getText());
		this.thread=new Thread(){
			long m1 = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
			@Override
			public void run(){

				double t_start= System.currentTimeMillis();
				
				Main.this.model.loadMesh(model.meshFilePath);
			
	
				prepare();
				
				model.loadData(model.dataFilePath);

		
				 if(model.magAnalysis) {
			
					runMag(); 

					}
		
		
			double t_end= System.currentTimeMillis();
			System.out.format("Total cpu time (s): %10.1f\n",(t_end-t_start)/1000.);	

				gui.writeLog(model.logFilePath);
				gui.Run.setBackground(Color.green);
				gui.Run.setEnabled(true);

			}
		};
		this.thread.start();		
	}

	public void runMag(){


			if(model.AC){
				 RunMagAC module=new RunMagAC();

				 module.runMag(model, this);
				
			}
			else if(model.nCLNstages>0){
				 RunCLN module=new RunCLN();

				 module.run(model, this);
			
				}
			else{
			 RunMag module=new RunMag();

			 module.runMag(model, this);
		
			}
		
}

		public void prepare(){
		Main.this.model.iterMax=Main.this.iterMax;
		Main.this.model.errCGmax=Main.this.errMax;
	
		
		
		try{

			PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(new File("").getAbsolutePath()+"\\_last_MagFEM_elected_path")));		
			pw.println(this.model.meshFilePath);
			pw.println(this.model.dataFilePath);

			pw.close();
		}
		catch(IOException e){}


		//this.model.fluxFilePath=System.getProperty("user.dir")+"\\flux.txt";
		//this.model.eddyFilePath=System.getProperty("user.dir")+"\\eddy.txt";
		this.gui.iccgArea.setText("");
		if(console){
		//Console.redirectOutput(this.gui.iccgArea);
		Console.redirectOutput(this.gui.dataArea);
		}
		this.gui.Run.setBackground(Color.gray);
		this.gui.Run.setEnabled(false);
	}



	
	public void loadMesh(){	

				System.gc();
				long m1 = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();

				Main.this.model.loadMesh(Main.this.model.meshFilePath);

				long m2 = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
				System.out.println("Used memory for  model setup: "+(m2-m1)/(1024*1024)+"MB");

					long m3 = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
					System.out.println("Used memory for drawing mesh: "+(m3-m2)/(1024*1024)+"MB");
					System.out.println();
					System.out.println(" Number of regions: "+Main.this.model.numberOfRegions);
					System.out.println(" Total number of elements: "+Main.this.model.numberOfElements);
					System.out.println(" Total number of nodes: "+Main.this.model.numberOfNodes);

					System.out.println();

	}		


	public String readFirst(String filePath){
		try{

			Scanner scr=new Scanner(new FileReader(filePath));
			String first= scr.next();
			scr.close();
			return first;
		}
		catch(IOException e){return null;}

	}


	

	@Override
	public void actionPerformed(ActionEvent e)
	{	

		 if(e.getSource()==this.gui.Run){
			this.gui.Run.setBackground(Color.gray);
			this.gui.Run.setEnabled(false);

			runGUI();

		}
		else if(e.getSource()==this.gui.bTerminate){
			this.model.solver.terminate(true);

		}

	
	}








}

