package components;

import static java.lang.Math.abs;

import java.awt.Color;

import javax.swing.JFrame;

import org.math.plot.Plot2DPanel;

import math.Vect;

public class DynamicPlotter {
	int width = GUI.screenWidth;
	int height = GUI.screenHeight;

	 public Plot2DPanel[] plot = new Plot2DPanel[4];
	 public  JFrame[] frame = new JFrame[4];
	
public void setPlotWinodws(){

	//******************************
	
	 
	 double[][] XY=new double[1][2];
	
	 int hh=40;
	 
	 for(int i=0;i<plot.length;i++){
		plot[i]=new Plot2DPanel();
	 
	 }

	 plot[0].addLinePlot("Ia", Color.red, XY);
	 plot[0].addLinePlot("Ib", Color.green, XY);
	 plot[0].addLinePlot("Ic", Color.blue, XY);
	 plot[0].addLinePlot("sumI", Color.black, XY);
		plot[0].setFixedBounds(0,0,50);
		plot[0].setFixedBounds(1,-10,10);

	   frame[0] = new JFrame("Currents");
	   frame[0].setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
	  frame[0].setSize(width/3,height/3);
	  frame[0].setContentPane(plot[0]);
	  frame[0].setLocation(width/3,height/3-hh);
	
	  
		plot[1] = new Plot2DPanel();

		plot[1].addLinePlot("Va", Color.red, XY);
		plot[1].addLinePlot("Vb", Color.green, XY);
		plot[1].addLinePlot("Vc", Color.blue, XY);
		plot[1].addLinePlot("Vfa", Color.cyan, XY);
		plot[1].addLinePlot("Vfb", Color.pink, XY);
		plot[1].addLinePlot("Vfc", Color.gray, XY);
			plot[1].setFixedBounds(0,0,50);
			plot[1].setFixedBounds(1,-100,100);

		 frame[1] = new JFrame("Voltages");
		  frame[1].setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
			 frame[1].setSize(width/3,height/3);
		 frame[1].setContentPane(plot[1]);
		 frame[1].setLocation(2*width/3,height/3-hh);


		  
		   plot[2] = new Plot2DPanel();
	
			 plot[2].addLinePlot("Torque", Color.red, XY);
			
				plot[2].setFixedBounds(0,0,50);
				plot[2].setFixedBounds(1,-50,50);
				

			   frame[2] = new JFrame("Torque");
			   frame[2].setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
			  frame[2].setSize(width/3,height/3);
			  frame[2].setContentPane(plot[2]);
			  frame[2].setLocation(width/3,2*height/3-hh);

	  
			   plot[3] = new Plot2DPanel();
				
				 plot[3].addLinePlot("Vg", Color.red, XY);
				
					plot[3].setFixedBounds(0,0,50);
					plot[3].setFixedBounds(1,-1.2,1.2);
					
				  frame[3] = new JFrame("Vg");
				  frame[3].setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
				 frame[3].setSize(width/3,height/3);
				 frame[3].setContentPane(plot[3]);
				 frame[3].setLocation(2*width/3,2*height/3-hh);

	 
}

public void setDynamicVisibility(){
	  frame[0].setVisible(true);
		 frame[1].setVisible(true);
		  frame[2].setVisible(true);

		 frame[3].setVisible(true);
}

public void updateData(int ix,Vect T, Vect Vn,double[][] Is,double[][] Vf,double[][] Vs){

	
	 double w1=-1e10,w2=1e-10;
	 
	double[][] XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Is[1][j];
		

		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[0].changePlotData(0, XY);

	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Is[2][j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[0].changePlotData(1, XY);

	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Is[3][j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[0].changePlotData(2, XY);

	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Is[4][j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[0].changePlotData(3, XY);
	
	plot[0].setFixedBounds(0,0,50*(1+ix/50));

	double mx=1.2*Math.max(w1,w2);
	
				
	//************************************************
	
	
	  w1=-1e10;
	  w2=1e-10;
	 
	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Vs[1][j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[1].changePlotData(0, XY);

	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Vs[2][j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[1].changePlotData(1, XY);

	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Vs[3][j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[1].changePlotData(2, XY);

	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Vf[1][j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[1].changePlotData(3, XY);

	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Vf[2][j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[1].changePlotData(4, XY);

	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Vf[3][j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[1].changePlotData(5, XY);
	
	plot[1].setFixedBounds(0,0,50*(1+ix/50));


	
	//************************************************
	
	
	  w1=-1e10;
	  w2=1e-10;
	 
	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=T.el[j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[2].changePlotData(0, XY);

	plot[2].setFixedBounds(0,0,50*(1+ix/50));
	 mx=1.2*Math.max(w1,w2);
	//plot[2].setFixedBounds(1,-mx,mx);
	
	//************************************************
	
	
	
	//************************************************
	
	
	  w1=-1e10;
	  w2=1e-10;
	 
	XY=new double[ix][2];
	for(int j=0;j<ix;j++){

		XY[j][0]=j;
		XY[j][1]=Vn.el[j];
		if(abs(XY[j][1])>w1) w1=abs(XY[j][1]);
	}
	plot[3].changePlotData(0, XY);
	
	plot[3].setFixedBounds(0,0,50*(1+ix/50));

	 mx=1.2*Math.max(w1,w2);
	 mx=400;
	//plot[3].setFixedBounds(1,-mx,mx);
	
	//************************************************
	

}
}
