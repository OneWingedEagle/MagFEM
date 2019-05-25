package triangulation;

import java.awt.Point;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.TreeSet;

import fem.Model;
import math.Vect;
import math.util;
import triangulation.DelaunayTriangulationold.GraphEdge;


public class DelaunayTriangulationold
{
	
	public static void main(String[] args) throws IOException {

		String bun=util.getFile(0);

		Model model=new Model(bun);
		
		String folder=new File(model.meshFilePath).getParentFile().getPath();

		int nNodes=model.numberOfNodes/1000;
	nNodes=5;
		Vect[] p = new Vect[nNodes];
		p[0]=new Vect(0,0);
		p[1]=new Vect(5,0);
		p[2]=new Vect(5,5);
		p[3]=new Vect(0,5);
		p[4]=new Vect(2,3);
		double[] x = new double[nNodes];
		double[] y = new double[nNodes];
		double[] z = new double[nNodes];
		for(int i=1;i<=nNodes;i++){
		//	p[i-1]=model.node[i].getCoord();

/*			x[i-1]=(int)(1000000*p[i-1].el[0]);
			y[i-1]=(int)(1000000*p[i-1].el[1]);*/
			x[i-1]=p[i-1].el[0];
			y[i-1]=p[i-1].el[1];
			//p[i-1].x=Math.random();
			//p[i-1].y=Math.random();
			//p[i-1]=model.node[i].getCoord();
			//p[i-1].x=model.node[i].getCoord(0);//(int)(1000*i*.001*cos(.1*i));
			//p[i-1].y=model.node[i].getCoord(1);
			//P[i-1].el[1]=Math.random();
			//P[i-1].hshow();
			//util.pr(p[i-1].x+"  "+p[i-1].y);
		}
	
		   DelaunayTriangulationold dt=new    DelaunayTriangulationold(nNodes);
		   
		   TreeSet<GraphEdge> treeSet=dt.getEdges(nNodes, x, y, z);
		   
		   util.pr(nNodes);
	
		    int[][] adjMat=dt.getAdj();
		   for (int i = 0; i < adjMat.length; i++) {
			   for (int j = i+1; j < adjMat[0].length; j++) {
				//	if (hull[i] != null)
			 //  GraphEdge edge=set.first();
						System.out.print(adjMat[i][j]+"  ");
						//if(adjMat[i][j]>0)
						//System.out.println("( "+p[i].el[0]+" , "+p[i].el[1]+") <-->  ( "+p[j].el[0]+" , "+p[j].el[1]+" )");
			   }
			   System.out.println();
				}
/*		Vect[] hull =getConvexHull(p);

		for (int i = 0; i < hull.length; i++) {
		//	if (hull[i] != null)
			//	System.out.print(hull[i]);
			
		}
		
		Vect center=new Vect(2);
		for (int i = 0; i < hull.length; i++) {
			center=center.add(hull[i]);
		}
		
		center.timesVoid(1./hull.length);

		Vect[] hull2 = new Vect[hull.length];
		for (int i = 0; i < hull.length; i++) {
			hull2[i]=hull[i].sub(center).times(1.01).add(center);

			//System.out.print(hull2[i]);
		}
*/
		/*int nRegions=1;
		Model hullMesh=new Model(nRegions,1,2*hull.length,"triangle");

		for(int i=1;i<=nRegions;i++){


			hullMesh.region[i].setMaterial("mat"+i);
			hullMesh.region[i].setName("reg"+i);

		}
	//	for(int i=1;i<=hullMesh.numberOfElements;i++)
		for(int i=0;i<hull.length;i++)
			hullMesh.node[i+1].setCoord(hull[i]);
		for(int i=hull.length;i<2*hull.length;i++)
			hullMesh.node[i+1].setCoord(hull2[i-hull.length]);

		for(int i=1;i<=hullMesh.numberOfElements;i++){
			int[] vn=new int[4];
			if(i<hullMesh.numberOfElements){
			vn[0]=i;
			vn[1]=i+1;
			vn[2]=hull.length+i+1;
			vn[3]=hull.length+i;
			}else{
				vn[0]=i;
				vn[1]=1;
				vn[2]=hull.length+1;
				vn[3]=hull.length+i;	
			}
			hullMesh.element[i].setVertNumb(vn);
		}
		
		hullMesh.region[1].setFirstEl(1);
		hullMesh.region[1].setLastEl(hullMesh.numberOfElements);

		
		hullMesh.writeMesh(folder+"\\hull.txt");*/
	}
	
	  class GraphEdge implements Comparable<GraphEdge>{
		  GraphPoint points[];  
		  GraphEdge(){
			  points=new GraphPoint[2];
		  }
		  GraphEdge(GraphPoint p1,GraphPoint p2){
			  points=new GraphPoint[2];
			  points[0]=p1;
			  points[1]=p2;
		  }
		@Override
		public int compareTo(GraphEdge o) {
			// TODO Auto-generated method stub
			return 0;
		}
	 }

	  class GraphPoint implements Comparable<GraphPoint>{
		//  class GraphPoint extends Point{
			  double x,y;
			  GraphPoint(double x, double y){
					this.x=x;
					this.y=y; 
			  }
			@Override
			public int compareTo(GraphPoint p) {
				// TODO Auto-generated method stub
				double r=0;
				if (this.x == p.x) {
					r=this.y- p.y;

				} else {
					r=this.x- p.x;
				}
				if(r>0) return 1;
				else if (r==0) return 0;
				else return -1;
			}
		  }



    int[][] adjMatrix;

   DelaunayTriangulationold(int size)
   {
     this.adjMatrix = new int[size][size];
   }
   public int[][] getAdj() {
     return this.adjMatrix;
   }

   public TreeSet<GraphEdge> getEdges(int n, double[] x, double[] y, double[] z)
   {
     TreeSet<GraphEdge> result = new TreeSet<GraphEdge>();

     if (n == 2)
     {
       this.adjMatrix[0][1] = 1;
       this.adjMatrix[1][0] = 1;
       result.add(new GraphEdge(new GraphPoint(x[0], y[0]), new GraphPoint(x[1], y[1])));

       return result;
     }

     for (int i = 0; i < n - 2; i++) {
       for (int j = i + 1; j < n; j++) {
    	//   int k=0;
         for (int k = i + 1; k < n; k++)
        {
          if (j == k) {
            continue;
         }
           double xn = (y[j] - y[i]) * (z[k] - z[i]) - (y[k] - y[i]) * (z[j] - z[i]);

           double yn = (x[k] - x[i]) * (z[j] - z[i]) - (x[j] - x[i]) * (z[k] - z[i]);

           double zn = (x[j] - x[i]) * (y[k] - y[i]) - (x[k] - x[i]) * (y[j] - y[i]);
           boolean flag;
           if (flag = (zn < 0 ? 1 : 0) != 0) {
             for (int m = 0; m < n; m++) {
               flag = (flag) && ((x[m] - x[i]) * xn + (y[m] - y[i]) * yn + (z[m] - z[i]) * zn <= 0);
             }

           }

           if (!flag)
           {
        	   //System.out.println(i+" "+j+" "+" continue ");
             continue;
           }

           GraphPoint p1=new GraphPoint(x[i], y[i]);
           GraphPoint p2=new GraphPoint(x[j], y[j]);
          // util.pr(p1.compareTo(p2));
           result.add(new GraphEdge(p1, p2));
           //System.out.println("----------");
           //System.out.println(x[i]+" "+ y[i] +"----"+x[j]+" "+y[j]);

          result.add(new GraphEdge(new GraphPoint(x[j], y[j]), new GraphPoint(x[k], y[k])));
          //System.out.println(x[j]+" "+ y[j] +"----"+x[k]+" "+y[k]);
          result.add(new GraphEdge(new GraphPoint(x[k], y[k]), new GraphPoint(x[i], y[i])));
           //System.out.println(x[k]+" "+ y[k] +"----"+x[i]+" "+y[i]);
           this.adjMatrix[i][j] = 1;
           this.adjMatrix[j][i] = 1;
           this.adjMatrix[k][i] = 1;
           this.adjMatrix[i][k] = 1;
           this.adjMatrix[j][k] = 1;
           this.adjMatrix[k][j] = 1;
           System.out.println(i+" "+j+" "+" -- "+this.adjMatrix[i][j]);

         }

       }

     }

     return result;
   }
}
